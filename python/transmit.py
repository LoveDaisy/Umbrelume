# ──────────────────────────────────────────────────────────────────────────────
# transmit.py —— 折射切向路径重构 + 沿程光学厚度
# ──────────────────────────────────────────────────────────────────────────────
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from typing import Dict

from atmo import Atmosphere1D, SpectralData, ciddor_n_air

# 建议：几何地球半径用于光路几何；USSA76 内部用的 6356.766 km 是“位势地球半径”，
# 这里采用常用几何半径 6371.0 km 即可（差异 << 1%）。
RE_KM = 6371.0
# 参考波长用于几何折射（同 α(h) 表的 λ0，保持一致）
LAMBDA0_NM = 550.0


@dataclass
class PathSegment:
    """
    折射切向路径的离散化表示（已经包含“入+出”两侧对称的总 ds）
    """
    z_km: np.ndarray        # 取每段的中点高度（km）
    ds_m: np.ndarray        # 该段对应的总路径长度（m），已包含两侧对称×2
    N_air_m3: np.ndarray    # 空气数密度 (m^-3) 在这些中点的插值
    O3_m3: np.ndarray       # 臭氧数密度 (m^-3) 在这些中点的插值


def _build_layer_edges_from_hstar(atmo: Atmosphere1D, h_star_km: float,
                                  z_top_km: float | None = None) -> np.ndarray:
    """
    用大气网格构造从 h* 到顶的“层边界”数组（km），第一段从 h* 开始。
    """
    z = np.asarray(atmo.z, float)
    if z_top_km is None:
        z_top_km = float(z.max())
    # 只取 >= h* 的段
    z_sel = z[(z >= h_star_km) & (z <= z_top_km)]
    if len(z_sel) == 0 or z_sel[0] > h_star_km + 1e-9:
        # 把 h* 插入为首个边界
        z_edges = np.concatenate([[h_star_km], z_sel])
    else:
        z_edges = z_sel.copy()
        z_edges[0] = h_star_km  # 保证首个边界正是 h*
    # 若顶端不是最后一个网格点，补上 z_top
    if z_edges[-1] < z_top_km - 1e-9:
        z_edges = np.concatenate([z_edges, [z_top_km]])
    return z_edges


def refracted_tangent_path(atmo: Atmosphere1D, h_star_km: float,
                           z_top_km: float | None = None,
                           n_lambda_nm: float = LAMBDA0_NM,
                           cospsi_floor: float = 1e-6) -> PathSegment:
    """
    生成通过切向高度 h* 的折射路径在分层大气中的离散段：
      - 在每个薄层 [z_i, z_{i+1}]，用折射不变量 μ = n(h*) * r_t 计算 cosψ(mid)，
        近似 ds_half = Δr / cosψ(mid)，对称两侧乘以 2 得 ds_total。
      - 返回各段的中点 z_mid、对应 ds、以及在中点插值的 N_air 与 O3。
    参考：
      • 折射不变量 n r sinψ 与非奇异积分法（Auer & Standish 2000）。 [oai_citation:5‡Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/2000AJ....119.2472A/abstract?utm_source=chatgpt.com) [oai_citation:6‡Green Flash Page](https://aty.sdsu.edu/explain/atmos_refr/models/Oriani.html?utm_source=chatgpt.com)
    """
    z_edges = _build_layer_edges_from_hstar(atmo, h_star_km, z_top_km)
    # 中点高度
    z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])
    # 半径（m）
    r_mid = (RE_KM + z_mid) * 1000.0
    r_t = (RE_KM + h_star_km) * 1000.0

    # 折射率：在参考 λ0 对整个剖面计算 n(z)
    n_star = float(ciddor_n_air(n_lambda_nm,
                                T=np.interp(h_star_km, atmo.z, atmo.T),
                                p=np.interp(h_star_km, atmo.z, atmo.p)))
    # 在中点的 n(z)
    n_mid = ciddor_n_air(n_lambda_nm,
                         T=np.interp(z_mid, atmo.z, atmo.T),
                         p=np.interp(z_mid, atmo.z, atmo.p))

    # 折射不变量 μ = n(h*) * r_t  (切向处 sinψ=1)
    mu = n_star * r_t

    # cosψ(mid) = sqrt(1 - (μ/(n*r))^2)；数值上做安全钳制
    ratio = mu / (n_mid * r_mid)
    ratio2 = np.minimum(ratio * ratio, 1.0 - 1e-15)
    cospsi = np.sqrt(1.0 - ratio2)
    cospsi = np.maximum(cospsi, cospsi_floor)  # 避免极端薄层处发散

    # Δr（m）
    r_edges = (RE_KM + z_edges) * 1000.0
    dr = np.diff(r_edges)

    # 单侧 ds（半程）≈ Δr / cosψ(mid)；两侧对称乘以 2
    ds_total = 2.0 * (dr / cospsi)

    # 在中点插值密度
    N_air = np.interp(z_mid, atmo.z, atmo.air_number_density)
    O3 = np.interp(z_mid, atmo.z, atmo.o3_number_density)

    return PathSegment(
        z_km=z_mid.astype(float),
        ds_m=ds_total.astype(float),
        N_air_m3=N_air.astype(float),
        O3_m3=O3.astype(float),
    )


def optical_depth_along_path(seg: PathSegment, spec: SpectralData) -> Dict[str, np.ndarray]:
    """
    计算该路径的多组分光学厚度：
      τ_R(λ)  = σ_R(λ)  * ∑_j N_air(z_j) * ds_j
      τ_O3(λ) = σ_O3(λ) * ∑_j n_O3(z_j)  * ds_j
    目前未含气溶胶项；后续可加入：τ_aer(λ) = ∑ β_ext(λ,z_j)*ds_j（或 SSA/相函数拆分）。
    """
    # 路径积分的“柱含量”
    S_air = float(np.sum(seg.N_air_m3 * seg.ds_m))   # [m^-2]
    S_o3 = float(np.sum(seg.O3_m3 * seg.ds_m))   # [m^-2]

    tau_R = spec.sigma_R * S_air
    tau_O3 = spec.sigma_O3 * S_o3
    # 可扩展：tau_aer = spec.sigma_aer * S_aer 或分光谱 β_ext(z,λ) 积分
    return {"R": tau_R, "O3": tau_O3}
