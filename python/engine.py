# ──────────────────────────────────────────────────────────────────────────────
# engine.py —— 月食颜色引擎: 直达折射 + 单次散射 (含降级近似)
# ──────────────────────────────────────────────────────────────────────────────
from __future__ import annotations
from dataclasses import dataclass
import os
import numpy as np

from atmo import Atmosphere1D, SpectralData
from raytrace import build_alpha_table, invert_h_from_alpha, default_solar_disk_samples, AlphaTable
from raytrace import invert_h_from_alpha  # ensure available
from color import spectrum_to_XYZ, XYZ_to_sRGB_D65
from transmit import refracted_tangent_path, optical_depth_along_path, PathSegment
from scatter import (build_sightline_segments, incoming_path_segment_to_cell, outgoing_path_segment_to_pixel,
                     rayleigh_phase, hg_phase)

# -------------------- 配置与数据类 --------------------
HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
OUT_DIR = os.path.join(HERE, "..", "out")


@dataclass
class EngineOptions:
    use_solar_disk: bool = False   # False=点太阳；True=有限太阳盘 (需要 K>1 的采样)
    lunar_albedo: float = 1.0      # 先当常数 (后续可接 ROLO)
    hg_g: float = 0.0              # 气溶胶 HG 参数；0 表示不用 HG（仅瑞利）
    use_aerosol: bool = False      # 先只做瑞利散射
    r_scale_arcmin: float = 57.0   # 固定的“物理刻度”：本影边界(~41′)+太阳半径(~16′)=~57′


# --------- r→α 的“反向映射”工具 ----------
def _alpha_required_for_r(r_rad: float, r_scale_rad: float, alpha_table) -> float:
    """
    把 r∈[0,r_scale] 单调映射到 α∈[α_max, α_min]：
    - r 越小 => α 越大（更深折射）；
    - r 超过 r_scale 时按 r_scale 夹紧。
    """
    a = np.asarray(alpha_table["alpha_rad"] if isinstance(alpha_table, dict) else alpha_table.alpha_rad)
    a_min, a_max = float(np.min(a)), float(np.max(a))
    t = 0.0 if r_scale_rad <= 0.0 else (r_rad / r_scale_rad)
    t = min(max(t, 0.0), 1.0)  # clip 到 [0,1]
    return a_max - t * (a_max - a_min)


# -------------------- 直达折射项 --------------------
def direct_term_spectrum(r_rad: float, atmo, spec, alpha_table, solar_samples, opt: EngineOptions) -> np.ndarray:
    """
    Ein_dir(λ) = ⟨ E⊙(λ) · T_dir(λ; h*(α(r))) ⟩_太阳盘
    这里 α(r) 的归一化只依赖 options.r_scale_arcmin，不依赖 r_arcmin 的最大值。
    """
    r_scale_rad = np.deg2rad(float(opt.r_scale_arcmin) / 60.0)
    acc = np.zeros_like(spec.lam, dtype=float)
    wsum = 0.0
    for u in solar_samples:
        # 依据固定刻度得到所需偏折角
        alpha_req = _alpha_required_for_r(r_rad, r_scale_rad, alpha_table)
        # 太阳盘样本的小修正（若 use_solar_disk=True）
        alpha_req -= getattr(u, "delta_alpha_rad", 0.0)

        h_star = invert_h_from_alpha(alpha_table, alpha_req)
        if h_star is None:
            continue

        seg = refracted_tangent_path(atmo, h_star)
        tau = optical_depth_along_path(seg, spec)
        T = np.exp(-sum(tau.values()))
        w = getattr(u, "weight", 1.0)
        acc += w * spec.Esun_AM0 * T
        wsum += w
    return acc / wsum if wsum > 0 else acc


# -------------------- 单次散射项 (瑞利起步)  --------------------
def scattered_term_spectrum(r_rad: float, atmo, spec,
                            solar_samples, phi_grid, z_grid,
                            opt) -> np.ndarray:
    """
    Ein_sca(λ) = ∫_LOS β_R(z,λ) · ⟨E⊙(λ)·T_in(z,λ)⟩ · P(θ≈r) · T_out(z→TOA,λ;r) · ds
    —— 沿“像素视线”做一维积分，不再对 φ 做环向求和。
    """
    # 1) 视线路径（近端半程）
    sl = build_sightline_segments(atmo, r_rad)
    if sl.ds_half_m.size == 0:
        return np.zeros_like(spec.lam, float)

    Nz = sl.z_mid_km.size
    lam = spec.lam
    # 2) 入射：J(z,λ) = 太阳盘平均 [E⊙ * T_in]（与 r 无关，可缓存）
    #    这里为不同 z_mid 逐点构造“顶→该层”的入射斜程
    Jz = np.zeros((Nz, lam.size), float)
    wsum = sum(getattr(u, "weight", 1.0) for u in solar_samples) or 1.0
    for k, u in enumerate(solar_samples):
        delta = getattr(u, "delta_alpha_rad", 0.0)
        for i, z in enumerate(sl.z_mid_km):
            seg_in = incoming_path_segment_to_cell(atmo, float(z), 0.0, sun_zenith_offset_rad=delta)
            tau_in = optical_depth_along_path(seg_in, spec)   # {'R','O3'}
            T_in = np.exp(-sum(tau_in.values()))
            Jz[i, :] += getattr(u, "weight", 1.0) * (spec.Esun_AM0 * T_in)
    Jz /= wsum

    # 3) 体散射系数 β_R(z,λ) 与 总消光 κ(z,λ)（用于出射透过）
    beta = (sl.N_air_m3[:, None]) * spec.sigma_R[None, :]             # (Nz,Nλ)
    kappa = (sl.N_air_m3[:, None]) * spec.sigma_R[None, :] + \
            (sl.O3_m3[:,  None]) * spec.sigma_O3[None, :]             # (Nz,Nλ)

    # 4) 出射：从每个层点到近端顶界的透过（沿该视线半程）
    #    τ_out[i,:] = sum_{j=i..end} κ[j,:] * ds_half[j]
    ds = sl.ds_half_m[:, None]                                        # (Nz,1)
    tau_cum = np.cumsum((kappa * ds)[::-1, :], axis=0)[::-1, :]       # (Nz,Nλ)
    T_out = np.exp(-tau_cum)

    # 5) 相函数（小角近似 θ≈r）
    theta = float(r_rad)
    P = hg_phase(theta, opt.hg_g) if opt.hg_g != 0.0 else rayleigh_phase(theta)

    # 6) 线积分
    integrand = beta * Jz * P * T_out * ds     # (Nz,Nλ)
    return np.sum(integrand, axis=0)


# -------------------- 顶层: 把谱 → 颜色条 --------------------


def compute_strip(r_arcmin: np.ndarray,
                  atmo, spec, alpha_table,
                  solar_samples=None,
                  options: EngineOptions | None = None,
                  phi_bins: int = 720,
                  z_min_km: float = 10.0, z_max_km: float = 40.0, dz_km: float = 1.0):
    """
    仍然保留 r_arcmin：它是要绘制的角半径采样。
    但内部 r→α 的归一化不再看 r_arcmin[-1]，而是用 options.r_scale_arcmin。
    """
    opt = options or EngineOptions()
    if solar_samples is None:
        # 点太阳: K=1；有限太阳盘: K=21（你也可以调整）
        solar_samples = default_solar_disk_samples(K=1 if not opt.use_solar_disk else 21)

    # φ 在当前几何近似下与散射项弱相关，这里用占位以兼容接口
    phi_grid = np.array([0.0])
    z_grid = np.arange(z_min_km, z_max_km + 1e-9, dz_km)

    rgb_list, XYZ_list, Lspec_list = [], [], []
    diag_rows = []  # 用于调试的诊断表格
    for r_am in r_arcmin:
        r_rad = np.deg2rad(r_am / 60.0)

        # 直达项
        Ein_dir = direct_term_spectrum(r_rad, atmo, spec, alpha_table, solar_samples, opt)

        # 散射项（你已有的矢量化实现）
        Ein_sca = scattered_term_spectrum(r_rad, atmo, spec, solar_samples, phi_grid, z_grid, opt)

        Llam = (Ein_dir + Ein_sca) * opt.lunar_albedo  # [W m^-2 nm^-1 sr^-1]（相对量）
        X, Y, Z = spectrum_to_XYZ(spec.lam, Llam, spec.cie_x, spec.cie_y, spec.cie_z)
        rgb = XYZ_to_sRGB_D65(X, Y, Z)

        rgb_list.append(rgb)
        XYZ_list.append([X, Y, Z])
        Lspec_list.append(Llam)

        Y_dir = spectrum_to_XYZ(spec.lam, Ein_dir, spec.cie_x, spec.cie_y, spec.cie_z)[1]
        Y_sca = spectrum_to_XYZ(spec.lam, Ein_sca, spec.cie_x, spec.cie_y, spec.cie_z)[1]
        diag_rows.append([r_am, Y_dir, Y_sca])

    np.savetxt(os.path.join(OUT_DIR, "diagnose_Y_dir_sca.csv"), np.array(
        diag_rows), delimiter=",", header="r_arcmin,Y_dir,Y_sca", comments="")

    rgb_arr = np.vstack(rgb_list)
    XYZ_arr = np.vstack(XYZ_list)
    Lspec = np.vstack(Lspec_list)
    return {"rgb": rgb_arr, "XYZ": XYZ_arr, "Lspec": Lspec}
