# ──────────────────────────────────────────────────────────────────────────────
# scatter.py —— 单次散射体积分所需: 相函数 + 入/出斜程路径几何
# ──────────────────────────────────────────────────────────────────────────────
from __future__ import annotations
import numpy as np
from dataclasses import dataclass

from atmo import Atmosphere1D, SpectralData
from transmit import PathSegment, optical_depth_along_path

# 地球几何半径 (用于斜程几何)
RE_KM = 6371.0


@dataclass
class Sightline:
    z_mid_km: np.ndarray   # 每个穿越层的中点高度
    ds_half_m: np.ndarray  # 该层在“靠近月球这一侧”的半程弦长（用于积分权重）
    N_air_m3: np.ndarray
    O3_m3: np.ndarray


def build_sightline_segments(atmo, r_rad: float, z_top_km: float | None = None) -> Sightline:
    """
    给定像素角半径 r（相对于地影轴），在直线光路近似下构造“从最近端顶界到最近端顶界”的视线。
    我们只取**靠近月球这一侧**的半程（closest approach -> 近端顶界），
    因为远侧到月球的透过要再穿过整条路径，贡献极小；这样也避免双计数。
    """
    z = np.asarray(atmo.z, float)
    if z_top_km is None:
        z_top_km = float(z.max())

    r_top = (RE_KM + z_top_km) * 1000.0
    # 直线光路的冲量参数 b：在顶界半径 r_top 处入射天顶角 ~ (90° - r)，sin(ζ_top)=cos r
    b = r_top * np.cos(float(r_rad))
    b = max(0.0, min(b, r_top*(1.0 - 1e-12)))

    # 只保留被视线穿过的壳层：r >= b
    z_edges = z[(RE_KM + z) * 1000.0 >= b]
    if z_edges[0] > 0.0:
        z_edges = np.concatenate([[0.0], z_edges])
    if z_edges[-1] < z_top_km - 1e-9:
        z_edges = np.concatenate([z_edges, [z_top_km]])

    r_edges = (RE_KM + z_edges) * 1000.0
    # s(r) = sqrt(r^2 - b^2)；弦长在半程上的层权重 = Δs
    s_edges = np.sqrt(np.maximum(r_edges*r_edges - b*b, 0.0))
    ds_half = np.diff(s_edges)                      # 只取“近端半程”：closest→近端顶界
    z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])

    N_air = np.interp(z_mid, atmo.z, atmo.air_number_density)
    O3 = np.interp(z_mid, atmo.z, atmo.o3_number_density)

    return Sightline(
        z_mid_km=z_mid.astype(float),
        ds_half_m=ds_half.astype(float),
        N_air_m3=N_air.astype(float),
        O3_m3=O3.astype(float),
    )

# -------------------- 相函数 --------------------


def rayleigh_phase(theta_rad: float) -> float:
    """
    瑞利相函数 (4π 归一) :
        P_R(θ) = (3 / 16π) * (1 + cos^2 θ)
    参考: NMSU PDS 大气教育资源. :contentReference[oaicite:4]{index=4}
    """
    return (3.0 / (16.0 * np.pi)) * (1.0 + np.cos(theta_rad) ** 2)


def hg_phase(theta_rad: float, g: float) -> float:
    """
    Henyey-Greenstein 相函数 (4π 归一) :
        P_HG(θ) = (1 - g^2) / [4π (1 + g^2 - 2g cosθ)^(3/2)]
    """
    ct = np.cos(theta_rad)
    return (1.0 - g * g) / (4.0 * np.pi * (1.0 + g * g - 2.0 * g * ct) ** 1.5)

# -------------------- 散射角 --------------------


def scattering_angle_to_pixel(z_km: float, phi_rad: float, r_rad: float) -> float:
    """
    体元 (z, φ) → 像素 (r) 的散射夹角 Θ. 
    远场小角近似 (r ≤ ~1.3°) 下, Θ ≈ r 已足够刻画颜色随 r 的主趋势. 
    若以后要纳入 φ 的几何细节, 可在此处替换为完整三角关系. 
    """
    return float(r_rad)

# -------------------- 斜程路径构造 (直线光路, 球壳弦长离散)  --------------------


def _slant_path_to_altitude(atmo: Atmosphere1D,
                            z_target_km: float,
                            zenith_rad_at_cell: float,
                            z_top_km: float | None = None) -> PathSegment:
    """
    以体元处的局地天顶角 (相对于径向外法线) 为边界条件, 在“直线光路”近似下, 
    从顶层积分到体元 (入射) , 或从体元到顶层 (出射) , 构造穿越各层的路径长度. 

    几何: 给定目标半径 r_t = R+z_t 与局地天顶角 ζ, 直线的冲量参数 b = r_t * sin ζ. 
    对任意半径 r ≥ r_t, 有在该直线上的弦参数 s(r) = sqrt(r^2 - b^2) (取正根, 朝外为正) . 
    穿越某层 [r_i, r_{i+1}] 的路径长度 |Δs| = | s(r_{i+1}) - s(r_i) |. 

    我们对每个薄层取高度中点 z_mid 做常系数近似: 把该层的数密度 *ds 累加, 
    交给 transmit.optical_depth_along_path 去做谱向量化的 τ(λ). 

    注意: 此函数既可用于“入射” (top→cell) , 也可用于“出射” (cell→top) , 
    只需给定不同的半径序列方向即可 (见下面两个封装) . 
    """
    z = np.asarray(atmo.z, float)
    if z_top_km is None:
        z_top_km = float(z.max())

    # 半径 (米)
    r_t = (RE_KM + float(z_target_km)) * 1000.0
    r_top = (RE_KM + float(z_top_km)) * 1000.0

    # 冲量参数 (直线)
    b = r_t * np.sin(float(zenith_rad_at_cell))
    # 防止数值误差导致 b>r_top
    b = min(b, r_top * (1.0 - 1e-12))

    # 构造半径边界: 从目标层到顶层
    z_edges = np.asarray([z_target_km] + list(z[(z > z_target_km) & (z <= z_top_km)]), float)
    if z_edges[-1] < z_top_km - 1e-9:
        z_edges = np.concatenate([z_edges, [z_top_km]])
    r_edges = (RE_KM + z_edges) * 1000.0
    # 中点
    z_mid = 0.5 * (z_edges[:-1] + z_edges[1:])
    r_mid = (RE_KM + z_mid) * 1000.0

    # s(r) = sqrt(r^2 - b^2)
    # 计算每个边界的 s 值, 然后差分得到各层 |Δs|
    s_edges = np.sqrt(np.maximum(r_edges * r_edges - b * b, 0.0))
    ds = np.diff(s_edges)  # 出射方向: s 递增；入射方向也取正长度
    ds_m = np.abs(ds)      # 路径段长度 (米)

    # 在中点插值数密度
    N_air = np.interp(z_mid, atmo.z, atmo.air_number_density)
    O3 = np.interp(z_mid, atmo.z, atmo.o3_number_density)

    return PathSegment(
        z_km=z_mid.astype(float),
        ds_m=ds_m.astype(float),
        N_air_m3=N_air.astype(float),
        O3_m3=O3.astype(float),
    )

# ---- 体元 ← 太阳盘 (入射)  ----


def incoming_path_segment_to_cell(atmo: Atmosphere1D,
                                  z_km: float,
                                  phi_rad: float,
                                  sun_zenith_offset_rad: float = 0.0) -> PathSegment:
    """
    入射路径: 太阳盘 → 体元 (z, φ). 
    月食几何下, 日光在周缘体元处的局地天顶角 ζ_in ≈ 90° (近地平穿越上层大气) . 
    允许给一个很小的 sun_zenith_offset_rad (默认 0) , 用于将来做太阳盘/方位修正. 
    """
    zenith = 0.5 * np.pi - float(sun_zenith_offset_rad)
    return _slant_path_to_altitude(atmo, z_target_km=float(z_km), zenith_rad_at_cell=zenith)

# ---- 体元 → 像素 (出射)  ----


def outgoing_path_segment_to_pixel(atmo: Atmosphere1D,
                                   z_km: float,
                                   phi_rad: float,
                                   r_rad: float) -> PathSegment:
    """
    出射路径: 体元 (z, φ) → 像素方向 (r) . 
    在月球远场小角近似下, 可把局地天顶角 ζ_out ≈ 90° - r. 
    """
    zenith = max(1e-6, 0.5 * np.pi - float(r_rad))
    return _slant_path_to_altitude(atmo, z_target_km=float(z_km), zenith_rad_at_cell=zenith)

# -------------------- 瑞利体散射系数 --------------------


def rayleigh_beta(air_number_density_m3: float, sigma_R_m2: np.ndarray) -> np.ndarray:
    """
    β_R(λ) = N_air * σ_R(λ)   [1/m]
    """
    return air_number_density_m3 * sigma_R_m2

# -------------------- 单次散射: 体元贡献 (给 engine 使用的细粒度函数)  --------------------


def single_cell_scatter_spectrum(atmo: Atmosphere1D, spec: SpectralData,
                                 z_km: float, phi_rad: float, r_rad: float,
                                 use_hg: bool = False, hg_g: float = 0.0,
                                 sun_zenith_offset_rad: float = 0.0) -> np.ndarray:
    """
    计算一个体元 (z, φ) 对像素 r 的光谱贡献: 
      β_R(λ) * J_λ(z,φ) * P(Θ) * T_out(λ; z→r)
    其中 J_λ 通过“太阳盘→体元”的入射斜程得到 (此处先按点太阳/极小偏移处理) . 
    """
    # 入射斜程 (top → cell)
    seg_in = incoming_path_segment_to_cell(atmo, z_km, phi_rad, sun_zenith_offset_rad)
    tau_in = optical_depth_along_path(seg_in, spec)    # {'R','O3'}
    T_in = np.exp(-sum(tau_in.values()))
    J = spec.Esun_AM0 * T_in                           # 点太阳；若做太阳盘平均, 请外层求均值

    # 体散射系数 (瑞利)
    Nz = float(np.interp(z_km, atmo.z, atmo.air_number_density))
    beta_R = rayleigh_beta(Nz, spec.sigma_R)

    # 相函数
    theta = scattering_angle_to_pixel(z_km, phi_rad, r_rad)
    P = (hg_phase(theta, hg_g) if use_hg and hg_g != 0.0 else rayleigh_phase(theta))

    # 出射斜程 (cell → top)
    seg_out = outgoing_path_segment_to_pixel(atmo, z_km, phi_rad, r_rad)
    tau_out = optical_depth_along_path(seg_out, spec)
    T_out = np.exp(-sum(tau_out.values()))

    return beta_R * J * P * T_out
