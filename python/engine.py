# ──────────────────────────────────────────────────────────────────────────────
# engine.py —— 月食颜色引擎: 直达折射 + 单次散射 (含降级近似)
# ──────────────────────────────────────────────────────────────────────────────
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from atmo import Atmosphere1D, SpectralData
from raytrace import build_alpha_table, invert_h_from_alpha, default_solar_disk_samples, AlphaTable
from raytrace import invert_h_from_alpha  # ensure available
from color import spectrum_to_XYZ, XYZ_to_sRGB_D65
from transmit import refracted_tangent_path, optical_depth_along_path, PathSegment
from scatter import single_cell_scatter_spectrum

# -------------------- 配置与数据类 --------------------


@dataclass
class EngineOptions:
    use_solar_disk: bool = False   # False=点太阳；True=有限太阳盘 (需要 K>1 的采样)
    lunar_albedo: float = 1.0      # 先当常数 (后续可接 ROLO)
    hg_g: float = 0.0              # 气溶胶 HG 参数；0 表示不用 HG（仅瑞利）
    use_aerosol: bool = False      # 先只做瑞利散射

# --------- r→α 的“反向映射”工具 ----------


def _alpha_required_for_r(r_rad: float, r_max_rad: float, alpha_table: AlphaTable) -> float:
    """
    把 r∈[0,r_max] 映射到 α∈[α_max, α_min] （中心最大偏折，外缘最小）。
    不依赖具体几何，只做单调映射以得到“中心暗→外侧亮”的正确趋势。
    """
    a = np.asarray(alpha_table["alpha_rad"] if isinstance(alpha_table, dict) else alpha_table.alpha_rad)
    a_min, a_max = float(np.min(a)), float(np.max(a))
    # 归一化 r∈[0, r_max] -> t∈[0,1]
    t = 0.0 if r_max_rad <= 0 else (r_rad / r_max_rad)
    t = min(max(t, 0.0), 1.0)
    return a_max - t * (a_max - a_min)


# -------------------- 直达折射项 --------------------


def direct_term_spectrum(r_rad: float, atmo: Atmosphere1D, spec: SpectralData,
                         alpha_table: AlphaTable,
                         solar_samples,  # list of SolarDiskSample
                         opt: EngineOptions) -> np.ndarray:
    """
    Ein_dir(λ) = ⟨ E⊙(λ) · T_dir(λ; h*(r,u)) ⟩_太阳盘
    使用 refracted_tangent_path + optical_depth_along_path 计算真实切向折射路径。
    """
    acc = np.zeros_like(spec.lam, dtype=float)
    wsum = 0.0
    # 需要 r_max_rad 来完成 r→α 映射；从样本的“最大 r”传入（见 compute_strip）
    r_max_rad = direct_term_spectrum._r_max_rad  # type: ignore[attr-defined]
    for u in solar_samples:
        alpha_req = _alpha_required_for_r(r_rad, r_max_rad, alpha_table)
        alpha_req -= getattr(u, "delta_alpha_rad", 0.0)  # 太阳盘修正（若启用）
        h_star = invert_h_from_alpha(alpha_table, alpha_req)
        if h_star is None:
            continue
        seg = refracted_tangent_path(atmo, h_star)
        tau = optical_depth_along_path(seg, spec)    # {'R','O3',('aer')}
        T = np.exp(-sum(tau.values()))
        acc += getattr(u, "weight", 1.0) * spec.Esun_AM0 * T
        wsum += getattr(u, "weight", 1.0)
    return acc / wsum if wsum > 0 else acc

# -------------------- 单次散射项 (瑞利起步)  --------------------


def scattered_term_spectrum(r_rad: float, atmo: Atmosphere1D, spec: SpectralData,
                            solar_samples,
                            phi_grid: np.ndarray, z_grid: np.ndarray,
                            opt: EngineOptions) -> np.ndarray:
    """
    Ein_sca(λ) = ∬ β_sca · Jλ · P(Θ) · T_out  dz dφ
    采用 scatter.single_cell_scatter_spectrum() 基于真实斜程路径的单元贡献，
    （当前版本按点太阳；如需太阳盘平均，可在此处对 solar_samples 做显式平均）
    """
    acc = np.zeros_like(spec.lam, dtype=float)
    dphi = float(phi_grid[1] - phi_grid[0]) if len(phi_grid) > 1 else 2*np.pi
    dz = float(z_grid[1] - z_grid[0]) if len(z_grid) > 1 else 1.0

    for phi in phi_grid:
        for z in z_grid:
            # 当前按点太阳：single_cell_scatter_spectrum 内部用入/出真实斜程
            contrib = single_cell_scatter_spectrum(
                atmo, spec,
                z_km=float(z), phi_rad=float(phi), r_rad=float(r_rad),
                use_hg=(opt.hg_g != 0.0), hg_g=opt.hg_g,
                sun_zenith_offset_rad=0.0
            )
            acc += contrib * dz * dphi
    return acc

# -------------------- 顶层: 把谱 → 颜色条 --------------------


def compute_strip(r_arcmin: np.ndarray,
                  atmo: Atmosphere1D, spec: SpectralData,
                  alpha_table: AlphaTable,
                  solar_samples=None,
                  options: EngineOptions | None = None,
                  phi_bins: int = 720,
                  z_min_km: float = 10.0, z_max_km: float = 40.0, dz_km: float = 1.0):
    """
    计算一条 r→颜色的径向条带. 
    返回 dict:
      {
        "r_arcmin": r_arcmin,
        "rgb": [(r,g,b), ...],
        "XYZ": [(X,Y,Z), ...],
        "L_spectra": [ Lλ(r_i), ... ]  # 可选: 用于画色度/亮度剖面
      }
    """
    opt = options or EngineOptions()
    if solar_samples is None:
        solar_samples = default_solar_disk_samples(K=1 if not opt.use_solar_disk else 21)

    phi_grid = np.linspace(0.0, 2*np.pi, phi_bins, endpoint=False)
    z_grid = np.arange(z_min_km, z_max_km + 1e-9, dz_km)

    # 把“最大 r（弧分）”换成弧度，供 r→α 映射使用
    r_max_rad = np.deg2rad(float(r_arcmin[-1]) / 60.0)
    direct_term_spectrum._r_max_rad = r_max_rad   # 以函数属性的方式传入

    rgb_list, XYZ_list, Lspec_list = [], [], []
    for r_am in r_arcmin:
        r_rad = np.deg2rad(r_am / 60.0)

        Ein_dir = direct_term_spectrum(r_rad, atmo, spec, alpha_table, solar_samples, opt)
        Ein_sca = scattered_term_spectrum(r_rad, atmo, spec, solar_samples, phi_grid, z_grid, opt)

        Llam = (Ein_dir + Ein_sca) * opt.lunar_albedo
        X, Y, Z = spectrum_to_XYZ(spec.lam, Llam, spec.cie_x, spec.cie_y, spec.cie_z)
        rgb = XYZ_to_sRGB_D65(X, Y, Z)

        rgb_list.append(rgb)
        XYZ_list.append((X, Y, Z))
        Lspec_list.append(Llam)

    return {
        "r_arcmin": np.asarray(r_arcmin, float),
        "rgb": np.array(rgb_list, float),
        "XYZ": np.array(XYZ_list, float),
        "L_spectra": Lspec_list,
    }
