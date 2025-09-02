# ──────────────────────────────────────────────────────────────────────────────
# engine.py —— 直达项 + 散射项（按 r 逐点）
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
from atmo import Atmosphere1D, SpectralData
from raytrace import build_alpha_table, invert_h_from_alpha, default_solar_disk_samples
from transmit import refracted_tangent_path, optical_depth_along_path, beer_transmittance
from scatter import (rayleigh_beta, hg_phase, rayleigh_phase,
                     scattering_angle_to_pixel, incoming_path_segment_to_cell,
                     outgoing_path_segment_to_pixel)
from color import spectrum_to_XYZ, XYZ_to_sRGB_D65

def direct_term_spectrum(r_rad: float, atmo: Atmosphere1D, spec: SpectralData,
                         alpha_table, solar_samples, use_aerosol=False) -> np.ndarray:
    """直达折射项：⟨E⊙·T_dir⟩_solar_disk  → 返回 Ein_dir(λ)"""
    acc = np.zeros_like(spec.lam, dtype=float)
    wsum = 0.0
    for u in solar_samples:
        r_eff = r_rad - u.delta_alpha_rad
        h_star = invert_h_from_alpha(alpha_table, r_eff)
        if h_star is None:
            continue
        seg = refracted_tangent_path(atmo, h_star)
        tau = optical_depth_along_path(seg, spec)         # dict of components
        T = beer_transmittance(tau)                       # np.ndarray over λ
        acc += u.weight * spec.Esun_AM0 * T
        wsum += u.weight
    return acc/wsum if wsum>0 else acc

def scattered_term_spectrum(r_rad: float, atmo: Atmosphere1D, spec: SpectralData,
                            solar_samples, phi_grid: np.ndarray, z_grid: np.ndarray,
                            use_aerosol=False) -> np.ndarray:
    """单次散射体积分：φ×z 上累积  → 返回 Ein_sca(λ)"""
    acc = np.zeros_like(spec.lam, dtype=float)
    dphi = float(phi_grid[1] - phi_grid[0])
    dz   = float(z_grid[1] - z_grid[0])
    for phi in phi_grid:
        for z in z_grid:
            # 入射源项 Jλ (太阳盘平均)
            J = np.zeros_like(spec.lam, dtype=float); wsum = 0.0
            for u in solar_samples:
                seg_in = incoming_path_segment_to_cell(z, phi, u)
                tau_in = optical_depth_along_path(seg_in, spec)
                T_in   = beer_transmittance(tau_in)
                J += u.weight * spec.Esun_AM0 * T_in
                wsum += u.weight
            if wsum == 0: 
                continue
            J /= wsum

            # 体散射系数 β 和相函数 P(θ)
            N_air = np.interp(z, atmo.z, atmo.air_number_density)
            beta_R = rayleigh_beta(N_air, spec.sigma_R)               # [1/m] vs λ
            theta  = scattering_angle_to_pixel(z, phi, r_rad)
            P = rayleigh_phase(theta)                                 # start with Rayleigh only

            # 出射消光
            seg_out = outgoing_path_segment_to_pixel(z, phi, r_rad)
            tau_out = optical_depth_along_path(seg_out, spec)
            T_out   = beer_transmittance(tau_out)

            acc += beta_R * J * P * T_out * dz * dphi
    return acc

def spectrum_to_srgb(Llam: np.ndarray, spec: SpectralData) -> tuple[float,float,float]:
    X, Y, Z = spectrum_to_XYZ(spec.lam, Llam, spec.cie_x, spec.cie_y, spec.cie_z)
    return XYZ_to_sRGB_D65(X, Y, Z)
