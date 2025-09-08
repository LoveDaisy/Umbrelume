# ──────────────────────────────────────────────────────────────────────────────
# scatter.py —— 单次散射体积分（φ × z）
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
from atmo import Atmosphere1D, SpectralData, AerosolProfile
from transmit import optical_depth_along_path, PathSegment


def rayleigh_beta(air_number_density_m3: float, sigma_R_m2: np.ndarray) -> np.ndarray:
    """β_R(λ) = N_air * σ_R(λ)  [1/m]"""
    return air_number_density_m3 * sigma_R_m2


def hg_phase(theta_rad: float, g: float) -> float:
    """Henyey-Greenstein phase function (normalized over 4π).
    P(θ)= (1 - g^2) / (4π * (1 + g^2 - 2g cosθ)^(3/2))"""
    ct = np.cos(theta_rad)
    return (1.0 - g*g) / (4.0*np.pi * (1 + g*g - 2*g*ct)**1.5)


def rayleigh_phase(theta_rad: float) -> float:
    """Rayleigh phase normalized over 4π: (3/16π)(1 + cos^2θ)."""
    return (3.0/(16.0*np.pi))*(1.0 + np.cos(theta_rad)**2)


def scattering_angle_to_pixel(z_km: float, phi_rad: float, r_rad: float) -> float:
    """Geometry to get Θ between (z,φ)→pixel direction.
    TODO: implement small-angle approximation or full geometry."""
    raise NotImplementedError


def incoming_path_segment_to_cell(z_km: float, phi_rad: float, sample_u) -> PathSegment:
    """Build path from Sun sample u to cell (z,φ). Can approximate by slant through limb."""
    raise NotImplementedError


def outgoing_path_segment_to_pixel(z_km: float, phi_rad: float, r_rad: float) -> PathSegment:
    """Build path from cell (z,φ) to pixel direction (r)."""
    raise NotImplementedError
