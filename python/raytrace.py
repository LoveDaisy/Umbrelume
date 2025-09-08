# ──────────────────────────────────────────────────────────────────────────────
# raytrace.py —— 折射几何：α(h) 表 + 反解 h*(r[,u])
# ──────────────────────────────────────────────────────────────────────────────
from dataclasses import dataclass
import numpy as np


@dataclass
class AlphaTable:
    h_km: np.ndarray        # monotonic increasing
    alpha_rad: np.ndarray   # corresponding total deflection at λ0


def build_alpha_table(atmo: Atmosphere1D, lam0_nm: float = 550.0) -> AlphaTable:
    """Compute α(h) using Auer-Standish style ray-bending in stratified n(z).
    Steps:
      1) compute n(z) at lam0 via ciddor_n_air
      2) integrate ray equation through layered atmosphere to get total bending α
    TODO: implement nonsingular integral / piecewise-constant layers."""
    raise NotImplementedError


def invert_h_from_alpha(alpha_table: AlphaTable, r_rad: float) -> float | None:
    """Binary-search / monotone spline invert α(h)=r.
    Return h* in km or None if no solution (e.g., outside range)."""
    a = alpha_table
    if r_rad < a.alpha_rad.min() or r_rad > a.alpha_rad.max():
        return None
    # TODO: replace with robust monotone cubic interpolation + invert
    return float(np.interp(r_rad, a.alpha_rad[::-1], a.h_km[::-1]))


@dataclass
class SolarDiskSample:
    delta_alpha_rad: float  # small angular offset relative to disk center
    weight: float           # includes limb-darkening if enabled


def default_solar_disk_samples(K: int = 1) -> list[SolarDiskSample]:
    """K=1 → point Sun (δ=0, w=1). For finite disk, stratified samples + weights."""
    if K == 1:
        return [SolarDiskSample(0.0, 1.0)]
    # TODO: build ring/quadrature with optional limb darkening
    raise NotImplementedError
