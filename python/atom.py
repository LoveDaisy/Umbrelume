# ──────────────────────────────────────────────────────────────────────────────
# atmo.py  —— 大气与常数、数据接口
# ──────────────────────────────────────────────────────────────────────────────
from dataclasses import dataclass
import numpy as np

@dataclass
class Atmosphere1D:
    z: np.ndarray         # [K] altitude grid, km
    T: np.ndarray         # [K] temperature
    p: np.ndarray         # [Pa] pressure
    air_number_density: np.ndarray  # [m^-3] from p, T
    o3_number_density: np.ndarray   # [m^-3] simple ozone profile (can be φ-invariant first)

@dataclass
class AerosolProfile:
    beta_ext: callable    # beta_ext(λ, z) -> [1/m]
    ssa: callable         # single-scattering albedo ssa(λ, z) -> [0..1]
    g: callable           # HG parameter g(λ, z) -> [0..1)

@dataclass
class SpectralData:
    lam: np.ndarray                 # [nm]
    Esun_AM0: np.ndarray            # [W/m^2/nm]
    sigma_R: np.ndarray             # [m^2] molecular Rayleigh xs (per molecule)
    sigma_O3: np.ndarray            # [m^2] ozone xs at given T (interpolate if needed)
    cie_x: np.ndarray               # CIE 1931 2° CMFs
    cie_y: np.ndarray
    cie_z: np.ndarray

def us_standard_atmosphere_1976(z_grid_km: np.ndarray) -> Atmosphere1D:
    """Build a simple 1D USSA76 profile on z_grid.
    TODO: implement hydrostatic T(z), p(z) → number densities. (USSA76 tables)"""
    raise NotImplementedError

def ciddor_n_air(wavelength_nm: float, T: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Refractive index of air n(λ, T, p) at standard CO2 / dry air.
    TODO: implement or call a known param (Ciddor/Edlén)."""
    raise NotImplementedError

def load_astm_e490() -> tuple[np.ndarray, np.ndarray]:
    """Return (lam_nm, AM0_solar_W_m2_nm).
    TODO: load from file / resource; unit consistency check."""
    raise NotImplementedError

def load_cie_1931_2deg() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return CIE 1931 2° CMFs sampled to your λ grid (interpolate if needed)."""
    raise NotImplementedError

def rayleigh_sigma_bodhaine(lam_nm: np.ndarray) -> np.ndarray:
    """Return molecular Rayleigh cross-section per molecule σ_R(λ).
    TODO: implement Bodhaine 1999 parameterization."""
    raise NotImplementedError

def ozone_sigma_serdyuchenko(lam_nm: np.ndarray, T_K: float=293.0) -> np.ndarray:
    """Return ozone absorption cross-section σ_O3(λ, T) on λ grid.
    TODO: load AMT 2014 dataset or subset for visible band."""
    raise NotImplementedError
