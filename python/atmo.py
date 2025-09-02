# ──────────────────────────────────────────────────────────────────────────────
# atmo.py  —— 大气与常数、数据接口（已实现 USSA76 与 Bodhaine σ_R）
# ──────────────────────────────────────────────────────────────────────────────
from dataclasses import dataclass
import numpy as np

@dataclass
class Atmosphere1D:
    z: np.ndarray         # [km] geometric altitude grid
    T: np.ndarray         # [K] temperature
    p: np.ndarray         # [Pa] pressure
    air_number_density: np.ndarray  # [m^-3] from p, T (dry air)
    o3_number_density: np.ndarray   # [m^-3] ozone profile (here: zeros placeholder)

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
    """
    Build COESA 'U.S. Standard Atmosphere 1976' (0–86 km) on a geometric altitude grid.

    Implements the layer-wise geopotential formulation:
      - Convert geometric altitude z -> geopotential h = Re * z / (Re + z)
      - For each layer b with base (h_b, T_b, p_b, lapse L_b), compute:
          If L_b != 0:
            T = T_b + L_b * (h - h_b)
            p = p_b * (T / T_b) ** (-g0 * M / (R* * L_b))
          If L_b == 0:
            T = T_b
            p = p_b * exp(-g0 * M * (h - h_b) / (R* * T_b))
      - Density from ideal gas: ρ = p / (R_specific * T);
        Number density N = p / (k_B * T)

    Units:
      z, h in km (internally converted to m where needed); L_b in K/m.

    Constants (USSA76 / SI):
      Re = 6356.766 km (Earth radius for geopotential)  [COESA 1976]
      g0 = 9.80665 m/s^2
      R* = 8.31432 J/(mol·K)
      M  = 28.9644e-3 kg/mol (dry air mean molecular weight)
      kB = 1.380649e-23 J/K

    Layer bases (geopotential, km), lapse rates (K/km), and base temperatures (K) per USSA76:
      h_b: [0, 11, 20, 32, 47, 51, 71, 84.852]
      L_b: [-6.5, 0, +1.0, +2.8, 0, -2.8, -2.0]  (K/km)
      T_b: [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946]
      p_b[0] = 101325 Pa; other p_b computed recursively from layer formulas.

    References:
      - COESA / NOAA “U.S. Standard Atmosphere, 1976” formulas & layer limits. :contentReference[oaicite:0]{index=0}
      - Quick lookup of layer table / lapse rates (cross-check). :contentReference[oaicite:1]{index=1}
    """
    z = np.asarray(z_grid_km, dtype=float)
    # constants
    Re_km = 6356.766
    g0 = 9.80665
    Rstar = 8.31432
    M = 28.9644e-3
    kB = 1.380649e-23

    # geometric -> geopotential (km)
    h = (Re_km * z) / (Re_km + z)
    h_m = h * 1000.0

    # Layer bases (geopotential, km)
    hb = np.array([0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852])
    # Lapse rates [K/km] -> [K/m]
    Lb_km = np.array([-6.5, 0.0, +1.0, +2.8, 0.0, -2.8, -2.0])
    Lb = Lb_km / 1000.0
    # Base temperatures [K]
    Tb = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946])
    # Base pressures [Pa] (p_b[0] known; others by recurrence)
    pb = np.zeros_like(Tb)
    pb[0] = 101325.0
    for i in range(len(hb) - 1):
        if Lb[i] == 0.0:
            dh = (hb[i + 1] - hb[i]) * 1000.0
            pb[i + 1] = pb[i] * np.exp(-g0 * M * dh / (Rstar * Tb[i]))
        else:
            Tb_next = Tb[i + 1]
            pb[i + 1] = pb[i] * (Tb_next / Tb[i]) ** (-g0 * M / (Rstar * Lb[i]))

    # For each altitude, pick its layer and compute T, p
    T = np.empty_like(h_m)
    p = np.empty_like(h_m)
    for idx, hm in enumerate(h_m):
        # clamp to model top
        if hm/1000.0 >= hb[-1]:
            ib = len(hb) - 2  # use top layer
            hm = min(hm, hb[-1] * 1000.0 - 1e-6)
        else:
            ib = np.searchsorted(hb * 1000.0, hm, side='right') - 1
            ib = np.clip(ib, 0, len(hb) - 2)
        h0 = hb[ib] * 1000.0
        L = Lb[ib]
        T0 = Tb[ib]
        p0 = pb[ib]
        if L == 0.0:
            T[idx] = T0
            p[idx] = p0 * np.exp(-g0 * M * (hm - h0) / (Rstar * T0))
        else:
            T[idx] = T0 + L * (hm - h0)
            p[idx] = p0 * (T[idx] / T0) ** (-g0 * M / (Rstar * L))

    # Number density (dry air)
    air_number_density = p / (kB * T)

    # Placeholder ozone number density (all zeros; plug real profile later)
    o3_number_density = np.zeros_like(air_number_density)

    return Atmosphere1D(
        z=z,
        T=T,
        p=p,
        air_number_density=air_number_density,
        o3_number_density=o3_number_density,
    )

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
    """
    Molecular Rayleigh scattering cross-section per molecule σ_R(λ) for *dry air*
    following Bodhaine et al. (1999) recommended approach:

      σ(λ) = (24 * π^3 / (N_s^2 * λ^4)) * [ (n_s(λ)^2 - 1)^2 / (n_s(λ)^2 + 2)^2 ] * F_air(λ)

    where:
      - N_s = molecular number density at STP (dry air, 1013.25 hPa, 288.15 K)
              ≈ 2.5469×10^25 m^-3 (2.5469×10^19 cm^-3).  [Colour/ CODATA]
      - n_s(λ): refractive index of standard dry air incl. CO2 correction (Bodhaine-style).
      - F_air(λ): King depolarization factor for air (mixture of N2, O2, Ar, CO2).

    Inputs:
      lam_nm : wavelength(s) in nanometers

    Return:
      σ_R(λ) in m^2 per molecule

    References:
      - Core σ formula / King factor discussion (Bodhaine 1999). :contentReference[oaicite:2]{index=2}
      - n_s(λ) parameterization with CO2 correction (λ in µm). :contentReference[oaicite:3]{index=3}
      - Mixture King factor using Bates (1984) components as compiled/used in later ATBDs. :contentReference[oaicite:4]{index=4}
      - Standard molecular density N_s at STP (288.15 K). :contentReference[oaicite:5]{index=5}
    """
    lam_nm = np.asarray(lam_nm, dtype=float)
    lam_um = lam_nm * 1e-3
    lam_m = lam_nm * 1e-9

    # 1) Refractive index of standard dry air (Bodhaine-style, λ in µm)
    #    ns = 1 + 1e-8 * [8060.51 + 2480990/(132.274 - 1/λ^2) + 17455.7/(39.32957 - 1/λ^2)]
    #    × [1 + 0.54(CCO2 - 0.0003)] ; here take CCO2=300 ppm → correction factor = 1
    inv_l2 = 1.0 / (lam_um ** 2)
    ns = 1.0 + 1e-8 * (
        8060.51
        + 2480990.0 / (132.274 - inv_l2)
        + 17455.7 / (39.32957 - inv_l2)
    )
    # (optional CO2 correction if you later expose it)
    # CCO2 = 0.0003  # 300 ppm
    # ns *= (1.0 + 0.54 * (CCO2 - 0.0003))

    # 2) King factor for air mixture F_air(λ) (λ in µm)
    #    Component King factors (Bates 1984; used e.g. in AerGOM / GOMOS ATBD):
    #      F_N2 = 1.034 + 3.17e-4 / λ^2
    #      F_O2 = 1.096 + 1.385e-3 / λ^2 + 1.448e-4 / λ^4
    #      F_Ar = 1.0
    #      F_CO2 = 1.15
    ln2 = 1.0 / (lam_um ** 2)
    ln4 = ln2 ** 2
    F_N2 = 1.034 + 3.17e-4 * ln2
    F_O2 = 1.096 + 1.385e-3 * ln2 + 1.448e-4 * ln4
    F_Ar = 1.0
    F_CO2 = 1.15

    # Volume fractions in % by volume (CO2 at 300 ppm → 0.03%)
    C_N2 = 78.084
    C_O2 = 20.946
    C_Ar = 0.934
    C_CO2 = 0.03

    F_air = (C_N2 * F_N2 + C_O2 * F_O2 + C_Ar * F_Ar + C_CO2 * F_CO2) / (C_N2 + C_O2 + C_Ar + C_CO2)

    # 3) Loschmidt number (STP, 288.15 K) in m^-3
    N_s = 2.546902e25  # consistent with 2.5469e19 cm^-3  (Colour docs / CODATA)

    # 4) Rayleigh σ(λ) [m^2] per molecule
    ns2 = ns ** 2
    term = ((ns2 - 1.0) / (ns2 + 2.0)) ** 2
    sigma = (24.0 * np.pi ** 3 / (N_s ** 2 * lam_m ** 4)) * term * F_air

    return sigma
