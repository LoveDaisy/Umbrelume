# ──────────────────────────────────────────────────────────────────────────────
# atmo.py  —— 大气与常数、数据接口（已实现 USSA76 与 Bodhaine σ_R）
# ──────────────────────────────────────────────────────────────────────────────
from dataclasses import dataclass
import numpy as np
import os, io, csv, math, pathlib
import pandas as pd

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

# --- Ciddor 1996 折射率（干空气；按密度线性缩放），λ 输入 nm ---
def ciddor_n_air(wavelength_nm: float | np.ndarray,
                 T: np.ndarray,  # [K]
                 p: np.ndarray,  # [Pa]
                 co2_ppm: float = 450.0) -> np.ndarray:
    """
    n(λ, T, p) for dry air using Ciddor (1996) standard-air dispersion with density scaling.
    - Step 1: standard dry air (15 °C, 101325 Pa, ~450 ppm CO2):
        n_std - 1 = 0.05792105/(238.0185 - 1/λ^2) + 0.00167917/(57.362 - 1/λ^2)
        (λ in µm; result is absolute refractivity, not ×1e-6)
    - Step 2: scale by number density ~ p/T (ideal-gas) to get n(λ,T,p):
        n - 1 ≈ (n_std - 1) * (p/101325) * (288.15/T)
      NOTE: This omits humidity/compressibility corrections; for our月食折射用途误差在 ppm 量级。
    References: Ciddor 1996; NIST EMT toolbox notes; refractiveindex.info Ciddor page. 
    """
    lam_nm = np.asarray(wavelength_nm, dtype=float)
    lam_um = lam_nm * 1e-3
    s2 = 1.0 / (lam_um**2)                      # [µm^-2]
    # Ciddor standard-air dispersion (dry air, 15 °C, 101325 Pa, ~450 ppm CO2)
    n_minus_1_std = (0.05792105 / (238.0185 - s2)) + (0.00167917 / (57.362 - s2))
    # density scaling to T,p
    T = np.asarray(T, dtype=float)
    p = np.asarray(p, dtype=float)
    scale = (p / 101325.0) * (288.15 / T)
    return 1.0 + n_minus_1_std * scale

# --- ASTM E490 AM0 太阳光谱加载（自动识别单位；必要时从µm→nm） ---
def load_astm_e490(path: str | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    返回 (λ_nm, E_AM0_W_m2_nm)。优先读取本地 CSV/TXT。
    若 path 为空，抛错并给出官方下载源（见 fetch 脚本）。
    支持常见格式：两列/多列，自动识别 λ 单位（nm 或 µm），自动挑选辐照度列。
    官方数据：NREL “2000 ASTM Standard Extraterrestrial Spectrum Reference E-490-00”.
    """
    if path is None or not os.path.exists(path):
        raise FileNotFoundError(
            "未找到 AM0 数据文件。请先下载 NREL 官方文件，再传入路径。\n"
            "下载源见 https://www.nrel.gov/grid/solar-resource/spectra-astm-e490 （Data Files 部分）。"
        )
    # 读入
    txt = pathlib.Path(path).read_text(encoding="utf-8", errors="ignore")
    # 去掉注释行
    lines = [ln for ln in txt.splitlines() if ln.strip() and not ln.strip().startswith(("#","!","//","%"))]
    # 用空白分隔或逗号都支持
    sample = lines[0].replace(",", " ").split()
    # 找出两列或多列：第一列波长，后一列中选一个“irradiance-like”的
    data = np.genfromtxt(io.StringIO("\n".join(lines)))
    if data.ndim == 1:
        data = data[None, :]
    # 选择列
    lam_col = 0
    # heuristic: 优先找单位匹配的列
    if data.shape[1] < 2:
        raise ValueError("AM0 文件列数不足：至少需要波长与辐照度两列。")
    # 辐照度列：取数值波动最大的非第一列
    ir_cols = [i for i in range(1, data.shape[1])]
    irr_col = max(ir_cols, key=lambda i: np.nanstd(data[:, i]))
    lam = data[:, lam_col].astype(float)
    E = data[:, irr_col].astype(float)

    # 单位处理：判断 λ 是否在 [0.1, 10] 范围（µm），若是则转换为 nm 并把 E 从 W/m^2/µm 转为 W/m^2/nm
    lam_min, lam_max = np.nanmin(lam), np.nanmax(lam)
    if 0.1 <= lam_min <= 10.0 and lam_max <= 1e3:
        # 认为是 µm
        lam_nm = lam * 1e3
        # 如果 E 看起来像 W/m^2/µm → 除以 1000 变成 /nm
        # （无法完全自动判断，按常见 NREL/ASTM 文件默认 /µm 处理）
        E_W_m2_nm = E / 1000.0
    else:
        # 认为是 nm
        lam_nm = lam
        E_W_m2_nm = E
    # 去 NaN/排序
    ok = np.isfinite(lam_nm) & np.isfinite(E_W_m2_nm)
    lam_nm = lam_nm[ok]; E_W_m2_nm = E_W_m2_nm[ok]
    idx = np.argsort(lam_nm)
    return lam_nm[idx], E_W_m2_nm[idx]

# --- CIE 1931 2° CMFs 加载（官方 CSV） ---
def load_cie_1931_2deg(path: str | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    返回 (λ_nm, x̄, ȳ, z̄)，来自 CIE 官方开放数据集 CSV。
    如果未提供本地路径，则报错并提示去 CIE 下载。
    官方数据集：CIE_xyz_1931_2deg.csv（1 nm 采样，CIE 018:2019 表 6）。
    """
    if path is None or not os.path.exists(path):
        raise FileNotFoundError(
            "未找到 CIE 1931 2° CMFs 文件。请从 CIE 官方开放数据集下载 CSV：\n"
            "https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer"
        )
    df = pd.read_csv(path)
    # 尝试匹配列名（不同分发版本可能不同）
    cols = [c.lower() for c in df.columns]
    # 找波长列
    w_idx = np.argmin([abs("wavelength" in c or "lambda" in c or c.strip()=="nm") == False for c in cols])  # fallback
    # 更稳妥：优先找包含 'wavelength' / 'nm'
    w_candidates = [i for i,c in enumerate(cols) if ("wavelength" in c) or (c.strip()=="nm")]
    w_col = w_candidates[0] if w_candidates else 0
    # 找三条 CMFs（取后三列或名字含 xbar/ybar/zbar）
    def find_col(keys, default):
        for i,c in enumerate(cols):
            if any(k in c for k in keys): return i
        return default
    x_col = find_col(["x_bar","xbar","x 2"], 1)
    y_col = find_col(["y_bar","ybar","y 2"], 2)
    z_col = find_col(["z_bar","zbar","z 2"], 3)

    lam_nm = df.iloc[:, w_col].to_numpy(float)
    xbar = df.iloc[:, x_col].to_numpy(float)
    ybar = df.iloc[:, y_col].to_numpy(float)
    zbar = df.iloc[:, z_col].to_numpy(float)

    # 清理与排序
    ok = np.isfinite(lam_nm) & np.isfinite(xbar) & np.isfinite(ybar) & np.isfinite(zbar)
    lam_nm = lam_nm[ok]; xbar=xbar[ok]; ybar=ybar[ok]; zbar=zbar[ok]
    idx = np.argsort(lam_nm)
    return lam_nm[idx], xbar[idx], ybar[idx], zbar[idx]

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
