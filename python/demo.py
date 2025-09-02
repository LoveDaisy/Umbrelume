# ──────────────────────────────────────────────────────────────────────────────
# demo.py —— 最小演示：扫一条 r 径向曲线
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
from atmo import (us_standard_atmosphere_1976, load_astm_e490, load_cie_1931_2deg,
                  rayleigh_sigma_bodhaine, ozone_sigma_serdyuchenko, SpectralData)
from raytrace import build_alpha_table, default_solar_disk_samples
from engine import direct_term_spectrum, scattered_term_spectrum, spectrum_to_srgb

def main():
    # λ网格与光谱数据
    lam, Esun = load_astm_e490()                     # ASTM E490 AM0  :contentReference[oaicite:0]{index=0}
    xbar, ybar, zbar = load_cie_1931_2deg()          # CIE 1931 CMFs   :contentReference[oaicite:1]{index=1}
    sigma_R = rayleigh_sigma_bodhaine(lam)           # Rayleigh σ_R     :contentReference[oaicite:2]{index=2}
    sigma_O3 = ozone_sigma_serdyuchenko(lam)         # O3 σ(λ)         :contentReference[oaicite:3]{index=3}
    spec = SpectralData(lam, Esun, sigma_R, sigma_O3, xbar, ybar, zbar)

    # 大气与折射表
    z_grid = np.arange(0.0, 80.5, 1.0)
    atmo = us_standard_atmosphere_1976(z_grid)       # USSA76           :contentReference[oaicite:4]{index=4}
    alpha_table = build_alpha_table(atmo, lam0_nm=550.0)  # Auer–Standish  :contentReference[oaicite:5]{index=5}
    solar_samples = default_solar_disk_samples(K=1)       # 点太阳；可改 K>1 加肢变暗

    # r 采样（弧分→弧度）
    r_arcmin = np.arange(0.0, 80.1, 1.0)
    r_rad = np.deg2rad(r_arcmin / 60.0)

    # 计算直达项 + 散射项（单次瑞利）
    phi_grid = np.linspace(0.0, 2*np.pi, 720, endpoint=False)
    z_scat   = np.arange(10.0, 40.1, 1.0)

    rgb_strip = []
    for r in r_rad:
        Ein_dir = direct_term_spectrum(r, atmo, spec, alpha_table, solar_samples)
        Ein_sca = scattered_term_spectrum(r, atmo, spec, solar_samples, phi_grid, z_scat)
        Llam = (Ein_dir + Ein_sca)  # * R_moon(λ)  # 先当常数 1
        rgb_strip.append(spectrum_to_srgb(Llam, spec))

    # TODO: 保存/显示 r→sRGB（例如保存为 PNG 条形图）
    # TODO: 输出 Y 亮度与 (x,y) 色度随 r 的曲线

if __name__ == "__main__":
    main()
