#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
demo.py — 端到端：USSA76 + AM0 + CIE -> 折射/散射 -> 颜色条 PNG/CSV

需要你项目里的模块：
- atmo.py: Atmosphere1D, SpectralData, us_standard_atmosphere_1976,
           rayleigh_sigma_bodhaine, load_astm_e490, load_cie_1931_2deg
- raytrace.py: build_alpha_table, invert_h_from_alpha（若未实现，会用本文件的 fallback）
- transmit.py: refracted_tangent_path, optical_depth_along_path
- scatter.py: single_cell_scatter_spectrum
- engine.py: compute_strip, EngineOptions

运行:
    python demo.py
会输出:
    out/strip_demo.png
    out/strip_demo_srgb.csv
"""
from __future__ import annotations
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ---- 项目内模块 ----
from color import XYZ_to_sRGB_D65
from atmo import (
    Atmosphere1D, SpectralData,
    us_standard_atmosphere_1976,
    rayleigh_sigma_bodhaine,
    load_astm_e490, load_cie_1931_2deg,
    ciddor_n_air,
)

from engine import compute_strip, EngineOptions

# 优先用你的 raytrace 实现；若不可用则回退
try:
    from raytrace import build_alpha_table as _build_alpha_table_native
    from raytrace import invert_h_from_alpha as _invert_h_from_alpha_native
    HAVE_NATIVE_RAYTRACE = True
except Exception:
    HAVE_NATIVE_RAYTRACE = False

# -------------------- I/O 目录 --------------------
HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
OUT_DIR = os.path.join(HERE, "..", "out")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# 你之前用 fetch_data.py 下载的默认文件名（可按需修改）
AM0_PATH = os.path.join(DATA_DIR, "e490_00a_amo.csv")
CIE_PATH = os.path.join(DATA_DIR, "CIE_xyz_1931_2deg.csv")

# -------------------- O3：轻量混合比与 Chappuis 近似（仅 demo 用） --------------------


def ozone_vmr_profile_ppmv(z_km: np.ndarray) -> np.ndarray:
    """简化的 O3 体积分数（ppmv）剖面，让外缘蓝绿能显现。"""
    base = 0.04
    peak = 8.0
    z0, sigma = 25.0, 7.0
    gauss = peak * np.exp(-0.5 * ((z_km - z0) / sigma) ** 2)
    upper = 0.6 * np.exp(-(z_km - 35.0) / 12.0)
    vmr = base + gauss + upper
    vmr[z_km < 8.0] *= 0.5
    return vmr


def sigma_O3_chappuis_fallback(lam_nm: np.ndarray) -> np.ndarray:
    """Chappuis 带极简双高斯近似（m^2/分子），仅为演示。"""
    lam = np.asarray(lam_nm, float)
    g1 = 4.0e-27 * np.exp(-0.5 * ((lam - 575.0) / 30.0) ** 2)
    g2 = 3.0e-27 * np.exp(-0.5 * ((lam - 603.0) / 35.0) ** 2)
    wing = 5.0e-29 * np.exp(-0.5 * ((lam - 700.0) / 80.0) ** 2)
    out = (g1 + g2 + wing)
    out[(lam < 430.0) | (lam > 780.0)] = 0.0
    return out

# -------------------- 回退 α(h) 表（若 raytrace 未就绪） --------------------


def _build_alpha_table_fallback(atmo: Atmosphere1D,
                                lam0_nm: float = 550.0,
                                zmin_km: float = 6.0,
                                zmax_km: float = 60.0,
                                nz: int = 200) -> dict:
    """
    用电波掩星常见公式从 n(r) 数值积分 α(a)，建立 h*→α 的查找表。
      α(a) = -2 a ∫_{r_t}^{∞} (d ln n / dr) / sqrt( (n r)^2 - a^2 )  dr,
      a = n(r_t) r_t,  r_t = R + h*
    """
    z_full = np.asarray(atmo.z, float)
    m = (z_full >= 0.0) & (z_full <= 80.0)
    z = z_full[m]
    r = (6371.0 + z) * 1000.0
    n_z = ciddor_n_air(lam0_nm,
                       T=np.interp(z, atmo.z, atmo.T),
                       p=np.interp(z, atmo.z, atmo.p))
    x = n_z * r
    ln_n = np.log(n_z)
    dlnn_dz = np.gradient(ln_n, z)   # per km
    dlnn_dr = dlnn_dz / 1000.0       # per m

    h_grid = np.linspace(zmin_km, zmax_km, nz)
    alpha = np.zeros_like(h_grid)
    for i, h in enumerate(h_grid):
        rt = (6371.0 + h) * 1000.0
        n_t = float(ciddor_n_air(lam0_nm,
                                 T=np.interp(h, atmo.z, atmo.T),
                                 p=np.interp(h, atmo.z, atmo.p)))
        a = n_t * rt
        mm = (r >= rt)
        r_seg = r[mm]
        dlnn_dr_seg = dlnn_dr[mm]
        denom = np.sqrt(np.maximum((n_z[mm]*r_seg)**2 - a*a, 0.0))
        denom[0] = max(denom[0], 1.0)
        integrand = dlnn_dr_seg / np.maximum(denom, 1e-30)
        alpha[i] = -2.0 * a * np.trapz(integrand, r_seg)
    return {"h_km": h_grid.astype(float), "alpha_rad": alpha.astype(float)}


def _invert_h_from_alpha_fallback(alpha_table: dict, r_eff_rad: float) -> float | None:
    h = alpha_table["h_km"]
    a = alpha_table["alpha_rad"]
    if not (np.all(np.diff(a) < 0) or np.all(np.diff(a) > 0)):
        idx = np.argsort(a)
        a = a[idx]
        h = h[idx]
    if (r_eff_rad < a.min()) or (r_eff_rad > a.max()):
        return None
    return float(np.interp(r_eff_rad, a[::-1], h[::-1]))

# -------------------- 主流程 --------------------


def main():
    # 1) 加载外部数据（只调用 atmo.py 的加载器）
    lam_am0, Esun_AM0 = load_astm_e490(AM0_PATH)  # NREL / ASTM E-490 AM0 数据
    lam_cie, xbar, ybar, zbar = load_cie_1931_2deg(CIE_PATH)  # CIE 1931 2° CMFs

    # 2) 构造目标波长轴（可见 420–760 nm, 2 nm 步长），并插值到该轴
    lam = np.arange(420.0, 760.0 + 1e-9, 2.0)
    Esun_i = np.interp(lam, lam_am0, Esun_AM0)
    x_i = np.interp(lam, lam_cie, xbar)
    y_i = np.interp(lam, lam_cie, ybar)
    z_i = np.interp(lam, lam_cie, zbar)

    # 3) 大气：USSA76（0–80 km, 0.5 km）
    z_km = np.arange(0.0, 80.0 + 1e-9, 0.5)
    atmo = us_standard_atmosphere_1976(z_km)

    # O3 数密度（简化剖面）
    vmr_ppmv = ozone_vmr_profile_ppmv(atmo.z)  # ppmv
    o3_m3 = (vmr_ppmv * 1e-6) * atmo.air_number_density
    atmo.o3_number_density = o3_m3

    # 4) 光学截面：Rayleigh 用 Bodhaine；O3 用 Chappuis 近似（可替换为正式表）
    sigma_R = rayleigh_sigma_bodhaine(lam)
    sigma_O3 = sigma_O3_chappuis_fallback(lam)

    spec = SpectralData(
        lam=lam.astype(float),
        Esun_AM0=Esun_i.astype(float),
        sigma_R=sigma_R.astype(float),
        sigma_O3=sigma_O3.astype(float),
        cie_x=x_i.astype(float),
        cie_y=y_i.astype(float),
        cie_z=z_i.astype(float),
    )

    # 5) 折射 α(h) 表
    if HAVE_NATIVE_RAYTRACE:
        try:
            alpha_table = _build_alpha_table_native(atmo, lam0_nm=550.0)
            invert_h = _invert_h_from_alpha_native
            print("[raytrace] 使用项目内置 build_alpha_table()。")
        except Exception:
            alpha_table = _build_alpha_table_fallback(atmo, lam0_nm=550.0)
            invert_h = _invert_h_from_alpha_fallback
            print("[raytrace] 内置实现失败，使用 fallback α(h) 。")
    else:
        alpha_table = _build_alpha_table_fallback(atmo, lam0_nm=550.0)
        invert_h = _invert_h_from_alpha_fallback
        print("[raytrace] 未检测到内置实现，使用 fallback α(h) 。")

    # 6) 径向采样与散射网格
    r_arcmin = np.linspace(0.0, 60.0, 241)
    phi_bins = 360              # φ 方向的积分分辨率
    z_min, z_max, dz = 8.0, 55.0, 1.0

    opt = EngineOptions(
        use_solar_disk=False,    # 点太阳；后续可打开有限太阳盘
        lunar_albedo=1.0,
        hg_g=0.0,                # 仅瑞利；要看前向偏置可设 0.6–0.8
        use_aerosol=False,
    )

    # 7) 计算颜色条（engine 返回 dict）
    result = compute_strip(
        r_arcmin=r_arcmin,
        atmo=atmo,
        spec=spec,
        alpha_table=alpha_table,
        solar_samples=None,             # None=自动：点太阳 / 或有限太阳盘样本
        options=opt,
        phi_bins=phi_bins,
        z_min_km=z_min, z_max_km=z_max, dz_km=dz,
    )
    # sRGB（绝对亮度编码）通常会 >1 被裁剪，为更直观的“色带”，可把亮度归一化再编码一次
    # 不影响物理量，仅用于 PNG 可视化
    rgb_raw = result["rgb"]     # (N,3), 这是未归一直接 sRGB 编码，可能会饱和
    RESCALE_FOR_DISPLAY = True
    if RESCALE_FOR_DISPLAY:
        XYZ = result["XYZ"]     # (N,3)
        Ymax = float(np.max(XYZ[:, 1]))
        if Ymax > 0:
            rgb = np.array([XYZ_to_sRGB_D65(*(XYZ[i]/Ymax)) for i in range(XYZ.shape[0])],
                           dtype=float)
        else:
            rgb = rgb_raw
    else:
        rgb = rgb_raw

    # 8) 输出 PNG/CSV
    H = 40
    strip = np.clip(np.broadcast_to(rgb[None, :, :], (H, rgb.shape[0], 3)), 0, 1)
    plt.figure(figsize=(10, 1))
    plt.imshow(strip, origin="lower",
               extent=[r_arcmin[0], r_arcmin[-1], 0, 1], aspect="auto")
    plt.yticks([])
    plt.xlabel("r (arcmin)")
    plt.tight_layout()
    png_path = os.path.join(OUT_DIR, "strip_demo.png")
    plt.savefig(png_path, dpi=200)
    plt.close()

    csv_path = os.path.join(OUT_DIR, "strip_demo_srgb.csv")
    np.savetxt(csv_path, np.column_stack([r_arcmin, rgb]),
               delimiter=",", header="r_arcmin,R,G,B", comments="")

    print(f"[OK] 写出：{png_path}")
    print(f"[OK] 写出：{csv_path}")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        # 当数据文件未找到时给出明确提示（仅加载器会抛）
        sys.stderr.write(
            "未找到数据文件，请先下载后放入 data/ 目录：\n"
            "  - AM0（ASTM E-490）：见 NREL 下载页（Data Files）\n"
            "  - CIE 1931 2° CMFs：见 CIE 官方数据集\n"
        )
        raise
