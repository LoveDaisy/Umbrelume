# selftest_min.py
import numpy as np
from atmo import us_standard_atmosphere_1976, rayleigh_sigma_bodhaine, load_astm_e490, load_cie_1931_2deg
import os

def approx(a, b, rtol=0.02, atol=0.0):
    return abs(a - b) <= (atol + rtol * abs(b))

def main():
    # ---- 1) USSA76: 基本层点值 ----
    z = np.array([0.0, 11.0, 20.0])  # geometric km
    atmo = us_standard_atmosphere_1976(z)

    T0, T11, T20 = atmo.T
    p0, p11, p20 = atmo.p
    N0 = atmo.air_number_density[0]

    # 对照：海平面 & 11 km（对照 ISA/USSA）
    # 资料：T0=288.15 K, p0=101325 Pa；T(11km)=216.65 K, p(11km)=22632 Pa
    assert approx(T0, 288.15, rtol=1e-4), f"T0 wrong: {T0}"
    assert approx(p0, 101325.0, rtol=1e-3), f"p0 wrong: {p0}"
    assert approx(T11, 216.65, rtol=1e-3), f"T11 wrong: {T11}"
    assert approx(p11, 22632.0, rtol=0.02), f"p11 off: {p11}"  # 2% 余量

    # 20 km：等温层延伸，T≈216.65 K；p≈(5.47~5.53)e3 Pa（不同表/舍入略有差异）
    assert approx(T20, 216.65, rtol=1e-3), f"T20 wrong: {T20}"
    assert 5.3e3 <= p20 <= 5.7e3, f"p20 out of range: {p20}"

    # 海平面分子数密度 ~ 2.547e25 m^-3（USSA76 表 2）
    assert 2.50e25 <= N0 <= 2.60e25, f"N0 out of range: {N0}"

    print("USSA76 checks passed.")
    print(f"  T(z km)  = [{T0:.2f}, {T11:.2f}, {T20:.2f}] K")
    print(f"  p(z km)  = [{p0:.1f}, {p11:.1f}, {p20:.1f}] Pa")
    print(f"  N0       = {N0:.4e} m^-3")

    # ---- 2) Rayleigh σ_R: 数量级 & λ^-4 标度 ----
    lam_nm = np.array([450.0, 550.0, 650.0])
    sig = rayleigh_sigma_bodhaine(lam_nm)

    s450, s550, s650 = sig
    # 550 nm 数量级 ~ 5e-31 m^2
    assert 3.0e-31 <= s550 <= 7.0e-31, f"sigma(550nm) unexpected: {s550:.3e}"

    # 近似 λ^-4：sigma(450)/sigma(650) ≈ (650/450)^4 ≈ 4.35，给10–15%的余量
    ratio = s450 / s650
    expect = (650.0/450.0)**4  # ≈ 4.353
    assert 0.85*expect <= ratio <= 1.15*expect, f"λ^-4 scaling off: {ratio:.3f} vs {expect:.3f}"

    print("Rayleigh sigma checks passed.")
    print(f"  σ_R(450, 550, 650 nm) = {s450:.3e}, {s550:.3e}, {s650:.3e} m^2")
    print(f"  ratio σ(450)/σ(650)   = {ratio:.3f} (expected ~{expect:.3f})")

    # ---- 3) ASTM E490 AM0 Solar Spectrum Loading ----
    data_path = "../data/e490_00a_amo.csv"
    if os.path.exists(data_path):
        try:
            lam_nm, E_AM0 = load_astm_e490(data_path)
            
            # Basic sanity checks on loaded data
            assert len(lam_nm) > 1000, f"Too few data points: {len(lam_nm)}"
            assert len(lam_nm) == len(E_AM0), "Wavelength and irradiance arrays length mismatch"
            
            # Wavelength range should cover UV to far IR (ASTM E490 covers 119.5 nm to 1000 µm)
            assert 100 <= np.min(lam_nm) <= 200, f"Minimum wavelength unexpected: {np.min(lam_nm):.1f} nm"
            assert 500000 <= np.max(lam_nm) <= 1500000, f"Maximum wavelength unexpected: {np.max(lam_nm):.1f} nm"
            
            # Solar constant integration check (rough): total should be ~1360-1370 W/m²
            total_irradiance = np.trapezoid(E_AM0, lam_nm)
            assert 1300 <= total_irradiance <= 1400, f"Solar constant off: {total_irradiance:.1f} W/m²"
            
            # Peak should be somewhere in visible range (solar spectrum can vary)
            vis_mask = (lam_nm >= 400) & (lam_nm <= 700)
            if np.any(vis_mask):
                peak_idx = np.argmax(E_AM0[vis_mask])
                peak_wavelength = lam_nm[vis_mask][peak_idx]
                assert 400 <= peak_wavelength <= 600, f"Visible peak at unexpected λ: {peak_wavelength:.1f} nm"
            
            print("ASTM E490 loading checks passed.")
            print(f"  Data points: {len(lam_nm)}")
            print(f"  λ range: {np.min(lam_nm):.1f} - {np.max(lam_nm):.1f} nm")
            print(f"  Solar constant: {total_irradiance:.1f} W/m²")
            print(f"  Visible peak: ~{peak_wavelength:.1f} nm")
            
        except Exception as e:
            print(f"ASTM E490 loading test failed: {e}")
    else:
        print("ASTM E490 data file not found, skipping solar spectrum test")

    # ---- 4) CIE 1931 2° Color Matching Functions Loading ----
    cie_path = "../data/CIE_xyz_1931_2deg.csv"
    if os.path.exists(cie_path):
        try:
            lam_nm, xbar, ybar, zbar = load_cie_1931_2deg(cie_path)
            
            # Basic sanity checks on loaded CIE data
            assert len(lam_nm) > 300, f"Too few CIE data points: {len(lam_nm)}"
            assert len(lam_nm) == len(xbar) == len(ybar) == len(zbar), "CIE CMF arrays length mismatch"
            
            # Wavelength range should cover visible spectrum (360-830 nm for CIE 1931)
            assert 350 <= np.min(lam_nm) <= 370, f"CIE minimum wavelength unexpected: {np.min(lam_nm):.1f} nm"
            assert 820 <= np.max(lam_nm) <= 840, f"CIE maximum wavelength unexpected: {np.max(lam_nm):.1f} nm"
            
            # Y CMF should peak around 555 nm (photopic peak)
            y_peak_idx = np.argmax(ybar)
            y_peak_wavelength = lam_nm[y_peak_idx]
            assert 550 <= y_peak_wavelength <= 560, f"Y peak at unexpected λ: {y_peak_wavelength:.1f} nm"
            
            # CMFs should be non-negative and have reasonable magnitudes  
            assert np.all(xbar >= 0) and np.all(ybar >= 0) and np.all(zbar >= 0), "CMFs should be non-negative"
            assert 1.0 <= np.max(ybar) <= 2.0, f"Y peak value unexpected: {np.max(ybar):.3f}"
            
            # X and Z should have peaks in blue/red regions
            x_peak_wavelength = lam_nm[np.argmax(xbar)]
            z_peak_wavelength = lam_nm[np.argmax(zbar)]
            assert 590 <= x_peak_wavelength <= 610, f"X peak at unexpected λ: {x_peak_wavelength:.1f} nm" 
            assert 440 <= z_peak_wavelength <= 460, f"Z peak at unexpected λ: {z_peak_wavelength:.1f} nm"
            
            print("CIE 1931 2° CMFs loading checks passed.")
            print(f"  Data points: {len(lam_nm)}")
            print(f"  λ range: {np.min(lam_nm):.1f} - {np.max(lam_nm):.1f} nm")
            print(f"  Y peak: {y_peak_wavelength:.1f} nm (value: {np.max(ybar):.3f})")
            print(f"  X peak: {x_peak_wavelength:.1f} nm, Z peak: {z_peak_wavelength:.1f} nm")
            
        except Exception as e:
            print(f"CIE 1931 2° CMFs loading test failed: {e}")
    else:
        print("CIE 1931 2° CMFs data file not found, skipping color matching functions test")

if __name__ == "__main__":
    main()
