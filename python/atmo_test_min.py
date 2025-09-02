# selftest_min.py
import numpy as np
from .atmo import us_standard_atmosphere_1976, rayleigh_sigma_bodhaine

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
    assert approx(T11, 216.65, rtol=1e-4), f"T11 wrong: {T11}"
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

if __name__ == "__main__":
    main()
