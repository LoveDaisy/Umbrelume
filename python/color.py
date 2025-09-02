# ──────────────────────────────────────────────────────────────────────────────
# color.py —— 光谱→XYZ→sRGB（D65）
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np

def spectrum_to_XYZ(lam_nm: np.ndarray, Llam: np.ndarray,
                    xbar: np.ndarray, ybar: np.ndarray, zbar: np.ndarray) -> tuple[float,float,float]:
    """Simple Riemann sum (nm grid): XYZ = ∑ Lλ * CMF * Δλ"""
    dlam = np.gradient(lam_nm)
    X = float(np.sum(Llam * xbar * dlam))
    Y = float(np.sum(Llam * ybar * dlam))
    Z = float(np.sum(Llam * zbar * dlam))
    return X, Y, Z

def XYZ_to_sRGB_D65(X: float, Y: float, Z: float) -> tuple[float,float,float]:
    """sRGB (IEC 61966-2-1) with D65; linear → gamma-encoded."""
    # XYZ to linear sRGB matrix (D65)
    M = np.array([[ 3.2406, -1.5372, -0.4986],
                  [-0.9689,  1.8758,  0.0415],
                  [ 0.0557, -0.2040,  1.0570]])
    rgb_lin = M @ np.array([X, Y, Z])
    # OETF
    def enc(c):
        return (1.055*(c**(1/2.4)) - 0.055) if c > 0.0031308 else 12.92*c
    r, g, b = (enc(float(max(0.0, v))) for v in rgb_lin)
    return r, g, b
