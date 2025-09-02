# ──────────────────────────────────────────────────────────────────────────────
# transmit.py —— 光学厚度与透过率（沿“折射路径”与“入/出散射路径”）
# ──────────────────────────────────────────────────────────────────────────────
import numpy as np
from atmo import Atmosphere1D, SpectralData

@dataclass
class PathSegment:
    z_km: np.ndarray      # sample heights along path
    ds_m: np.ndarray      # segment lengths in meters

def refracted_tangent_path(atmo: Atmosphere1D, h_star_km: float) -> PathSegment:
    """Return symmetric tangent path around h*: heights & segment lengths.
    TODO: either export segments from raytrace integration or reconstruct with geometry."""
    raise NotImplementedError

def optical_depth_along_path(seg: PathSegment, spec: SpectralData) -> dict[str, np.ndarray]:
    """Compute τ fields on λ-grid: {'R': τ_R, 'O3': τ_O3, 'aer': τ_aer (optional)}."""
    # TODO: accumulate τ = ∑ k(λ,z) * ds using per-layer number densities & cross-sections
    raise NotImplementedError

def beer_transmittance(tau_components: dict[str, np.ndarray]) -> np.ndarray:
    """Return T(λ) = exp(-∑ τ_i)."""
    return np.exp(-sum(tau_components.values()))
