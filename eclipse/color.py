"""CIE 1931 colorimetry and an explicit, shared display transform.

Spectral integration returns unscaled XYZ. Display values use the uneclipsed
solar spectrum as one scene-wide white and luminance reference. Bradford
adaptation maps that white to D65; it is never re-estimated per pixel.

The sRGB matrix/transfer function and Bradford matrix follow the sample color
conversion code in https://www.w3.org/TR/css-color-4/#color-conversion-code .
This is a photopic colorimetric visualization, not a mesopic appearance model.
"""
from __future__ import annotations

import numpy as np

from .spectra import cie_xyz


XYZ_TO_SRGB = np.array([
    [3.2409699419045226, -1.537383177570094, -0.4986107602930034],
    [-0.9692436362808796, 1.8759675015077202, 0.04155505740717559],
    [0.05563007969699366, -0.20397695888897652, 1.0569715142428786],
])
BRADFORD = np.array([
    [0.8951, 0.2664, -0.1614],
    [-0.7502, 1.7135, 0.0367],
    [0.0389, -0.0685, 1.0296],
])
D65 = np.array([0.9504559270516716, 1.0, 1.0890577507598784])


def _spectral_inputs(spectra, wavelength_nm) -> tuple[np.ndarray, np.ndarray]:
    wl = np.asarray(wavelength_nm, dtype=float)
    values = np.asarray(spectra, dtype=float)
    if wl.ndim != 1 or wl.size < 2 or not np.all(np.isfinite(wl)):
        raise ValueError("wavelength_nm must be a finite one-dimensional grid of length >= 2")
    if not np.all(np.diff(wl) > 0):
        raise ValueError("wavelength_nm must be strictly increasing")
    if values.ndim == 0 or values.shape[-1] != wl.size:
        raise ValueError("the last spectra axis must match wavelength_nm")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("spectra must be finite and nonnegative")
    return values, wl


def spectrum_to_xyz(spectra, wavelength_nm) -> np.ndarray:
    """Integrate spectra (..., W) against CIE CMFs, returning (..., 3).

    With spectral irradiance in W m^-2 nm^-1, Y is the photopic weighted
    irradiance divided by 683 lm/W; no factor of 683 is applied here. Relative
    Y ratios are consequently also relative photopic irradiance ratios.
    """
    values, wl = _spectral_inputs(spectra, wavelength_nm)
    cmf = cie_xyz(wl)
    return np.trapezoid(values[..., :, None] * cmf, x=wl, axis=-2)


def xyz_to_linear_srgb(xyz) -> np.ndarray:
    """Convert D65-referred XYZ (..., 3) to unclipped linear sRGB."""
    values = np.asarray(xyz, dtype=float)
    if values.ndim == 0 or values.shape[-1] != 3 or not np.all(np.isfinite(values)):
        raise ValueError("xyz must be finite with final axis length 3")
    return values @ XYZ_TO_SRGB.T


def encode_srgb(linear_rgb) -> np.ndarray:
    """Encode linear RGB using the sRGB transfer function, clipped to [0, 1]."""
    values = np.asarray(linear_rgb, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("linear_rgb must be finite")
    values = np.clip(values, 0.0, 1.0)
    return np.where(values <= 0.0031308, 12.92 * values,
                    1.055 * values ** (1.0 / 2.4) - 0.055)


def linear_display_rgb(spectra, wavelength_nm, reference_spectrum, *,
                       normalize_luminance: bool = False) -> np.ndarray:
    """Return solar-white-adapted linear RGB before exposure and gamut handling.

    Normalization sets the unadapted CIE Y of every nonzero spectrum to one;
    that option is an explicitly nonphotometric chromaticity diagnostic.
    Ordinary output retains irradiance relative to the solar reference.
    """
    values, wl = _spectral_inputs(spectra, wavelength_nm)
    reference, _ = _spectral_inputs(reference_spectrum, wl)
    if reference.ndim != 1:
        raise ValueError("reference_spectrum must be one-dimensional")
    ref_xyz = spectrum_to_xyz(reference, wl)
    if ref_xyz[1] <= 0:
        raise ValueError("reference_spectrum must have positive photopic irradiance")
    source_white = ref_xyz / ref_xyz[1]
    source_cones = BRADFORD @ source_white
    if np.any(source_cones <= 0):
        raise ValueError("reference_spectrum does not define a valid Bradford white")
    adaptation = np.linalg.solve(
        BRADFORD, np.diag((BRADFORD @ D65) / source_cones) @ BRADFORD
    )
    xyz = spectrum_to_xyz(values, wl) / ref_xyz[1]
    if normalize_luminance:
        luminance = xyz[..., 1:2]
        xyz = np.divide(xyz, luminance, out=np.zeros_like(xyz), where=luminance > 0)
    return xyz_to_linear_srgb(xyz @ adaptation.T)


def finish_display(linear_rgb, exposure_stops: float = 0.0) -> np.ndarray:
    """Apply global exposure, discard negative RGB and compress bright colors.

    Above display range, all three channels at a pixel are divided by the
    same maximum channel. This retains RGB ratios in highlights but compresses
    their luminance. Below that limit, exposure is strictly linear. Negative
    channels from colors outside sRGB are clipped; no color appearance or
    perceptual gamut mapping is claimed.
    """
    if not np.isfinite(exposure_stops) or abs(exposure_stops) > 100:
        raise ValueError("exposure_stops must be finite and between -100 and 100")
    values = np.asarray(linear_rgb, dtype=float)
    if values.ndim == 0 or values.shape[-1] != 3 or not np.all(np.isfinite(values)):
        raise ValueError("linear_rgb must be finite with final axis length 3")
    values = np.maximum(values, 0.0) * 2.0 ** exposure_stops
    values /= np.maximum(1.0, np.max(values, axis=-1, keepdims=True))
    return encode_srgb(values)


def display_rgb(spectra, wavelength_nm, reference_spectrum, exposure_stops: float = 0.0,
                normalize_luminance: bool = False) -> np.ndarray:
    """Convert spectra to display RGB using one solar white and global exposure."""
    return finish_display(
        linear_display_rgb(spectra, wavelength_nm, reference_spectrum,
                           normalize_luminance=normalize_luminance),
        exposure_stops,
    )
