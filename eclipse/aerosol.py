"""Fixed, independent OPAC free-tropospheric and stratospheric-background priors.

The default describes a climatological/model component, not a measured profile
for the June 2011 eclipse. Extinction includes scattering out of the direct
beam and absorption. No scattering-into-beam or aerosol refractivity is added.
The production default composes two layers: the OPAC free troposphere (2-12 km)
and the OPAC stratospheric background (12-35 km).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Provenance of the pinned prior rows. The numeric values themselves are
# anchored to these rows by pytest (json-anchor test); the SHA identifies the
# OPAC source paper recorded in sources.json alongside the parameter rows.
PRIOR_PARAMETERS_JSON = (
    "validation/results/extinction_gap_diagnosis/prior_parameters.json")
OPAC_PDF_SHA256 = ("68a20ba9488db1e0f65a9976d4ee46e6c879f89054fcff583ba67b"
                   "b5d68d425c")
FREE_TROPOSPHERE_SOURCE = {
    "candidate": "free_troposphere_reference",
    "rows": ("aod550_2_to_12km_reference", "angstrom_0.35_to_0.5um",
             "angstrom_0.5_to_0.8um", "scale_height_m"),
}
STRATOSPHERE_BACKGROUND_SOURCE = {
    "candidate": "stratosphere_background",
    "rows": ("aod550", "angstrom_0.35_to_0.5um", "angstrom_0.5_to_0.8um",
             "layer_bottom_m", "layer_top_m", "scale_height_m"),
}
# Layer order of ``AerosolColumn.opac_default()`` paired with each layer's
# pinned prior-row provenance. CLI metadata and validation products consume
# this single pairing, so the layer order and its source rows are written in
# exactly one place and cannot drift apart when the default column evolves.
OPAC_DEFAULT_LAYER_SOURCES = (FREE_TROPOSPHERE_SOURCE, STRATOSPHERE_BACKGROUND_SOURCE)


@dataclass(frozen=True)
class AerosolProfile:
    """Finite exponential layer with a fixed 550 nm vertical optical depth.

    Distances are geometric metres and extinction is m^-1. The defaults come
    from OPAC Tables 3 and 5 (Hess, Koepke & Schult, 1998): free troposphere
    2–12 km, scale height 8 km, AOD550 0.013, and spectral exponents 1.58
    over 500–800 nm and 1.21 over 350–500 nm.
    """

    aod550: float = 0.013
    bottom_m: float = 2_000.0
    top_m: float = 12_000.0
    scale_height_m: float = 8_000.0
    angstrom_visible: float = 1.58
    angstrom_blue: float = 1.21
    break_nm: float = 500.0

    def __post_init__(self):
        values = (self.aod550, self.bottom_m, self.top_m, self.scale_height_m,
                  self.angstrom_visible, self.angstrom_blue, self.break_nm)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("Aerosol parameters must be finite.")
        if self.aod550 < 0:
            raise ValueError("AOD550 must be nonnegative.")
        if self.bottom_m < 0 or self.top_m <= self.bottom_m:
            raise ValueError("Aerosol layer requires 0 <= bottom_m < top_m.")
        if self.scale_height_m <= 0 or self.break_nm <= 0:
            raise ValueError("Scale height and spectral break must be positive.")
        if self.angstrom_visible < 0 or self.angstrom_blue < 0:
            raise ValueError("Angstrom exponents must be nonnegative.")

    @classmethod
    def opac_free_troposphere(cls, aod550: float = 0.013) -> AerosolProfile:
        """Published OPAC shape and spectral law, with explicit column control.

        Parameters are pinned to ``prior_parameters.json`` rows
        ``free_troposphere_reference`` (aod550_2_to_12km_reference 0.013,
        Angstrom 1.21 over 350-500 nm and 1.58 over 500-800 nm, scale height
        8 km; the 2 km bottom and 12 km top are stated in that row's source
        location). The OPAC PDF sha256 is the module constant
        ``OPAC_PDF_SHA256``.
        """
        return cls(aod550=aod550)

    @classmethod
    def opac_stratospheric_background(cls) -> AerosolProfile:
        """OPAC stratospheric-background layer: 12-35 km, near homogeneous.

        Parameters are pinned to ``prior_parameters.json`` rows
        ``stratosphere_background`` (aod550 0.005, Angstrom 0.74 over
        350-500 nm and 1.14 over 500-800 nm, layer 12-35 km, scale height
        99 km). The 99 km scale height over the 23 km thickness is the OPAC
        representation of a homogeneous layer (exp(-23/99) ~ 0.79). OPAC
        fixes this layer independent of volcanic activity.
        """
        return cls(aod550=0.005, bottom_m=12_000.0, top_m=35_000.0,
                   scale_height_m=99_000.0, angstrom_visible=1.14,
                   angstrom_blue=0.74)

    @classmethod
    def none(cls) -> AerosolProfile:
        """Zero extinction, retaining the same profile API for counterfactuals."""
        return cls(aod550=0.0)

    def extinction550(self, altitude_m) -> np.ndarray:
        """Extinction coefficient at 550 nm [m^-1], zero outside the layer.

        The shifted exponential and expm1 denominator preserve the exact
        requested vertical AOD for layers far from sea level or very thin
        compared with their scale height. Layer endpoints have measure zero;
        quadrature must split at both boundaries.
        """
        altitude = np.asarray(altitude_m, dtype=float)
        if np.any(~np.isfinite(altitude)) or np.any(altitude < 0):
            raise ValueError("Geometric altitude must be finite and nonnegative [m].")
        thickness = self.top_m - self.bottom_m
        denominator = self.scale_height_m * -np.expm1(-thickness/self.scale_height_m)
        # Clipping the shifted coordinate avoids evaluating an exponentially
        # large unused value below a high-altitude layer.
        shifted = np.maximum(altitude-self.bottom_m, 0)
        value = self.aod550/denominator * np.exp(-shifted/self.scale_height_m)
        return np.where((altitude >= self.bottom_m) & (altitude <= self.top_m), value, 0.0)

    def spectral_factor(self, wavelength_nm) -> np.ndarray:
        """Continuous broken power law beta(lambda)/beta(550 nm).

        Normalization is tied to the physical 550 nm reference, independent
        of whether 550 nm is present in the caller's wavelength array.
        """
        wavelength = np.asarray(wavelength_nm, dtype=float)
        if np.any(~np.isfinite(wavelength)) or np.any(wavelength <= 0):
            raise ValueError("Wavelengths must be positive and finite [nm].")
        exponent = np.where(wavelength >= self.break_nm,
                            self.angstrom_visible, self.angstrom_blue)
        raw = (wavelength/self.break_nm)**(-exponent)
        reference_exponent = (self.angstrom_visible if 550.0 >= self.break_nm
                              else self.angstrom_blue)
        reference = (550.0/self.break_nm)**(-reference_exponent)
        return raw/reference


@dataclass(frozen=True)
class AerosolColumn:
    """Ordered composition of finite exponential aerosol layers.

    A column's extinction is the sum of its layers' coefficients; layer
    supports may overlap only on measure-zero boundaries, so the vertical
    AOD of a column is exactly the sum of its layers' AOD550 (documented
    invariant, asserted in tests). Each layer carries its own geometry, AOD
    and broken-power-law spectral exponents; no merged spectral factor exists,
    the only legal expansion is the per-layer (tau, profile) pair.
    """

    layers: tuple[AerosolProfile, ...]

    def __post_init__(self):
        if not self.layers:
            raise ValueError("AerosolColumn requires at least one layer.")
        if not all(isinstance(layer, AerosolProfile) for layer in self.layers):
            raise ValueError("AerosolColumn layers must all be AerosolProfile instances.")

    @classmethod
    def opac_default(cls, ft_aod550: float = 0.013) -> AerosolColumn:
        """Production default: OPAC free troposphere + stratospheric background.

        ``ft_aod550`` overrides only the free-troposphere layer (the CLI
        ``--aod550`` semantics); the stratospheric background stays at its
        published 0.005.
        """
        return cls((AerosolProfile.opac_free_troposphere(aod550=ft_aod550),
                    AerosolProfile.opac_stratospheric_background()))


# The complete type face of an aerosol prior argument: single layer, layer
# column, or no aerosol at all.
AerosolPrior = AerosolProfile | AerosolColumn | None


def active_aerosol_layers(aerosol: AerosolPrior) -> tuple[AerosolProfile, ...]:
    """Single normalization point for every aerosol prior argument.

    Returns the layers that carry extinction, in order. Zero-AOD layers are
    dropped: they change no physical quantity, and ray tables store tau tables
    only for layers that have one. A column mixing zero and nonzero layers
    therefore normalizes to its active subset, which is exactly what the
    layer-isolation counterfactuals compare.
    """
    if aerosol is None:
        return ()
    if isinstance(aerosol, AerosolColumn):
        layers = aerosol.layers
    elif isinstance(aerosol, AerosolProfile):
        layers = (aerosol,)
    else:
        raise TypeError(f"Unsupported aerosol prior: {type(aerosol).__name__}")
    return tuple(layer for layer in layers if layer.aod550 > 0)
