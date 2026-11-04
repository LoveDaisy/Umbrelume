"""A fixed, independent OPAC free-tropospheric extinction prior.

The default describes a climatological/model component, not a measured profile
for the June 2011 eclipse. Extinction includes scattering out of the direct
beam and absorption. No scattering-into-beam or aerosol refractivity is added.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


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
        """Published OPAC shape and spectral law, with explicit column control."""
        return cls(aod550=aod550)

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
