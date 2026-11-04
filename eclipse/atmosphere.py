"""A spherical dry atmosphere and an explicit idealised ozone profile.

US Standard Atmosphere 1976 geopotential layers are used below approximately
86 km geometric altitude. An isothermal, constant-composition hydrostatic tail
extends them to the model top. This tail is an optical convenience, not the
upper-atmosphere portion of US1976. All public altitudes are geometric metres.
"""

from dataclasses import dataclass
from math import erf, sqrt

import numpy as np


EARTH_RADIUS_M = 6_371_000.0
REFERENCE_NUMBER_DENSITY = 101325.0 / (1.380649e-23 * 288.15)
DOBSON_UNIT = 2.687e20  # molecules per m^2

# The standard uses an effective geopotential radius, distinct from the
# physical radius of the spherical Earth used for ray geometry.
_GEOPOTENTIAL_RADIUS_M = 6_356_766.0
_BOLTZMANN = 1.380649e-23
_G0 = 9.80665
_R_AIR = 8314.32 / 28.9644  # US1976 convention [J kg^-1 K^-1]
_H_BASE = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0,
                    51000.0, 71000.0, 84852.0])
_LAPSE = np.array([-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0])
LAYER_BOUNDARIES_M = _GEOPOTENTIAL_RADIUS_M * _H_BASE / (_GEOPOTENTIAL_RADIUS_M - _H_BASE)
LAYER_BOUNDARIES_M.setflags(write=False)


def _layer_bases():
    temperature = [288.15]
    pressure = [101325.0]
    for index, dh in enumerate(np.diff(_H_BASE)):
        lower_t = temperature[-1]
        upper_t = lower_t + _LAPSE[index] * dh
        if _LAPSE[index] == 0.0:
            upper_p = pressure[-1] * np.exp(-_G0 * dh / (_R_AIR * lower_t))
        else:
            upper_p = pressure[-1] * (upper_t / lower_t) ** (
                -_G0 / (_R_AIR * _LAPSE[index])
            )
        temperature.append(upper_t)
        pressure.append(upper_p)
    return np.asarray(temperature), np.asarray(pressure)


_T_BASE, _P_BASE = _layer_bases()


def _altitudes(altitude_m):
    altitude = np.asarray(altitude_m, dtype=float)
    if np.any(~np.isfinite(altitude)) or np.any(altitude < 0.0):
        raise ValueError("Geometric altitude must be finite and at least zero [m].")
    return altitude


def _state(altitude_m):
    """Untruncated hydrostatic temperature, pressure, and d(log N)/dz."""
    altitude = _altitudes(altitude_m)
    geopotential = _GEOPOTENTIAL_RADIUS_M * altitude / (_GEOPOTENTIAL_RADIUS_M + altitude)
    # Select layers in the public geometric coordinate. Converting an exact
    # layer boundary to geopotential and back can round just below it, giving
    # the wrong one-sided density gradient at the ray's turning point.
    layer = np.minimum(np.searchsorted(LAYER_BOUNDARIES_M, altitude, side="right") - 1, 7)
    dh = geopotential - _H_BASE[layer]
    lapse = _LAPSE[layer]
    temperature = _T_BASE[layer] + lapse * dh
    # Use a harmless nonzero denominator for the isothermal entries, then
    # select their exponential form explicitly.
    safe_lapse = np.where(lapse == 0.0, 1.0, lapse)
    pressure = _P_BASE[layer] * np.where(
        lapse == 0.0,
        np.exp(-_G0 * dh / (_R_AIR * _T_BASE[layer])),
        (temperature / _T_BASE[layer]) ** (-_G0 / (_R_AIR * safe_lapse)),
    )
    d_geopotential_dz = (_GEOPOTENTIAL_RADIUS_M / (_GEOPOTENTIAL_RADIUS_M + altitude)) ** 2
    log_density_gradient = -(_G0 / _R_AIR + lapse) / temperature * d_geopotential_dz
    return temperature, pressure, log_density_gradient


@dataclass(frozen=True)
class Atmosphere:
    """Dry molecular atmosphere plus an independently normalised ozone column.

    ``density_scale`` scales molecular density and refractivity only. Ozone is
    controlled by its own column ``ozone_du``. Ozone's Gaussian width is one
    standard deviation. Density is zero strictly above ``top_m``; the ray
    tracer treats the top as the interface to vacuum. The tiny residual air
    density there makes the finite-top approximation explicit.
    """

    ozone_du: float = 300.0
    density_scale: float = 1.0
    top_m: float = 120_000.0
    ozone_center_m: float = 25_000.0
    ozone_width_m: float = 8_000.0

    def __post_init__(self):
        values = (self.ozone_du, self.density_scale, self.top_m,
                  self.ozone_center_m, self.ozone_width_m)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("Atmosphere parameters must be finite.")
        if self.ozone_du < 0.0 or self.density_scale < 0.0:
            raise ValueError("Ozone column and molecular density scale must be nonnegative.")
        if self.top_m <= 0.0 or self.ozone_width_m <= 0.0:
            raise ValueError("Atmosphere top and ozone width must be positive [m].")
        if not 0.0 <= self.ozone_center_m <= self.top_m:
            raise ValueError("Ozone centre must lie between the ground and atmosphere top.")

    def number_density(self, altitude_m):
        """Molecular number density [molecules m^-3], zero above the top."""
        altitude = _altitudes(altitude_m)
        temperature, pressure, _ = _state(altitude)
        return np.where(
            altitude <= self.top_m,
            self.density_scale * pressure / (_BOLTZMANN * temperature),
            0.0,
        )

    def density_gradient(self, altitude_m):
        """Derivative dN/dz [molecules m^-4] within the atmosphere.

        At a US1976 layer boundary the upper-layer one-sided derivative is
        returned. The top cutoff's distributional jump is not included.
        """
        altitude = _altitudes(altitude_m)
        _, _, log_gradient = _state(altitude)
        return self.number_density(altitude) * log_gradient

    def ozone_density(self, altitude_m):
        """Ozone number density [molecules m^-3], column-normalised Gaussian."""
        altitude = _altitudes(altitude_m)
        width = self.ozone_width_m
        centre = self.ozone_center_m
        normalization_m = width * sqrt(np.pi / 2.0) * (
            erf((self.top_m - centre) / (sqrt(2.0) * width))
            + erf(centre / (sqrt(2.0) * width))
        )
        value = self.ozone_du * DOBSON_UNIT / normalization_m * np.exp(
            -0.5 * ((altitude - centre) / width) ** 2
        )
        return np.where(altitude <= self.top_m, value, 0.0)
