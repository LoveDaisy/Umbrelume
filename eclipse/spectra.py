"""Solar and molecular spectra, with wavelengths in nm and outputs in SI.

The packaged tables are independent copies of the source data documented in
``data/SOURCES.md``. Nothing is imported from the sky-colour project.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np


_DATA_DIR = Path(__file__).resolve().parent / "data"
BOLTZMANN_CONSTANT = 1.380649e-23
# Molecular number density of the US1976 sea-level standard state; the single
# definition, re-exported by eclipse.atmosphere for its own consumers.
REFERENCE_NUMBER_DENSITY = 101325.0 / (BOLTZMANN_CONSTANT * 288.15)


def _wavelengths(wavelength_nm):
    values = np.asarray(wavelength_nm, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("Wavelengths must be positive and finite [nm].")
    return values


@lru_cache(maxsize=None)
def _table(filename):
    table = np.loadtxt(_DATA_DIR / filename, delimiter=",", skiprows=1)
    if not np.all(np.isfinite(table)) or np.any(np.diff(table[:, 0]) <= 0):
        raise ValueError(f"Invalid spectral table: {filename}")
    table.setflags(write=False)
    return table


def _resample(filename, wavelength_nm):
    wl = _wavelengths(wavelength_nm)
    table = _table(filename)
    if np.any(wl < table[0, 0]) or np.any(wl > table[-1, 0]):
        raise ValueError(
            f"{filename} covers {table[0, 0]:g}–{table[-1, 0]:g} nm; "
            "extrapolation is not supported."
        )
    return np.stack(
        [np.interp(wl, table[:, 0], column) for column in table[:, 1:].T], axis=-1
    )


def solar_irradiance(wavelength_nm):
    """Extraterrestrial solar irradiance at 1 au [W m^-2 nm^-1]."""
    return _resample("solar_am0.csv", wavelength_nm)[..., 0]


def cie_xyz(wavelength_nm):
    """CIE 1931 2-degree colour matching functions, shape ``(..., 3)``."""
    return _resample("cie1931_2deg.csv", wavelength_nm)


def ozone_cross_section(wavelength_nm):
    """Ozone absorption cross section [m^2 molecule^-1], JPL-2006 table.

    Source values are in cm^2 per molecule, hence the factor of 1e-4.
    The visible-band table has no resolved temperature dependence.
    """
    return _resample("o3_cross_section.csv", wavelength_nm)[..., 0] * 1e-4


def refractive_index_air(wavelength_nm):
    """Edlén (1966) standard dry-air index at 288.15 K, 101325 Pa.

    This optical formula is used over 300–1700 nm only. Its standard air
    composition is historical; humidity and varying CO2 are not modelled.
    """
    wl = _wavelengths(wavelength_nm)
    if np.any(wl < 300.0) or np.any(wl > 1700.0):
        raise ValueError("The dry-air refractivity formula is restricted to 300–1700 nm.")
    inverse_um_squared = (1000.0 / wl) ** 2
    return 1.0 + 1e-8 * (
        8342.13
        + 2406030.0 / (130.0 - inverse_um_squared)
        + 15997.0 / (38.9 - inverse_um_squared)
    )


def rayleigh_cross_section(wavelength_nm):
    """Dry-air Rayleigh scattering cross section [m^2 molecule^-1].

    The first-principles expression discussed by Bodhaine et al. (1999) is
    used with Edlén refractivity and a constant depolarisation ratio 0.0279.
    This is a visible-band baseline, not the full species-dependent King
    correction recommended for precision atmospheric retrievals.
    """
    wl = _wavelengths(wavelength_nm)
    refractive_index = refractive_index_air(wl)
    n_squared = refractive_index**2
    depolarization = 0.0279
    king_factor = (6.0 + 3.0 * depolarization) / (6.0 - 7.0 * depolarization)
    return (
        24.0 * np.pi**3 / (REFERENCE_NUMBER_DENSITY**2 * (wl * 1e-9) ** 4)
        * ((n_squared - 1.0) / (n_squared + 2.0)) ** 2
        * king_factor
    )
