"""Checks against independently stated physical units and reference values."""

import numpy as np
import pytest

from eclipse.spectra import (
    cie_xyz,
    ozone_cross_section,
    rayleigh_cross_section,
    refractive_index_air,
    solar_irradiance,
)


def test_solar_visible_energy_and_absolute_scale():
    # ASTM G173 extraterrestrial 400–700 nm power is approximately 530 W/m2.
    # This detects nm/um or W/m2/relative-spectrum mistakes.
    wl = np.arange(400.0, 701.0)
    assert np.trapezoid(solar_irradiance(wl), wl) == pytest.approx(530.0, rel=0.001)
    assert 1.5 < solar_irradiance(550.0) < 2.1


def test_cie_photopic_peak_and_axes():
    xyz = cie_xyz(np.array([450.0, 555.0, 650.0]))
    assert xyz.shape == (3, 3)
    assert xyz[1, 1] == pytest.approx(1.0, abs=1e-8)
    assert xyz[0, 2] > xyz[0, 0]
    assert xyz[2, 0] > xyz[2, 2]


def test_ozone_units_and_chappuis_shape():
    # The source table is in cm2; 600 nm is 5.13e-21 cm2/molecule.
    sigma = ozone_cross_section(np.array([400.0, 600.0, 780.0]))
    assert sigma[1] == pytest.approx(5.13e-25, rel=1e-12)
    assert sigma[1] > 20 * sigma[0]
    assert sigma[1] > 10 * sigma[2]


def test_air_reference_index_and_rayleigh_cross_section():
    assert refractive_index_air(550.0) == pytest.approx(1.00027782, abs=1e-8)
    sigma = rayleigh_cross_section(np.array([400.0, 550.0, 700.0]))
    assert sigma[1] == pytest.approx(4.51e-31, rel=0.005)
    # Rayleigh lambda^-4 behaviour, including the small refractivity dispersion.
    assert 9.0 < sigma[0] / sigma[2] < 11.0
    assert np.all(np.diff(sigma) < 0)


def test_input_shapes_and_no_silent_extrapolation():
    assert cie_xyz(555.0).shape == (3,)
    assert solar_irradiance(np.array([[500.0, 550.0]])).shape == (1, 2)
    for function in (solar_irradiance, cie_xyz, ozone_cross_section):
        for bad in ([299.0], [831.0], [float("nan")]):
            with pytest.raises(ValueError):
                function(bad)
    for bad in ([0.0], [299.0], [1701.0], [float("inf")]):
        with pytest.raises(ValueError):
            refractive_index_air(bad)
