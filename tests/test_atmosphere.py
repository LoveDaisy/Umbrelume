"""US1976 anchors, hydrostatic consistency, and explicit ozone column checks."""

import numpy as np
import pytest
from scipy.integrate import quad

from eclipse.atmosphere import (
    Atmosphere,
    DOBSON_UNIT,
    LAYER_BOUNDARIES_M,
    REFERENCE_NUMBER_DENSITY,
)
from eclipse.spectra import rayleigh_cross_section


def test_us1976_reference_number_densities():
    atmosphere = Atmosphere()
    assert atmosphere.number_density(0.0) == pytest.approx(REFERENCE_NUMBER_DENSITY)
    # Independent US1976 geometric-altitude table anchors (T, P).
    # At 20 km: 216.65 K and 5529.3 Pa. At 30 km: 226.509 K, 1197.0 Pa.
    expected = np.array([5529.3 / (1.380649e-23 * 216.65),
                         1197.0 / (1.380649e-23 * 226.509)])
    np.testing.assert_allclose(atmosphere.number_density([20000.0, 30000.0]),
                               expected, rtol=5e-5)


def test_gradient_matches_independent_finite_difference():
    atmosphere = Atmosphere()
    altitude = np.array([0.5, 6000.0, 15000.0, 25000.0, 40000.0, 60000.0, 100000.0])
    finite_difference = (
        atmosphere.number_density(altitude + 0.1)
        - atmosphere.number_density(altitude - 0.1)
    ) / 0.2
    np.testing.assert_allclose(atmosphere.density_gradient(altitude),
                               finite_difference, rtol=1e-8)


def test_layer_density_continuity_and_monotonicity():
    atmosphere = Atmosphere()
    boundaries = LAYER_BOUNDARIES_M[1:]
    np.testing.assert_allclose(atmosphere.number_density(boundaries - 1e-4),
                               atmosphere.number_density(boundaries + 1e-4), rtol=5e-8)
    assert np.all(np.diff(atmosphere.number_density(np.linspace(0, 120000, 1000))) < 0)


def test_exact_layer_boundaries_use_upper_side_density_gradient():
    atmosphere = Atmosphere()
    boundaries = LAYER_BOUNDARIES_M[1:]
    np.testing.assert_allclose(atmosphere.density_gradient(boundaries),
                               atmosphere.density_gradient(boundaries + 1e-4), rtol=5e-8)


@pytest.mark.parametrize("column,centre,width", [(300.0,25000.0,8000.0),
                                               (420.0,10000.0,15000.0),
                                               (0.0,25000.0,8000.0)])
def test_ozone_profile_has_requested_column(column, centre, width):
    atmosphere = Atmosphere(ozone_du=column, ozone_center_m=centre, ozone_width_m=width)
    integral = quad(atmosphere.ozone_density, 0.0, atmosphere.top_m, epsabs=1.0)[0]
    assert integral == pytest.approx(column * DOBSON_UNIT, rel=1e-10)


def test_density_and_ozone_controls_are_independent_and_top_is_explicit():
    standard = Atmosphere()
    half_air = Atmosphere(density_scale=0.5)
    altitude = np.array([0.0, 25000.0, 120001.0])
    np.testing.assert_allclose(half_air.number_density(altitude),
                               standard.number_density(altitude) * 0.5)
    np.testing.assert_array_equal(half_air.ozone_density(altitude),
                                  standard.ozone_density(altitude))
    assert standard.number_density(120001.0) == 0.0
    assert standard.density_gradient(120001.0) == 0.0
    assert standard.ozone_density(120001.0) == 0.0


def test_vertical_rayleigh_optical_depth_has_physical_scale():
    atmosphere = Atmosphere()
    column = quad(atmosphere.number_density, 0.0, atmosphere.top_m,
                  points=LAYER_BOUNDARIES_M[1:].tolist(), epsabs=1e12, epsrel=1e-9)[0]
    tau_550 = column * rayleigh_cross_section(550.0)
    assert 0.095 < tau_550 < 0.100


@pytest.mark.parametrize("kwargs", [{"ozone_du":-1}, {"density_scale":-1},
                                    {"top_m":0}, {"ozone_width_m":0},
                                    {"ozone_center_m":-1}, {"ozone_center_m":200000},
                                    {"density_scale":float("nan")}])
def test_invalid_atmosphere_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        Atmosphere(**kwargs)


def test_invalid_altitudes_rejected():
    for bad in (-0.1, float("inf"), float("nan")):
        for method in (Atmosphere().number_density, Atmosphere().density_gradient,
                       Atmosphere().ozone_density):
            with pytest.raises(ValueError):
                method(bad)
