"""Geometry limits and finite-source, multi-branch transport checks."""
import numpy as np
import pytest
from scipy.integrate import quad

from eclipse.atmosphere import Atmosphere
from eclipse.rays import RayTable, build_ray_table
from eclipse.transport import (
    Geometry, circle_overlap, geometric_transmission, ring_fraction,
    signed_landing_radius, transmission_profile,
)


def test_circle_overlap_known_intersection_and_limits():
    np.testing.assert_allclose(circle_overlap([0, 1, 2, 3], 1, 1),
                               [np.pi, 2*np.pi/3-np.sqrt(3)/2, 0, 0], atol=1e-14)
    assert circle_overlap(0, 2, 1) == pytest.approx(np.pi)


def test_ring_azimuth_fraction_against_direct_angular_count():
    theta = np.linspace(0, 2*np.pi, 100_000, endpoint=False)
    for r, s, a in [(3., 2., 3.), (0., 1., 2.), (5., 1., 2.), (1., 1., 3.)]:
        fraction = np.mean((s*np.cos(theta)-r)**2 + (s*np.sin(theta))**2 < a*a)
        assert ring_fraction(r, s, a) == pytest.approx(fraction, abs=2e-5)


def test_vacuum_annuli_recover_analytic_geometric_eclipse():
    table = build_ray_table(Atmosphere(density_scale=0, ozone_du=0), height_step_m=500)
    radii = np.linspace(0, 9e6, 101)
    actual = transmission_profile(table, radii, [450, 600, 750])
    expected = geometric_transmission(radii)
    np.testing.assert_allclose(actual, np.broadcast_to(expected[:, None], actual.shape), atol=2e-7)
    np.testing.assert_array_equal(table.bending_rad, 0)
    np.testing.assert_array_equal(table.air_column_m2, 0)


def test_ozone_counterfactual_only_removes_extinction():
    table = build_ray_table(Atmosphere(), height_step_m=200)
    r = np.array([0, 3e6, 4.5e6, 5e6, 9e6])
    wl = np.array([450., 600., 750.])
    on = transmission_profile(table, r, wl)
    off = transmission_profile(table, r, wl, ozone=False)
    assert np.all(on <= off + 1e-15)
    assert on[2, 1]/off[2, 1] < on[2, 0]/off[2, 0]
    np.testing.assert_allclose(on[-1], 1, atol=1e-14)


def test_nonmonotone_landing_accumulates_every_branch():
    geometry = Geometry()
    b = np.array([6.370e6, 6.380e6, 6.390e6, 6.400e6])
    landing = np.array([-2e6, 3e6, 1e6, 6e6])
    alpha = (geometry.magnification*b-landing)/geometry.moon_distance_m
    table = RayTable(b-b[0], b, alpha, np.zeros(4), np.zeros(4), 550., True)
    for r in [0, 1e6, 4e6]:
        a, m = geometry.solar_blur_m, geometry.magnification
        expected = 1-circle_overlap(r, a, m*b[-1])/(np.pi*a*a)
        for lo, hi in zip(b[:-1], b[1:]):
            expected += quad(lambda bb: 2*m*m*bb/a**2*float(ring_fraction(
                r, abs(np.interp(bb, b, landing)), a)), lo, hi,
                epsabs=1e-9, limit=500)[0]
        actual = transmission_profile(table, [r], [550], extinction=False, integration_order=80)
        assert actual[0, 0] == pytest.approx(expected, abs=2e-7)


def test_bending_and_columns_have_physical_sign_and_scale():
    table = build_ray_table(Atmosphere(), height_step_m=100)
    assert 1.0 < np.degrees(table.bending_rad[0]) < 1.2
    assert np.all(table.bending_rad >= 0)
    assert np.all(table.air_column_m2 >= 0)
    assert np.all(table.ozone_column_m2 >= 0)
    assert table.air_column_m2[-1] == 0
    assert table.bending_rad[-1] == 0
    # Idealized standard-atmosphere layer boundaries can turn the map back;
    # transport explicitly keeps these branches instead of inverting globally.
    assert np.any(np.diff(signed_landing_radius(table, Geometry())) < 0)


@pytest.mark.parametrize("kw", [{"height_step_m":0}, {"height_step_m":float("nan")},
                              {"quadrature_order":2}, {"reference_wavelength_nm":-1}])
def test_bad_table_settings_rejected(kw):
    with pytest.raises(ValueError):
        build_ray_table(Atmosphere(), **kw)


def test_low_top_with_non_negligible_interface_is_rejected():
    with pytest.raises(ValueError, match="top refractivity"):
        build_ray_table(Atmosphere(top_m=1000, ozone_center_m=1000))


@pytest.mark.parametrize("radius", [float("nan"), float("inf"), -1, 0])
def test_invalid_circle_radii_rejected(radius):
    with pytest.raises(ValueError):
        circle_overlap([0, 1], radius, 1)
