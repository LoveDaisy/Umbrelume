"""Physical column, source anchors, and wavelength-independent normalization."""
from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from scipy.integrate import quad

from eclipse.aerosol import AerosolProfile


def test_published_opac_free_troposphere_parameters():
    profile = AerosolProfile.opac_free_troposphere()
    assert (profile.aod550, profile.bottom_m, profile.top_m, profile.scale_height_m) == (
        0.013, 2000.0, 12000.0, 8000.0)
    assert (profile.angstrom_visible, profile.angstrom_blue, profile.break_nm) == (
        1.58, 1.21, 500.0)
    # Independent numerical value from the published optical-depth/scale-height
    # relation: about 7.39e-4 km^-1, not 7.39e-4 m^-1.
    assert profile.extinction550(11000.0)*1000 == pytest.approx(7.394e-4, rel=1e-3)


@pytest.mark.parametrize("profile", [
    AerosolProfile(), AerosolProfile(aod550=0.04), AerosolProfile.none(),
    AerosolProfile(aod550=.02, bottom_m=50000, top_m=50001, scale_height_m=100000),
])
def test_vertical_integral_is_requested_aod(profile):
    column = quad(profile.extinction550, profile.bottom_m, profile.top_m,
                  epsabs=1e-12, epsrel=1e-12)[0]
    assert column == pytest.approx(profile.aod550, abs=1e-12)


def test_support_bounds_and_zero_counterfactual():
    profile = AerosolProfile()
    altitude = np.array([[0, 1999.9, 2000], [11000, 12000, 12000.1]])
    extinction = profile.extinction550(altitude)
    assert extinction.shape == altitude.shape
    assert np.all(extinction[np.array([[True, True, False], [False, False, True]])] == 0)
    assert np.all(extinction[np.array([[False, False, True], [True, True, False]])] > 0)
    np.testing.assert_array_equal(AerosolProfile.none().extinction550(altitude), 0)


def test_spectral_reference_and_slicing_are_independent_of_query_grid():
    profile = AerosolProfile()
    assert profile.spectral_factor(550) == 1
    grid = np.array([380, 503, 550, 677, 780])
    values = profile.spectral_factor(grid)
    assert profile.spectral_factor([677])[0] == values[3]
    assert profile.spectral_factor(503) == values[1]
    np.testing.assert_allclose(values[1]/values[3], (677/503)**1.58, rtol=1e-14)


def test_broken_spectrum_is_continuous_and_uses_both_published_slopes():
    profile = AerosolProfile()
    np.testing.assert_allclose(profile.spectral_factor(380)/profile.spectral_factor(500),
                               (500/380)**1.21, rtol=1e-14)
    np.testing.assert_allclose(profile.spectral_factor(500-1e-6),
                               profile.spectral_factor(500+1e-6), rtol=1e-8)
    assert np.all(np.diff(profile.spectral_factor(np.arange(380, 781))) < 0)


def test_aod_control_preserves_profile_and_spectral_shape():
    low = AerosolProfile.opac_free_troposphere(.013)
    high = AerosolProfile.opac_free_troposphere(.039)
    np.testing.assert_allclose(high.extinction550([3000, 8000, 11000]),
                               3*low.extinction550([3000, 8000, 11000]))
    np.testing.assert_array_equal(high.spectral_factor([380, 503, 677]),
                                  low.spectral_factor([380, 503, 677]))


def test_custom_break_still_has_exact_physical_550nm_reference():
    assert AerosolProfile(break_nm=600).spectral_factor(550) == 1


def test_frozen_profile():
    with pytest.raises(FrozenInstanceError):
        AerosolProfile().aod550 = .1


@pytest.mark.parametrize("kwargs", [
    {"aod550": -1}, {"aod550": float("nan")}, {"bottom_m": -1},
    {"top_m": 2000}, {"scale_height_m": 0}, {"scale_height_m": float("inf")},
    {"angstrom_visible": -1}, {"angstrom_blue": -1}, {"break_nm": 0},
])
def test_invalid_profile_parameters(kwargs):
    with pytest.raises(ValueError):
        AerosolProfile(**kwargs)


@pytest.mark.parametrize("altitude", [-1, float("nan"), float("inf")])
def test_invalid_altitudes(altitude):
    with pytest.raises(ValueError):
        AerosolProfile().extinction550(altitude)


@pytest.mark.parametrize("wavelength", [0, -1, float("nan"), float("inf")])
def test_invalid_wavelengths(wavelength):
    with pytest.raises(ValueError):
        AerosolProfile().spectral_factor(wavelength)
