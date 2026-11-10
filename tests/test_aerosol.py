"""Physical column, source anchors, and wavelength-independent normalization."""
import json
from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from eclipse.aerosol import (
    FREE_TROPOSPHERE_SOURCE, OPAC_PDF_SHA256, PRIOR_PARAMETERS_JSON,
    STRATOSPHERE_BACKGROUND_SOURCE, AerosolColumn, AerosolProfile,
)

ROOT = Path(__file__).resolve().parents[1]
PRIOR_PATH = ROOT/PRIOR_PARAMETERS_JSON


def prior_row(candidate, name):
    """One pinned row of prior_parameters.json; a missing file fails loud."""
    rows = json.loads(PRIOR_PATH.read_text())["parameters"]
    return next(row for row in rows
                if row["candidate"] == candidate and row["name"] == name)


def test_published_opac_free_troposphere_parameters():
    profile = AerosolProfile.opac_free_troposphere()
    assert (profile.aod550, profile.bottom_m, profile.top_m, profile.scale_height_m) == (
        0.013, 2000.0, 12000.0, 8000.0)
    assert (profile.angstrom_visible, profile.angstrom_blue, profile.break_nm) == (
        1.58, 1.21, 500.0)
    # Independent numerical value from the published optical-depth/scale-height
    # relation: about 7.39e-4 km^-1, not 7.39e-4 m^-1.
    assert profile.extinction550(11000.0)*1000 == pytest.approx(7.394e-4, rel=1e-3)


def test_published_opac_stratospheric_background_parameters():
    profile = AerosolProfile.opac_stratospheric_background()
    assert (profile.aod550, profile.bottom_m, profile.top_m, profile.scale_height_m) == (
        0.005, 12000.0, 35000.0, 99000.0)
    assert (profile.angstrom_visible, profile.angstrom_blue, profile.break_nm) == (
        1.14, 0.74, 500.0)
    # Homogeneous-layer representation: the 23 km thickness against the 99 km
    # scale height retains about 79% of the bottom extinction at the layer top.
    ratio = profile.extinction550(35_000.0)/profile.extinction550(12_000.0)
    assert ratio == pytest.approx(np.exp(-23_000.0/99_000.0), rel=1e-14)
    # The vertical integral is the published 0.005.
    column = quad(profile.extinction550, profile.bottom_m, profile.top_m,
                  epsabs=1e-12, epsrel=1e-12)[0]
    assert column == pytest.approx(0.005, abs=1e-12)


def test_priors_are_anchored_to_pinned_prior_parameter_rows():
    """Parameter hard constraint: factories take only the pinned json rows.

    Every factory value must equal its prior_parameters.json row value and
    every cited row name must exist there; the OPAC sha constant must equal
    the sha recorded in sources.json. Reading the files fails loud if they
    move (a01: independence of the anchor, not docstring self-agreement).
    """
    strat = AerosolProfile.opac_stratospheric_background()
    rows = STRATOSPHERE_BACKGROUND_SOURCE["rows"]
    assert STRATOSPHERE_BACKGROUND_SOURCE["candidate"] == "stratosphere_background"
    pinned = {name: prior_row("stratosphere_background", name)["value"] for name in rows}
    assert strat.aod550 == pinned["aod550"]
    assert strat.angstrom_blue == pinned["angstrom_0.35_to_0.5um"]
    assert strat.angstrom_visible == pinned["angstrom_0.5_to_0.8um"]
    assert strat.bottom_m == pinned["layer_bottom_m"]
    assert strat.top_m == pinned["layer_top_m"]
    assert strat.scale_height_m == pinned["scale_height_m"]

    ft = AerosolProfile.opac_free_troposphere()
    ft_rows = FREE_TROPOSPHERE_SOURCE["rows"]
    assert FREE_TROPOSPHERE_SOURCE["candidate"] == "free_troposphere_reference"
    ft_pinned = {name: prior_row("free_troposphere_reference", name)["value"]
                 for name in ft_rows}
    assert ft.aod550 == ft_pinned["aod550_2_to_12km_reference"]
    assert ft.angstrom_blue == ft_pinned["angstrom_0.35_to_0.5um"]
    assert ft.angstrom_visible == ft_pinned["angstrom_0.5_to_0.8um"]
    assert ft.scale_height_m == ft_pinned["scale_height_m"]

    sources = json.loads((PRIOR_PATH.parent/"sources.json").read_text())
    opac = next(item for item in sources if item["file_name"] == "opac1998.pdf")
    assert opac["sha256"] == OPAC_PDF_SHA256


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


def test_column_vertical_aod_is_the_sum_of_layer_aods():
    column = AerosolColumn.opac_default()
    total = quad(lambda h: sum(layer.extinction550(h) for layer in column.layers),
                 column.layers[0].bottom_m, column.layers[1].top_m,
                 epsabs=1e-12, epsrel=1e-12)[0]
    assert total == pytest.approx(sum(layer.aod550 for layer in column.layers), abs=1e-12)


def test_column_support_and_zero_outside_all_layers():
    strat = AerosolProfile.opac_stratospheric_background()
    column = AerosolColumn.opac_default()
    altitude = np.array([0.0, 1999.9, 12_000.1, 34_999.9, 35_000.1, 50_000.0])
    total = sum(layer.extinction550(altitude) for layer in column.layers)
    # 0-2 km is bare (both layers start at 2/12 km); 2-12 km carries the free
    # troposphere; 12-35 km the stratosphere; above 35 km exactly zero.
    np.testing.assert_array_equal(total > 0, np.array(
        [False, False, True, True, False, False]))
    np.testing.assert_array_equal(strat.extinction550(np.array([11_999.9, 35_000.1])),
                                  [0.0, 0.0])


def test_column_spectral_factors_stay_layer_independent():
    ft, strat = AerosolProfile.opac_free_troposphere(), \
        AerosolProfile.opac_stratospheric_background()
    grid = np.array([380., 450., 503., 550., 677., 780.])
    column = AerosolColumn((ft, strat))
    for layer, reference in zip(column.layers, (ft, strat)):
        np.testing.assert_array_equal(layer.spectral_factor(grid),
                                      reference.spectral_factor(grid))
    # The two layers' published exponents differ; no merged factor exists.
    assert not np.array_equal(ft.spectral_factor(grid), strat.spectral_factor(grid))


def test_default_column_has_the_two_production_layers():
    column = AerosolColumn.opac_default()
    assert isinstance(column, AerosolColumn) and len(column.layers) == 2
    assert column.layers[0] == AerosolProfile.opac_free_troposphere()
    assert column.layers[1] == AerosolProfile.opac_stratospheric_background()
    ft_override = AerosolColumn.opac_default(ft_aod550=0.5)
    assert ft_override.layers[0].aod550 == 0.5
    assert ft_override.layers[1] == AerosolProfile.opac_stratospheric_background()


def test_column_rejects_empty_or_non_profile_layers():
    with pytest.raises(ValueError):
        AerosolColumn(())
    with pytest.raises(ValueError):
        AerosolColumn((AerosolProfile(), "not a profile"))
    with pytest.raises(FrozenInstanceError):
        AerosolColumn.opac_default().layers = ()
