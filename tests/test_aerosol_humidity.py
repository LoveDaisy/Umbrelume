"""Humidity-axis tests (task-humidity-axis, plan Step 1 / M1).

Gate scope: the pinned WASO knots, the q_rest derivation, the bit-exact 50%
reference, the domain guard, the single growth-application path, the
geometry-premise gate (J6), and the transcription product's numeric contract
against humidity_growth.json. All humidity tests live in this one file so the
component stays removable: deleting the file together with the humidity
section of eclipse/aerosol.py restores the previous state exactly.
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from eclipse.aerosol import (
    FT_AOD550_PUBLISHED, FT_BOTTOM_KM, FT_GEOMETRY_FACTOR_KM, FT_SCALE_HEIGHT_KM,
    FT_TOP_KM, HUMIDITY_AXIS_PRIORS_JSON, HUMIDITY_GROWTH_JSON, OPAC_PDF_SHA256,
    Q_REST, RH_MAX, RH_MIN, WASO_EXT550_PER_PARTICLE_KM, WASO_EXT550_SHA256,
    WASO_NUMBER_MIXING_RATIO, FT_SEA_LEVEL_NUMBER_DENSITY_PER_CM3,
    AerosolColumn, AerosolProfile, ft_growth_factor, ft_humidity_provenance,
    _waso_column_aod550,
)

ROOT = Path(__file__).resolve().parents[1]
GROWTH = json.loads((ROOT/HUMIDITY_GROWTH_JSON).read_text())
PRIORS = json.loads((ROOT/HUMIDITY_AXIS_PRIORS_JSON).read_text())


def _waso_row(rh_percent: int) -> dict:
    component = next(item for item in GROWTH["components"]
                     if item["component"] == "WASO")
    return component["per_rh"][str(rh_percent)]


# --- Pinned knots and transcription anchors -------------------------------------


def test_waso_knots_anchored_to_humidity_growth_product():
    """Every pinned knot value and sha256 must equal the committed survey
    product's row (a01: the anchor is the independent product, not the
    module's own constants)."""
    for rh, value in WASO_EXT550_PER_PARTICLE_KM.items():
        row = _waso_row(int(round(rh*100)))
        assert value == row["ext_coef_550nm_per_particle_1km"], rh
        assert WASO_EXT550_SHA256[rh] == row["source"]["sha256"], rh
        assert row["source"]["file_name"] == f"waso{int(round(rh*100))}"


def test_humidity_section_reuses_the_opac_paper_constant():
    """The humidity provenance must reference the module's single OPAC PDF
    sha256 constant (whose equality with sources.json is anchored in
    test_aerosol.py) - the humidity section reuses, never redefines (r2
    Suggestion 1)."""
    provenance = ft_humidity_provenance(0.7, FT_AOD550_PUBLISHED)
    assert provenance is not None
    assert "opac1998_sha256" not in provenance  # added by the metadata merge
    assert PRIORS["opac_pdf"]["sha256"] == OPAC_PDF_SHA256
    assert WASO_EXT550_SHA256[0.5] == _waso_row(50)["source"]["sha256"]


def test_ft_geometry_matches_repo_layer():
    """J6 geometry-premise gate: the paper's FT geometry used by the growth
    math must equal the repository free-troposphere layer's pinned geometry."""
    ft = AerosolProfile.opac_free_troposphere()
    assert FT_BOTTOM_KM == ft.bottom_m/1000.0
    assert FT_TOP_KM == ft.top_m/1000.0
    assert FT_SCALE_HEIGHT_KM == ft.scale_height_m/1000.0
    assert FT_BOTTOM_KM == PRIORS["table3_rh_semantics"]["ft_row"]["layer_km"][0]
    assert FT_TOP_KM == PRIORS["table3_rh_semantics"]["ft_row"]["layer_km"][1]


# --- q_rest and growth factor -----------------------------------------------------


def _recompute_e_waso(rh: float) -> float:
    """Independent recompute of E_waso from the stated formula and knots."""
    knots = (0.5, 0.7, 0.99)
    values = tuple(WASO_EXT550_PER_PARTICLE_KM[knot] for knot in knots)
    ext550 = values[0]
    for lo, hi, lo_value, hi_value in zip(knots, knots[1:], values, values[1:]):
        if rh <= hi:
            weight = (rh-lo)/(hi-lo)
            ext550 = lo_value + weight*(hi_value-lo_value)
            break
    geometry = 8.0*(math.exp(-2.0/8.0) - math.exp(-12.0/8.0))
    return WASO_NUMBER_MIXING_RATIO*FT_SEA_LEVEL_NUMBER_DENSITY_PER_CM3*ext550*geometry


def test_q_rest_gate_and_derivation_recompute():
    """The q_rest gate and its derivation must hold and recompute from the
    pinned knots (import-time gate mirrored here; plan risk 1)."""
    assert 0.0 < Q_REST < FT_AOD550_PUBLISHED
    e50 = _recompute_e_waso(0.5)
    assert Q_REST == FT_AOD550_PUBLISHED-e50
    assert Q_REST == PRIORS["q_rest_derivation"]["q_rest"]
    assert e50 == PRIORS["q_rest_derivation"]["e_waso_50"]
    assert FT_GEOMETRY_FACTOR_KM == PRIORS["q_rest_derivation"]["geometry_factor_km"]
    assert Q_REST/FT_AOD550_PUBLISHED == \
        PRIORS["q_rest_derivation"]["q_rest_fraction_of_reference"]
    # the composition constants must equal the preregistered product values
    assert WASO_NUMBER_MIXING_RATIO == \
        PRIORS["ft_composition"]["number_mixing_ratios"]["waso"]
    assert FT_SEA_LEVEL_NUMBER_DENSITY_PER_CM3 == \
        PRIORS["ft_composition"]["sea_level_number_density_per_cm3"]


def test_ft_growth_factor_reference_is_bit_exact():
    """The default bit-exactness contract: f(0.5) is exactly 1.0 and the
    factories' reference-state output equals the historical values."""
    assert ft_growth_factor(0.5) == 1.0
    assert ft_growth_factor(RH_MIN) == PRIORS["ft_growth_factor_pinned"]["0.5"]
    column = AerosolColumn.opac_default()
    assert column.layers[0].aod550 == 0.013
    assert column.layers[0] == AerosolProfile.opac_free_troposphere()
    assert AerosolProfile.opac_free_troposphere(rh=0.5).aod550 == 0.013
    assert AerosolColumn.opac_default(rh=0.5).layers[0].aod550 == 0.013


def test_ft_growth_factor_pinned_values_match_product():
    """The pinned growth values must equal the transcription product's
    records (read-back gate; the product carries the exact floats)."""
    for rh in (0.5, 0.7, 0.99):
        assert ft_growth_factor(rh) == PRIORS["ft_growth_factor_pinned"][f"{rh:g}"]
    assert ft_growth_factor(0.7) == pytest.approx(1.2317, rel=1e-4)
    assert ft_growth_factor(0.99) == pytest.approx(7.1459, rel=1e-4)


def test_ft_growth_factor_monotone_on_pinned_domain():
    grid = [RH_MIN + (RH_MAX-RH_MIN)*i/98 for i in range(99)]
    factors = [ft_growth_factor(rh) for rh in grid]
    assert all(later >= earlier for earlier, later in zip(factors, factors[1:]))
    assert factors[0] == 1.0 and factors[-1] == ft_growth_factor(RH_MAX)
    # strictly increasing between the knots (the growth is real everywhere)
    assert all(later > earlier for earlier, later in zip(factors, factors[1:]))


@pytest.mark.parametrize("rh", [0.0, 0.1, 0.49, 0.499, 1.0, 1.5, float("nan"),
                                float("inf"), -0.5])
def test_ft_growth_factor_domain_guard(rh):
    with pytest.raises(ValueError):
        ft_growth_factor(rh)
    with pytest.raises(ValueError):
        AerosolProfile.opac_free_troposphere(rh=rh)
    with pytest.raises(ValueError):
        AerosolColumn.opac_default(rh=rh)


def test_growth_changes_only_the_column_aod():
    """J1: the humidity growth rescales only the column AOD - vertical shape
    and spectral law are untouched (frozen law is a preregistered invariant)."""
    reference = AerosolProfile.opac_free_troposphere()
    humid = AerosolProfile.opac_free_troposphere(rh=0.99)
    factor = ft_growth_factor(0.99)
    ratio = humid.extinction550(np.array([3000.0, 8000.0, 11000.0])) \
        / reference.extinction550(np.array([3000.0, 8000.0, 11000.0]))
    # the per-altitude ratio reproduces the growth factor to the last ulps
    # (the profile evaluates aod550/denominator*exp, whose rounding sequence
    # differs from the factor's own chain by <= 1 ulp)
    np.testing.assert_allclose(ratio, np.full(3, factor), rtol=1e-14)
    np.testing.assert_array_equal(humid.spectral_factor(np.array([380., 503., 677.])),
                                  reference.spectral_factor(np.array([380., 503., 677.])))


def test_single_growth_application_path():
    """J2 composition rule: with an explicit ``ft_aod550`` base and ``rh`` the
    effective FT AOD is exactly ``ft_aod550*ft_growth_factor(rh)`` - the growth
    is applied once, at the factory (the caller never pre-multiplies)."""
    base, rh = 0.5, 0.7
    column = AerosolColumn.opac_default(ft_aod550=base, rh=rh)
    assert column.layers[0].aod550 == base*ft_growth_factor(rh)
    # the stratosphere stays at its published state and the FT spectral law
    # is unchanged
    assert column.layers[1] == AerosolProfile.opac_stratospheric_background()
    assert (column.layers[0].angstrom_visible, column.layers[0].angstrom_blue) \
        == (1.58, 1.21)
    ft = AerosolProfile.opac_free_troposphere(aod550=base, rh=rh)
    assert ft.aod550 == base*ft_growth_factor(rh)


def test_humidified_ray_table_structural_identity():
    """The kernel-path identity the A/B re-inversion relies on: the humidified
    default column's ray table keeps the production grid bitwise, leaves the
    stratosphere tau bitwise unchanged, and scales the FT tau linearly with
    the growth factor (well within 1e-12)."""
    from eclipse.atmosphere import Atmosphere
    from eclipse.rays import build_ray_table
    production = build_ray_table(Atmosphere(),
                                 aerosol_profile=AerosolColumn.opac_default())
    for rh in (0.5, 0.7, 0.99):
        table = build_ray_table(
            Atmosphere(), aerosol_profile=AerosolColumn.opac_default(rh=rh))
        np.testing.assert_array_equal(table.impact_m, production.impact_m)
        np.testing.assert_array_equal(table.altitude_m, production.altitude_m)
        np.testing.assert_array_equal(table.air_column_m2, production.air_column_m2)
        np.testing.assert_array_equal(table.ozone_column_m2, production.ozone_column_m2)
        np.testing.assert_array_equal(table.aerosol_layer_taus[1],
                                      production.aerosol_layer_taus[1])
        factor = ft_growth_factor(rh)
        deviation = np.abs(table.aerosol_layer_taus[0]
                           - factor*production.aerosol_layer_taus[0])
        scale = np.maximum(np.abs(production.aerosol_layer_taus[0]), 1e-300)
        assert float(np.max(deviation/scale)) < 1e-12, rh
        if rh == 0.5:
            # full bitwise identity at the reference state
            np.testing.assert_array_equal(table.aerosol_layer_taus[0],
                                          production.aerosol_layer_taus[0])


def test_ft_humidity_provenance_shape():
    provenance = ft_humidity_provenance(0.99, 0.013)
    assert provenance["rh"] == 0.99
    assert provenance["ft_growth_factor"] == ft_growth_factor(0.99)
    assert provenance["ft_aod550_reference_base"] == 0.013
    assert provenance["humidity_growth_json"] == HUMIDITY_GROWTH_JSON
    assert provenance["humidity_axis_priors_json"] == HUMIDITY_AXIS_PRIORS_JSON
    for rh, sha256 in WASO_EXT550_SHA256.items():
        assert provenance["waso_ext550_sha256"][f"{rh:g}"] == sha256


def test_priors_product_transcription_contract():
    """The transcription product's quotes and pinned rows must agree with the
    committed survey products (values + sha256 per knot; a01 independence)."""
    priors_knots = PRIORS["waso_growth_knots"]["per_rh"]
    for rh_percent, row in (("50", "0.5"), ("70", "0.7"), ("99", "0.99")):
        survey_row = _waso_row(int(rh_percent))
        assert priors_knots[row]["ext_coef_550nm_per_particle_1km"] \
            == survey_row["ext_coef_550nm_per_particle_1km"]
        assert priors_knots[row]["sha256"] == survey_row["source"]["sha256"]
    assert PRIORS["ft_composition"]["quote_sec5a_p842"].startswith(
        "The aerosol in the free troposphere is modeled")
    assert "730 particles cm-3" in PRIORS["ft_composition"]["quote_sec5a_p842"]
    assert PRIORS["table3_rh_semantics"]["quote_sec5a_p842"].startswith(
        "They are valid for a relative humidity of 80%")
    assert PRIORS["ft_vertical_profile_evidence"]["sec5b_upper_boundary"][
        "quote"].startswith("For the proposed aerosol types")
    assert PRIORS["opac_pdf"]["sha256"] == OPAC_PDF_SHA256


# --- CLI --rh (plan Step 2) -------------------------------------------------------


CLI_COMMON = ["--no-plots", "--radial-step-km", "1000", "--wavelength-step-nm",
              "50", "--height-step-m", "500"]


def _run_cli(tmp_path, *extra):
    from eclipse.cli import main
    out = tmp_path/"out"
    main(["simulate", "--out", str(out), *CLI_COMMON, *extra])
    metadata = json.loads((out/"metadata.json").read_text())
    return metadata, np.load(out/"rays.npz")


def test_cli_default_metadata_carries_reference_humidity(tmp_path):
    """The default run's FT source record carries the humidity provenance at
    the 50% reference: rh=0.5, factor exactly 1.0, effective AOD unchanged -
    schema extension only, numerical sections untouched."""
    metadata, rays = _run_cli(tmp_path)
    ft_source = metadata["aerosol"]["layers"][0]["source"]
    assert ft_source["rh"] == 0.5
    assert ft_source["ft_growth_factor"] == 1.0
    assert ft_source["ft_aod550_reference_base"] == 0.013
    assert ft_source["humidity_growth_json"] == HUMIDITY_GROWTH_JSON
    assert ft_source["humidity_axis_priors_json"] == HUMIDITY_AXIS_PRIORS_JSON
    for rh, sha256 in WASO_EXT550_SHA256.items():
        assert ft_source["waso_ext550_sha256"][f"{rh:g}"] == sha256
    assert metadata["aerosol"]["layers"][0]["aod550"] == 0.013
    assert metadata["aerosol"]["layers"][1]["aod550"] == 0.005
    # the stratosphere source record gains no humidity keys
    assert "rh" not in metadata["aerosol"]["layers"][1]["source"]
    assert "humidity" not in metadata["approximations"][3]
    assert rays["aerosol_layer_taus"].shape[0] == 2


def test_cli_explicit_rh_grows_free_troposphere_only(tmp_path):
    factor = ft_growth_factor(0.7)
    metadata, rays = _run_cli(tmp_path, "--rh", "0.7")
    layers = metadata["aerosol"]["layers"]
    assert layers[0]["aod550"] == 0.013*factor
    assert layers[0]["source"]["rh"] == 0.7
    assert layers[0]["source"]["ft_growth_factor"] == factor
    assert layers[1]["aod550"] == 0.005
    assert "rh" not in layers[1]["source"]
    assert "humidity growth at rh=0.7" in metadata["approximations"][3]
    # the FT layer tau carries the growth, the stratosphere tau is unchanged:
    # the table shape and landing map stay identical
    assert rays["aerosol_layer_taus"].shape[0] == 2


def test_cli_opac_ft_accepts_rh(tmp_path):
    factor = ft_growth_factor(0.9)
    metadata, rays = _run_cli(tmp_path, "--aerosol", "opac-ft", "--rh", "0.9")
    layers = metadata["aerosol"]["layers"]
    assert len(layers) == 1
    assert layers[0]["aod550"] == 0.013*factor
    assert layers[0]["source"]["rh"] == 0.9
    assert rays["aerosol_layer_taus"].shape[0] == 1


@pytest.mark.parametrize("argv", [
    ("--aerosol", "calibrated:june2011", "--rh", "0.7"),
    ("--aerosol", "calibrated:june2011", "--rh", "0.5"),  # explicit reference too
    ("--aerosol", "none", "--rh", "0.7"),
    ("--aerosol", "opac", "--aod550", "0", "--rh", "0.7"),
    ("--aerosol", "opac", "--rh", "0.49"),
    ("--aerosol", "opac", "--rh", "1.0"),
    ("--aerosol", "opac", "--rh", "nan"),
])
def test_cli_rh_illegal_combinations_exit_loud(tmp_path, argv):
    from eclipse.cli import main
    with pytest.raises(SystemExit):
        main(["simulate", "--out", str(tmp_path/"out"), "--no-plots", *argv,
              *CLI_COMMON])
