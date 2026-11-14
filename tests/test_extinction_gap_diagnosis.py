"""Offline tests for the extinction-gap diagnosis machine (plan Step 1 / M1).

The synthetic-table zeroth-moment gate and the June 2011 full-scale anchor are
the two non-negotiable pytest residents (plan M1); the full-scale four-aperture
gates stay in the gate too (measured cost ~10 s, budget 60 s).
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere, EARTH_RADIUS_M
from eclipse.rays import RayTable, build_ray_table
from eclipse.transport import Geometry, transmission_profile

from validation import extinction_gap_diagnosis as egd
from docs_references import resolve_product_fields

ROOT = Path(__file__).resolve().parents[1]
PRIOR_PATH = ROOT/"validation/results/extinction_gap_diagnosis/prior_parameters.json"


def synthetic_ray_table():
    """Small monotone table whose landing map reaches umbral receivers.

    bending_rad maps the landing radius to ``b - R + 50 km`` so the default
    geometry's receivers see table rays; the aerosol column exercises the
    aerosol branch of kernel_weights against transmission_profile.
    """
    altitude = np.linspace(0., 120_000., 129)
    impact = EARTH_RADIUS_M + altitude
    geometry = Geometry()
    bending = (geometry.magnification*impact - (impact - EARTH_RADIUS_M + 50_000.)) \
        / geometry.moon_distance_m
    return RayTable(altitude, impact, bending,
                    8e28*np.exp(-altitude/8_000.), 3e22*np.exp(-altitude/20_000.),
                    550.0, False, (0.05*np.exp(-altitude/10_000.),),
                    (egd.DEFAULT_OPAC,))


def test_synthetic_kernel_reproduces_transmission_profile():
    """Gate logic resident: sum(weights) + vacuum == transmission_profile."""
    table = synthetic_ray_table()
    geometry = Geometry()
    radius = 1_000_000.0
    wavelength = 550.0
    bb_quad, area = egd.receiver_quadrature(table, radius, geometry)
    weights = egd.kernel_weights(table, bb_quad, area, wavelength)
    moment = float(np.sum(weights)) + egd.vacuum_branch(table, radius, geometry)
    reference = float(transmission_profile(table, [radius], [wavelength],
                                           geometry=geometry)[0, 0])
    assert abs(moment - reference) <= egd.MOMENT_TOLERANCE*abs(reference)
    stats = egd.altitude_statistics(table, [bb_quad], [weights])
    assert stats["p10_p50_p90_height_m"][0] <= stats["p10_p50_p90_height_m"][1] \
        <= stats["p10_p50_p90_height_m"][2]
    fractions = [stats["band_fraction_below_2km"], stats["band_fraction_2_to_12km"],
                 stats["band_fraction_over_12km"]]
    assert all(value >= 0 for value in fractions)
    assert sum(fractions) == pytest.approx(1., abs=1e-12)


def test_extra_aerosol_zero_restores_kernel_exactly():
    """Counterfactual interface: a zero-AOD extra column changes nothing."""
    table = synthetic_ray_table()
    bb_quad, area = egd.receiver_quadrature(table, 1_000_000.0, Geometry())
    zero_profile = AerosolProfile(aod550=0.0)
    extra = (table.impact_m, np.zeros_like(table.impact_m), zero_profile)
    base = egd.kernel_weights(table, bb_quad, area, np.array([450., 677.]))
    assert np.array_equal(base, egd.kernel_weights(table, bb_quad, area,
                                                   np.array([450., 677.]), extra))


def test_extra_aerosol_darkens_everywhere():
    """Sign discipline: more low-layer extinction must dim the model."""
    table = synthetic_ray_table()
    bb_quad, area = egd.receiver_quadrature(table, 1_000_000.0, Geometry())
    probe_profile = AerosolProfile(aod550=0.01, bottom_m=0., top_m=2_000.,
                                   scale_height_m=2_000.)
    probe_table = build_ray_table(Atmosphere(), aerosol_profile=probe_profile)
    extra = (probe_table.impact_m, probe_table.aerosol_layer_taus[0], probe_profile)
    wavelengths = np.array([450., 677.])
    base = egd.kernel_weights(table, bb_quad, area, wavelengths)
    probed = egd.kernel_weights(table, bb_quad, area, wavelengths, extra)
    # The 0-2 km extra column cannot touch tangent heights above the layer.
    assert np.all(probed <= base)
    assert np.any(probed < base)
    darkened = np.log(probed[probed < base]/base[probed < base])
    assert np.all(darkened < 0)


@pytest.mark.parametrize("index", range(len(egd.june_records())))
def test_june_anchor_matches_round2_committed_quantiles(index):
    """Full-scale anchor resident: exact equality with the committed values."""
    record = egd.june_records()[index]
    geometry = egd.june_geometry()
    radius = geometry.moon_distance_m*np.tan(np.deg2rad(record["shadow_angle_deg"]))
    computed = egd.june_anchor_quantiles(radius, record["wavelength_nm"])
    assert computed == record["molecular_tangent_height_p10_p50_p90_m"]


def test_full_scale_zeroth_moment_gates_and_vacuum_branches():
    """All four apertures x two variants reproduce the committed predictions."""
    stats = egd.build_height_stats()
    assert len(stats["records"]) == 18*2 + 2*2 + 7*2 + 2
    for record in stats["records"]:
        assert record["relative_gap"] <= egd.MOMENT_TOLERANCE
        # Physically closed vacuum branch carries only float noise (the
        # geometric branch opens at order-unity overlap, never gradually).
        assert record["vacuum_branch_fraction"] <= 1e-9
    for record in stats["records"]:
        if "per_wavelength_max_relative_gap" in record:
            assert record["per_wavelength_max_relative_gap"] <= egd.MOMENT_TOLERANCE
        assert record["p10_p50_p90_height_m"][0] <= record["p10_p50_p90_height_m"][1] \
            <= record["p10_p50_p90_height_m"][2]
        fractions = [record["band_fraction_below_2km"], record["band_fraction_2_to_12km"],
                     record["band_fraction_over_12km"]]
        assert all(value >= 0 for value in fractions)
        assert sum(fractions) == pytest.approx(1., abs=1e-12)


def test_correlation_event_level_and_hst_cross_check():
    """Event-level aggregation and the committed-median cross-check (plan Step 3)."""
    stats = egd.build_height_stats()
    correlation = egd.build_correlation(stats)
    assert set(correlation["event_level"]) == set(egd.EVENT_ORDER)
    for variant in egd.VARIANTS:
        for event in egd.EVENT_ORDER:
            row = correlation["event_level"][event][variant]
            assert row["residual_readout_count"] >= 1
            assert row["p50_height_m_median"] > 0
        assert correlation["event_level"]["event_1990_02_09"][variant]["shape_axis"] is None
        assert correlation["event_level"]["june_2011"][variant]["shape_axis"]["paired_angle_count"] == 5
        pairwise = correlation["cross_event"][variant]["pairwise_direction"]
        assert pairwise["pairs_total"] == 6
    for check in correlation["hst_consistency_check"].values():
        assert check["pooled_window_median"] == check["committed_model_to_observed_median"]
    # Result pins frozen in committed products: the molecular variant's depth
    # ordering is perfectly rank-consistent (deeper paths carry larger
    # residuals); any change routes as a finding, not a silent pass.
    assert correlation["cross_event"]["molecular"]["spearman_rho_residual_ln_vs_p50_median"] == -1.0
    assert correlation["cross_event"]["molecular"]["pairwise_direction"]["deeper_has_larger_residual_count"] == 6


def test_prior_parameters_schema_every_numeric_row_sourced():
    """Literature gate: no number without a pinned source (plan Step 4).

    Delegates to the command-side ``validate_literature_products`` (single
    authority): the CLI and this test enforce the identical schema walk.
    """
    checked = egd.validate_literature_products()
    assert checked["prior_parameters.json"] >= 40
    prior = json.loads(PRIOR_PATH.read_text())
    rows = prior["parameters"]
    for row in rows:
        if isinstance(row["value"], (int, float)) and not isinstance(row["value"], bool):
            assert row["source"], f"{row['candidate']}/{row['name']}: numeric row without source"
    sources = json.loads((PRIOR_PATH.parent/"sources.json").read_text())
    opac = next(item for item in sources if item["file_name"] == "opac1998.pdf")
    assert opac["sha256"] == "68a20ba9488db1e0f65a9976d4ee46e6c879f89054fcff583ba67bb5d68d425c"


def test_prior_boundary_layer_aod_derivation_consistency():
    """Recorded derived BL AOD reproduces the documented subtraction formula."""
    prior = json.loads(PRIOR_PATH.read_text())
    rows = prior["parameters"]

    def get(candidate, name):
        return next(row for row in rows
                    if row["candidate"] == candidate and row["name"] == name)

    ft_tau = get("free_troposphere_reference", "aod550_2_to_12km_reference")["value"]
    z_ft_km = get("free_troposphere_reference", "scale_height_m")["value"]/1000.
    strat = get("stratosphere_background", "aod550")["value"]

    def ft_part(h_bl_km):
        return ft_tau*(math.exp(-h_bl_km/z_ft_km) - math.exp(-12/z_ft_km)) \
            / (math.exp(-2/z_ft_km) - math.exp(-12/z_ft_km))

    boundary_types = sorted({row["candidate"] for row in rows
                             if row["name"] == "aod550_boundary_layer_only"})
    assert len(boundary_types) == 10
    for candidate in boundary_types:
        tau = get(candidate, "tau_0.55um_type_total")["value"]
        top_km = get(candidate, "layer_top_m")["value"]/1000.
        recorded = get(candidate, "aod550_boundary_layer_only")["value"]
        assert recorded == pytest.approx(tau - ft_part(top_km) - strat, abs=1e-6)


def test_discrimination_product_schema_and_signs():
    """Probe product structure and sign discipline (plan Step 5)."""
    product = json.loads((ROOT/"validation/results/extinction_gap_diagnosis"
                          "/discrimination.json").read_text())
    expected = {row["candidate"] for row in
                json.loads(PRIOR_PATH.read_text())["parameters"]}
    expected.discard("free_troposphere_reference")
    assert set(product["candidates"]) == expected
    verdicts = {"compatible", "tension", "excluded_leaning"}
    for item in product["candidates"].values():
        for event, row in item["events"].items():
            brightness = row["brightness"]
            # Median-of-sums float noise is ~1e-16; real deltas are orders larger.
            assert brightness["probe_delta_ln_median"] <= 1e-9
            assert brightness["verdict"] in verdicts
            assert brightness["max_readout_gap_closure_fraction"] >= 0
            shape = row["shape"]
            if event == "june_2011":
                assert shape["verdict"] in verdicts
                for metric in shape["by_angle"]:
                    assert metric["verdict"] in verdicts and metric["new_ratio"] > 0
            elif event == "hst_2019":
                for dataset_metrics in shape.values():
                    for metric in dataset_metrics.values():
                        assert metric["verdict"] in verdicts and metric["new_slope"] > 0
            elif event == "dec_2011":
                for metric in shape.values():
                    assert metric["verdict"] in verdicts and metric["new_ratio"] > 0
            else:
                assert shape == {}
        assert sum(item["brightness_verdict_counts"].values()) == 4
    assert product["ranking"] == [name for name, item in sorted(
        product["candidates"].items(),
        key=lambda kv: -kv[1]["min_gap_closure_fraction"])]
    assert product["not_probed"][0]["candidate"] == "free_troposphere_humidified"


def test_docs_field_paths_resolve_into_products():
    """Existing and mainline citations share one offline field resolver."""
    paths = []
    for name in ["records/OBSERVATIONAL_VALIDATION.md", "records/MULTI_EVENT_OBSERVATIONS.md",
                 "records/ROUND2.md", "ROADMAP.zh.md", "MAINLINE.zh.md"]:
        found = resolve_product_fields((ROOT/"docs"/name).read_text(),
                                       ROOT/"validation/results")
        if name == "MAINLINE.zh.md":
            assert found, "mainline field coverage must not be empty"
        paths += found
    assert len(paths) >= 12


def test_hst_slope_baseline_matches_committed():
    """Shape baseline anchor: recomputed slopes == committed fields exactly.

    The committed slope aggregates band medians over the FULL window grid;
    the anchored baseline (fail loud inside hst_slope_baseline) is what the
    probe's HST shape branch must use (review round 1 Major#1).
    """
    baselines = egd.hst_slope_baseline()
    for variant in egd.VARIANTS:
        for dataset, slopes in baselines[variant].items():
            committed = egd.hst_comparison()["metrics"][variant][dataset]["slope_model_vs_observed"]
            assert slopes == committed


def test_hst_probe_slope_zero_aod_restores_baseline():
    """AOD=0 probe restores the slopes: exact vs the fresh base kernel,
    committed-congruent within the reproduction tolerance vs the anchored
    CSV baseline (the r1 drift was 1.2-6.6%, far above this budget)."""
    molecular = egd.ray_table("molecular")
    zero_extra = (molecular.impact_m, np.zeros_like(molecular.impact_m),
                  AerosolProfile(aod550=0.0))
    baselines = egd.hst_slope_baseline()["fixed_opac"]
    restored = egd.hst_probe_slopes(zero_extra)
    exact = egd.hst_probe_slopes(None)
    for dataset in restored:
        assert restored[dataset] == exact[dataset]
        for band, slope in restored[dataset].items():
            assert slope == pytest.approx(baselines[dataset][band], rel=1e-9)


def test_broken_power_law_spectral_mirrors_aerosol_profile():
    """Mirror anchor: the validation-side spectral factor equals the
    production form on nonnegative exponents (review round 1 Minor#2)."""
    grid = np.array([300., 350., 430., 499.9, 500., 550., 570., 750., 800.])
    for blue, visible in [(1.10, 1.42), (0.0, 0.0), (1.21, 1.58), (2.0, 0.5)]:
        mirror = egd.BrokenPowerLawSpectral(blue, visible)
        profile = AerosolProfile(aod550=1.0, angstrom_blue=blue, angstrom_visible=visible)
        assert np.array_equal(mirror.spectral_factor(grid), profile.spectral_factor(grid))
    # Negative-exponent branch (mineral transported): 550 nm normalization
    # holds and extinction rises toward the red.
    dust = egd.BrokenPowerLawSpectral(-0.10, -0.13).spectral_factor(np.array([350., 550., 750.]))
    assert dust[1] == 1.0
    assert dust[0] < 1.0 < dust[2]


def test_docs_height_stats_index_citations_semantically_pinned():
    """Position-index citations in docs are pinned to the readout semantics
    they describe; any record reordering fails loud instead of silently
    re-pointing the citation (review round 2 Minor#1/#2)."""
    stats = json.loads((ROOT/"validation/results/extinction_gap_diagnosis"
                        "/height_stats.json").read_text())
    pinned = {41: ("dec_2011", "450nm", "fixed_opac"),
              53: ("dec_2011", "750nm", "fixed_opac")}
    for index, (event, readout, variant) in pinned.items():
        record = stats["records"][index]
        assert (record["event"], record["readout"], record["variant"]) \
            == (event, readout, variant)


def test_full_scale_probe_zero_aod_restores_committed_moments():
    """Adapter-level zero-value restoration on every aperture's own grid."""
    zero_profile = AerosolProfile(aod550=0.0)
    molecular = egd.ray_table("molecular")
    zero_extra = (molecular.impact_m, np.zeros_like(molecular.impact_m), zero_profile)
    record = egd.june_records()[0]
    geometry = egd.june_geometry()
    radius = geometry.moon_distance_m*np.tan(np.deg2rad(record["shadow_angle_deg"]))
    _, _, _, base = egd.june_kernel("fixed_opac", radius, record["wavelength_nm"])
    _, _, _, probed = egd.june_kernel("fixed_opac", radius, record["wavelength_nm"],
                                      zero_extra)
    assert base == probed
    dataset = egd.hst_comparison()["grouping"]["umbral_in_eclipse"][0]
    wavelengths, _ = egd.hst_window_data()
    nodes_b, nodes_w, _, moment_b = egd.hst_kernel("fixed_opac", dataset, wavelengths)
    nodes_p, weights_p, _, moment_p = egd.hst_kernel("fixed_opac", dataset, wavelengths,
                                                     zero_extra)
    assert moment_b == moment_p
    assert all(np.array_equal(a, b) for a, b in zip(nodes_b, nodes_p))
    assert all(np.array_equal(a, b) for a, b in zip(nodes_w, weights_p))
    wavelength = egd.dec_comparison()["records"][0]["wavelength_nm"]
    _, _, _, base = egd.dec_kernel("fixed_opac", wavelength)
    _, _, _, probed = egd.dec_kernel("fixed_opac", wavelength, zero_extra)
    assert base == probed
    _, _, _, base = egd.mv1990_kernel("fixed_opac")
    _, _, _, probed = egd.mv1990_kernel("fixed_opac", zero_extra)
    assert base == probed


def test_strat_rerun_reproduces_committed_and_fulfills_preregistration():
    """Live three-variant rerun gate and preregistered check (plan Step 6).

    The production variant walks the real AerosolColumn chain; the committed
    variants must reproduce their committed predictions within the moment
    budget, and the preregistered probe-vs-production comparison must hold.
    The preregistered reading itself (brightness tension on all four events,
    June shape compatible) is pinned via the committed product's verdicts.
    """
    predictions = egd.build_strat_predictions(egd.build_height_stats())
    for event, worst in predictions["moment_gate"][
            "committed_variant_worst_relative_gap"].items():
        for variant, gap in worst.items():
            assert gap <= egd.MOMENT_TOLERANCE, f"{event}/{variant}: {gap}"
    check = egd.build_preregistered_check(predictions)
    assert check["verdict"] == "within_tolerance"
    for row in check["rows"]:
        assert row["verdict"] in ("within_tolerance", "systematic_deviation")
        assert row["tolerance"] == max(egd.PREREGISTERED_RELATIVE_TOLERANCE
                                       * abs(row["probe_value"]),
                                       egd.PREREGISTERED_FLOOR_ABSOLUTE)
        assert row["abs_diff"] <= row["tolerance"]
        assert row["probe_source_path"].startswith(
            "validation/results/extinction_gap_diagnosis/discrimination.json:")
        assert row["production_source_path"].startswith(
            "validation/results/stratosphere_background/predictions.json:")
    committed = json.loads((ROOT/"validation/results/extinction_gap_diagnosis"
                            "/discrimination.json").read_text())
    candidate = committed["candidates"][egd.PREREGISTERED_CANDIDATE]
    # Probe single source: every brightness probe value is the committed
    # product's field, never a retyped number.
    for row, event in zip(check["rows"][:4], egd.EVENT_ORDER):
        assert row["axis"] == "brightness"
        assert row["probe_value"] == candidate["events"][event]["brightness"][
            "gap_closure_fraction"]
    # Preregistered reading: the brightness residual stays in the tension band
    # (not closed by this component) and the June shape axis is compatible.
    assert all(candidate["events"][event]["brightness"]["verdict"] == "tension"
               for event in egd.EVENT_ORDER)
    assert candidate["events"]["june_2011"]["shape"]["verdict"] == "compatible"
