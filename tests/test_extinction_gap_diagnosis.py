"""Offline tests for the extinction-gap diagnosis machine (plan Step 1 / M1).

The synthetic-table zeroth-moment gate and the June 2011 full-scale anchor are
the two non-negotiable pytest residents (plan M1); the full-scale four-aperture
gates stay in the gate too (measured cost ~10 s, budget 60 s).
"""
import numpy as np
import pytest

from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere, EARTH_RADIUS_M
from eclipse.rays import RayTable, build_ray_table
from eclipse.transport import Geometry, transmission_profile

from validation import extinction_gap_diagnosis as egd


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
                    550.0, False, 0.05*np.exp(-altitude/10_000.),
                    egd.DEFAULT_OPAC)


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
    extra = (probe_table.impact_m, probe_table.aerosol_tau550, probe_profile)
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
