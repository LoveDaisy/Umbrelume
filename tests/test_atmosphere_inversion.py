"""Offline tests for the cross-event atmosphere inversion (plan Step 1 / M1).

Gate scope (plan J7): mathematical units and mechanism tests only - the
four-event ledger is a command-side product, gated through its committed
artifact by the product-schema tests (M2). The anchoring gates here are the
inversion baseline's independent evidence: the production two-layer P(0)
reproduces the committed predictions bitwise and the committed-variant chain
identity holds within the moment tolerance.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table

from validation import atmosphere_inversion as ai
from validation import extinction_gap_diagnosis as egd

ROOT = Path(__file__).resolve().parents[1]


# --- Preregistered morphology columns -------------------------------------------


def test_unit_column_spectral_normalization_and_forms():
    """Unit columns exist for every preregistered form; 550 nm factor is 1."""
    for form in (ai.MAIN_FORM,)+ai.SENSITIVITY_FORMS:
        b, tau_unit, spectral = ai.unit_column(form)
        assert tau_unit.shape == b.shape
        assert spectral.spectral_factor(550.0) == 1.0
        # The unit column's vertical AOD is 1.0 by construction (AerosolProfile
        # invariant); the tau table must be nonnegative everywhere.
        assert np.all(tau_unit >= 0)


def test_unit_column_rows_match_pinned_prior_rows():
    """The preregistered morphology literals equal the committed pinned rows."""
    for form, rows in ai.MORPHOLOGY_ROWS.items():
        for name, value in rows.items():
            assert ai.prior_row(form, name)["value"] == value, f"{form}/{name}"
    # The anchor spectrum bounds are the min/max of the pinned rows.
    values = sorted(row["value"] for row in ai.prior_rows()
                    if row["name"] == ai.BL_AOD_ROW)
    assert values[0] == ai.PHYSICAL_ANCHOR_MIN and values[-1] == ai.PHYSICAL_ANCHOR_MAX
    assert len(values) == 10


def test_linear_scaling_identity_two_times_unit_tau():
    """2x tau table == 2 x unit tau (tight allclose), and the kernel-level
    scaled-extra evaluation matches the doubled column (plan performance
    design point; the linear-scaling contract of the unit column)."""
    form = ai.MAIN_FORM
    b_unit, tau_unit, spectral = ai.unit_column(form)
    doubled = build_ray_table(Atmosphere(),
                              aerosol_profile=AerosolProfile(
                                  aod550=2.0, bottom_m=ai.MORPHOLOGY_ROWS[form]["layer_bottom_m"],
                                  top_m=ai.MORPHOLOGY_ROWS[form]["layer_top_m"],
                                  scale_height_m=ai.MORPHOLOGY_ROWS[form]["scale_height_m"]))
    np.testing.assert_allclose(doubled.aerosol_layer_taus[0], 2.0*tau_unit,
                               rtol=1e-12, atol=0)
    # Kernel-level: weights with the doubled-table column == weights with the
    # scaled unit column, tight allclose.
    kernel = ai.event_kernels()["june_2011"]
    readout = kernel.scalars[0]
    member = readout.members[0]
    table = egd.ray_table(ai.PRODUCTION_VARIANT)
    base = egd.kernel_weights(table, member.bb, member.area, readout.wavelength_nm)
    scaled = egd.kernel_weights(table, member.bb, member.area, readout.wavelength_nm,
                                (b_unit, 2.0*tau_unit, spectral))
    doubled_extra = egd.kernel_weights(table, member.bb, member.area,
                                       readout.wavelength_nm,
                                       (doubled.impact_m, doubled.aerosol_layer_taus[0],
                                        spectral))
    np.testing.assert_allclose(scaled, doubled_extra, rtol=1e-12, atol=0)
    np.testing.assert_allclose(scaled, base, rtol=2.0, atol=0)  # sanity: positive


def test_zero_aod_extra_restores_kernel_bitwise():
    """A zero-AOD extra column changes no kernel value bitwise (counterfactual
    interface); the full P(0) bitwise reproduction is the anchoring gate."""
    table = egd.ray_table(ai.PRODUCTION_VARIANT)
    kernels = ai.event_kernels()
    # June scalar readout
    readout = kernels["june_2011"].scalars[0]
    member = readout.members[0]
    base = egd.kernel_weights(table, member.bb, member.area, readout.wavelength_nm)
    zero = egd.kernel_weights(table, member.bb, member.area, readout.wavelength_nm,
                              ai.extra_at(ai.MAIN_FORM, 0.0))
    assert np.array_equal(base, zero)
    # HST dataset row sums on the full grid
    row_kernel = kernels["hst_2019"].rows[0]
    _, grid_rows = egd.hst_grid_data()
    wavelengths = np.array([float(row["wavelength_nm"]) for row in grid_rows])
    member = row_kernel.members[0]
    base = egd.kernel_weights(table, member.bb, member.area, wavelengths)
    zero = egd.kernel_weights(table, member.bb, member.area, wavelengths,
                              ai.extra_at(ai.MAIN_FORM, 0.0))
    assert np.array_equal(base, zero)
    # December scalar readout
    readout = kernels["dec_2011"].scalars[0]
    member = readout.members[0]
    base = egd.kernel_weights(table, member.bb, member.area, readout.wavelength_nm)
    zero = egd.kernel_weights(table, member.bb, member.area, readout.wavelength_nm,
                              ai.extra_at(ai.MAIN_FORM, 0.0))
    assert np.array_equal(base, zero)


def test_anchoring_gates_green_and_bitwise_reproduction():
    """The inversion baseline's independent evidence (plan Step 1 anchoring).

    The reproduction gate is asserted inside run_anchoring_gates (fail loud);
    here the recorded worst relative diffs must additionally sit at moment-
    tolerance level (the structural-reuse claim), not merely inside the looser
    preregistered reproduction tolerance.
    """
    gates = ai.run_anchoring_gates()
    assert gates["reproduction"]["verdict"] == "within_tolerance"
    for event, worst in gates["reproduction"]["worst_relative_diff"].items():
        assert worst <= egd.MOMENT_TOLERANCE, f"{event}: {worst}"
    assert gates["chain_identity"]["verdict"] == "within_tolerance"
    for event, worst in gates["chain_identity"]["kernel_leg_worst_relative_dev"].items():
        assert worst <= egd.MOMENT_TOLERANCE, f"{event}: {worst}"
    expected = {"june_2011": (ai.JUNE_FIT_COUNT, ai.JUNE_HOLDOUT_COUNT),
                "dec_2011": (ai.DEC_FIT_COUNT, ai.DEC_HOLDOUT_COUNT),
                "hst_2019": (ai.HST_FIT_ROW_COUNT*ai.HST_UMBRAL_DATASET_COUNT,
                             ai.HST_HOLDOUT_ROW_COUNT*ai.HST_UMBRAL_DATASET_COUNT),
                "event_1990_02_09": (ai.MV1990_FIT_COUNT, 0)}
    for event, (fit, holdout) in expected.items():
        counts = gates["anchored_readout_counts"][event]
        assert (counts["fit"], counts["holdout"]) == (fit, holdout)


def test_g_monotone_on_june_fit_subset():
    """g is nonincreasing in a on the real June fit subset (preregistered
    monotonicity, exercised at the readout level)."""
    fit_kernel = ai.fit_restricted_kernels()["june_2011"]
    values = [ai._fit_g(fit_kernel, a, ai.MAIN_FORM)
              for a in (0.0, 0.01, 0.03, 0.06, 0.12)]
    assert all(value is not None for value in values)
    for lower, higher in zip(values, values[1:]):
        assert higher <= lower


def test_root_finder_deterministic_on_june():
    """Two inversions of the same event agree bitwise (deterministic contract)."""
    kernels = ai.event_kernels()
    fit_kernels = ai.fit_restricted_kernels()
    first = ai.invert_event(kernels["june_2011"], fit_kernels["june_2011"], ai.MAIN_FORM)
    second = ai.invert_event(kernels["june_2011"], fit_kernels["june_2011"], ai.MAIN_FORM)
    assert first["aod550"] == second["aod550"]
    assert first["g0"] == second["g0"]
    assert first["status"] == second["status"] == ai.STATUS_INVERTED
    # The root closes the fit median: |g(root)| is root-finder noise level.
    residual = ai._fit_g(fit_kernels["june_2011"], first["aod550"], ai.MAIN_FORM)
    assert abs(residual) <= 1e-8


def test_dark_bracket_top_halving_branch():
    """The preregistered branch_dark_top: a dark bracket top routes through
    geometric halving to a finite negative bracket top; the root is unchanged
    by monotonicity (June's top is dark at the preregistered 2.0)."""
    fit_kernel = ai.fit_restricted_kernels()["june_2011"]
    assert ai._fit_g(fit_kernel, ai.FIT_BRACKET_HIGH, ai.MAIN_FORM) is None
    kernels = ai.event_kernels()
    result = ai.invert_event(kernels["june_2011"], fit_kernel, ai.MAIN_FORM)
    assert result["bracket_high_used"] < ai.FIT_BRACKET_HIGH
    assert result["bracket_high_used"] > 0
    # Every a below the used top is finite (monotone underflow cliff).
    assert ai._fit_g(fit_kernel, result["bracket_high_used"]/2, ai.MAIN_FORM) is not None


def test_physical_anchor_brackets_and_out_of_range():
    """Anchor positions against the 10 pinned boundary-layer rows."""
    anchor = ai.physical_anchor(0.05)
    assert anchor["within_pinned_range"] is True
    assert anchor["bracketed_by"] == ["continental_clean", "antarctic"]
    assert anchor["ratio_to_spectrum_min"] == pytest.approx(0.05/0.038)
    above = ai.physical_anchor(0.7)
    assert above["within_pinned_range"] is False
    assert above["bracketed_by"] == ["urban", None]
    below = ai.physical_anchor(0.02)
    assert below["within_pinned_range"] is False
    assert below["bracketed_by"] == [None, "maritime_tropical"]


def test_consistency_verdict_bands_synthetic():
    """The three preregistered bands plus the degenerate branch (synthetic)."""
    def events(aod550s):
        out = {}
        for event, value in zip(ai.EVENT_ORDER, aod550s):
            if value == "unclosed":
                out[event] = {"status": ai.STATUS_UNCLOSED, "aod550": None}
            elif value == "zero":
                out[event] = {"status": ai.STATUS_ZERO_INCREMENT, "aod550": 0.0}
            else:
                out[event] = {"status": ai.STATUS_INVERTED, "aod550": value}
        return out
    bands = ai.consistency_verdict(events([0.05, 0.08, 0.06, 0.07]))
    assert bands["verdict"] == "systematic_gap" and bands["R"] == pytest.approx(0.08/0.05)
    tension = ai.consistency_verdict(events([0.05, 0.20, 0.06, 0.35]))
    assert tension["verdict"] == "partially_consistent"
    specific = ai.consistency_verdict(events([0.05, 0.60, 0.06, 0.07]))
    assert specific["verdict"] == "event_specific"
    degenerate = ai.consistency_verdict(events([0.05, "unclosed", 0.06, 0.07]))
    assert degenerate["verdict"] == "not_applicable"
    assert degenerate["degenerate_events"] == {"hst_2019": ai.STATUS_UNCLOSED}
    zero = ai.consistency_verdict(events([0.05, "zero", 0.06, 0.07]))
    assert zero["verdict"] == "not_applicable"
    assert zero["degenerate_events"] == {"hst_2019": ai.STATUS_ZERO_INCREMENT}


def test_retrieval_spec_rules_synthetic():
    """Rule -> field: range pool, zero-increment min, unclosed exclusion and
    the out-of-range branch rule (review round 2 item 3)."""
    findings = []

    def events(aod550s):
        out = {}
        for event, value in zip(ai.EVENT_ORDER, aod550s):
            if value == "unclosed":
                out[event] = {"status": ai.STATUS_UNCLOSED, "aod550": None}
            elif value == "zero":
                out[event] = {"status": ai.STATUS_ZERO_INCREMENT, "aod550": 0.0}
            else:
                out[event] = {"status": ai.STATUS_INVERTED, "aod550": value}
        return out

    spec = ai.retrieval_spec(events([0.05, 0.08, 0.06, 0.07]), findings)
    assert spec["aod550_range"] == [0.005, 0.8]
    assert spec["branch_rule_fired"] is False
    assert findings == []

    zero_min = ai.retrieval_spec(events([0.05, "zero", 0.06, 0.07]), findings)
    assert zero_min["aod550_range"][0] == 0.0
    assert zero_min["aod550_pool"]["hst_2019"] == 0.0

    excluded = ai.retrieval_spec(events([0.05, "unclosed", 0.06, 0.07]), findings)
    assert excluded["aod550_range"] == [0.05/10.0, 0.07*10.0]
    assert excluded["excluded_events"] == {"hst_2019": ai.STATUS_UNCLOSED}
    assert excluded["aod550_pool"]["hst_2019"] is None

    branch = ai.retrieval_spec(events([0.05, 0.80, 0.06, 0.07]), findings)
    assert branch["branch_rule_fired"] is True
    assert branch["source_types"][-1] == "volcanic / stratospheric-perturbation types"

    dead = ai.retrieval_spec(events(["unclosed", "unclosed", "unclosed", "unclosed"]),
                             findings)
    assert dead["aod550_range"] is None
    assert any(item["code"] == "retrieval_range_degenerate" for item in findings)
    # The alpha envelope covers the pinned forms (union with the main form).
    assert spec["alpha_blue_range"] == [0.2, 1.11]
    assert spec["alpha_visible_range"] == [0.17, 1.42]


def test_hst_shape_machinery_matches_committed_baseline():
    """Structure anchor: the ledger's band-median slope machinery evaluated on
    the committed fixed_opac variant reproduces the anchored committed slope
    baselines within the moment tolerance (the aggregation structure is the
    committed compare_hst2019 one; the two-layer baseline has no committed
    counterpart, so this is its only structural witness)."""
    kernel = ai.event_kernels()["hst_2019"]
    predictions = ai.event_predictions(kernel, 0.0, ai.MAIN_FORM,
                                       table=egd.ray_table("fixed_opac"))
    baselines = egd.hst_slope_baseline()["fixed_opac"]
    for row_kernel in kernel.rows:
        dataset = row_kernel.readout_id
        mine = ai._hst_band_slopes(predictions[dataset], row_kernel.observed,
                                   row_kernel.wavelengths)
        for band, slope in mine.items():
            assert slope == pytest.approx(baselines[dataset][band], rel=1e-12)


def test_fit_restricted_hst_kernel_matches_full_grid_bitwise():
    """HST wavelength subsetting: the fit-restricted evaluation reproduces the
    full-grid fit values bitwise (wavelength columns are independent)."""
    kernel = ai.event_kernels()["hst_2019"]
    fit_kernel = ai.fit_restricted_kernels()["hst_2019"]
    full = ai.event_predictions(kernel, 0.05, ai.MAIN_FORM)
    restricted = ai.event_predictions(fit_kernel, 0.05, ai.MAIN_FORM)
    for row_kernel in kernel.rows:
        np.testing.assert_array_equal(
            full[row_kernel.readout_id][row_kernel.fit_mask],
            restricted[row_kernel.readout_id])
