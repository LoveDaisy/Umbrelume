"""Offline tests for the cross-event atmosphere inversion (plan Step 1 / M1).

Gate scope (plan J7): mathematical units and mechanism tests only - the
four-event ledger is a command-side product, gated through its committed
artifact by the product-schema tests (M2). The anchoring gates here are the
inversion baseline's independent evidence: the production two-layer P(0)
reproduces the committed predictions bitwise and the committed-variant chain
identity holds within the moment tolerance.
"""
import csv
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


# --- Committed product schema (M2: Step 2+3 command products) -------------------

PRODUCT_PATH = ROOT/"validation/results/atmosphere_inversion/inversion.json"
LEDGER_PATH = ROOT/"validation/results/atmosphere_inversion/ledger.csv"

EVENT_FIELDS = {"morphology", "observation_source", "readout_counts", "g0",
                "g_at_bracket_high", "bracket_high_used", "status", "aod550",
                "fit", "holdout", "baseline_contrast", "shape", "physical_anchor"}


def test_committed_product_preregistered_verbatim():
    """The product carries the preregistered block verbatim (plan M2)."""
    product = json.loads(PRODUCT_PATH.read_text())
    assert product["preregistered"] == ai.PREREGISTERED


def test_committed_product_schema_and_anchoring():
    """Product structure, event fields and anchoring verdicts."""
    product = json.loads(PRODUCT_PATH.read_text())
    assert set(product) == {"scope", "preregistered", "anchoring", "events",
                            "cross_event", "form_sensitivity", "retrieval_spec",
                            "findings"}
    assert set(product["events"]) == set(ai.EVENT_ORDER)
    for event, record in product["events"].items():
        assert EVENT_FIELDS <= set(record)
        assert record["status"] in ai.PREREGISTERED["objective"]["status_vocabulary"]
        assert record["morphology"] == ai.MAIN_FORM
        if record["status"] == ai.STATUS_INVERTED:
            assert record["aod550"] > 0
            assert record["physical_anchor"]["inverted_aod550"] == record["aod550"]
        assert record["fit"]["readout_count"] == record["readout_counts"]["fit"]
    anchoring = product["anchoring"]
    assert anchoring["reproduction"]["verdict"] == "within_tolerance"
    assert anchoring["chain_identity"]["verdict"] == "within_tolerance"


def test_committed_cross_event_verdict_is_mechanical():
    """The cross-event verdict re-derives mechanically from the product's own
    event values (never hand-written); the matrix diagonal is root-level."""
    product = json.loads(PRODUCT_PATH.read_text())
    events = {event: {"status": record["status"], "aod550": record["aod550"]}
              for event, record in product["events"].items()}
    recomputed = ai.consistency_verdict(events)
    cross = product["cross_event"]
    assert cross["verdict"] == recomputed["verdict"]
    assert cross["R"] == recomputed["R"]
    assert cross["inverted_aod550"] == recomputed["inverted_aod550"]
    assert cross["degenerate_events"] == recomputed["degenerate_events"]
    for event, row in cross["prediction_matrix"].items():
        assert abs(row[event]) <= 1e-9
        for value in row.values():
            if value is not None:
                assert np.isfinite(value)


def test_committed_shape_and_holdout_declarations():
    """1990 has no shape metric and no event-internal holdout; the other
    events declare the holdout-band involvement of their shape metrics."""
    product = json.loads(PRODUCT_PATH.read_text())
    events = product["events"]
    assert events["event_1990_02_09"]["shape"] is None
    assert events["event_1990_02_09"]["holdout"] is None
    for event in ("june_2011", "dec_2011", "hst_2019"):
        shape = events[event]["shape"]
        if shape is not None:
            assert shape["involves_holdout_bands"] is True


def test_committed_findings_are_mechanical_consequences():
    """Findings re-derive from the product's own values: the main-form
    out-of-range anchor obligation and the sensitivity bracket branches."""
    product = json.loads(PRODUCT_PATH.read_text())
    findings = product["findings"]
    known = {"physical_range_exceeded", "zero_increment", "unclosed",
             "consistency_degenerate", "retrieval_range_degenerate",
             "sensitivity_unclosed", "sensitivity_zero_increment"}
    assert all(item["code"] in known for item in findings)
    # Anchor obligation: every out-of-range main-form value carries a finding.
    out_of_range = {event for event, record in product["events"].items()
                    if record["physical_anchor"] is not None
                    and record["physical_anchor"]["within_pinned_range"] is False}
    recorded = {item["event"] for item in findings
                if item["code"] == "physical_range_exceeded"}
    assert recorded == out_of_range
    # Sensitivity: every unclosed sensitivity cell carries a finding.
    unclosed_cells = sum(
        1 for form in product["form_sensitivity"]["forms"].values()
        for cell in form.values() if cell["status"] == ai.STATUS_UNCLOSED)
    assert sum(1 for item in findings
               if item["code"] == "sensitivity_unclosed") == unclosed_cells


def test_committed_ledger_rows_match_events():
    """Ledger rows: one per readout (HST per window row), subsets partition."""
    product = json.loads(PRODUCT_PATH.read_text())
    with LEDGER_PATH.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1046  # 18 + 2*510 + 7 + 1
    per_event = {}
    for row in rows:
        per_event.setdefault(row["event"], []).append(row)
    assert {event: len(items) for event, items in per_event.items()} == {
        "june_2011": 18, "hst_2019": 1020, "dec_2011": 7,
        "event_1990_02_09": 1}
    for event, items in per_event.items():
        fit = sum(1 for item in items if item["subset"] == "fit")
        holdout = sum(1 for item in items if item["subset"] == "holdout")
        assert (fit, holdout) == (
            product["events"][event]["readout_counts"]["fit"],
            product["events"][event]["readout_counts"]["holdout"])
        for item in items:
            assert item["subset"] in ("fit", "holdout")


def test_calibrated_registry_anchors_to_committed_product():
    """Registry values anchor mechanically to the committed inversion product
    (plan Step 4): preset AOD == the product's event value, morphology rows ==
    the pinned main form, provenance sha == the committed file's actual sha;
    any regeneration that breaks the anchor fails loud instead of silently
    drifting. The preset-name-to-event mapping lives solely in the registry's
    'event' field (review round 2 item 4)."""
    import hashlib

    import eclipse.cli as cli_module

    registry = json.loads(cli_module.CALIBRATED_PRESETS_JSON.read_text())
    product = json.loads(PRODUCT_PATH.read_text())
    sha = hashlib.sha256(PRODUCT_PATH.read_bytes()).hexdigest()
    assert set(registry["presets"]) == {"june2011", "dec2011", "hst2019", "feb1990"}
    rows = ai.MORPHOLOGY_ROWS[ai.MAIN_FORM]
    for name, entry in registry["presets"].items():
        provenance = entry["provenance"]
        assert provenance["kind"] == "calibrated"
        assert provenance["calibrated_not_physical_prior"] is True
        assert provenance["inversion_product"] \
            == "validation/results/atmosphere_inversion/inversion.json"
        assert provenance["sha256"] == sha
        event = entry["event"]
        assert event in ai.EVENT_ORDER
        layer = entry["increment_layer"]
        assert layer["aod550"] == product["events"][event]["aod550"]
        assert layer["bottom_m"] == rows["layer_bottom_m"]
        assert layer["top_m"] == rows["layer_top_m"]
        assert layer["scale_height_m"] == rows["scale_height_m"]
        assert layer["angstrom_blue"] == rows["angstrom_0.35_to_0.5um"]
        assert layer["angstrom_visible"] == rows["angstrom_0.5_to_0.8um"]
        assert layer["break_nm"] == 500.0
        assert provenance["fit_holdout_split"] \
            == ai.PREREGISTERED["fit_holdout_split"][event]
