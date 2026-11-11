"""Schema and anchor gates for the humidity-axis A/B re-inversion product
(task-humidity-axis, plan Step 3 / M3).

Light-weight contract checks over the committed artifact only (no kernel
evaluation): the B(50%) bitwise pool reproduction anchor (independent chain:
file-vs-file, a01), the preregistered verbatim copy, the R_RH band-edge
consistency, and the preregistered share/participation rule applications.
The full inversion is a validation script; this file must stay removable
without touching the data products.
"""
import json
from pathlib import Path

import pytest

from eclipse.aerosol import FT_AOD550_PUBLISHED, ft_growth_factor
from validation.extinction_gap_diagnosis import VERDICT_COMPATIBLE, VERDICT_TENSION
from validation import humidity_axis_inversion as hai

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = ROOT/"validation"/"results"/"humidity_axis"/"inversion_ab.json"
TASK5 = ROOT/"validation"/"results"/"atmosphere_inversion"/"inversion.json"
PRIORS = ROOT/"validation"/"results"/"prior_source_survey"/"humidity_axis_priors.json"
EVENT_ORDER = ("june_2011", "hst_2019", "dec_2011", "event_1990_02_09")


@pytest.fixture(scope="module")
def product():
    return json.loads(PRODUCT.read_text())


@pytest.fixture(scope="module")
def task5():
    return json.loads(TASK5.read_text())


def test_committed_product_preregistered_verbatim(product):
    """The product carries the preregistered block verbatim (task-5 pattern)."""
    assert product["preregistered"] == hai.PREREGISTERED


def test_b50_pool_bitwise_reproduction_anchor(product, task5):
    """B(50%) must equal the committed task-5 pool event by event, BITWISE.

    Independent chain (r2 Minor 2, a01): this reads both committed files and
    compares values - it shares no code path with the inversion machine, so a
    machine regression cannot silently self-certify here."""
    pool = task5["retrieval_spec"]["aod550_pool"]
    ladder_b50 = product["variant_b"]["ladder"]["0.5"]
    for event in EVENT_ORDER:
        assert ladder_b50[event]["aod550"] == pool[event], event
        assert ladder_b50[event]["status"] == "inverted", event
    assert product["b50_reproduction_gate"]["verdict"] == "bitwise_equal"
    # the ratio over the reproduced pool must be the committed raw R
    assert product["r_rh"]["0.5"]["R_RH"] == task5["cross_event"]["R"]


def test_r_rh_band_edges_match_mechanical_constants(product):
    """The descriptive band edges in the priors product must equal the
    imported mechanical edges (edge-value comparison, r3 Minor 2), and every
    R_RH verdict must be the mechanical band application."""
    priors = json.loads(PRIORS.read_text())
    band = priors["r_rh_band"]
    assert band["compatible_upper"] == VERDICT_COMPATIBLE[1]
    assert band["partially_consistent_upper"] == VERDICT_TENSION[1]
    for label, account in product["r_rh"].items():
        assert "n_participating" in account, label
        ratio, verdict = account["R_RH"], account["verdict"]
        if verdict == "degenerate":
            assert account["n_participating"] < 2 and ratio is None, label
        else:
            if ratio <= VERDICT_COMPATIBLE[1]:
                expected = "compatible"
            elif ratio <= VERDICT_TENSION[1]:
                expected = "partially_consistent"
            else:
                expected = "tension"
            assert verdict == expected, label


def test_share_table_follows_preregistered_definition(product, task5):
    """share = 1 - a_res/a_pool for inverted knots; capped at 1.0 for a zero
    residual; None (listed) for over_closed/unclosed."""
    pool = task5["retrieval_spec"]["aod550_pool"]
    for label, knot in product["variant_b"]["ladder"].items():
        shares = product["share_table"][label]
        for event in EVENT_ORDER:
            record, result = shares[event], knot[event]
            if result["status"] == "inverted":
                assert record["share"] == pytest.approx(
                    1.0 - result["aod550"]/pool[event], abs=1e-15), (label, event)
            elif result["status"] == "zero_increment":
                assert record["share"] == 1.0, (label, event)
            else:
                assert record["share"] is None, (label, event)
                assert record["status"] == result["status"], (label, event)


def test_humidity_sufficiency_recomputes(product):
    """The preregistered proxy values recompute independently; the model
    verdict follows variant A's closure status; disagreements carry the
    measured effectiveness ratio (r1 Major 1/4)."""
    f_top = ft_growth_factor(0.99)
    suff = product["humidity_sufficiency"]
    proxy = suff["preregistered_proxy"]
    assert "falsified" in proxy["grade"]
    model = suff["model_verdict"]
    for event, record in proxy["per_event"].items():
        assert record["f_ft_at_domain_top"] == f_top
        f_req = (record["a_pool"] + FT_AOD550_PUBLISHED)/FT_AOD550_PUBLISHED
        assert record["f_required"] == pytest.approx(f_req, rel=1e-15)
        expected = "humidity_sufficient" if f_req <= f_top else "humidity_insufficient"
        assert record["verdict"] == expected, event
    for event, record in model["per_event"].items():
        closed = product["variant_a"][event]["status"] in ("rh_inverted",
                                                           "closed_at_reference")
        assert record["verdict"] == ("humidity_sufficient" if closed
                                     else "humidity_insufficient_for_closure"), event
        assert record["basis"]["rh_star"] == product["variant_a"][event]["rh_star"]
    for event, record in suff["proxy_vs_model"].items():
        assert record["agree"] == (record["proxy"] == record["model"]), event
        if not record["agree"]:
            ratio = record["effectiveness_ratio_at_0_7_knot"]
            assert ratio is not None and abs(ratio-1.0) > 0.1, event
    # the effectiveness account itself: absorbed per unit extra FT AOD550
    extra_ft = FT_AOD550_PUBLISHED*(ft_growth_factor(0.7)-1.0)
    for event, record in suff["effectiveness"]["per_event"].items():
        assert record["extra_ft_aod550"] == pytest.approx(extra_ft, rel=1e-15)
        if record["ratio"] is not None:
            residual = product["variant_b"]["ladder"]["0.7"][event]["aod550"]
            pool = record["absorbed_increment_aod550"] + residual
            task5 = json.loads((ROOT/"validation"/"results"/"atmosphere_inversion"/
                                "inversion.json").read_text())
            assert pool == task5["retrieval_spec"]["aod550_pool"][event]


def test_variant_a_account_contrast_semantics(product):
    """The A accounts measure what humidity changed: before == the production
    reference baseline (equal to B(50%)'s reference accounts bitwise), after
    == the rh* closure; the contrast is non-vacuous (r1 Major 2/3)."""
    b50_accounts = product["variant_b"]["accounts"]["0.5"]
    for event, account in product["variant_a"].items():
        if account["status"] != "rh_inverted":
            continue
        fit = account["fit"]
        assert fit["ln_residual_median_before"] \
            == account["g_at_reference"], event
        assert fit["ln_residual_median_after"] == account["g_at_root"], event
        assert fit["ln_residual_median_before"] > 0, event
        assert abs(fit["ln_residual_median_after"]) < 1e-6, event
        holdout = account["holdout"]
        if holdout is not None:
            # 1990 has no event-internal holdout (its holdout account is the
            # cross-event matrix row, preregistered None)
            assert holdout["ln_rms_before"] \
                == b50_accounts[event]["holdout"]["ln_rms_before"], event
    # the reference-state shape baselines equal B(50%)'s bitwise, and the
    # humidification actually moved the shape axis (factors != 1)
    june = product["variant_a"]["june_2011"]["shape"]
    b50_june = b50_accounts["june_2011"]["shape"]
    assert june["median_base"] == b50_june["median_base"]
    assert june["median_new"] != june["median_base"]
    assert all(angle["shape_factor"] != 1.0 for angle in june["by_angle"])
    assert june["disagreement_closed_fraction"] != 0.0
    for ds, bands in product["variant_a"]["hst_2019"]["shape"][
            "metrics_by_dataset"].items():
        for band, record in bands.items():
            reference = b50_accounts["hst_2019"]["shape"]["metrics_by_dataset"][
                ds][band]
            assert record["base_slope"] == reference["base_slope"], (ds, band)
            assert record["new_slope"] != record["base_slope"], (ds, band)


def test_variant_a_statuses_and_domain(product):
    """Variant A statuses stay inside the preregistered vocabulary and any
    closed RH* lies in the pinned domain."""
    vocabulary = set(hai.PREREGISTERED["variant_a"]["status_vocabulary"])
    for event, account in product["variant_a"].items():
        assert account["status"] in vocabulary, event
        if account["status"] == "rh_inverted":
            assert 0.5 <= account["rh_star"] <= 0.99, event


def test_variant_b_statuses_and_gate_before_results(product):
    """Variant B statuses stay inside the preregistered vocabulary; the
    bitwise gate precedes the results (its verdict is recorded)."""
    vocabulary = set(hai.PREREGISTERED["variant_b"]["status_vocabulary"])
    for label, knot in product["variant_b"]["ladder"].items():
        for event, record in knot.items():
            assert record["status"] in vocabulary, (label, event)
    assert product["b50_reproduction_gate"]["verdict"] == "bitwise_equal"


def test_humidified_path_identity_records(product):
    """Every recorded humidified table walked a verified path: production grid
    bitwise, stratosphere tau bitwise, FT tau linear within the budget."""
    for label, record in product["anchoring"]["humidified_path_identity"].items():
        assert record["grid_bitwise"] is True
        assert record["stratosphere_tau_bitwise"] is True
        assert record["ft_tau_max_rel_dev_from_linear"] \
            <= record["ft_tau_linear_rtol_budget"], label
        rh = float(label)
        assert record["bitwise_production_table"] == (rh == 0.5), label


def test_log_rh_bracket_statuses_match_linear(product):
    """Knot values are identical under both interpolation semantics, so the
    closure statuses cannot differ; closed events carry a two-point bracket."""
    for event, record in product["sensitivity_appendix"][
            "log_rh_interpolation_bracket"]["per_event"].items():
        linear = product["variant_a"][event]
        assert record["status"] == linear["status"], event
        if record["status"] == "rh_inverted":
            lo, hi = record["bracket"]
            assert min(record["rh_star_log"], record["rh_star_linear"]) == lo
            assert max(record["rh_star_log"], record["rh_star_linear"]) == hi
            assert 0.5 <= lo <= hi <= 0.99, event
