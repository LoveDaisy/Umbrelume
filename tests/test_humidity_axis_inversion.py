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
    """f_req = (a_pool + 0.013)/0.013 vs f_FT(0.99), recomputed independently."""
    f_top = ft_growth_factor(0.99)
    for event, record in product["humidity_sufficiency"]["per_event"].items():
        assert record["f_ft_at_domain_top"] == f_top
        f_req = (record["a_pool"] + FT_AOD550_PUBLISHED)/FT_AOD550_PUBLISHED
        assert record["f_required"] == pytest.approx(f_req, rel=1e-15)
        expected = "humidity_sufficient" if f_req <= f_top else "humidity_insufficient"
        assert record["verdict"] == expected, event


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
