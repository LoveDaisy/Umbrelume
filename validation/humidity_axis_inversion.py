"""A/B preregistered humidity-axis re-inversion of the four events
(calibration lane; task-humidity-axis).

Question: how much of the R=7.110 cross-event spread of the task-5 low-layer
increment pool can the pinned OPAC humidity growth absorb? Two preregistered
variants re-invert all four events on top of the production two-layer default
whose free-troposphere layer carries the pinned humidity growth
(``eclipse/aerosol.py`` ``ft_growth_factor``; f(50%) is exactly 1.0, so the
default output is bit-exact unchanged):

- Variant A (RH-only, one free parameter per event): zero increment; the FT
  humidity is the only free parameter, root-found on the pinned domain
  [0.5, 0.99] (monotone closure, brentq) or recorded unclosed.
- Variant B (RH + residual increment): ladder over the pinned RH knots
  {50, 70, 99}% - at each knot the task-5 inversion machine runs verbatim
  (same objective, fit/holdout splits, antarctic morphology, bracket, root
  finder) through the ``table`` override, root-finding the residual increment
  AOD550. The joint (rh, a_res) family is carried by the per-knot scalar
  roots: the task-5 closure target reduces per event to a single scalar, so
  a 2-parameter joint root is not recoverable and the ladder is the
  preregistered realization.

B(50%) must reproduce the committed task-5 pool BITWISE (mechanical gate,
run before any other result). Lane discipline: nothing here enters the
default model; every value is relative to the humidified production baseline
and is calibrated, not a physical prior. The preregistered block below is
written before results and copied verbatim into the product (pytest-pinned);
the R_RH verdict band is generated from the committed transcription product
``humidity_axis_priors.json`` (single source, never hand-copied) with
band-edge value asserts against the mechanical edges imported from
``extinction_gap_diagnosis``.

Per-table identity gates: every humidified ray table is verified at
construction (production grid bitwise, stratosphere tau bitwise, FT tau
linear in the growth factor within 1e-12; full bitwise equality at the
reference state), so every evaluation walks a verified kernel path.
Command products (``inversion_ab.json``, ``ledger_ab.csv``) regenerate
byte-deterministically: no timestamps, deterministic iteration order, and
the R_RH verdicts are the mechanical application of the preregistered bands.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from eclipse.aerosol import (
    FT_AOD550_PUBLISHED, HUMIDITY_AXIS_PRIORS_JSON, Q_REST, RH_MAX, RH_MIN,
    WASO_EXT550_PER_PARTICLE_KM, AerosolColumn, AerosolProfile, ft_growth_factor,
)
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table

from . import atmosphere_inversion as ai
from .extinction_gap_diagnosis import (
    EVENT_ORDER, PRODUCTION_VARIANT, ROOT, VERDICT_COMPATIBLE, VERDICT_TENSION,
    ray_table,
)

# --- Preregistered constants (written before any inversion result) ------------

LADDER_RH = (0.5, 0.7, 0.99)
STATUS_RH_INVERTED = "rh_inverted"
STATUS_CLOSED_AT_REFERENCE = "closed_at_reference"
STATUS_OVER_CLOSED = "over_closed"
# The humidified FT tau must reproduce the linear reparameterization of the
# production FT tau (same quadrature grid, extinction linear in aod550) well
# within this budget; the residual is last-ulp rounding of the product chain.
FT_TAU_LINEAR_RTOL = 1e-12
COMMITTED_INVERSION_JSON = (
    "validation/results/atmosphere_inversion/inversion.json")

RH_LABEL = {rh: f"{rh:g}" for rh in LADDER_RH}


def _load_r_rh_band() -> dict:
    """The R_RH verdict band, generated from the committed transcription
    product (single source; never hand-copied into this module).

    Band-edge consistency assert (plan review r3 Minor 2): the descriptive
    edges are compared EDGE-VALUE against the mechanical edges imported from
    ``extinction_gap_diagnosis`` - a plain equality of two numbers, not an
    interval-set comparison (which would contradict the committed
    R=7.110 -> partially_consistent fact and fail spuriously).
    """
    priors = json.loads((ROOT/HUMIDITY_AXIS_PRIORS_JSON).read_text())
    band = priors["r_rh_band"]
    if band["compatible_upper"] != VERDICT_COMPATIBLE[1]:
        raise RuntimeError(f"R_RH band edge drift: compatible_upper "
                           f"{band['compatible_upper']!r} != VERDICT_COMPATIBLE[1] "
                           f"{VERDICT_COMPATIBLE[1]!r}")
    if band["partially_consistent_upper"] != VERDICT_TENSION[1]:
        raise RuntimeError(f"R_RH band edge drift: partially_consistent_upper "
                           f"{band['partially_consistent_upper']!r} != "
                           f"VERDICT_TENSION[1] {VERDICT_TENSION[1]!r}")
    return band


@lru_cache(maxsize=1)
def committed_pool() -> dict:
    """The committed task-5 retrieval pool (single source for a_pool)."""
    product = json.loads((ROOT/COMMITTED_INVERSION_JSON).read_text())
    pool = product["retrieval_spec"]["aod550_pool"]
    if set(pool) != set(EVENT_ORDER):
        raise RuntimeError(f"committed pool events {sorted(pool)} != EVENT_ORDER")
    return {event: pool[event] for event in EVENT_ORDER}


@lru_cache(maxsize=1)
def committed_raw_R() -> float:
    product = json.loads((ROOT/COMMITTED_INVERSION_JSON).read_text())
    return product["cross_event"]["R"]


PREREGISTERED = {
    "baseline": "every value is relative to the PRODUCTION two-layer default "
                "(fixed_opac_strat) whose free-troposphere layer carries the "
                "pinned humidity growth (eclipse/aerosol.py ft_growth_factor; "
                "f(0.5) exactly 1.0, q_rest gate enforced at import); the "
                "increment (variant B) is the task-5 antarctic-form extra "
                "column; a zero increment restores the humidified baseline "
                "bitwise and at rh=0.5 the humidified table is bitwise the "
                "production table",
    "ladder": {
        "knots_rh": list(LADDER_RH),
        "knot_source": "pinned WASO ext550 knots (humidity_axis_priors.json "
                       "waso_growth_knots; sha256 per knot)",
        "semantics": "variant B evaluates the task-5 inversion machine "
                     "verbatim at each knot (table override; same objective, "
                     "splits, morphology, bracket, root finder); the ladder "
                     "carries the joint (rh, a_res) family through per-knot "
                     "scalar roots - the task-5 closure target reduces per "
                     "event to a single scalar (brentq on one AOD), so a "
                     "2-parameter joint root is not recoverable and this "
                     "realization is preregistered",
        "interpolation": "between knots the FT growth is the pinned "
                         "linear-in-RH interpolation (eclipse/aerosol.py); "
                         "the ladder itself never interpolates",
    },
    "variant_a": {
        "free_parameter": "rh in [0.5, 0.99] (pinned domain), zero increment",
        "objective": "g(rh) = fit-subset median ln(P(rh)/O) with P through the "
                     "humidified table; monotone nonincreasing in rh (more "
                     "humidity -> more extinction -> every prediction "
                     "decreases)",
        "root_finder": "scipy.optimize.brentq on [0.5, 0.99], xtol=1e-12",
        "branch_unclosed": "g(0.99) > 0 -> the pinned growth range cannot close "
                           "the event; status unclosed, not forced; its "
                           "accounts evaluate the before-state at the reference "
                           "and the after-state at the domain top (the closest "
                           "achievable humidity-only state)",
        "branch_closed_at_reference": "g(0.5) <= 0 -> the reference state "
                                      "already closes the event; status "
                                      "closed_at_reference (not expected: all "
                                      "four task-5 pools are positive)",
        "branch_dark": "a dark g evaluation inside the domain fails loud (the "
                       "pinned growth range is far from the underflow cliff)",
        "status_vocabulary": [STATUS_RH_INVERTED, STATUS_CLOSED_AT_REFERENCE,
                              "unclosed"],
    },
    "variant_b": {
        "free_parameters": "(rh, a_res) jointly, realized as the ladder: rh "
                           "takes the pinned knot values, a_res is the "
                           "per-knot scalar root of the task-5 objective",
        "branch_over_closed": "g(0) < 0 at a knot: the humidified baseline is "
                              "already darker than the fit subset, so the "
                              "closure equation has no nonnegative root; the "
                              "residual is 0 (nothing applied) and the status "
                              "is over_closed - the humidity axis alone "
                              "overshoots this event at this knot",
        "branch_zero_increment": "g(0) == 0 at a knot: status zero_increment "
                                 "(task-5 convention)",
        "status_vocabulary": [ai.STATUS_INVERTED, ai.STATUS_ZERO_INCREMENT,
                              ai.STATUS_UNCLOSED, STATUS_OVER_CLOSED],
    },
    "b50_reproduction_gate": "B(50%) must reproduce the committed task-5 pool "
                             "(validation/results/atmosphere_inversion/"
                             "inversion.json retrieval_spec.aod550_pool) "
                             "BITWISE for every event, before any other "
                             "result is produced; f(0.5)==1.0 exactly makes "
                             "the humidified table bitwise the production "
                             "table, so the gate proves the machine "
                             "reproduction, not an assumption",
    "r_rh_band": _load_r_rh_band(),
    "share_table": {
        "definition": "share(rh, event) = 1 - a_res(rh, event)/a_pool(event) "
                      "for inverted knots; a_res == 0 -> share capped at 1.0; "
                      "over_closed and unclosed knots carry no share (listed "
                      "with their status)",
        "reading": "share is a measured absorption share under the "
                   "preregistered semantics, not a causal claim that humidity "
                   "explains the event; the residual event-specificity stays "
                   "a separate account",
    },
    "humidity_sufficiency": {
        "criterion": "f_req = (a_pool + 0.013)/0.013 <= f_FT(0.99) -> "
                     "humidity-sufficient; otherwise humidity-insufficient",
        "pool_source": "committed task-5 retrieval_spec.aod550_pool",
    },
    "holdout_accounts": "per variant per ladder point per event: fit (g0 -> g "
                        "at the root), holdout (before/after median and RMS), "
                        "shape ledger (committed aggregation structure); the "
                        "before-state is the production reference baseline "
                        "(rh=0.5, bitwise the task-5 baseline) and the "
                        "after-state is the variant's terminal state - variant "
                        "A: rh* (or the domain top when unclosed), variant B: "
                        "the knot's root state (over-closed knots evaluate at "
                        "the zero increment, the applied model state)",
    "anchoring": "production-baseline anchoring gates "
                 "(atmosphere_inversion.run_anchoring_gates) plus per-table "
                 "humidified path identity (production grid bitwise, "
                 "stratosphere tau bitwise, FT tau linear within 1e-12, full "
                 "bitwise equality at rh=0.5), enforced inside the table "
                 "factory so every evaluation walks a verified path",
}


# --- Humidified tables -----------------------------------------------------------


@lru_cache(maxsize=None)
def humidified_table(rh: float):
    """Production ray table with the FT layer humidified to ``rh`` (cached).

    The FT layer's effective AOD is ``0.013*ft_growth_factor(rh)`` through the
    factory's single growth-consumption point. Identity gates run at
    construction (fail loud): the production grid bitwise, the stratosphere
    tau bitwise, and the FT tau linear in the growth factor within
    ``FT_TAU_LINEAR_RTOL``; at the reference state the whole table must be
    bitwise the production table.
    """
    table = build_ray_table(
        Atmosphere(), aerosol_profile=AerosolColumn.opac_default(
            ft_aod550=FT_AOD550_PUBLISHED, rh=rh))
    production = ray_table(PRODUCTION_VARIANT)
    if not np.array_equal(table.impact_m, production.impact_m) \
            or not np.array_equal(table.altitude_m, production.altitude_m):
        raise RuntimeError(f"humidified table grid drifted at rh={rh!r}: the "
                           f"layer geometry is rh-independent, so the impact "
                           f"grid must be bitwise the production grid")
    if not np.array_equal(table.aerosol_layer_taus[1],
                          production.aerosol_layer_taus[1]):
        raise RuntimeError(f"humidified table stratosphere tau drifted at "
                           f"rh={rh!r}: the stratospheric layer carries no "
                           f"humidity parameter (preregistered ruling)")
    factor = ft_growth_factor(rh)
    deviation = np.abs(table.aerosol_layer_taus[0]
                       - factor*production.aerosol_layer_taus[0])
    scale = np.maximum(np.abs(production.aerosol_layer_taus[0]), 1e-300)
    worst = float(np.max(deviation/scale))
    if worst > FT_TAU_LINEAR_RTOL:
        raise RuntimeError(f"humidified FT tau deviates from the linear "
                           f"reparameterization by {worst:.3e} (> "
                           f"{FT_TAU_LINEAR_RTOL:.0e}) at rh={rh!r}")
    if rh == RH_MIN:
        bitwise = (np.array_equal(table.aerosol_layer_taus[0],
                                  production.aerosol_layer_taus[0])
                   and np.array_equal(table.air_column_m2, production.air_column_m2)
                   and np.array_equal(table.ozone_column_m2, production.ozone_column_m2)
                   and np.array_equal(table.bending_rad, production.bending_rad))
        if not bitwise:
            raise RuntimeError("humidified table at the reference rh must be "
                               "bitwise the production table")
    return table


def humidified_identity_record(rh: float) -> dict:
    """The identity metrics of one humidified table (for the product)."""
    table = humidified_table(rh)
    production = ray_table(PRODUCTION_VARIANT)
    deviation = np.abs(table.aerosol_layer_taus[0]
                       - ft_growth_factor(rh)*production.aerosol_layer_taus[0])
    scale = np.maximum(np.abs(production.aerosol_layer_taus[0]), 1e-300)
    return {"rh": rh,
            "grid_bitwise": True,
            "stratosphere_tau_bitwise": True,
            "ft_tau_max_rel_dev_from_linear": float(np.max(deviation/scale)),
            "ft_tau_linear_rtol_budget": FT_TAU_LINEAR_RTOL,
            "bitwise_production_table": bool(
                rh == RH_MIN
                and np.array_equal(table.aerosol_layer_taus[0],
                                   production.aerosol_layer_taus[0]))}


# --- Variant A: the FT humidity is the only free parameter ------------------------


def event_predictions_at(kernel, rh: float) -> dict:
    """Full-kernel predictions of the humidified baseline (zero increment)."""
    return ai.event_predictions(kernel, 0.0, ai.MAIN_FORM,
                                table=humidified_table(rh))


def _g_rh(fit_kernel, rh: float) -> float | None:
    """g(rh) of variant A: fit-subset median ln residual of the humidified
    baseline (zero increment) at rh."""
    return ai._fit_g(fit_kernel, 0.0, ai.MAIN_FORM, table=humidified_table(rh))


def invert_event_rh(kernel, fit_kernel) -> dict:
    """Variant A per event: bracket checks on the pinned domain, then brentq."""
    g_reference = _g_rh(fit_kernel, RH_MIN)
    if g_reference is None:
        raise RuntimeError(f"{kernel.event}: g(0.5) dark - impossible: the "
                           f"production baseline predictions are positive")
    g_top = _g_rh(fit_kernel, RH_MAX)
    if g_top is None:
        raise RuntimeError(f"{kernel.event}: g(0.99) dark - the pinned growth "
                           f"range (f(0.99)={ft_growth_factor(RH_MAX):.6g}) is "
                           f"far from the underflow cliff; investigate")
    if g_top > 0.0:
        return {"status": "unclosed", "rh_star": None,
                "g_at_reference": g_reference, "g_at_domain_top": g_top,
                "p0": event_predictions_at(kernel, RH_MIN),
                "preds_root": event_predictions_at(kernel, RH_MAX)}
    if g_reference <= 0.0:
        return {"status": STATUS_CLOSED_AT_REFERENCE, "rh_star": RH_MIN,
                "g_at_reference": g_reference, "g_at_domain_top": g_top,
                "p0": event_predictions_at(kernel, RH_MIN),
                "preds_root": event_predictions_at(kernel, RH_MIN)}
    rh_star = float(brentq(lambda rh: _g_rh(fit_kernel, rh),
                           RH_MIN, RH_MAX, xtol=ai.BRENTQ_XTOL))
    # The three accounts measure what humidity CHANGED: the before-state is the
    # production reference baseline (rh=0.5, bitwise the task-5 baseline), the
    # after-state is the rh* closure. Binding both ends to the closure state
    # would make the fit/holdout/shape contrasts vacuous (code-review r1
    # Major 2/3).
    p0 = event_predictions_at(kernel, RH_MIN)
    preds_root = event_predictions_at(kernel, rh_star)
    g_root = ai.fit_median(kernel, preds_root)
    if abs(g_root) > 1e-9:
        raise RuntimeError(f"{kernel.event}: variant A root g(rh*)={g_root!r} "
                           f"not closed (|g| > 1e-9) - monotonicity or bracket "
                           f"assumption violated")
    return {"status": STATUS_RH_INVERTED, "rh_star": rh_star,
            "g_at_reference": g_reference, "g_at_domain_top": g_top,
            "p0": p0, "preds_root": preds_root, "g_at_root": g_root}


# --- Variant B: ladder over the pinned knots --------------------------------------


def ladder_point(kernel, fit_kernel, rh: float) -> dict:
    """One task-5-style inversion at a pinned RH knot (verbatim machine).

    The preregistered branch structure adds over_closed to the task-5
    statuses: a negative g(0) means the humidified baseline alone overshoots
    the fit subset, so no nonnegative increment root exists and the residual
    is 0 (the applied model state is the humidified baseline itself).
    """
    table = humidified_table(rh)
    g0 = ai._fit_g(fit_kernel, 0.0, ai.MAIN_FORM, table=table)
    if g0 is None:
        raise RuntimeError(f"{kernel.event}: g(0) dark at rh={rh!r} - the "
                           f"humidified baseline is far from the underflow "
                           f"cliff; investigate")
    if g0 < 0.0:
        p0 = event_predictions_at(kernel, rh)
        return {"status": STATUS_OVER_CLOSED, "aod550": 0.0, "g0": g0,
                "p0": p0, "preds_root": p0, "bracket_high_used": None,
                "g_bracket_high": None, "rh": rh}
    if g0 == 0.0:
        p0 = event_predictions_at(kernel, rh)
        return {"status": ai.STATUS_ZERO_INCREMENT, "aod550": 0.0, "g0": g0,
                "p0": p0, "preds_root": p0, "bracket_high_used": 0.0,
                "g_bracket_high": None, "rh": rh}
    result = ai.invert_event(kernel, fit_kernel, ai.MAIN_FORM, table=table)
    result["rh"] = rh
    return result


# --- Cross-knot accounts ------------------------------------------------------------


def r_rh_account(knot_results: dict) -> dict:
    """Mechanical application of the preregistered R_RH band at one knot."""
    band = PREREGISTERED["r_rh_band"]
    participating = {event: result["aod550"]
                     for event, result in knot_results.items()
                     if result["status"] == ai.STATUS_INVERTED}
    listed = {event: result["status"] for event, result in knot_results.items()
              if result["status"] != ai.STATUS_INVERTED}
    account = {"n_participating": len(participating),
               "participating_residuals": participating,
               "listed_separately": listed,
               "structural_domain": band["structural_domain"],
               "participation_rule": band["participation_rule"]}
    if len(participating) < 2:
        account.update({"R_RH": None, "verdict": "degenerate",
                        "evidence_grade_note": band["evidence_grade_note"]})
        return account
    values = np.array([participating[event] for event in EVENT_ORDER
                       if event in participating])
    ratio = float(values.max()/values.min())
    if ratio <= VERDICT_COMPATIBLE[1]:
        verdict = "compatible"
    elif ratio <= VERDICT_TENSION[1]:
        verdict = "partially_consistent"
    else:
        verdict = "tension"
    n = len(participating)
    account.update({"R_RH": ratio, "verdict": verdict,
                    "band_edges": {"compatible_upper": VERDICT_COMPATIBLE[1],
                                   "partially_consistent_upper": VERDICT_TENSION[1]},
                    "evidence_grade_note": f"ratio over {n} participating "
                                           f"events; a {n}-event ratio is not a "
                                           f"4-event R (weaker evidence grade)",
                    "interpretation_note": band["interpretation_note"],
                    "raw_task5_R_reference": committed_raw_R()})
    return account


def share_account(knot_results: dict, pool: dict) -> dict:
    """(B-A) absorption shares at one knot (preregistered definition)."""
    shares = {}
    for event, result in knot_results.items():
        if result["status"] == ai.STATUS_INVERTED:
            share = 1.0 - result["aod550"]/pool[event]
        elif result["status"] == ai.STATUS_ZERO_INCREMENT:
            share = 1.0  # capped: nothing left for the increment to absorb
        else:
            share = None  # over_closed / unclosed: listed with status
        shares[event] = {"status": result["status"], "a_res": result["aod550"],
                         "share": share}
    return shares


def effectiveness_account(ladder: dict, pool: dict) -> dict:
    """Measured per-event absorption effectiveness of the FT humidification
    against the antarctic-form increment, at the 0.7 knot (machine-derived).

    A ratio far from 1 falsifies the 1:1 AOD interchangeability premise of the
    preregistered budget proxy - the root cause of the proxy's closure-status
    disagreement with the machine results (code-review r1 Major 1/4).
    """
    extra_ft = FT_AOD550_PUBLISHED*(ft_growth_factor(0.7)-1.0)
    per_event = {}
    for event in EVENT_ORDER:
        result = ladder["0.7"][event]
        if result["status"] == ai.STATUS_INVERTED:
            absorbed = pool[event]-result["aod550"]
            per_event[event] = {"absorbed_increment_aod550": absorbed,
                                "extra_ft_aod550": extra_ft,
                                "ratio": absorbed/extra_ft}
        else:
            per_event[event] = {"status": result["status"], "ratio": None}
    return {"definition": "absorbed increment AOD550 (a_pool - residual at the "
                          "0.7 knot) per unit extra FT AOD550 (0.013*(f(0.7)-1)); "
                          "measured on this product's own machine results",
            "per_event": per_event}


def humidity_sufficiency_account(pool: dict, variant_a_summary: dict,
                                 effectiveness: dict) -> dict:
    """The preregistered linear-AOD-budget proxy (retained for audit, graded)
    plus the model-level sufficiency verdict from variant A's machine results.

    Code-review r1 Major 1/4: the preregistered criterion assumes the
    antarctic-form increment and the FT humidification are interchangeable
    1:1 in AOD550 - a premise this product's own results falsified (the
    measured effectiveness ratios differ per event), so the proxy values must
    not be published as the sufficiency verdict. The model-level verdict uses
    variant A's closure status and domain-endpoint residuals; the preregistered
    formula and its values stay untouched for audit.
    """
    f_top = ft_growth_factor(RH_MAX)
    proxy_per_event = {}
    model_per_event = {}
    comparison = {}
    for event in EVENT_ORDER:
        a_pool = pool[event]
        if a_pool is None:
            proxy_per_event[event] = {"a_pool": None, "f_required": None,
                                      "f_ft_at_domain_top": f_top,
                                      "verdict": "not_applicable"}
        else:
            f_req = (a_pool + FT_AOD550_PUBLISHED)/FT_AOD550_PUBLISHED
            proxy_per_event[event] = {"a_pool": a_pool, "f_required": f_req,
                                      "f_ft_at_domain_top": f_top,
                                      "verdict": "humidity_sufficient" if f_req <= f_top
                                      else "humidity_insufficient"}
        summary = variant_a_summary[event]
        closed_in_domain = summary["status"] in (STATUS_RH_INVERTED,
                                                 STATUS_CLOSED_AT_REFERENCE)
        model_per_event[event] = {
            "verdict": "humidity_sufficient" if closed_in_domain
            else "humidity_insufficient_for_closure",
            "basis": {"status": summary["status"], "rh_star": summary["rh_star"],
                      "g_at_reference": summary["g_at_reference"],
                      "g_at_domain_top": summary["g_at_domain_top"]}}
        proxy_verdict = proxy_per_event[event]["verdict"]
        model_verdict = model_per_event[event]["verdict"]
        comparison[event] = {"proxy": proxy_verdict, "model": model_verdict,
                             "agree": proxy_verdict == model_verdict,
                             "effectiveness_ratio_at_0_7_knot":
                                 effectiveness["per_event"][event]["ratio"]}
    return {
        "preregistered_proxy": {
            "grade": "linear_aod_budget_proxy: assumes the antarctic-form "
                     "increment and the FT humidification are interchangeable "
                     "1:1 in AOD550; this product's own machine results "
                     "falsified that premise (see effectiveness and "
                     "proxy_vs_model), so these values are retained for audit "
                     "only and are NOT the sufficiency verdict",
            "criterion": PREREGISTERED["humidity_sufficiency"]["criterion"],
            "per_event": proxy_per_event},
        "model_verdict": {
            "basis": "variant A machine results: the humidified baseline alone "
                     "closes the fit subset within the pinned domain (rh*), or "
                     "the domain-top residual stays positive (unclosed)",
            "per_event": model_per_event},
        "proxy_vs_model": comparison,
        "effectiveness": effectiveness,
    }


def b50_reproduction_gate(b50: dict) -> dict:
    """B(50%) bitwise pool reproduction (runs before any other result)."""
    pool = committed_pool()
    per_event = {}
    for event in EVENT_ORDER:
        result = b50[event]
        expected = pool[event]
        equal = result["aod550"] == expected
        per_event[event] = {"ladder_b50_aod550": result["aod550"],
                            "ladder_b50_status": result["status"],
                            "committed_pool": expected,
                            "bitwise_equal": equal}
        if not equal:
            raise RuntimeError(f"B(50%) reproduction gate failed for {event}: "
                               f"ladder {result['aod550']!r} != committed pool "
                               f"{expected!r} (the humidified table at the "
                               f"reference state must be bitwise the "
                               f"production table)")
    return {"rule": PREREGISTERED["b50_reproduction_gate"],
            "per_event": per_event, "verdict": "bitwise_equal"}


# --- Accounts per variant and knot ---------------------------------------------------


def fit_account(kernel, p0: dict, preds_root: dict | None, readout_count: int) -> dict:
    before = ai.fit_median(kernel, p0)
    return {"ln_residual_median_before": before,
            "ln_residual_median_after": (ai.fit_median(kernel, preds_root)
                                         if preds_root is not None else None),
            "readout_count": readout_count}


def event_account(kernel, p0: dict, preds_root: dict | None) -> dict:
    """Fit / holdout / shape accounts of one knot state (committed structure)."""
    counts = ai.run_anchoring_gates()["anchored_readout_counts"][kernel.event]
    return {"fit": fit_account(kernel, p0, preds_root, counts["fit"]),
            "holdout": ai.holdout_account(kernel, p0, preds_root),
            "shape": ai.shape_ledger(kernel, p0, preds_root)}


# --- Sensitivity appendix -------------------------------------------------------------


_LOG_KNOTS = (0.5, 0.7, 0.99)


def _log_rh_growth(rh: float) -> float:
    """Sensitivity-only alternative: the pinned ext550 knots interpolated
    linearly in log(RH), composed with the same q_rest arithmetic.

    This is a validation-side alternative semantics for the appendix bracket;
    the production axis stays the preregistered linear-in-RH interpolation
    (eclipse/aerosol.py). Exact knot queries return the pinned knot value, so
    f at every knot is identical under both semantics and only the
    between-knot curve differs.
    """
    if not RH_MIN <= rh <= RH_MAX:
        raise ValueError(f"rh {rh!r} outside the pinned domain")
    if rh in _LOG_KNOTS:
        ext550 = WASO_EXT550_PER_PARTICLE_KM[rh]
    else:
        values = tuple(WASO_EXT550_PER_PARTICLE_KM[knot] for knot in _LOG_KNOTS)
        log_knots = tuple(math.log(knot) for knot in _LOG_KNOTS)
        ext550 = values[0]
        for x_hi, lo_log, hi_log, lo_value, hi_value in zip(
                _LOG_KNOTS[1:], log_knots, log_knots[1:], values, values[1:]):
            if rh <= x_hi:
                weight = (math.log(rh)-lo_log)/(hi_log-lo_log)
                ext550 = lo_value + weight*(hi_value-lo_value)
                break
    geometry = 8.0*(math.exp(-2.0/8.0) - math.exp(-12.0/8.0))
    e_waso = 0.6*730.0*ext550*geometry
    return (e_waso+Q_REST)/FT_AOD550_PUBLISHED


@lru_cache(maxsize=None)
def _log_rh_table(rh: float):
    """Ray table for the log-RH sensitivity alternative.

    The sensitivity column sets the FT AOD to ``0.013*f_log(rh)`` explicitly
    (the production growth path is not consumed; the log alternative is
    validation-side only). The FT geometry and spectral law are unchanged, so
    the structural identity of the grid and the stratosphere holds as for the
    production humidified tables.
    """
    column = AerosolColumn((
        AerosolProfile.opac_free_troposphere(
            aod550=FT_AOD550_PUBLISHED*_log_rh_growth(rh)),
        AerosolProfile.opac_stratospheric_background()))
    return build_ray_table(Atmosphere(), aerosol_profile=column)


def invert_event_rh_log(kernel, fit_kernel, linear_result: dict) -> dict:
    """Variant A under the log-RH alternative (appendix bracket only)."""
    if linear_result["status"] != STATUS_RH_INVERTED:
        # f at every knot is identical under both semantics, so the closure
        # status at the domain endpoints cannot change; record the statement.
        return {"status": linear_result["status"], "rh_star_log": None,
                "note": "closure status is knot-determined and identical under "
                        "both interpolation semantics"}
    def g_log(rh: float) -> float:
        predictions = ai.event_predictions(
            fit_kernel, 0.0, ai.MAIN_FORM, table=_log_rh_table(rh))
        return ai.fit_median(fit_kernel, predictions)
    rh_star_log = float(brentq(g_log, RH_MIN, RH_MAX, xtol=ai.BRENTQ_XTOL))
    return {"status": STATUS_RH_INVERTED, "rh_star_log": rh_star_log,
            "rh_star_linear": linear_result["rh_star"],
            "bracket": sorted([linear_result["rh_star"], rh_star_log])}


def sensitivity_appendix(variant_a_results: dict, kernels: dict,
                         fit_kernels: dict) -> dict:
    """Preregistered appendix contents (plan section 3)."""
    priors = json.loads((ROOT/HUMIDITY_AXIS_PRIORS_JSON).read_text())
    log_bracket = {event: invert_event_rh_log(
        kernels[event], fit_kernels[event], result)
        for event, result in variant_a_results.items()}
    return {
        "log_rh_interpolation_bracket": {
            "semantics": "the pinned ext550 knots interpolated linearly in "
                         "log(RH) (validation-side alternative; the production "
                         "axis stays linear-in-RH); exact knot queries return "
                         "the pinned knot value, so closure statuses at the "
                         "domain endpoints are identical under both semantics",
            "per_event": log_bracket},
        "j3_dormant_ruling": priors["layer_scope_rulings"]["increment_layers"],
        "j4_stratosphere_ruling":
            priors["layer_scope_rulings"]["stratosphere_background"],
        "scalar_axis_limitation": "the humidity axis is a 550 nm scalar "
                                  "column-AOD factor: real humidification also "
                                  "enlarges particles and flattens the spectral "
                                  "slope (Angstrom decrease), which would change "
                                  "the shape axis beyond this model's fixed "
                                  "per-layer power laws; that wavelength-"
                                  "resolved response is not simulated (it would "
                                  "entangle the pinned spectral rows) and the "
                                  "shape ledgers above carry the resulting cost "
                                  "as measured, not modeled",
        "rh_star_interpretation_grade": "a closed RH* is the arithmetic amount "
                                        "of humidity the closure equation "
                                        "requires under the preregistered "
                                        "semantics; until an independent "
                                        "meteorological check (e.g. ERA5 along "
                                        "the limb path) is available it must "
                                        "not be cited as an observed humidity",
    }


# --- Product assembly -------------------------------------------------------------------

SCOPE_TEXT = (
    "A/B preregistered humidity-axis re-inversion (task-humidity-axis): how "
    "much of the R=7.110 cross-event spread of the task-5 low-layer increment "
    "pool can the pinned OPAC humidity growth absorb? Variant A (RH-only) "
    "root-finds the free-troposphere humidity per event on the pinned domain "
    "[0.5, 0.99]; variant B ladders the task-5 inversion machine verbatim over "
    "the pinned RH knots {50, 70, 99}% for the residual antarctic-form "
    "increment. Calibration lane: every value is relative to the humidified "
    "production baseline and is calibrated, not a physical prior; nothing here "
    "enters the default model. The preregistered block (verbatim module "
    "constant, pytest-pinned; the R_RH band generated from "
    "humidity_axis_priors.json) is written before results; B(50%) reproduces "
    "the committed task-5 pool bitwise before any other result; the R_RH "
    "verdicts are the mechanical application of the preregistered bands with "
    "the preregistered participation rule.")

LEDGER_FIELDS = ["event", "variant", "rh", "readout", "subset", "wavelength_nm",
                 "observed", "prediction_baseline", "prediction_at_root",
                 "ln_residual_baseline", "ln_residual_at_root",
                 "model_to_observed_baseline", "model_to_observed_at_root"]


def _ledger_rows_ab(kernel, p0: dict, preds_root: dict | None, variant: str,
                    rh_label: str) -> list:
    rows = ai._ledger_rows(kernel, p0, preds_root)
    return [[row[0], variant, rh_label, *row[1:]] for row in rows]


def build_product() -> tuple:
    """Build the product and ledger rows (deterministic; no timestamps)."""
    findings = []
    kernels = ai.event_kernels()
    fit_kernels = ai.fit_restricted_kernels()
    pool = committed_pool()
    ledger = []

    # Variant B, knot 0.5 first: the bitwise reproduction gate precedes every
    # other result (preregistered).
    b50 = {event: ladder_point(kernels[event], fit_kernels[event], RH_MIN)
           for event in EVENT_ORDER}
    gate = b50_reproduction_gate(b50)

    ladder = {RH_LABEL[RH_MIN]: b50}
    for rh in LADDER_RH[1:]:
        ladder[RH_LABEL[rh]] = {
            event: ladder_point(kernels[event], fit_kernels[event], rh)
            for event in EVENT_ORDER}

    # Variant B per-knot accounts
    variant_b = {"ladder": {}, "accounts": {}}
    for label, knot in ladder.items():
        per_event = {}
        for event, result in knot.items():
            kernel = kernels[event]
            account = event_account(kernel, result["p0"], result["preds_root"])
            account.update({"status": result["status"], "aod550": result["aod550"],
                            "g0": result["g0"], "rh": result["rh"],
                            "physical_anchor":
                                ai.physical_anchor(result["aod550"])
                                if result["status"] == ai.STATUS_INVERTED else None})
            if result["status"] == STATUS_OVER_CLOSED:
                findings.append({
                    "code": "over_closed", "event": event, "rh": result["rh"],
                    "detail": f"g(0) = {result['g0']:.6g} < 0 at this knot: the "
                              f"humidified baseline alone is darker than the fit "
                              f"subset; no nonnegative increment root exists and "
                              f"the residual is 0 (preregistered over_closed)"})
            per_event[event] = account
            ledger.extend(_ledger_rows_ab(kernel, result["p0"],
                                          result["preds_root"], "B", label))
        variant_b["ladder"][label] = {
            event: {"status": result["status"], "aod550": result["aod550"],
                    "g0": result["g0"]}
            for event, result in knot.items()}
        variant_b["accounts"][label] = per_event

    # Variant A
    a_raw = {}
    variant_a = {}
    for event in EVENT_ORDER:
        a_raw[event] = invert_event_rh(kernels[event], fit_kernels[event])
        result = a_raw[event]
        account = event_account(kernels[event], result["p0"], result["preds_root"])
        account.update({"status": result["status"], "rh_star": result["rh_star"],
                        "g_at_reference": result["g_at_reference"],
                        "g_at_domain_top": result["g_at_domain_top"]})
        if "g_at_root" in result:
            account["g_at_root"] = result["g_at_root"]
        if result["status"] == "unclosed":
            findings.append({
                "code": "variant_a_unclosed", "event": event,
                "detail": f"g(0.99) = {result['g_at_domain_top']:.6g} > 0: the "
                          f"pinned growth range cannot close this event without "
                          f"an increment; recorded unclosed, not forced"})
        variant_a[event] = account
        ledger.extend(_ledger_rows_ab(
            kernels[event], result["p0"], result["preds_root"], "A",
            "unclosed" if result["rh_star"] is None else repr(result["rh_star"])))

    r_rh = {label: r_rh_account(knot) for label, knot in ladder.items()}
    for label, account in r_rh.items():
        if account["verdict"] == "degenerate":
            findings.append({"code": "r_rh_degenerate", "rh_label": label,
                             "detail": f"fewer than 2 participating events "
                                       f"({account['n_participating']}); no band "
                                       f"verdict (preregistered branch)"})
    shares = {label: share_account(knot, pool) for label, knot in ladder.items()}
    effectiveness = effectiveness_account(ladder, pool)
    sufficiency = humidity_sufficiency_account(
        pool,
        {event: {"status": a_raw[event]["status"],
                 "rh_star": a_raw[event]["rh_star"],
                 "g_at_reference": a_raw[event]["g_at_reference"],
                 "g_at_domain_top": a_raw[event]["g_at_domain_top"]}
         for event in EVENT_ORDER},
        effectiveness)
    appendix = sensitivity_appendix(a_raw, kernels, fit_kernels)

    identity_rhs = list(LADDER_RH) + [a_raw[event]["rh_star"]
                                      for event in EVENT_ORDER
                                      if a_raw[event]["status"] == STATUS_RH_INVERTED]
    identity = {repr(rh): humidified_identity_record(rh) for rh in identity_rhs}

    product = {
        "scope": SCOPE_TEXT,
        "preregistered": PREREGISTERED,
        "anchoring": {
            "production_baseline": ai.run_anchoring_gates(),
            "humidified_path_identity": identity,
        },
        "variant_a": variant_a,
        "variant_b": variant_b,
        "b50_reproduction_gate": gate,
        "share_table": shares,
        "r_rh": r_rh,
        "humidity_sufficiency": sufficiency,
        "sensitivity_appendix": appendix,
        "findings": findings,
    }
    return product, ledger


def write_products(out: Path) -> dict:
    """Write ``inversion_ab.json`` + ``ledger_ab.csv`` into ``out``."""
    product, ledger = build_product()
    out.mkdir(parents=True, exist_ok=True)
    (out/"inversion_ab.json").write_text(json.dumps(product, indent=2)+"\n")
    with (out/"ledger_ab.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(LEDGER_FIELDS)
        for row in ledger:
            writer.writerow([json.dumps(value) if isinstance(value, (list, dict))
                             else ("" if value is None else value) for value in row])
    return {
        "b50_gate": product["b50_reproduction_gate"]["verdict"],
        "variant_a": {event: {"status": account["status"],
                              "rh_star": account["rh_star"]}
                      for event, account in product["variant_a"].items()},
        "variant_b": {label: {event: {"status": record["status"],
                                      "aod550": record["aod550"]}
                              for event, record in knot.items()}
                      for label, knot in product["variant_b"]["ladder"].items()},
        "r_rh": {label: {"R_RH": account["R_RH"], "verdict": account["verdict"],
                         "n_participating": account["n_participating"]}
                 for label, account in product["r_rh"].items()},
        "sufficiency_model": {event: record["verdict"] for event, record
                              in product["humidity_sufficiency"]["model_verdict"][
                                  "per_event"].items()},
        "sufficiency_proxy_audit": {event: record["verdict"] for event, record
                                    in product["humidity_sufficiency"][
                                        "preregistered_proxy"]["per_event"].items()},
        "findings": product["findings"],
        "ledger_rows": len(ledger),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("outputs/humidity_axis"))
    args = parser.parse_args()
    summary = write_products(args.out)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
