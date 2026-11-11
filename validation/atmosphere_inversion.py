"""Cross-event atmospheric inversion of the required low-layer aerosol
increment (calibration lane).

For each of the four committed observational apertures (June 2011 curve
readouts, 2011-12-10 seven-point spectra, 2019-01-21 HST/STIS umbral spectra,
1990-02-09 whole-disk photometry) this module root-finds the extra low-layer
AOD550 that closes the event's fit-subset brightness residual on top of the
PRODUCTION two-layer default (``fixed_opac_strat``), with layer morphology
and spectral law pinned to ``prior_parameters.json`` rows and only AOD550
free.

Lane discipline: the inverted values and the ``calibrated:`` presets built
from them are calibration-lane products. They never enter the default model,
never touch a physical parameter, and every consumer must carry the
"calibrated, not a physical prior" marker. The preregistered block below is
the written-before-results contract (splits, objective, criterion, spec
rules) and is copied verbatim into the product; the pytest gate asserts the
copy is literal.

Prediction machinery: P(a) is evaluated through the same kernel adapters and
aggregation orders as the committed predictions (``extinction_gap_diagnosis``
quadrature and ``kernel_weights``, reused by import); the increment enters as
an ``extra_aerosol`` column of the pinned morphology at AOD550 = a, so a zero
increment restores the baseline bitwise and the unit column scales linearly
in its tau table (both pytest-pinned). Per-readout receiver quadrature is
cached once and every evaluation walks the fused ``kernel_weights`` path, so
P(0) equals the adapter prediction bitwise by construction. Root-finding
evaluations restrict HST to the fit-window wavelength subset (wavelength
columns of ``kernel_weights`` are independent, pytest-pinned bitwise); all
product-side accounts evaluate the full grid.

Anchoring gates reproduce the committed production predictions from
``validation/results/stratosphere_background/predictions.json`` per readout,
and the committed-variant ratio chain, before any inversion value is
produced. Command products (``inversion.json``, ``ledger.csv``) regenerate
byte-identically: no timestamps, deterministic iteration order, and the
cross-event verdict is the mechanical application of the preregistered bands.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.transport import Geometry

from .compare_hst2019 import (
    G430L_MAX_NM, SLOPE_BANDS_NM, UMBRAL_QUANTITATIVE_MIN_NM, band_median,
)
from .extinction_gap_diagnosis import (
    EVENT_ORDER, MOMENT_TOLERANCE, MV1990_APERTURE_ORDER, MV1990_EVENT,
    MV1990_WAVELENGTH_STEP_NM, PREREGISTERED_FLOOR_ABSOLUTE,
    PREREGISTERED_RELATIVE_TOLERANCE, PRODUCTION_VARIANT, ROOT,
    VERDICT_COMPATIBLE, VERDICT_TENSION, build_height_stats,
    closure_verdict, dec_comparison, dec_time_weights, hst_comparison,
    hst_grid_data, hst_members, june_geometry, june_records, kernel_weights,
    mv1990_analysis, mv1990_record, prior_row, prior_rows, ray_table,
    receiver_quadrature, shape_closed_fraction, vacuum_branch,
)
from .multiple_events import aperture_quadrature, photopic_response

# --- Preregistered constants (written before any inversion result) ------------

MAIN_FORM = "antarctic"
SENSITIVITY_FORMS = ("desert", "continental_average")
BL_AOD_ROW = "aod550_boundary_layer_only"
FIT_BRACKET_HIGH = 2.0
BRENTQ_XTOL = 1e-12
STATUS_INVERTED = "inverted"
STATUS_ZERO_INCREMENT = "zero_increment"
STATUS_UNCLOSED = "unclosed"

# Fit/holdout splits, pinned per event (review round 2 item 1: HST window-row
# counts are hard numbers, asserted against the committed grid).
JUNE_FIT_WL_NM = 677.0
JUNE_HOLDOUT_WL_NM = 503.0
JUNE_FIT_COUNT = 10
JUNE_HOLDOUT_COUNT = 8
DEC_FIT_WL_NM = (550.0, 600.0, 650.0, 700.0, 750.0)
DEC_HOLDOUT_WL_NM = (450.0, 500.0)
DEC_FIT_COUNT = 5
DEC_HOLDOUT_COUNT = 2
HST_FIT_WINDOW_NM = (500.0, 570.0)
HST_HOLDOUT_WINDOW_NM = (430.0, 500.0)  # upper bound exclusive
HST_FIT_ROW_COUNT = 255
HST_HOLDOUT_ROW_COUNT = 255
HST_UMBRAL_DATASET_COUNT = 2
MV1990_FIT_COUNT = 1

# Physical-plausibility anchor spectrum: the 10 pinned boundary-layer rows.
PHYSICAL_ANCHOR_MIN = 0.038
PHYSICAL_ANCHOR_MAX = 0.625

COMMITTED_STRAT_PREDICTIONS = (
    "validation/results/stratosphere_background/predictions.json")

# Pinned morphology rows (literal preregistration; ``unit_column`` asserts
# equality with the committed prior_parameters.json rows, fail loud on drift).
MORPHOLOGY_ROWS = {
    "antarctic": {"layer_bottom_m": 0.0, "layer_top_m": 10000.0,
                  "scale_height_m": 8000.0,
                  "angstrom_0.35_to_0.5um": 0.34, "angstrom_0.5_to_0.8um": 0.73},
    "desert": {"layer_bottom_m": 0.0, "layer_top_m": 6000.0,
               "scale_height_m": 2000.0,
               "angstrom_0.35_to_0.5um": 0.2, "angstrom_0.5_to_0.8um": 0.17},
    "continental_average": {"layer_bottom_m": 0.0, "layer_top_m": 2000.0,
                            "scale_height_m": 8000.0,
                            "angstrom_0.35_to_0.5um": 1.11,
                            "angstrom_0.5_to_0.8um": 1.42},
}


def _morphology_block() -> dict:
    return {
        "main": {"form": MAIN_FORM, "rows": dict(MORPHOLOGY_ROWS[MAIN_FORM])},
        "sensitivity": [{"form": form, "rows": dict(MORPHOLOGY_ROWS[form])}
                        for form in SENSITIVITY_FORMS],
        "discipline": "layer geometry and spectral law are pinned "
                      "prior_parameters.json rows (nonnegative Angstrom "
                      "exponents, direct AerosolProfile); AOD550 is the only "
                      "free parameter, applied as a unit-AOD column scaled "
                      "linearly; the main conclusion is carried by the main "
                      "form only and the sensitivity forms enter the appendix "
                      "table only",
    }


MORPHOLOGIES = _morphology_block()

PREREGISTERED = {
    "baseline": "every value in this product is relative to the PRODUCTION "
                "two-layer default (fixed_opac_strat: AerosolColumn.opac_default()) "
                "with zero increment; the committed fixed_opac variant appears only "
                "as a monitoring contrast, never as the inversion baseline",
    "fit_holdout_split": {
        "june_2011": {"fit": {"wavelength_nm": JUNE_FIT_WL_NM,
                              "readout_count": JUNE_FIT_COUNT},
                      "holdout": {"wavelength_nm": JUNE_HOLDOUT_WL_NM,
                                  "readout_count": JUNE_HOLDOUT_COUNT},
                      "source": "ROUND2 calibration split reuse (677 fit / 503 withheld)"},
        "dec_2011": {"fit": {"wavelength_nm": list(DEC_FIT_WL_NM),
                             "readout_count": DEC_FIT_COUNT},
                     "holdout": {"wavelength_nm": list(DEC_HOLDOUT_WL_NM),
                                 "readout_count": DEC_HOLDOUT_COUNT}},
        "hst_2019": {"fit": {"window_nm": list(HST_FIT_WINDOW_NM),
                             "window_row_count": HST_FIT_ROW_COUNT},
                     "holdout": {"window_nm": list(HST_HOLDOUT_WINDOW_NM),
                                 "upper_exclusive": True,
                                 "window_row_count": HST_HOLDOUT_ROW_COUNT},
                     "pool": {"datasets": "grouping.umbral_in_eclipse (committed)",
                              "umbral_dataset_count": HST_UMBRAL_DATASET_COUNT,
                              "row_pools_jointly": True}},
        "event_1990_02_09": {"fit": {"readout_count": MV1990_FIT_COUNT,
                                     "definition": "single whole-disk broadband readout"},
                             "holdout": "cross-event prediction matrix rows are its "
                                        "holdout account (no event-internal holdout)"},
    },
    "objective": {
        "definition": "g(a) = median over the fit subset of ln(P(a)/O); P(a) is the "
                      "production two-layer baseline kernel prediction with an extra "
                      "increment column of the pinned morphology at AOD550 = a, built "
                      "through exactly the committed kernel adapters, solar scaling "
                      "and observed-ratio benchmarks",
        "monotonicity": "g is nonincreasing in a (each ln ratio decreases with more "
                        "extinction; the median preserves monotonicity)",
        "root_finder": "scipy.optimize.brentq on [0.0, 2.0], xtol=1e-12",
        "branch_zero_increment": "g(0) <= 0 -> required increment is 0; recorded as "
                                 "status zero_increment with a finding, not an error",
        "branch_unclosed": "g(2.0) > 0 -> the preregistered bracket cannot close the "
                           "event; recorded as status unclosed with a finding; no "
                           "forced closure at any larger AOD (lane discipline)",
        "branch_dark_top": "if any fit-subset P(a) underflows to <= 0 at the bracket "
                           "top (the model transmits no measurable light; the vacuum "
                           "branch carries ~1e-16 float noise), the top is 'dark' and "
                           "the event is closable within the preregistered bracket; "
                           "the finite bracket top is found by geometric halving, and "
                           "the first finite point is guaranteed negative: tau*a there "
                           "is within a factor 2 of the underflow cliff (~745), hence "
                           "P << O; brentq then runs on [0, a_finite] and monotonicity "
                           "leaves the root unchanged (every a <= a_finite has P > 0 "
                           "because the underflow cliff is monotone in a)",
        "status_vocabulary": [STATUS_INVERTED, STATUS_ZERO_INCREMENT, STATUS_UNCLOSED],
    },
    "observation_sources": {
        "june_2011": {"product": "validation/results/round2/analysis.json",
                      "field": "records[*].observed_curve_brightness"},
        "dec_2011": {"product": "validation/results/dec2011/comparison.json",
                     "field": "records[*].observed_curve_ratio"},
        "hst_2019": {"product": "validation/results/hst2019/comparison.csv",
                     "field": "rows[*].{dataset}_obs_ratio"},
        "event_1990_02_09": {"product": "validation/results/multiple_events/analysis.json",
                             "field": "observed = prediction / model_to_observed "
                                      "(variant=fixed_opac, exact same-record recovery)"},
    },
    "consistency_criterion": {
        "statistic": "R = max/min over the four events' main-form inverted AOD550",
        "systematic_gap": f"R <= {VERDICT_COMPATIBLE[1]}",
        "partially_consistent": f"{VERDICT_COMPATIBLE[1]} < R <= {VERDICT_TENSION[1]}",
        "event_specific": f"R > {VERDICT_TENSION[1]}",
        "band_edge_source": "VERDICT_COMPATIBLE[1] / VERDICT_TENSION[1] imported from "
                            "extinction_gap_diagnosis (single verdict vocabulary)",
        "degenerate_branch": "any event without a valid inversion value (zero increment "
                             "or unclosed) -> not_applicable + findings[] record (the "
                             "max/min ratio degenerates)",
        "auxiliary": "geometric mean and population geometric standard deviation "
                     "describe the spread only and never carry the verdict",
    },
    "morphologies": MORPHOLOGIES,
    "prediction_matrix": "event i re-evaluates its fit-subset median ln residual with "
                         "event j's (j != i) inverted AOD550 on the main form; the "
                         "diagonal is the event's own root (approximately zero by "
                         "construction); the 1990 row is its holdout account",
    "shape_ledger": {
        "june_2011": "503/677 model-over-observed transmission ratios at the paired "
                     "angles carrying both bands",
        "dec_2011": "450|500 over 750 transmission ratio model-over-observed",
        "hst_2019": "the two committed slope bands, full-window-grid band-median "
                    "aggregation (same committed aggregation structure)",
        "event_1990_02_09": None,
        "declaration": "the shape axis is a cost ledger, not a fit object; its metrics "
                       "naturally involve holdout bands and this involvement is "
                       "declared, not hidden; closure semantics reuse "
                       "shape_closed_fraction (negative = worsening); the baseline is "
                       "the production default at zero increment",
    },
    "physical_anchor": {
        "reference": "the 10 pinned aod550_boundary_layer_only rows of "
                     "prior_parameters.json",
        "spectrum_min": PHYSICAL_ANCHOR_MIN,
        "spectrum_max": PHYSICAL_ANCHOR_MAX,
        "obligation": "report the inverted value's ratio to the spectrum bounds and "
                      "its bracketing pinned rows; out-of-range values are recorded "
                      "in findings[], never hidden (report obligation, no pass/fail "
                      "threshold)",
    },
    "retrieval_spec_rules": {
        "aod550_range": "[min/10, max*10] over events with valid inversion values; a "
                        "zero-increment event enters min as 0; unclosed events are "
                        "excluded from the range and recorded in findings[]",
        "spectral_coverage": "sensitivity-morphology Angstrom envelope "
                             "(alpha_blue in [0.2, 1.11], alpha_visible in "
                             "[0.17, 1.42]) union the main form's values",
        "layer_geometry": "0-10 km class exponential layer, scale height of the "
                          "pinned row in use (main form: antarctic/scale_height_m)",
        "source_types": "boundary-layer / lower-tropospheric aerosol optical types "
                        "(OPAC-like) with humidity variants (the RH rows' unavailable "
                        "status in the pinned source is the known gap)",
        "branch_rule": "any inversion value above the pinned boundary-layer upper "
                       "bound (0.625) appends volcanic / stratospheric-perturbation "
                       "source types",
    },
    "dimensional_form": {
        "statement": "P(a) is built through exactly the kernel adapters, solar "
                     "scaling and observed-ratio benchmarks that produced the "
                     "committed predictions; per readout residual = ln(P(a)/O) with "
                     "O recovered from committed fields only",
        "mechanical_criterion": "P(0)/O reproduces the committed *_model_to_observed "
                                "per readout within MOMENT_TOLERANCE (the kernel "
                                "re-sums its terms, so exact bitwise equality holds "
                                "only at the committed-value level, which is asserted "
                                "separately); a zero increment restores the baseline "
                                "prediction bitwise",
    },
}


# --- Pinned morphology columns -------------------------------------------------


@lru_cache(maxsize=None)
def unit_column(form: str) -> tuple:
    """Unit-AOD extra column of one pinned morphology.

    Mirrors ``egd.candidate_probe_extra``'s parameter reading: layer geometry
    from the pinned prior rows, AOD normalized to 1.0 so the inversion scales
    the column linearly; the spectral factor comes from ``AerosolProfile``
    itself (all preregistered forms have nonnegative Angstrom rows - a
    negative row fails loud here instead of silently switching representation).
    """
    rows = MORPHOLOGY_ROWS[form]
    committed = {name: prior_row(form, name)["value"] for name in rows}
    if committed != rows:
        raise RuntimeError(f"Preregistered morphology rows drifted from the pinned "
                           f"prior_parameters.json rows for {form}: "
                           f"preregistered {rows} != committed {committed}")
    if rows["angstrom_0.35_to_0.5um"] < 0 or rows["angstrom_0.5_to_0.8um"] < 0:
        raise ValueError(f"{form}: preregistered forms require nonnegative rows")
    spectral = AerosolProfile(aod550=1.0, bottom_m=rows["layer_bottom_m"],
                              top_m=rows["layer_top_m"],
                              scale_height_m=rows["scale_height_m"],
                              angstrom_blue=rows["angstrom_0.35_to_0.5um"],
                              angstrom_visible=rows["angstrom_0.5_to_0.8um"])
    column = AerosolProfile(aod550=1.0, bottom_m=rows["layer_bottom_m"],
                            top_m=rows["layer_top_m"],
                            scale_height_m=rows["scale_height_m"])
    table = build_ray_table(Atmosphere(), aerosol_profile=column)
    return table.impact_m, table.aerosol_layer_taus[0], spectral


def extra_at(form: str, aod550: float) -> tuple:
    """The ``extra_aerosol`` tuple at AOD550 = aod550 (linear column scaling).

    ``aod550 = 0.0`` yields an exactly-zero tau column, so ``kernel_weights``
    adds exactly zero optical depth and the baseline is restored bitwise.
    """
    b, tau_unit, spectral = unit_column(form)
    return b, aod550*tau_unit, spectral


# --- Cached event kernels -------------------------------------------------------


@dataclass(frozen=True)
class Member:
    """One cached receiver-quadrature contribution with its shares."""
    bb: np.ndarray
    area: np.ndarray
    share: float          # multiplies the member's summed weights
    vacuum_share: float   # a-independent outer-branch contribution


@dataclass(frozen=True)
class ScalarReadout:
    """One scalar prediction readout (June / December / 1990)."""
    readout_id: str
    subset: str           # "fit" | "holdout"
    wavelength_nm: float | None
    observed: float
    members: tuple[Member, ...]
    shadow_angle_deg: float | None = None


@dataclass(frozen=True)
class RowKernel:
    """One per-wavelength-row prediction block (HST dataset)."""

    readout_id: str
    members: tuple[Member, ...]
    wavelengths: np.ndarray
    observed: np.ndarray
    scaling: float
    fit_mask: np.ndarray
    holdout_mask: np.ndarray

    def fit_only(self) -> "RowKernel":
        """The same members restricted to the fit-window wavelength subset.

        ``kernel_weights`` treats each wavelength column independently, so the
        restricted evaluation reproduces the full-grid fit columns bitwise
        (pytest-pinned through the fit-median self-check in ``invert_event``).
        """
        return RowKernel(self.readout_id, self.members,
                         self.wavelengths[self.fit_mask],
                         self.observed[self.fit_mask], self.scaling,
                         np.ones(int(self.fit_mask.sum()), dtype=bool),
                         np.zeros(int(self.fit_mask.sum()), dtype=bool))


@dataclass(frozen=True)
class EventKernel:
    event: str
    scalars: tuple[ScalarReadout, ...]
    rows: tuple[RowKernel, ...]
    photopic_grid: tuple | None = None  # (wavelengths, response, norm) for 1990


@lru_cache(maxsize=None)
def june_event_kernel() -> EventKernel:
    table = ray_table(PRODUCTION_VARIANT)
    geometry = june_geometry()
    scalars = []
    for record in june_records():
        radius = geometry.moon_distance_m*np.tan(np.deg2rad(record["shadow_angle_deg"]))
        bb, area = receiver_quadrature(table, float(radius), geometry)
        vacuum = vacuum_branch(table, float(radius), geometry)
        wl = float(record["wavelength_nm"])
        subset = "fit" if wl == JUNE_FIT_WL_NM else "holdout"
        scalars.append(ScalarReadout(
            readout_id=f'{wl:g}nm@{record["shadow_angle_deg"]:g}deg', subset=subset,
            wavelength_nm=wl, observed=float(record["observed_curve_brightness"]),
            members=(Member(bb, area, 1.0, vacuum),),
            shadow_angle_deg=float(record["shadow_angle_deg"])))
    return EventKernel("june_2011", tuple(scalars), ())


@lru_cache(maxsize=None)
def dec_event_kernel() -> EventKernel:
    table = ray_table(PRODUCTION_VARIANT)
    geom_records = dec_comparison()["geometry"]["records"]
    time_weights = dec_time_weights()
    scalars = []
    for record in dec_comparison()["records"]:
        members = []
        for rec, time_weight in zip(geom_records, time_weights):
            geometry = Geometry(moon_distance_m=rec["receiver_axial_distance_m"],
                                sun_distance_m=rec["sun_earth_distance_m"])
            radius = float(rec["shadow_radius_m"])
            bb, area = receiver_quadrature(table, radius, geometry)
            outer = vacuum_branch(table, radius, geometry)
            share = float(time_weight*rec["solar_flux_relative_to_b2"])
            members.append(Member(bb, area, share, outer*share))
        wl = float(record["wavelength_nm"])
        subset = "fit" if wl in DEC_FIT_WL_NM else "holdout"
        scalars.append(ScalarReadout(f"{wl:g}nm", subset, wl,
                                     float(record["observed_curve_ratio"]),
                                     tuple(members)))
    return EventKernel("dec_2011", tuple(scalars), ())


@lru_cache(maxsize=None)
def hst_event_kernel() -> EventKernel:
    table = ray_table(PRODUCTION_VARIANT)
    comparison = hst_comparison()
    wavelengths, rows = hst_grid_data()
    solar_scale = comparison["solar_scale_in_over_out"]
    datasets = comparison["grouping"]["umbral_in_eclipse"]
    if len(datasets) != HST_UMBRAL_DATASET_COUNT:
        raise RuntimeError(f"HST umbral dataset count changed: {len(datasets)} != "
                           f"{HST_UMBRAL_DATASET_COUNT} (preregistered)")
    row_kernels = []
    for dataset in datasets:
        members_list = []
        members = hst_members(dataset)
        total = sum(member["texp_s"] for member in members)
        for member in members:
            radius = float(member["shadow_axis_distance_m"])
            geometry = Geometry(moon_distance_m=member["receiver_axial_distance_m"],
                                sun_distance_m=member["sun_site_distance_m"])
            bb, area = receiver_quadrature(table, radius, geometry)
            outer = vacuum_branch(table, radius, geometry)
            share = member["texp_s"]/total
            members_list.append(Member(bb, area, share, outer*share))
        observed = np.array([float(row[f"{dataset}_obs_ratio"]) for row in rows])
        window = ((wavelengths >= UMBRAL_QUANTITATIVE_MIN_NM)
                  & (wavelengths <= G430L_MAX_NM))
        fit = (window & (wavelengths >= HST_FIT_WINDOW_NM[0])
               & (wavelengths <= HST_FIT_WINDOW_NM[1]))
        holdout = window & (wavelengths < HST_HOLDOUT_WINDOW_NM[1])
        if (int(fit.sum()) != HST_FIT_ROW_COUNT
                or int(holdout.sum()) != HST_HOLDOUT_ROW_COUNT):
            raise RuntimeError(f"HST window split counts drifted: fit {int(fit.sum())} "
                               f"!= {HST_FIT_ROW_COUNT} or holdout "
                               f"{int(holdout.sum())} != {HST_HOLDOUT_ROW_COUNT}")
        if int((fit & holdout).sum()) or not np.array_equal(fit | holdout, window):
            raise RuntimeError("HST fit/holdout masks do not partition the "
                               "quantitative window")
        row_kernels.append(RowKernel(dataset, tuple(members_list),
                                     np.asarray(wavelengths, dtype=float), observed,
                                     float(solar_scale[dataset]), fit, holdout))
    return EventKernel("hst_2019", (), tuple(row_kernels))


@lru_cache(maxsize=None)
def mv1990_event_kernel() -> EventKernel:
    analysis = mv1990_analysis()
    record = mv1990_record("fixed_opac")
    event = next(item for item in analysis["geometry"]["events"]
                 if item["event_date"] == MV1990_EVENT)
    geometry = Geometry(moon_distance_m=event["moon_distance_m"],
                        sun_distance_m=event["sun_distance_m"])
    wavelengths = np.arange(380., 780.+MV1990_WAVELENGTH_STEP_NM/2, MV1990_WAVELENGTH_STEP_NM)
    response = photopic_response(wavelengths)
    response_norm = np.trapezoid(response, wavelengths)
    samples, area_samples = aperture_quadrature(record["shadow_center_offset_m"],
                                                record["aperture_radius_m"],
                                                MV1990_APERTURE_ORDER)
    table = ray_table(PRODUCTION_VARIANT)
    members = []
    for sample_radius, sample_area in zip(samples, area_samples):
        radius = float(sample_radius)
        bb, area = receiver_quadrature(table, radius, geometry)
        outer = vacuum_branch(table, radius, geometry)
        members.append(Member(bb, area, float(sample_area), outer*float(sample_area)))
    observed = record["prediction"]/record["model_to_observed"]
    readout = ScalarReadout("whole_disk_nominal_v", "fit", None, float(observed),
                            tuple(members))
    return EventKernel("event_1990_02_09", (readout,), (),
                       photopic_grid=(wavelengths, response, response_norm))


@lru_cache(maxsize=None)
def event_kernels() -> dict:
    """All four event kernels, keyed in the preregistered event order."""
    builders = {"june_2011": june_event_kernel, "hst_2019": hst_event_kernel,
                "dec_2011": dec_event_kernel,
                "event_1990_02_09": mv1990_event_kernel}
    if set(builders) != set(EVENT_ORDER):
        raise RuntimeError("Event registry drifted from EVENT_ORDER")
    return {event: builders[event]() for event in EVENT_ORDER}


@lru_cache(maxsize=None)
def fit_restricted_kernels() -> dict:
    """Per event, the kernel restricted to its fit subset (root-finding use).

    Values are bitwise-identical to the corresponding full-kernel fit values
    (scalar subsets drop readouts; HST slices the wavelength grid, whose
    columns are independent); ``invert_event`` asserts the equality.
    """
    out = {}
    for event, kernel in event_kernels().items():
        if kernel.rows:
            out[event] = EventKernel(kernel.event, (),
                                     tuple(row.fit_only() for row in kernel.rows))
        else:
            out[event] = EventKernel(
                kernel.event,
                tuple(r for r in kernel.scalars if r.subset == "fit"), (),
                kernel.photopic_grid)
    return out


# --- Prediction evaluation (fused kernel_weights path, adapter orders) ----------


def _scalar_prediction(table, readout: ScalarReadout, wavelength, extra) -> float:
    """June/December moment aggregation (``egd.dec_kernel`` accumulation order;
    June's single member mirrors ``egd.june_kernel``)."""
    moment = 0.0
    for member in readout.members:
        weights = kernel_weights(table, member.bb, member.area, wavelength, extra)
        moment += float(np.sum(weights))*member.share + member.vacuum_share
    return moment


def _photopic_prediction(table, readout: ScalarReadout, grid, extra) -> float:
    """Whole-disk photopic aggregation (``egd.mv1990_kernel`` order)."""
    wavelengths, response, response_norm = grid
    moment = 0.0
    for member in readout.members:
        weights = kernel_weights(table, member.bb, member.area, wavelengths, extra)
        photopic = np.trapezoid(weights*response[None, :], wavelengths,
                                axis=1)/response_norm
        moment += float(photopic.sum())*member.share + member.vacuum_share
    return moment


def _row_predictions(table, row_kernel: RowKernel, extra) -> np.ndarray:
    """Per-wavelength predictions of one HST dataset (``egd.hst_kernel`` order)."""
    contributions = []
    for member in row_kernel.members:
        weights = kernel_weights(table, member.bb, member.area,
                                 row_kernel.wavelengths, extra)
        contributions.append(np.sum(weights, axis=0)*member.share
                             + member.vacuum_share)
    return np.sum(contributions, axis=0)*row_kernel.scaling


def event_predictions(kernel: EventKernel, aod550: float, form: str,
                      table=None) -> dict:
    """Per-readout predictions of one event at increment AOD550 on ``form``.

    ``table`` overrides the production two-layer table (the anchoring chain
    gate evaluates the committed fixed_opac variant through the same cached
    quadrature; impact parameters and landing maps are aerosol-independent,
    so the quadrature is bitwise identical across variants).
    ``aod550 = 0.0`` restores the table's own baseline bitwise.
    """
    if table is None:
        table = ray_table(PRODUCTION_VARIANT)
    extra = extra_at(form, aod550)
    predictions = {}
    for readout in kernel.scalars:
        if kernel.photopic_grid is not None:
            predictions[readout.readout_id] = _photopic_prediction(
                table, readout, kernel.photopic_grid, extra)
        else:
            predictions[readout.readout_id] = _scalar_prediction(
                table, readout, readout.wavelength_nm, extra)
    for row_kernel in kernel.rows:
        predictions[row_kernel.readout_id] = _row_predictions(table, row_kernel, extra)
    return predictions


def ln_residual_values(kernel: EventKernel, predictions: dict, subset: str) -> list:
    """Per-readout ln(P/O) values of one subset (HST rows pooled jointly)."""
    values = []
    for readout in kernel.scalars:
        if readout.subset == subset:
            observed = readout.observed
            predicted = predictions[readout.readout_id]
            if not observed > 0:
                raise RuntimeError(f"{kernel.event}/{readout.readout_id}: nonpositive "
                                   f"observed {observed}")
            if not predicted > 0:
                raise RuntimeError(f"{kernel.event}/{readout.readout_id}: nonpositive "
                                   f"prediction {predicted} in the {subset} account")
            values.append(math.log(predicted/observed))
    for row_kernel in kernel.rows:
        mask = row_kernel.fit_mask if subset == "fit" else row_kernel.holdout_mask
        predicted = predictions[row_kernel.readout_id][mask]
        observed = row_kernel.observed[mask]
        if not np.all(observed > 0):
            raise RuntimeError(f"{kernel.event}/{row_kernel.readout_id}: nonpositive "
                               f"observed rows")
        if not np.all(predicted > 0):
            raise RuntimeError(f"{kernel.event}/{row_kernel.readout_id}: nonpositive "
                               f"prediction rows in the {subset} account")
        values.extend(np.log(predicted/observed).tolist())
    return values


def fit_median(kernel: EventKernel, predictions: dict) -> float:
    return float(np.median(ln_residual_values(kernel, predictions, "fit")))


# --- Anchoring gates -------------------------------------------------------------


@lru_cache(maxsize=1)
def committed_strat_predictions() -> dict:
    return json.loads((ROOT/COMMITTED_STRAT_PREDICTIONS).read_text())["predictions"]


def _committed_chain_rows(kernel: EventKernel) -> tuple[list, str]:
    """Per-readout (O, committed prediction, committed model_to_observed) of the
    committed fixed_opac variant, flattened in the kernel's readout order.

    Returns the rows and the exact-leg note: for June/December/HST the
    observed field is a committed product field and the committed ratio was
    computed as prediction/observed, so the identity must hold bitwise; for
    1990 the observed is DEFINED as prediction/model_to_observed, so the
    exact leg is the recovery of the record's own observed field.
    """
    rows = []
    if kernel.event == "june_2011":
        for record in june_records():
            rows.append((float(record["observed_curve_brightness"]),
                         float(record["fixed_opac_prediction"]),
                         float(record["fixed_opac_model_to_observed"])))
        note = "committed_prediction/observed == committed model_to_observed bitwise"
    elif kernel.event == "dec_2011":
        for record in dec_comparison()["records"]:
            rows.append((float(record["observed_curve_ratio"]),
                         float(record["fixed_opac_prediction"]),
                         float(record["fixed_opac_model_to_observed"])))
        note = "committed_prediction/observed == committed model_to_observed bitwise"
    elif kernel.event == "hst_2019":
        _, grid_rows = hst_grid_data()
        for row_kernel in kernel.rows:
            dataset = row_kernel.readout_id
            for row in grid_rows:
                wl = float(row["wavelength_nm"])
                if not (UMBRAL_QUANTITATIVE_MIN_NM <= wl <= G430L_MAX_NM):
                    continue
                rows.append((float(row[f"{dataset}_obs_ratio"]),
                             float(row[f"{dataset}_fixed_opac_prediction"]),
                             float(row[f"{dataset}_fixed_opac_model_to_observed"])))
        note = "committed_prediction/observed == committed model_to_observed bitwise"
    else:
        record = mv1990_record("fixed_opac")
        observed = record["prediction"]/record["model_to_observed"]
        if observed != record["observed_minimum_to_reference"]:
            raise RuntimeError("1990 observed recovery drifted from the committed "
                               "record's observed_minimum_to_reference")
        rows.append((float(observed), float(record["prediction"]),
                     float(record["model_to_observed"])))
        note = "observed = prediction/model_to_observed recovers the record's own " \
               "observed_minimum_to_reference bitwise"
    return rows, note


def _kernel_chain_predictions(kernel: EventKernel, table) -> list:
    """Per-readout (O, kernel P(0)) pairs, flattened in the same order as
    ``_committed_chain_rows``."""
    predictions = event_predictions(kernel, 0.0, MAIN_FORM, table=table)
    rows = []
    for readout in kernel.scalars:
        rows.append((readout.observed, predictions[readout.readout_id]))
    for row_kernel in kernel.rows:
        window = row_kernel.fit_mask | row_kernel.holdout_mask
        observed = row_kernel.observed[window]
        mine = predictions[row_kernel.readout_id][window]
        for o, m in zip(observed, mine):
            rows.append((float(o), float(m)))
    return rows


def _committed_subset_flags(kernel: EventKernel) -> list:
    """Fit/holdout flags flattened in the same order as the chain rows.

    Row kernels expand per WINDOW row (``_committed_chain_rows`` slices the
    window via ``fit_mask | holdout_mask``), not per full-grid row - the two
    expansions must agree row for row or the committed-contrast pairing
    silently misaligns (code-review round 2 Major 1). The window partition is
    asserted in the kernel builder.
    """
    flags = [readout.subset for readout in kernel.scalars]
    for row_kernel in kernel.rows:
        window = row_kernel.fit_mask | row_kernel.holdout_mask
        flags.extend("fit" if m else "holdout"
                     for m in row_kernel.fit_mask[window])
    return flags


@lru_cache(maxsize=1)
def run_anchoring_gates() -> dict:
    """Anchoring gates for the inversion baseline (plan Step 1; fail loud).

    Legs: (a) the production two-layer P(0) reproduces the committed
    ``fixed_opac_strat`` predictions per readout within the preregistered
    reproduction tolerances imported from ``extinction_gap_diagnosis``;
    (b) the committed-variant chain identity - the committed prediction over
    the recovered observed equals the committed ``model_to_observed`` bitwise,
    and the kernel's own fixed_opac P(0)/O reproduces it within the moment
    tolerance (the mechanical criterion of the preregistered dimensional
    form); (c) anchored readout counts equal the preregistered counts;
    (d) the zeroth-moment gate of the committed variants (single-source
    ``build_height_stats``).
    """
    kernels = event_kernels()
    build_height_stats()  # moment gate; raises on any failure
    strat = committed_strat_predictions()
    reproduction = {}
    chain = {}
    counts = {}
    for event in EVENT_ORDER:
        kernel = kernels[event]
        fit_scalars = sum(1 for r in kernel.scalars if r.subset == "fit")
        holdout_scalars = sum(1 for r in kernel.scalars if r.subset == "holdout")
        fit_rows = sum(int(r.fit_mask.sum()) for r in kernel.rows)
        holdout_rows = sum(int(r.holdout_mask.sum()) for r in kernel.rows)
        counts[event] = {"fit": fit_scalars + fit_rows,
                         "holdout": holdout_scalars + holdout_rows}
        # (a) reproduction of the committed production predictions
        p0 = event_predictions(kernel, 0.0, MAIN_FORM)
        strat_records = strat[event]
        worst = 0.0
        within = True
        if event == "hst_2019":
            if len(strat_records) != len(kernel.rows):
                raise RuntimeError("hst_2019: strat record count mismatch")
            for row_kernel, strat_record in zip(kernel.rows, strat_records):
                window = row_kernel.fit_mask | row_kernel.holdout_mask
                if not np.array_equal(np.array(strat_record["window_wavelength_nm"]),
                                      row_kernel.wavelengths[window]):
                    raise RuntimeError(f"hst_2019/{row_kernel.readout_id}: window "
                                       f"wavelength grid drifted vs committed rerun")
                committed = np.array(strat_record["prediction"][PRODUCTION_VARIANT])
                mine = p0[row_kernel.readout_id][window]
                if len(committed) != len(mine):
                    raise RuntimeError(f"hst_2019/{row_kernel.readout_id}: window "
                                       f"length mismatch vs committed predictions")
                diffs = np.abs(mine-committed)
                tolerances = np.maximum(PREREGISTERED_RELATIVE_TOLERANCE*np.abs(committed),
                                        PREREGISTERED_FLOOR_ABSOLUTE)
                worst = max(worst, float(np.max(diffs/np.abs(committed))))
                within = within and bool(np.all(diffs <= tolerances))
        else:
            if len(strat_records) != len(kernel.scalars):
                raise RuntimeError(f"{event}: strat record count mismatch")
            for readout, strat_record in zip(kernel.scalars, strat_records):
                committed = strat_record["prediction"][PRODUCTION_VARIANT]
                mine = p0[readout.readout_id]
                diff = abs(mine-committed)
                tolerance = max(PREREGISTERED_RELATIVE_TOLERANCE*abs(committed),
                                PREREGISTERED_FLOOR_ABSOLUTE)
                worst = max(worst, diff/abs(committed))
                within = within and diff <= tolerance
        if not within:
            raise RuntimeError(f"{event}: production reproduction outside the "
                               f"preregistered tolerance (worst relative {worst:.3e})")
        reproduction[event] = worst
        # (b) committed-variant chain identity
        committed_rows, exact_note = _committed_chain_rows(kernel)
        kernel_rows = _kernel_chain_predictions(kernel, ray_table("fixed_opac"))
        if len(committed_rows) != len(kernel_rows):
            raise RuntimeError(f"{event}: chain row count mismatch")
        kernel_worst = 0.0
        for (observed, committed_prediction, committed_ratio), (observed_k, mine) \
                in zip(committed_rows, kernel_rows):
            if committed_prediction/observed != committed_ratio:
                raise RuntimeError(f"{event}: exact chain leg failed ({exact_note})")
            kernel_worst = max(kernel_worst,
                               abs(mine/observed-committed_ratio)/abs(committed_ratio))
        if kernel_worst > MOMENT_TOLERANCE:
            raise RuntimeError(f"{event}: kernel chain identity {kernel_worst:.3e} "
                               f"exceeds moment tolerance {MOMENT_TOLERANCE:.0e}")
        chain[event] = kernel_worst
    # preregistered count assertions (review round 2 item 1 + plan M2)
    expected = {"june_2011": (JUNE_FIT_COUNT, JUNE_HOLDOUT_COUNT),
                "dec_2011": (DEC_FIT_COUNT, DEC_HOLDOUT_COUNT),
                "hst_2019": (HST_FIT_ROW_COUNT*HST_UMBRAL_DATASET_COUNT,
                             HST_HOLDOUT_ROW_COUNT*HST_UMBRAL_DATASET_COUNT),
                "event_1990_02_09": (MV1990_FIT_COUNT, 0)}
    for event, (fit, holdout) in expected.items():
        if (counts[event]["fit"], counts[event]["holdout"]) != (fit, holdout):
            raise RuntimeError(f"{event}: anchored readout counts {counts[event]} "
                               f"!= preregistered {(fit, holdout)}")
    return {
        "reproduction": {
            "rule": f"abs_diff <= max({PREREGISTERED_RELATIVE_TOLERANCE:g} * "
                    f"|committed|, {PREREGISTERED_FLOOR_ABSOLUTE:g}) per readout "
                    f"(tolerances imported from extinction_gap_diagnosis)",
            "worst_relative_diff": reproduction,
            "verdict": "within_tolerance"},
        "chain_identity": {
            "exact_leg": _committed_chain_rows(kernels["june_2011"])[1],
            "exact_leg_1990": _committed_chain_rows(
                kernels["event_1990_02_09"])[1],
            "kernel_leg_worst_relative_dev": chain,
            "budget": MOMENT_TOLERANCE,
            "verdict": "within_tolerance"},
        "anchored_readout_counts": counts,
        "moment_gate": "zeroth-moment gate over the committed variants passed "
                       "(build_height_stats, single source)",
    }


# --- Inversion -------------------------------------------------------------------


def _fit_g(fit_kernel: EventKernel, a: float, form: str) -> float | None:
    """g(a) of the preregistered objective, or None when the model is dark.

    None means at least one fit-subset prediction underflowed to <= 0 at this
    increment: the model transmits no measurable light, so the residual is
    below every finite value (the vacuum branch's ~1e-16 float noise is the
    only thing that can push the summed weights negative). The preregistered
    branch_dark_top rule consumes this.
    """
    predictions = event_predictions(fit_kernel, a, form)
    for readout in fit_kernel.scalars:
        if not predictions[readout.readout_id] > 0:
            return None
    for row_kernel in fit_kernel.rows:
        if not np.all(predictions[row_kernel.readout_id] > 0):
            return None
    return fit_median(fit_kernel, predictions)


def invert_event(kernel: EventKernel, fit_kernel: EventKernel, form: str) -> dict:
    """One event's preregistered inversion: bracket checks, then brentq.

    The g evaluations walk the fit-restricted kernel (bitwise-equal to the
    full kernel's fit values, asserted against the full evaluation of g(0));
    the product-side predictions come from the full kernel.
    """
    p0 = event_predictions(kernel, 0.0, form)
    g0_full = fit_median(kernel, p0)
    g0 = _fit_g(fit_kernel, 0.0, form)
    if g0 is None:
        raise RuntimeError(f"{kernel.event}: g(0) dark - impossible: the baseline "
                           f"predictions are positive")
    if g0 != g0_full:
        raise RuntimeError(f"{kernel.event}: fit-restricted g(0) {g0!r} != full-grid "
                           f"g(0) {g0_full!r} (wavelength subsetting must be bitwise)")
    if g0 <= 0.0:
        return {"p0": p0, "preds_root": p0, "status": STATUS_ZERO_INCREMENT,
                "aod550": 0.0, "g0": g0, "g_bracket_high": None,
                "bracket_high_used": 0.0}

    g_high = _fit_g(fit_kernel, FIT_BRACKET_HIGH, form)
    if g_high is not None and g_high > 0.0:
        return {"p0": p0, "preds_root": None, "status": STATUS_UNCLOSED,
                "aod550": None, "g0": g0, "g_bracket_high": g_high,
                "bracket_high_used": FIT_BRACKET_HIGH}
    bracket_high_used = FIT_BRACKET_HIGH
    if g_high is None:
        # Preregistered branch_dark_top: geometric halving to the first finite
        # point, guaranteed negative (within a factor 2 of the underflow cliff).
        probe = FIT_BRACKET_HIGH
        g_bracket = None
        for _ in range(60):
            probe /= 2.0
            g_bracket = _fit_g(fit_kernel, probe, form)
            if g_bracket is not None:
                break
        if g_bracket is None:
            raise RuntimeError(f"{kernel.event}: no finite g below the dark bracket "
                               f"top after 60 halvings")
        if g_bracket > 0.0:
            raise RuntimeError(f"{kernel.event}: first finite point after a dark "
                               f"bracket top is positive ({g_bracket!r}) - violates "
                               f"the underflow-cliff argument; investigate")
        bracket_high_used = probe
        g_high = g_bracket
    root = float(brentq(lambda a: _fit_g(fit_kernel, a, form),
                        0.0, bracket_high_used, xtol=BRENTQ_XTOL))
    preds_root = event_predictions(kernel, root, form)
    return {"p0": p0, "preds_root": preds_root, "status": STATUS_INVERTED,
            "aod550": root, "g0": g0, "g_bracket_high": g_high,
            "bracket_high_used": bracket_high_used}


def holdout_account(kernel: EventKernel, p0: dict, preds_root: dict | None) -> dict | None:
    """Holdout-subset ln-residual median and RMS, before -> after.

    Returns None when the event has no event-internal holdout (1990, whose
    holdout account is the cross-event prediction matrix row).
    """
    before = ln_residual_values(kernel, p0, "holdout")
    if not before:
        return None
    after = (ln_residual_values(kernel, preds_root, "holdout")
             if preds_root is not None else None)
    return {"readout_count": len(before),
            "ln_residual_median_before": float(np.median(before)),
            "ln_residual_median_after": (float(np.median(after))
                                         if after is not None else None),
            "ln_rms_before": float(np.sqrt(np.mean(np.square(before)))),
            "ln_rms_after": (float(np.sqrt(np.mean(np.square(after))))
                             if after is not None else None)}


# --- Shape ledgers ----------------------------------------------------------------


def june_shape_ledger(kernel: EventKernel, p0: dict, preds_root: dict) -> dict:
    """503/677 shape cost at the inverted increment (committed aggregation)."""
    by_angle = {}
    for readout in kernel.scalars:
        by_angle.setdefault(readout.shadow_angle_deg, {})[readout.wavelength_nm] = readout
    metrics = []
    for angle in sorted({r.shadow_angle_deg for r in kernel.scalars}):
        bands = by_angle[angle]
        if not {503.0, 677.0} <= bands.keys():
            continue
        blue, red = bands[503.0], bands[677.0]
        observed_ratio = blue.observed/red.observed
        base = (p0[blue.readout_id]/p0[red.readout_id])/observed_ratio
        factor = ((preds_root[blue.readout_id]/p0[blue.readout_id])
                  /(preds_root[red.readout_id]/p0[red.readout_id]))
        new = base*factor
        closed = shape_closed_fraction(base, new)
        metrics.append({"shadow_angle_deg": angle, "base_ratio": base,
                        "new_ratio": new, "shape_factor": factor,
                        "disagreement_closed_fraction": closed,
                        "verdict": closure_verdict(closed)})
    median_base = float(np.median([m["base_ratio"] for m in metrics]))
    median_new = float(np.median([m["new_ratio"] for m in metrics]))
    median_closed = shape_closed_fraction(median_base, median_new)
    return {"metric": "503_over_677_transmission_ratio_model_over_observed",
            "median_base": median_base, "median_new": median_new,
            "disagreement_closed_fraction": median_closed,
            "verdict": closure_verdict(median_closed),
            "involves_holdout_bands": True, "by_angle": metrics}


def dec_shape_ledger(kernel: EventKernel, p0: dict, preds_root: dict) -> dict:
    """450|500 over 750 shape cost at the inverted increment."""
    index = {readout.wavelength_nm: readout for readout in kernel.scalars}
    metrics = {}
    for key, blue_wl in [("450_to750_transmission_ratio_model_over_observed", 450.0),
                         ("500_to750_transmission_ratio_model_over_observed", 500.0)]:
        blue, red = index[blue_wl], index[750.0]
        base = (p0[blue.readout_id]/p0[red.readout_id])/(blue.observed/red.observed)
        factor = ((preds_root[blue.readout_id]/p0[blue.readout_id])
                  /(preds_root[red.readout_id]/p0[red.readout_id]))
        new = base*factor
        closed = shape_closed_fraction(base, new)
        metrics[key] = {"base_ratio": base, "new_ratio": new, "shape_factor": factor,
                        "disagreement_closed_fraction": closed,
                        "verdict": closure_verdict(closed)}
    return {"metrics": metrics, "involves_holdout_bands": True}


def _hst_band_slopes(prediction: np.ndarray, observed: np.ndarray,
                     wavelengths: np.ndarray) -> dict:
    """The two committed slope bands of one dataset (band-median ratios over
    the full window grid, the committed ``compare_hst2019`` structure)."""
    return {f"{num:g}_over_{den:g}":
            (band_median(prediction, wavelengths, num)
             /band_median(prediction, wavelengths, den))
            /(band_median(observed, wavelengths, num)
              /band_median(observed, wavelengths, den))
            for num, den in SLOPE_BANDS_NM}


def hst_shape_ledger(kernel: EventKernel, p0: dict, preds_root: dict) -> dict:
    """Two committed slope bands on the full window grid (committed structure)."""
    metrics = {}
    for row_kernel in kernel.rows:
        dataset = row_kernel.readout_id
        base_slopes = _hst_band_slopes(p0[dataset], row_kernel.observed,
                                       row_kernel.wavelengths)
        new_slopes = _hst_band_slopes(preds_root[dataset], row_kernel.observed,
                                      row_kernel.wavelengths)
        bands = {}
        for band, base_slope in base_slopes.items():
            new_slope = new_slopes[band]
            closed = shape_closed_fraction(base_slope, new_slope)
            bands[band] = {
                "base_slope": base_slope, "new_slope": new_slope,
                "disagreement_closed_fraction": closed,
                "verdict": closure_verdict(closed)}
        metrics[dataset] = bands
    return {"metrics_by_dataset": metrics, "involves_holdout_bands": True,
            "aggregation": "committed band-median structure over the full window grid"}


def shape_ledger(kernel: EventKernel, p0: dict, preds_root: dict | None):
    if preds_root is None:
        return None
    if kernel.event == "june_2011":
        return june_shape_ledger(kernel, p0, preds_root)
    if kernel.event == "dec_2011":
        return dec_shape_ledger(kernel, p0, preds_root)
    if kernel.event == "hst_2019":
        return hst_shape_ledger(kernel, p0, preds_root)
    return None  # 1990: no shape metric (preregistered)


# --- Physical anchor ----------------------------------------------------------------


def physical_anchor(value: float) -> dict:
    """Position of one inverted value in the pinned boundary-layer spectrum."""
    rows = sorted(((row["candidate"], float(row["value"])) for row in prior_rows()
                   if row["name"] == BL_AOD_ROW), key=lambda item: item[1])
    if len(rows) != 10:
        raise RuntimeError(f"Expected 10 pinned boundary-layer rows, found {len(rows)}")
    lo, hi = rows[0][1], rows[-1][1]
    if lo != PHYSICAL_ANCHOR_MIN or hi != PHYSICAL_ANCHOR_MAX:
        raise RuntimeError(f"Pinned boundary-layer spectrum drifted: [{lo}, {hi}] != "
                           f"preregistered [{PHYSICAL_ANCHOR_MIN}, {PHYSICAL_ANCHOR_MAX}]")
    below = [candidate for candidate, row_value in rows if row_value <= value]
    above = [candidate for candidate, row_value in rows if row_value >= value]
    return {"inverted_aod550": value,
            "pinned_rows": [{"candidate": candidate, "aod550_boundary_layer_only": row_value}
                            for candidate, row_value in rows],
            "ratio_to_spectrum_min": value/lo,
            "ratio_to_spectrum_max": value/hi,
            "within_pinned_range": lo <= value <= hi,
            "bracketed_by": [below[-1] if below else None,
                             above[0] if above else None]}


# --- Cross-event accounts ------------------------------------------------------------


def committed_fit_median(kernel: EventKernel) -> float:
    """Fit-subset median ln of the committed fixed_opac model_to_observed."""
    ratios = _committed_chain_rows(kernel)[0]
    flags = _committed_subset_flags(kernel)
    if len(flags) != len(ratios):
        raise RuntimeError(f"{kernel.event}: committed fit/holdout flags "
                           f"({len(flags)}) and chain rows ({len(ratios)}) "
                           f"misaligned - the contrast would silently mispair")
    values = [math.log(ratio) for (_, _, ratio), flag in zip(ratios, flags)
              if flag == "fit"]
    return float(np.median(values))


def consistency_verdict(events: dict) -> dict:
    """Mechanical application of the preregistered consistency criterion."""
    values = {}
    degenerate = {}
    for event, result in events.items():
        if result["status"] == STATUS_INVERTED:
            values[event] = result["aod550"]
        else:
            degenerate[event] = result["status"]
    if degenerate:
        return {"verdict": "not_applicable", "R": None,
                "inverted_aod550": values, "degenerate_events": degenerate,
                "rule": PREREGISTERED["consistency_criterion"]["degenerate_branch"]}
    aods = np.array([values[event] for event in EVENT_ORDER])
    r_value = float(aods.max()/aods.min())
    if r_value <= VERDICT_COMPATIBLE[1]:
        verdict = "systematic_gap"
    elif r_value <= VERDICT_TENSION[1]:
        verdict = "partially_consistent"
    else:
        verdict = "event_specific"
    logs = np.log(aods)
    return {"verdict": verdict, "R": r_value, "inverted_aod550": values,
            "degenerate_events": {},
            "rule": PREREGISTERED["consistency_criterion"]["statistic"],
            "geometric_mean": float(np.exp(np.mean(logs))),
            "geometric_std_population": float(np.exp(np.std(logs))),
            "band_edges": {"systematic_gap_upper": VERDICT_COMPATIBLE[1],
                           "partially_consistent_upper": VERDICT_TENSION[1]}}


def prediction_matrix(kernels: dict, fit_kernels: dict, events: dict) -> dict:
    """Fit-subset median ln residual of event i under event j's inverted AOD."""
    matrix = {}
    for event_i in EVENT_ORDER:
        fit_kernel = fit_kernels[event_i]
        row = {}
        for event_j in EVENT_ORDER:
            result = events[event_j]
            if result["status"] != STATUS_INVERTED:
                row[event_j] = None
            else:
                preds = event_predictions(fit_kernel, result["aod550"], MAIN_FORM)
                row[event_j] = fit_median(fit_kernel, preds)
        matrix[event_i] = row
    return matrix


def retrieval_spec(events: dict, findings: list) -> dict:
    """Retrieval spec for the follow-up prior-source survey (preregistered rules)."""
    valid = {event: result["aod550"] for event, result in events.items()
             if result["status"] == STATUS_INVERTED}
    zeros = [event for event, result in events.items()
             if result["status"] == STATUS_ZERO_INCREMENT]
    unclosed = [event for event, result in events.items()
                if result["status"] == STATUS_UNCLOSED]
    pool = dict(valid)
    for event in zeros:
        pool[event] = 0.0  # a zero-increment event enters min as 0 (review round 2)
    forms = (MAIN_FORM,)+SENSITIVITY_FORMS
    blue = [MORPHOLOGY_ROWS[form]["angstrom_0.35_to_0.5um"] for form in forms]
    visible = [MORPHOLOGY_ROWS[form]["angstrom_0.5_to_0.8um"] for form in forms]
    spec = {"aod550_pool": {event: pool.get(event) for event in EVENT_ORDER},
            "excluded_events": {event: STATUS_UNCLOSED for event in unclosed},
            "alpha_blue_range": [min(blue), max(blue)],
            "alpha_visible_range": [min(visible), max(visible)],
            "alpha_envelope_source": "pinned morphology rows (main + sensitivity forms)",
            "layer_geometry": {"class": "0-10 km exponential layer",
                               "scale_height_m": MORPHOLOGY_ROWS[MAIN_FORM]["scale_height_m"],
                               "source_row": f"{MAIN_FORM}/scale_height_m"},
            "source_types": ["boundary-layer / lower-tropospheric aerosol optical "
                             "types (OPAC-like)",
                             "humidity variants (RH rows unavailable in the pinned "
                             "source; known gap)"],
            "branch_rule_fired": None}
    if not pool:
        findings.append({"code": "retrieval_range_degenerate",
                         "detail": "no event carries a valid inversion value or a "
                                   "zero increment; the AOD range rule has no input"})
        spec["aod550_range"] = None
        spec["branch_rule_fired"] = False
        return spec
    lo, hi = min(pool.values()), max(pool.values())
    spec["aod550_range"] = [lo/10.0, hi*10.0]
    branch = hi > PHYSICAL_ANCHOR_MAX
    if branch:
        spec["source_types"].append("volcanic / stratospheric-perturbation types")
    spec["branch_rule_fired"] = branch
    return spec


# --- Product assembly ----------------------------------------------------------------

SCOPE_TEXT = (
    "Cross-event inversion of the low-layer aerosol increment required to close "
    "each event's fit-subset brightness residual on top of the production two-layer "
    "default (fixed_opac_strat). Calibration lane: every value is relative to the "
    "production default baseline and is calibrated, not a physical prior; nothing "
    "here enters the default model. Splits, objective, consistency criterion, "
    "morphologies and retrieval rules are preregistered in the 'preregistered' "
    "block (verbatim module constant, pytest-pinned); the cross-event verdict is "
    "the mechanical application of the preregistered bands.")

LEDGER_FIELDS = ["event", "readout", "subset", "wavelength_nm", "observed",
                 "prediction_baseline", "prediction_at_root",
                 "ln_residual_baseline", "ln_residual_at_root",
                 "model_to_observed_baseline", "model_to_observed_at_root"]


def _ledger_rows(kernel: EventKernel, p0: dict, preds_root: dict | None) -> list:
    rows = []
    for readout in kernel.scalars:
        baseline = p0[readout.readout_id]
        after = preds_root[readout.readout_id] if preds_root is not None else None
        rows.append([kernel.event, readout.readout_id, readout.subset,
                     readout.wavelength_nm, readout.observed,
                     baseline, after,
                     math.log(baseline/readout.observed),
                     math.log(after/readout.observed) if after is not None else None,
                     baseline/readout.observed,
                     after/readout.observed if after is not None else None])
    for row_kernel in kernel.rows:
        for index in np.flatnonzero(row_kernel.fit_mask | row_kernel.holdout_mask):
            wl = float(row_kernel.wavelengths[index])
            subset = ("fit" if row_kernel.fit_mask[index]
                      else "holdout" if row_kernel.holdout_mask[index] else "outside")
            baseline = float(p0[row_kernel.readout_id][index])
            after = (float(preds_root[row_kernel.readout_id][index])
                     if preds_root is not None else None)
            observed = float(row_kernel.observed[index])
            rows.append([kernel.event, f"{row_kernel.readout_id}@{wl:g}nm", subset,
                         wl, observed, baseline, after,
                         math.log(baseline/observed),
                         math.log(after/observed) if after is not None else None,
                         baseline/observed,
                         after/observed if after is not None else None])
    return rows


def build_inversion() -> tuple:
    """Build the product and its ledger rows (deterministic; no timestamps)."""
    findings = []
    kernels = event_kernels()
    anchoring = run_anchoring_gates()
    fit_kernels = fit_restricted_kernels()
    events = {}
    ledger = []
    for event in EVENT_ORDER:
        kernel = kernels[event]
        result = invert_event(kernel, fit_kernels[event], MAIN_FORM)
        p0, preds_root = result["p0"], result["preds_root"]
        if result["status"] == STATUS_ZERO_INCREMENT:
            findings.append({"code": "zero_increment", "event": event,
                             "detail": "g(0) <= 0: the production baseline already "
                                       "closes the fit-subset residual; the required "
                                       "increment is 0"})
        if result["status"] == STATUS_UNCLOSED:
            findings.append({"code": "unclosed", "event": event,
                             "detail": f"g({FIT_BRACKET_HIGH:g}) > 0: the preregistered "
                                       f"bracket cannot close this event; recorded "
                                       f"unclosed, not forced (lane discipline)"})
        fit_after = (fit_median(kernel, preds_root)
                     if preds_root is not None else None)
        shape = shape_ledger(kernel, p0, preds_root)
        anchor = (physical_anchor(result["aod550"])
                  if result["aod550"] is not None else None)
        if anchor is not None and not anchor["within_pinned_range"]:
            findings.append({"code": "physical_range_exceeded", "event": event,
                             "detail": f"inverted AOD550 {result['aod550']:.6g} lies "
                                       f"outside the pinned boundary-layer spectrum "
                                       f"[{PHYSICAL_ANCHOR_MIN}, {PHYSICAL_ANCHOR_MAX}]"})
        events[event] = {
            "morphology": MAIN_FORM,
            "observation_source": PREREGISTERED["observation_sources"][event],
            "readout_counts": dict(anchoring["anchored_readout_counts"][event]),
            "g0": result["g0"],
            "g_at_bracket_high": result["g_bracket_high"],
            "bracket_high_used": result["bracket_high_used"],
            "status": result["status"],
            "aod550": result["aod550"],
            "fit": {"ln_residual_median_before": result["g0"],
                    "ln_residual_median_after": fit_after,
                    "readout_count": anchoring["anchored_readout_counts"][event]["fit"]},
            "holdout": holdout_account(kernel, p0, preds_root),
            "baseline_contrast": {
                "fit_ln_residual_median_fixed_opac_committed":
                    committed_fit_median(kernel),
                "fit_ln_residual_median_production_two_layer": result["g0"],
                "note": "monitoring contrast only (plan risk 1); the inversion "
                        "baseline is the production two-layer default"},
            "shape": shape,
            "physical_anchor": anchor,
        }
        ledger.extend(_ledger_rows(kernel, p0, preds_root))
    consistency = consistency_verdict(events)
    if consistency["verdict"] == "not_applicable":
        findings.append({"code": "consistency_degenerate",
                         "detail": f"events without a valid inversion value: "
                                   f"{consistency['degenerate_events']}; the max/min "
                                   f"ratio degenerates (preregistered branch)"})
    cross_event = dict(consistency)
    cross_event["prediction_matrix"] = prediction_matrix(kernels, fit_kernels, events)
    cross_event["matrix_note"] = ("event i rows: fit-subset median ln residual under "
                                  "event j's inverted AOD550; the diagonal is the "
                                  "event's own root; the 1990 row is its holdout "
                                  "account")
    sensitivity = {}
    for form in SENSITIVITY_FORMS:
        per_event = {}
        for event in EVENT_ORDER:
            result = invert_event(kernels[event], fit_kernels[event], form)
            if result["status"] == STATUS_UNCLOSED:
                findings.append({"code": "sensitivity_unclosed", "event": event,
                                 "detail": f"sensitivity form {form}: g("
                                           f"{FIT_BRACKET_HIGH:g}) > 0 - the pinned "
                                           f"layer geometry cannot close this event "
                                           f"within the preregistered bracket "
                                           f"(appendix only; does not affect the "
                                           f"main-form conclusion)"})
            elif result["status"] == STATUS_ZERO_INCREMENT:
                findings.append({"code": "sensitivity_zero_increment", "event": event,
                                 "detail": f"sensitivity form {form}: g(0) <= 0 - the "
                                           f"production baseline already closes the "
                                           f"fit subset under this morphology"})
            holdout = holdout_account(kernels[event], result["p0"], result["preds_root"])
            per_event[event] = {
                "status": result["status"], "aod550": result["aod550"],
                "holdout_ln_rms_before": (holdout["ln_rms_before"]
                                          if holdout is not None else None),
                "holdout_ln_rms_after": (holdout["ln_rms_after"]
                                         if holdout is not None
                                         and result["preds_root"] is not None else None)}
        sensitivity[form] = per_event
    form_sensitivity = {"forms": sensitivity,
                        "scope": "appendix only: the main conclusion is carried by "
                                 "the antarctic main form; objective, splits and root "
                                 "finder identical to the preregistered contract"}
    retrieval = retrieval_spec(events, findings)
    product = {
        "scope": SCOPE_TEXT,
        "preregistered": PREREGISTERED,
        "anchoring": anchoring,
        "events": events,
        "cross_event": cross_event,
        "form_sensitivity": form_sensitivity,
        "retrieval_spec": retrieval,
        "findings": findings,
    }
    return product, ledger


def build_inversion_product() -> dict:
    """Product only (ledger rows discarded)."""
    return build_inversion()[0]


def write_products(out: Path) -> dict:
    """Write ``inversion.json`` + ``ledger.csv`` into ``out`` and summarize."""
    product, ledger = build_inversion()
    out.mkdir(parents=True, exist_ok=True)
    (out/"inversion.json").write_text(json.dumps(product, indent=2)+"\n")
    with (out/"ledger.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(LEDGER_FIELDS)
        for row in ledger:
            writer.writerow([json.dumps(value) if isinstance(value, (list, dict))
                             else ("" if value is None else value) for value in row])
    return {
        "anchoring": {name: product["anchoring"][name]["verdict"]
                      for name in ("reproduction", "chain_identity")},
        "events": {event: {"status": result["status"], "aod550": result["aod550"]}
                   for event, result in product["events"].items()},
        "cross_event": {"R": product["cross_event"]["R"],
                        "verdict": product["cross_event"]["verdict"]},
        "retrieval_spec": {"aod550_range": product["retrieval_spec"]["aod550_range"]},
        "findings": product["findings"],
        "ledger_rows": len(ledger),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("outputs/atmosphere_inversion"))
    args = parser.parse_args()
    summary = write_products(args.out)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
