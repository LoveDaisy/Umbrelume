"""Contribution-path altitude statistics, cross-event correlation and
missing-extinction discrimination across the four committed observational
comparisons (June 2011 curve readouts, 1990-02-09 whole-disk photometry,
2011-12-10 seven-point spectra, 2019-01-21 HST/STIS umbral spectra).

This module is the single authoritative implementation of contribution-kernel
altitude statistics for all four apertures. ``round2_analysis.altitude_quantiles``
is a historical product generator only: never extend it; its committed June
2011 values are the independent anchor this module reproduces exactly.

Kernel identity (zeroth-moment gate): every aperture kernel reuses the same
grid (``RayTable.impact_m``) and the same summation path with which
``eclipse.transport.transmission_profile`` produced the committed predictions,
so ``sum(weights) + vacuum branch == committed prediction`` holds by
construction up to a documented float-summation tolerance, checked in pytest.
The vacuum branch is the ``transmission_profile`` ``outer`` term: rays tangent
above the table top that traverse no atmosphere. It enters the moment
identity but never the atmospheric-branch statistics, because no candidate
component acts on it.

Counterfactual probes (missing-extinction candidates) evaluate through these
same adapters with an extra extinction column or a replaced aerosol column,
never a second evaluation path; a zero-AOD probe restores the fixed OPAC
baseline exactly.

Command products: ``height_stats``/``correlation``/``discrimination``
regenerate byte-identically; ``sources``/``prior_parameters`` are one-shot
literature products that the command only schema-checks when present.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.spectra import ozone_cross_section, rayleigh_cross_section
from eclipse.transport import Geometry, circle_overlap, receiver_quadrature

from .compare_dec2011 import exposure_average
from .compare_hst2019 import BAND_HALFWIDTH_NM, G430L_MAX_NM, SLOPE_BANDS_NM, UMBRAL_QUANTITATIVE_MIN_NM
from .multiple_events import aperture_quadrature, photopic_response

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"validation/results"
OBSERVATIONS = ROOT/"validation/observations"
VARIANTS = ("molecular", "fixed_opac")
QUANTILE_LEVELS = (.1, .5, .9)
# Component boundaries of the current default physics, taken from the default
# profile: 0-2 km bare, 2-12 km OPAC free troposphere, >12 km no component.
DEFAULT_OPAC = AerosolProfile.opac_free_troposphere()
BAND_BELOW_TOP_M = DEFAULT_OPAC.bottom_m
BAND_OVER_TOP_M = DEFAULT_OPAC.top_m
# Structural reuse makes the zeroth-moment equality hold by construction; the
# budget covers float summation-order differences only (plan design point 1).
MOMENT_TOLERANCE = 1e-12
# The 1990 branch reuses the committed refined photopic quadrature exactly.
MV1990_EVENT = "1990-02-09"
MV1990_RESPONSE = "nominal_v"
MV1990_APERTURE_ORDER = 256
MV1990_WAVELENGTH_STEP_NM = 2.5


@lru_cache(maxsize=None)
def ray_table(variant: str):
    """The two committed model variants' tables (molecular / fixed OPAC)."""
    if variant == "molecular":
        return build_ray_table(Atmosphere())
    if variant == "fixed_opac":
        return build_ray_table(Atmosphere(), aerosol_profile=DEFAULT_OPAC)
    raise ValueError(f"Unknown variant {variant!r}")


def vacuum_branch(table, radius_m: float, geometry: Geometry) -> float:
    """The ``transmission_profile`` outer term, verbatim.

    Rays whose tangent height lies above the table top contribute transmission
    one and are untouchable by any candidate component; they enter the
    zeroth-moment identity but never the atmospheric-branch statistics.
    """
    a = geometry.solar_blur_m
    outer = 1 - circle_overlap(np.array([radius_m]), a,
                               geometry.magnification*table.impact_m[-1])/(np.pi*a*a)
    return float(outer[0])


def kernel_weights(table, bb, area, wavelength_nm, extra_aerosol=None) -> np.ndarray:
    """Per-impact weights ``area*exp(-tau)`` on the transmission_profile path.

    tau mirrors ``transmission_profile`` term by term: Rayleigh then ozone,
    then the table's aerosol column times its spectral factor, then the
    optional counterfactual extra column evaluated on its own grid.
    ``wavelength_nm`` may be a scalar (June / December readouts) or a vector
    (HST window, 1990 photopic grid). With ``extra_aerosol=None`` and the
    committed tables this reproduces the committed prediction kernel exactly.
    """
    wl = np.asarray(wavelength_nm, dtype=float)
    scalar = wl.ndim == 0
    if scalar:
        wl = wl[None]
    b = table.impact_m
    air = np.interp(bb, b, table.air_column_m2)[:, None]
    o3 = np.interp(bb, b, table.ozone_column_m2)[:, None]
    tau = air*rayleigh_cross_section(wl)[None, :] + o3*ozone_cross_section(wl)[None, :]
    if table.aerosol_tau550 is not None:
        tau = tau + (np.interp(bb, b, table.aerosol_tau550)[:, None]
                     * table.aerosol_profile.spectral_factor(wl)[None, :])
    if extra_aerosol is not None:
        extra_b, extra_tau550, extra_profile = extra_aerosol
        tau = tau + (np.interp(bb, extra_b, extra_tau550)[:, None]
                     * extra_profile.spectral_factor(wl)[None, :])
    weights = area[:, None]*np.exp(-tau)
    return weights[:, 0] if scalar else weights


def altitude_statistics(table, node_bbs, node_weights, vacuum_weight=0.0) -> dict:
    """P10/P50/P90 tangent heights and three-band fractions, atmospheric branch.

    Quantiles follow the round2 ``altitude_quantiles`` algorithm (CDF linear
    interpolation onto tangent heights) so the committed June 2011 anchor is
    exactly reproducible. Band fractions normalise on the atmospheric branch
    (quadrature nodes only) and always sum to one; the vacuum branch is
    reported separately because no candidate component acts on it.
    """
    bb = np.concatenate(node_bbs)
    weights = np.concatenate(node_weights)
    order = np.argsort(bb, kind="stable")
    heights = np.interp(bb[order], table.impact_m, table.altitude_m)
    weights = weights[order]
    total = float(weights.sum())
    if not total > 0:
        raise ValueError("Atmospheric kernel weight must be positive")
    cdf = np.cumsum(weights)/total
    quantiles = np.interp(list(QUANTILE_LEVELS), cdf, heights).tolist()
    below = weights[heights < BAND_BELOW_TOP_M]
    over = weights[heights > BAND_OVER_TOP_M]
    vacuum = float(vacuum_weight/(total + vacuum_weight)) if vacuum_weight > 0 else 0.0
    return {"p10_p50_p90_height_m": quantiles,
            "band_fraction_below_2km": float(below.sum()/total),
            "band_fraction_2_to_12km": float((total - below.sum() - over.sum())/total),
            "band_fraction_over_12km": float(over.sum()/total),
            "vacuum_branch_fraction": vacuum}


def moment_gap(moment: float, committed: float) -> dict:
    """Zeroth-moment gate record for one readout."""
    gap = abs(moment - committed)
    return {"kernel_zeroth_moment": moment, "committed_prediction": committed,
            "relative_gap": gap/abs(committed) if committed else None}


def june_geometry() -> Geometry:
    data = json.loads((OBSERVATIONS/"ugolnikov2011/geometry.json").read_text())
    return Geometry(moon_distance_m=data["moon_distance_m"],
                    sun_distance_m=data["sun_distance_m"])


@lru_cache(maxsize=None)
def june_records() -> tuple:
    """Committed June 2011 readouts and predictions (round2 analysis.json)."""
    return tuple(json.loads((RESULTS/"round2/analysis.json").read_text())["records"])


def june_kernel(variant: str, radius_m: float, wavelength_nm: float,
                extra_aerosol=None) -> tuple:
    """One June 2011 readout kernel: (bb, weights, vacuum, zeroth moment)."""
    table = ray_table(variant)
    bb, area = receiver_quadrature(table, float(radius_m), june_geometry())
    weights = kernel_weights(table, bb, area, wavelength_nm, extra_aerosol)
    vacuum = vacuum_branch(table, float(radius_m), june_geometry())
    return bb, weights, vacuum, float(np.sum(weights)) + vacuum


def june_anchor_quantiles(radius_m: float, wavelength_nm: float) -> list:
    """Molecular height quantiles of one June readout (round2 anchor path)."""
    bb, weights, _, _ = june_kernel("molecular", radius_m, wavelength_nm)
    return altitude_statistics(ray_table("molecular"), [bb], [weights])["p10_p50_p90_height_m"]


def june_height_records() -> list:
    geometry = june_geometry()
    records = []
    for record in june_records():
        radius = geometry.moon_distance_m*np.tan(np.deg2rad(record["shadow_angle_deg"]))
        for variant in VARIANTS:
            bb, weights, vacuum, moment = june_kernel(variant, radius, record["wavelength_nm"])
            records.append({"event": "june_2011",
                            "readout": f'{record["wavelength_nm"]:g}nm@{record["shadow_angle_deg"]:g}deg',
                            "wavelength_nm": record["wavelength_nm"],
                            "shadow_angle_deg": record["shadow_angle_deg"],
                            "variant": variant,
                            **altitude_statistics(ray_table(variant), [bb], [weights], vacuum),
                            **moment_gap(moment, record[f"{variant}_prediction"])})
    return records


@lru_cache(maxsize=None)
def hst_comparison() -> dict:
    return json.loads((RESULTS/"hst2019/comparison.json").read_text())


@lru_cache(maxsize=None)
def hst_grid_data() -> tuple:
    """Committed full-window wavelengths and comparison.csv rows (300-570 nm).

    The committed ``slope_model_vs_observed`` aggregates band medians on this
    full grid (its 430 nm band includes rows in [425, 430) below the
    quantitative-window cutoff), so shape-axis probes must reuse it.
    """
    with (RESULTS/"hst2019/comparison.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    wavelengths = np.array([float(row["wavelength_nm"]) for row in rows])
    return wavelengths, rows


@lru_cache(maxsize=None)
def hst_window_data() -> tuple:
    """Committed quantitative-window wavelengths and comparison.csv rows."""
    wavelengths, rows = hst_grid_data()
    window = (wavelengths >= UMBRAL_QUANTITATIVE_MIN_NM) & (wavelengths <= G430L_MAX_NM)
    return wavelengths[window], [row for row, keep in zip(rows, window) if keep]


def band_median(values: np.ndarray, wavelength_nm: np.ndarray, center_nm: float) -> float:
    """The compare_hst2019 band-median aggregation (full window grid)."""
    mask = (wavelength_nm >= center_nm-BAND_HALFWIDTH_NM) & (wavelength_nm <= center_nm+BAND_HALFWIDTH_NM)
    return float(np.median(values[mask]))


@lru_cache(maxsize=None)
def hst_slope_baseline() -> dict:
    """Committed-aggregation HST slope baselines, anchored to the committed
    ``slope_model_vs_observed`` fields.

    The slopes are recomputed from the committed per-wavelength predictions
    and observed ratios with the exact ``compare_hst2019`` structure (band
    medians over the full window grid), then asserted equal to the committed
    metric fields; any drift fails loud instead of silently diverging (the
    asymmetry that let the shape baseline drift in review round 1).
    """
    wavelengths, rows = hst_grid_data()
    baselines = {}
    for variant in VARIANTS:
        baselines[variant] = {}
        for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]:
            observed = np.array([float(row[f"{dataset}_obs_ratio"]) for row in rows])
            committed = np.array([float(row[f"{dataset}_{variant}_prediction"]) for row in rows])
            slopes = {}
            for num, den in SLOPE_BANDS_NM:
                slopes[f"{num:g}_over_{den:g}"] = (
                    band_median(committed, wavelengths, num)
                    /band_median(committed, wavelengths, den)
                    /(band_median(observed, wavelengths, num)
                      /band_median(observed, wavelengths, den)))
            committed_slopes = hst_comparison()["metrics"][variant][dataset]["slope_model_vs_observed"]
            for key, value in slopes.items():
                if value != committed_slopes[key]:
                    raise RuntimeError(f"HST slope consistency check failed for "
                                       f"{dataset}/{variant}/{key}: recomputed {value!r} != "
                                       f"committed {committed_slopes[key]!r}")
            baselines[variant][dataset] = slopes
    return baselines


def hst_members(dataset: str) -> list:
    return [record for record in json.loads(
        (OBSERVATIONS/"hst2019/geometry.json").read_text())["records"]
        if record["dataset_id"] == dataset]


def hst_kernel(variant: str, dataset: str, wavelengths: np.ndarray,
               extra_aerosol=None) -> tuple:
    """Exposure-weighted joint kernel of one umbral dataset.

    Each member contributes texp-weighted nodes on its own quadrature grid;
    the joint kernel is their concatenation. Returns (node bbs, aggregated
    node weights, per-wavelength kernel sums, vacuum weight).
    """
    table = ray_table(variant)
    members = hst_members(dataset)
    total = sum(member["texp_s"] for member in members)
    node_bbs, node_weights, sums, vacuum = [], [], [], 0.0
    for member in members:
        radius = float(member["shadow_axis_distance_m"])
        geometry = Geometry(moon_distance_m=member["receiver_axial_distance_m"],
                            sun_distance_m=member["sun_site_distance_m"])
        bb, area = receiver_quadrature(table, radius, geometry)
        weights = kernel_weights(table, bb, area, wavelengths, extra_aerosol)
        outer = vacuum_branch(table, radius, geometry)
        share = member["texp_s"]/total
        node_bbs.append(bb)
        # Equal-weight mean over the quantitative window grid (equidistant
        # integral convention; window spacing is uniform to <=0.026 nm).
        node_weights.append(weights.mean(axis=1)*share)
        sums.append(np.sum(weights, axis=0)*share + outer*share)
        vacuum += outer*share
    return node_bbs, node_weights, np.sum(sums, axis=0), vacuum


def hst_height_records() -> list:
    comparison = hst_comparison()
    solar_scale = comparison["solar_scale_in_over_out"]
    wavelengths, window_rows = hst_window_data()
    records = []
    for dataset in comparison["grouping"]["umbral_in_eclipse"]:
        for variant in VARIANTS:
            node_bbs, node_weights, sums, vacuum = hst_kernel(variant, dataset, wavelengths)
            scaled = sums*solar_scale[dataset]
            committed = np.array([float(row[f"{dataset}_{variant}_prediction"])
                                  for row in window_rows])
            moment = float(scaled.mean())
            records.append({"event": "hst_2019", "readout": dataset,
                            "wavelength_window_nm": [UMBRAL_QUANTITATIVE_MIN_NM, G430L_MAX_NM],
                            "variant": variant,
                            **altitude_statistics(ray_table(variant), node_bbs, node_weights, vacuum),
                            **moment_gap(moment, float(committed.mean())),
                            "per_wavelength_max_relative_gap": float(np.max(np.abs(scaled-committed)/committed))})
    return records


@lru_cache(maxsize=None)
def dec_comparison() -> dict:
    return json.loads((RESULTS/"dec2011/comparison.json").read_text())


@lru_cache(maxsize=None)
def dec_time_weights() -> np.ndarray:
    """Simpson exposure-average weights, reusing the committed aggregator."""
    records = dec_comparison()["geometry"]["records"]
    seconds = np.array([record["seconds_since_start"] for record in records], dtype=float)
    return exposure_average(np.eye(len(records)), seconds)


def dec_kernel(variant: str, wavelength_nm: float, extra_aerosol=None) -> tuple:
    """One December 2011 readout kernel over the 31-sample U1 exposure."""
    records = dec_comparison()["geometry"]["records"]
    time_weights = dec_time_weights()
    table = ray_table(variant)
    node_bbs, node_weights, vacuum, moment = [], [], 0.0, 0.0
    for record, time_weight in zip(records, time_weights):
        geometry = Geometry(moon_distance_m=record["receiver_axial_distance_m"],
                            sun_distance_m=record["sun_earth_distance_m"])
        radius = float(record["shadow_radius_m"])
        bb, area = receiver_quadrature(table, radius, geometry)
        weights = kernel_weights(table, bb, area, wavelength_nm, extra_aerosol)
        outer = vacuum_branch(table, radius, geometry)
        share = float(time_weight*record["solar_flux_relative_to_b2"])
        node_bbs.append(bb)
        node_weights.append(weights*share)
        vacuum += outer*share
        moment += float(np.sum(weights))*share + outer*share
    return node_bbs, node_weights, vacuum, moment


def dec_height_records() -> list:
    records = []
    for record in dec_comparison()["records"]:
        for variant in VARIANTS:
            node_bbs, node_weights, vacuum, moment = dec_kernel(variant, record["wavelength_nm"])
            records.append({"event": "dec_2011", "readout": f'{record["wavelength_nm"]:g}nm',
                            "wavelength_nm": record["wavelength_nm"],
                            "variant": variant,
                            **altitude_statistics(ray_table(variant), node_bbs, node_weights, vacuum),
                            **moment_gap(moment, record[f"{variant}_prediction"])})
    return records


@lru_cache(maxsize=None)
def mv1990_analysis() -> dict:
    return json.loads((RESULTS/"multiple_events/analysis.json").read_text())


class BrokenPowerLawSpectral:
    """Validation-side broken power law allowing negative Angstrom exponents.

    ``eclipse.aerosol.AerosolProfile`` validates its exponents nonnegative;
    the OPAC mineral-transported type has negative tabulated exponents
    (extinction decreasing with wavelength). This class mirrors that class's
    ``spectral_factor`` functional form (500 nm break, 550 nm normalization)
    without touching the production module.
    """

    def __init__(self, angstrom_blue: float, angstrom_visible: float,
                 break_nm: float = 500.0):
        self.angstrom_blue = float(angstrom_blue)
        self.angstrom_visible = float(angstrom_visible)
        self.break_nm = float(break_nm)

    def spectral_factor(self, wavelength_nm) -> np.ndarray:
        wavelength = np.asarray(wavelength_nm, dtype=float)
        exponent = np.where(wavelength >= self.break_nm,
                            self.angstrom_visible, self.angstrom_blue)
        raw = (wavelength/self.break_nm)**(-exponent)
        reference_exponent = (self.angstrom_visible if 550.0 >= self.break_nm
                              else self.angstrom_blue)
        reference = (550.0/self.break_nm)**(-reference_exponent)
        return raw/reference


def mv1990_record(variant: str) -> dict:
    return next(record for record in mv1990_analysis()["records"]
                if record["event_date"] == MV1990_EVENT
                and record["response"] == MV1990_RESPONSE
                and record["variant"] == variant)


def mv1990_kernel(variant: str, extra_aerosol=None) -> tuple:
    """Whole-lunar-disk photopic kernel (nominal_v, order 256 / 2.5 nm grid).

    The committed refined prediction's aperture quadrature, wavelength grid
    and response are reused exactly; the kernel decomposes the disk-weighted
    mean spectrum into per-impact nodes with photopic wavelength aggregation.
    """
    analysis = mv1990_analysis()
    record = mv1990_record(variant)
    event = next(item for item in analysis["geometry"]["events"]
                 if item["event_date"] == MV1990_EVENT)
    geometry = Geometry(moon_distance_m=event["moon_distance_m"],
                        sun_distance_m=event["sun_distance_m"])
    wavelengths = np.arange(380., 780.+MV1990_WAVELENGTH_STEP_NM/2, MV1990_WAVELENGTH_STEP_NM)
    response = photopic_response(wavelengths)
    response_norm = np.trapezoid(response, wavelengths)
    samples, area_samples = aperture_quadrature(record["shadow_center_offset_m"],
                                                record["aperture_radius_m"], MV1990_APERTURE_ORDER)
    table = ray_table(variant)
    node_bbs, node_weights, vacuum, moment = [], [], 0.0, 0.0
    for sample_radius, sample_area in zip(samples, area_samples):
        radius = float(sample_radius)
        bb, area = receiver_quadrature(table, radius, geometry)
        weights = kernel_weights(table, bb, area, wavelengths, extra_aerosol)
        photopic = np.trapezoid(weights*response[None, :], wavelengths, axis=1)/response_norm
        outer = vacuum_branch(table, radius, geometry)
        node_bbs.append(bb)
        node_weights.append(photopic*sample_area)
        vacuum += outer*sample_area
        moment += float(photopic.sum())*sample_area + outer*sample_area
    return node_bbs, node_weights, vacuum, moment


def mv1990_height_records() -> list:
    records = []
    for variant in VARIANTS:
        node_bbs, node_weights, vacuum, moment = mv1990_kernel(variant)
        committed = mv1990_record(variant)["prediction"]
        records.append({"event": "event_1990_02_09", "readout": "whole_disk_nominal_v",
                        "variant": variant,
                        **altitude_statistics(ray_table(variant), node_bbs, node_weights, vacuum),
                        **moment_gap(moment, committed)})
    return records


def height_records() -> list:
    """All four apertures' per-readout altitude statistics and moment gates."""
    return (june_height_records() + hst_height_records()
            + dec_height_records() + mv1990_height_records())


EVENT_ORDER = ("june_2011", "hst_2019", "dec_2011", "event_1990_02_09")


def committed_residuals(variant: str) -> dict:
    """Committed per-readout model/observed values of one variant, per event.

    Read from the committed products only; no model recomputation (plan
    Step 3). Calibrated variants (red_calibrated etc.) never enter.
    """
    return {
        "june_2011": [record[f"{variant}_model_to_observed"] for record in june_records()],
        "hst_2019": {dataset: [float(row[f"{dataset}_{variant}_model_to_observed"])
                               for row in hst_window_data()[1]]
                     for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]},
        "dec_2011": [record[f"{variant}_model_to_observed"]
                     for record in dec_comparison()["records"]],
        "event_1990_02_09": [mv1990_record(variant)["model_to_observed"]],
    }


@lru_cache(maxsize=None)
def june_shape_ratios() -> dict:
    """503/677 model-over-observed transmission ratios per paired angle.

    June 2011 shape-axis metric (plan design point 4), computed from the
    committed round2 records at the five angles that carry both bands.
    """
    ratios = {variant: [] for variant in VARIANTS}
    for angle in sorted({record["shadow_angle_deg"] for record in june_records()}):
        bands = {record["wavelength_nm"]: record for record in june_records()
                 if record["shadow_angle_deg"] == angle}
        if not {503.0, 677.0} <= bands.keys():
            continue
        blue, red = bands[503.0], bands[677.0]
        observed_ratio = blue["observed_curve_brightness"]/red["observed_curve_brightness"]
        for variant in VARIANTS:
            ratios[variant].append({"shadow_angle_deg": angle,
                                    "ratio": (blue[f"{variant}_prediction"]
                                              /red[f"{variant}_prediction"])/observed_ratio})
    return ratios


def shape_axis(event: str, variant: str):
    """Committed shape metrics of one event x variant (plan design point 4)."""
    if event == "june_2011":
        ratios = june_shape_ratios()[variant]
        return {"metric": "503_over_677_transmission_ratio_model_over_observed",
                "median_over_paired_angles": float(np.median(
                    [item["ratio"] for item in ratios])),
                "paired_angle_count": len(ratios), "values_by_angle": ratios}
    if event == "hst_2019":
        metrics = hst_comparison()["metrics"][variant]
        return {dataset: metrics[dataset]["slope_model_vs_observed"]
                for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]}
    if event == "dec_2011":
        diagnostic = dec_comparison()["descriptive_diagnostics"][variant]
        return {key: diagnostic[key] for key in
                ["450_to750_transmission_ratio_model_over_observed",
                 "500_to750_transmission_ratio_model_over_observed"]}
    return None


def depth_statistics(height_stats: dict) -> dict:
    """Per event x variant medians of the per-readout path-depth statistics."""
    grouped = {}
    for record in height_stats["records"]:
        grouped.setdefault((record["event"], record["variant"]), []).append(record)
    result = {}
    for (event, variant), rows in grouped.items():
        result.setdefault(event, {})[variant] = {
            "p50_height_m_median": float(np.median(
                [row["p10_p50_p90_height_m"][1] for row in rows])),
            "band_fraction_below_2km_median": float(np.median(
                [row["band_fraction_below_2km"] for row in rows])),
            "band_fraction_2_to_12km_median": float(np.median(
                [row["band_fraction_2_to_12km"] for row in rows])),
            "band_fraction_over_12km_median": float(np.median(
                [row["band_fraction_over_12km"] for row in rows])),
            "depth_readout_count": len(rows),
        }
    return result


def pairwise_direction(depth_p50, depth_below2km, residual_ln):
    """Pairwise direction consistency of the four events (no significance claims)."""
    events = EVENT_ORDER
    pairs = []
    for i in range(len(events)):
        for j in range(i+1, len(events)):
            deeper = events[i] if depth_p50[i] < depth_p50[j] else events[j]
            larger = events[i] if residual_ln[i] > residual_ln[j] else events[j]
            pairs.append({"pair": [events[i], events[j]],
                          "p50_height_m_median": [depth_p50[i], depth_p50[j]],
                          "residual_ln_ratio_median": [residual_ln[i], residual_ln[j]],
                          "deeper_by_p50": deeper, "larger_residual": larger,
                          "deeper_has_larger_residual": bool(deeper == larger),
                          "more_below2km_has_larger_residual":
                              bool((depth_below2km[i] < depth_below2km[j])
                                   == (residual_ln[i] < residual_ln[j]))})
    return {"pairs_total": len(pairs), "pairs": pairs,
            "deeper_has_larger_residual_count":
                int(sum(pair["deeper_has_larger_residual"] for pair in pairs)),
            "more_below2km_has_larger_residual_count":
                int(sum(pair["more_below2km_has_larger_residual"] for pair in pairs))}


def build_correlation(height_stats: dict) -> dict:
    """Cross-event correlation of committed residuals vs contribution depth."""
    residuals = {variant: committed_residuals(variant) for variant in VARIANTS}
    depths = depth_statistics(height_stats)

    # HST cross-check: pooled per-dataset window medians must reproduce the
    # committed metric field exactly (plan Step 3; fail loud on any drift).
    hst_check = {}
    for variant in VARIANTS:
        for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]:
            pooled = float(np.median(residuals[variant]["hst_2019"][dataset]))
            committed = hst_comparison()["metrics"][variant][dataset]["model_to_observed_median"]
            if pooled != committed:
                raise RuntimeError(f"HST consistency check failed for {dataset}/{variant}: "
                                   f"pooled window median {pooled!r} != committed "
                                   f"{committed!r}")
            hst_check[f"{dataset}|{variant}"] = {
                "pooled_window_median": pooled,
                "committed_model_to_observed_median": committed}

    def flat(values):
        if isinstance(values, dict):
            return [item for sublist in values.values() for item in sublist]
        return list(values)

    event_level = {}
    for event in EVENT_ORDER:
        event_level[event] = {}
        for variant in VARIANTS:
            ratios = np.asarray(flat(residuals[variant][event]), dtype=float)
            event_level[event][variant] = {
                "residual_ln_ratio_median": float(np.median(np.log(ratios))),
                "residual_ratio_median": float(np.median(ratios)),
                "residual_readout_count": len(ratios),
                **depths[event][variant],
                "shape_axis": shape_axis(event, variant),
            }

    cross_event = {}
    for variant in VARIANTS:
        p50 = np.array([event_level[event][variant]["p50_height_m_median"]
                        for event in EVENT_ORDER])
        below = np.array([event_level[event][variant]["band_fraction_below_2km_median"]
                          for event in EVENT_ORDER])
        residual = np.array([event_level[event][variant]["residual_ln_ratio_median"]
                             for event in EVENT_ORDER])
        cross_event[variant] = {
            "spearman_rho_residual_ln_vs_p50_median": float(spearmanr(p50, residual)[0]),
            "spearman_rho_residual_ln_vs_below2km_median": float(spearmanr(below, residual)[0]),
            "pairwise_direction": pairwise_direction(p50, below, residual),
        }

    return {
        "scope": "Cross-event correlation of committed residual magnitudes against "
                 "contribution-path depth statistics (height_stats.json records). Brightness "
                 "and shape axes are kept in separate fields. Residual = event-level median of "
                 "ln(model/observed) over the aperture's committed quantitative readouts, with "
                 "ratio medians recorded alongside; only the two committed model variants enter "
                 "(calibrated variants excluded, discipline). N=4 events per variant: rank "
                 "orderings only, no significance claims. Negative rho(p50, residual) means "
                 "deeper contribution paths carry larger model/observed residuals (the issue "
                 "hypothesis direction); positive rho(below-2 km fraction, residual) is the "
                 "boundary-layer-relevant depth reading.",
        "residual_sources": {str(path.relative_to(ROOT)): {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in [RESULTS/"round2/analysis.json", RESULTS/"hst2019/comparison.json",
                         RESULTS/"hst2019/comparison.csv", RESULTS/"dec2011/comparison.json",
                         RESULTS/"multiple_events/analysis.json"]},
        "hst_consistency_check": hst_check,
        "event_level": event_level,
        "cross_event": cross_event,
    }


# --- Discrimination probes (plan Step 5) -------------------------------------

AOD_ROW_BY_CANDIDATE_KIND = {
    "boundary_layer": "aod550_boundary_layer_only",
    "layer_total": "tau_0.55um_layer",
    "background": "aod550",
}
# How each candidate kind stacks on the committed fixed_opac baseline.
STACKING_SEMANTICS = {
    "boundary_layer": "Adds the candidate's 0-2 km (desert 0-6 km, antarctic "
                      "0-10 km) layer on top of the committed fixed_opac free "
                      "troposphere; the 2 km types touch only kernel mass below 2 km.",
    "layer_total": "Adds the 2-3.5 km layer on top of the committed fixed_opac "
                   "free troposphere, whose 2-3.5 km band is therefore double "
                   "counted; this matches the OPAC construction where mineral "
                   "transported sits above the boundary layer.",
    "background": "Adds the 12-35 km layer on top of the committed fixed_opac "
                  "free troposphere; no overlap with existing components.",
}
# Prior-driven counterfactual probes: the AOD is the literature value, never
# tuned; verdicts compare the prior-order effect with the committed residual.
VERDICT_COMPATIBLE = (0.5, 2.0)
VERDICT_TENSION = (0.1, 10.0)


@lru_cache(maxsize=None)
def prior_rows() -> tuple:
    return tuple(json.loads(
        (RESULTS/"extinction_gap_diagnosis/prior_parameters.json").read_text())["parameters"])


def prior_row(candidate: str, name: str) -> dict:
    return next(row for row in prior_rows()
                if row["candidate"] == candidate and row["name"] == name)


def candidate_kind(candidate: str) -> str:
    if candidate == "stratosphere_background":
        return "background"
    if candidate == "mineral_transported":
        return "layer_total"
    return "boundary_layer"


@lru_cache(maxsize=None)
def candidate_probe_extra(candidate: str) -> tuple:
    """Extra extinction column of one literature candidate, on its own grid.

    The tau550 column depends only on the layer geometry (aod, bounds, scale
    height), so it is traced by the production ``build_ray_table`` with a
    profile carrying the candidate's geometry; the spectral factor carries
    the candidate's Angstrom exponents, via ``AerosolProfile`` itself when
    they are nonnegative and via the validation-side broken power law when
    they are negative (mineral transported). A zero-AOD candidate restores
    the baseline exactly (tested).
    """
    kind = candidate_kind(candidate)
    aod = prior_row(candidate, AOD_ROW_BY_CANDIDATE_KIND[kind])["value"]
    bottom = prior_row(candidate, "layer_bottom_m")["value"]
    top = prior_row(candidate, "layer_top_m")["value"]
    scale_height = prior_row(candidate, "scale_height_m")["value"]
    blue = prior_row(candidate, "angstrom_0.35_to_0.5um")["value"]
    visible = prior_row(candidate, "angstrom_0.5_to_0.8um")["value"]
    if blue >= 0 and visible >= 0:
        profile = AerosolProfile(aod550=aod, bottom_m=bottom, top_m=top,
                                 scale_height_m=scale_height,
                                 angstrom_blue=blue, angstrom_visible=visible)
    else:
        profile = BrokenPowerLawSpectral(blue, visible)
    column_profile = AerosolProfile(aod550=aod, bottom_m=bottom, top_m=top,
                                    scale_height_m=scale_height)
    table = build_ray_table(Atmosphere(), aerosol_profile=column_profile)
    return table.impact_m, table.aerosol_tau550, profile


def closure_verdict(closure: float) -> str:
    """Prior-order compatibility: the probe closes `closure` of the log-gap."""
    if VERDICT_COMPATIBLE[0] <= closure <= VERDICT_COMPATIBLE[1]:
        return "compatible"
    if VERDICT_TENSION[0] <= closure <= VERDICT_TENSION[1]:
        return "tension"
    return "excluded_leaning"


def shape_closed_fraction(base: float, new: float) -> float:
    """Fraction of the shape disagreement |base-1| closed by the probe."""
    disagreement = abs(base - 1.)
    if disagreement <= 0:
        return 0.0
    return (disagreement - abs(new - 1.))/disagreement


def probe_brightness_factor(candidate: str) -> dict:
    """Per-readout transmission factors of one candidate against fixed_opac.

    The baseline is the committed prediction itself (zeroth-moment-gated);
    the factor is probe kernel sum over committed prediction per readout.
    """
    extra = candidate_probe_extra(candidate)
    geometry = june_geometry()
    june = []
    for record in june_records():
        radius = geometry.moon_distance_m*np.tan(np.deg2rad(record["shadow_angle_deg"]))
        _, _, _, probed = june_kernel("fixed_opac", radius, record["wavelength_nm"], extra)
        june.append(probed/record["fixed_opac_prediction"])
    wavelengths, window_rows = hst_window_data()
    hst = {}
    for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]:
        committed = np.array([float(row[f"{dataset}_fixed_opac_prediction"])
                              for row in window_rows])
        _, _, sums, _ = hst_kernel("fixed_opac", dataset, wavelengths, extra)
        hst[dataset] = (sums*hst_comparison()["solar_scale_in_over_out"][dataset])/committed
    dec = []
    for record in dec_comparison()["records"]:
        _, _, _, probed = dec_kernel("fixed_opac", record["wavelength_nm"], extra)
        dec.append(probed/record["fixed_opac_prediction"])
    _, _, _, probed = mv1990_kernel("fixed_opac", extra)
    mv1990 = probed/mv1990_record("fixed_opac")["prediction"]
    return {"june_2011": june, "hst_2019": hst, "dec_2011": dec,
            "event_1990_02_09": mv1990}


def probe_event_brightness(candidate: str, event: str, factors: dict,
                           depths: dict) -> dict:
    """Brightness-axis probe result for one candidate x event."""
    base_rows = {variant: depths[event][variant] for variant in VARIANTS}
    fixed = base_rows["fixed_opac"]
    if event == "june_2011":
        ratios = np.array([record["fixed_opac_model_to_observed"]
                           for record in june_records()])
        factor = np.array(factors["june_2011"])
    elif event == "hst_2019":
        ratios = np.concatenate([
            np.array([float(row[f"{dataset}_fixed_opac_model_to_observed"])
                      for row in hst_window_data()[1]])
            for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]])
        factor = np.concatenate([np.asarray(values) for values
                                 in factors["hst_2019"].values()])
    elif event == "dec_2011":
        ratios = np.array([record["fixed_opac_model_to_observed"]
                           for record in dec_comparison()["records"]])
        factor = np.array(factors["dec_2011"])
    else:
        ratios = np.array([mv1990_record("fixed_opac")["model_to_observed"]])
        factor = np.array([factors["event_1990_02_09"]])
    new_ratios = ratios*factor
    base_ln = float(np.median(np.log(ratios)))
    new_ln = float(np.median(np.log(new_ratios)))
    delta_ln = new_ln - base_ln
    closure = (base_ln - new_ln)/base_ln if base_ln else 0.0
    per_readout = -np.log(factor)/np.log(ratios)
    return {"base_residual_ln_median": base_ln,
            "probe_delta_ln_median": delta_ln,
            "new_residual_ln_median": new_ln,
            "gap_closure_fraction": closure,
            "verdict": closure_verdict(closure),
            "max_readout_gap_closure_fraction": float(np.max(per_readout)),
            "base_p50_height_m_median": fixed["p50_height_m_median"],
            "base_band_fraction_below_2km_median": fixed["band_fraction_below_2km_median"],
            "base_band_fraction_over_12km_median": fixed["band_fraction_over_12km_median"]}


def hst_probe_slopes(extra_aerosol=None) -> dict:
    """HST slopes of the probe prediction on the full window grid.

    The aggregation structure is the committed ``compare_hst2019`` one (band
    medians over the full window grid, ratio of ratios); with
    ``extra_aerosol=None`` it reproduces the anchored baselines exactly.
    """
    wavelengths, rows = hst_grid_data()
    slopes = {}
    for dataset in hst_comparison()["grouping"]["umbral_in_eclipse"]:
        observed = np.array([float(row[f"{dataset}_obs_ratio"]) for row in rows])
        _, _, sums, _ = hst_kernel("fixed_opac", dataset, wavelengths, extra_aerosol)
        new_pred = sums*hst_comparison()["solar_scale_in_over_out"][dataset]
        slopes[dataset] = {f"{num:g}_over_{den:g}":
                           (band_median(new_pred, wavelengths, num)
                            /band_median(new_pred, wavelengths, den))
                           /(band_median(observed, wavelengths, num)
                             /band_median(observed, wavelengths, den))
                           for num, den in SLOPE_BANDS_NM}
    return slopes


def probe_event_shape(candidate: str, event: str, factors: dict) -> dict:
    """Shape-axis probe deltas for one candidate x event (fixed_opac).

    A probe closes shape disagreement when the new metric moves toward unity
    (perfect observed shape agreement); the closed fraction can go negative
    when the probe pushes the metric further away.
    """
    if event == "june_2011":
        factor_by_key = {(record["wavelength_nm"], record["shadow_angle_deg"]): factor
                         for record, factor in zip(june_records(), factors["june_2011"])}
        metrics = []
        for item in shape_axis("june_2011", "fixed_opac")["values_by_angle"]:
            angle = item["shadow_angle_deg"]
            shape_factor = factor_by_key[503.0, angle]/factor_by_key[677.0, angle]
            new = item["ratio"]*shape_factor
            closed = shape_closed_fraction(item["ratio"], new)
            metrics.append({"shadow_angle_deg": angle, "base_ratio": item["ratio"],
                            "shape_factor": shape_factor, "new_ratio": new,
                            "disagreement_closed_fraction": closed,
                            "verdict": closure_verdict(closed)})
        median_base = float(np.median([metric["base_ratio"] for metric in metrics]))
        median_new = float(np.median([metric["new_ratio"] for metric in metrics]))
        median_closed = shape_closed_fraction(median_base, median_new)
        return {"metric": "503_over_677_transmission_ratio_model_over_observed",
                "median_base": median_base, "median_new": median_new,
                "disagreement_closed_fraction": median_closed,
                "verdict": closure_verdict(median_closed), "by_angle": metrics}
    if event == "hst_2019":
        # The committed slope aggregates band medians over the FULL window
        # grid (the 430 nm band includes [425, 430) rows), so the probe
        # evaluates its kernel on that same grid and takes its baseline from
        # the anchored hst_slope_baseline() - not from the quantitative
        # window subset (review round 1 Major#1).
        baselines = hst_slope_baseline()["fixed_opac"]
        probe_slopes = hst_probe_slopes(candidate_probe_extra(candidate))
        metrics = {}
        for dataset, slope_values in probe_slopes.items():
            metrics[dataset] = {}
            for band, new_slope in slope_values.items():
                base_slope = baselines[dataset][band]
                closed = shape_closed_fraction(base_slope, new_slope)
                metrics[dataset][band] = {
                    "base_slope": base_slope, "new_slope": new_slope,
                    "disagreement_closed_fraction": closed,
                    "verdict": closure_verdict(closed)}
        return metrics
    if event == "dec_2011":
        records = dec_comparison()["records"]
        index = {record["wavelength_nm"]: position for position, record in enumerate(records)}
        base_axis = shape_axis("dec_2011", "fixed_opac")
        metrics = {}
        for key, blue_index in [("450_to750_transmission_ratio_model_over_observed", index[450.0]),
                                ("500_to750_transmission_ratio_model_over_observed", index[500.0])]:
            shape_factor = factors["dec_2011"][blue_index]/factors["dec_2011"][index[750.0]]
            base = base_axis[key]
            new = base*shape_factor
            closed = shape_closed_fraction(base, new)
            metrics[key] = {"base_ratio": base, "shape_factor": shape_factor,
                            "new_ratio": new,
                            "disagreement_closed_fraction": closed,
                            "verdict": closure_verdict(closed)}
        return metrics
    return {}


def build_discrimination(height_stats: dict) -> dict:
    """Prior-driven counterfactual discrimination across the four apertures."""
    depths = depth_statistics(height_stats)
    band_fractions = {event: {
        "band_fraction_below_2km_median": depths[event]["fixed_opac"]["band_fraction_below_2km_median"],
        "band_fraction_2_to_12km_median": depths[event]["fixed_opac"]["band_fraction_2_to_12km_median"],
        "band_fraction_over_12km_median": depths[event]["fixed_opac"]["band_fraction_over_12km_median"]}
        for event in EVENT_ORDER}
    order = []
    for row in prior_rows():
        if row["candidate"] not in order:
            order.append(row["candidate"])
    candidates = []
    for candidate in order:
        if candidate == "free_troposphere_reference":
            continue
        kind = candidate_kind(candidate)
        aod_row = prior_row(candidate, AOD_ROW_BY_CANDIDATE_KIND[kind])
        bottom = prior_row(candidate, "layer_bottom_m")["value"]
        top = prior_row(candidate, "layer_top_m")["value"]
        factors = probe_brightness_factor(candidate)
        events = {}
        for event in EVENT_ORDER:
            events[event] = {
                "brightness": probe_event_brightness(candidate, event, factors, depths),
                "shape": probe_event_shape(candidate, event, factors)}
        closures = [events[event]["brightness"]["gap_closure_fraction"] for event in EVENT_ORDER]
        verdict_counts = {}
        for event in EVENT_ORDER:
            verdict = events[event]["brightness"]["verdict"]
            verdict_counts[verdict] = verdict_counts.get(verdict, 0) + 1
        candidate_record = {
            "candidate": candidate,
            "action_band_m": [bottom, top],
            "probe_aod550": aod_row["value"],
            "probe_aod550_source_location": aod_row["source"]["location"],
            "stacking_semantics": STACKING_SEMANTICS[kind],
            "events": events,
            "brightness_verdict_counts": verdict_counts,
            "min_gap_closure_fraction": min(closures),
            "max_readout_gap_closure_fraction": max(
                events[event]["brightness"]["max_readout_gap_closure_fraction"]
                for event in EVENT_ORDER),
            # Zero-effect candidates close ~1e-16 (float noise around an
            # exactly-zero effect); clamped to 0 so the geometric mean of a
            # null effect is 0, not nan.
            "geometric_mean_gap_closure_fraction": float(np.exp(np.mean(
                np.log([max(closure, 0.) for closure in closures])))),
        }
        if candidate == "antarctic":
            # Interpretation pinned in the product (a08), not only in docs:
            # the antarctic type's 0-10 km layer double counts 2-10 km
            # extinction already carried by the fixed_opac free troposphere.
            candidate_record["interpretation"] = (
                "The antarctic type (0-10 km, weak spectral dependence, AOD "
                "about 5x the free-troposphere prior) reads on this baseline "
                "as a proxy for the 'more same-type extinction across 2-10 km' "
                "hypothesis: it double counts the 2-10 km band already carried "
                "by the fixed_opac free troposphere, so its brightness "
                "compatibility transfers to any low-layer scaled increment, "
                "and so does its shape cost.")
        candidates.append(candidate_record)
    ranked = sorted(candidates, key=lambda item: -item["min_gap_closure_fraction"])
    return {
        "scope": "Prior-driven counterfactual discrimination probes for the missing-extinction "
                 "candidates (plan Step 5). Each candidate's literature AOD550 enters as an extra "
                 "extinction column (or, for negative-Angstrom mineral transported, a "
                 "validation-side spectral factor) on top of the committed fixed_opac default, "
                 "evaluated through the same kernel adapters as the height statistics; a zero-AOD "
                 "probe restores the baseline exactly (pytest). Fixed priors only: nothing is "
                 "fitted, nothing becomes a default, and no causal attribution is made - the "
                 "product ranks prior-order compatibility only, per docs/ROUND2.md discipline.",        "variant_scope": "Probes run on the fixed_opac variant only: the candidates are missing "
                         "components on top of the committed default physics; probing the "
                         "molecular variant would conflate the free-troposphere candidate with "
                         "the boundary-layer ones.",
        "verdict_bands": {
            "brightness_semantics": "gap_closure_fraction = fraction of the event's committed "
                                    "log-residual closed by the probe (1 = exactly closes). "
                                    "Premise: every committed event residual is positive (model "
                                    "brighter than observed across all four apertures, verified in "
                                    "correlation.json); a future negative residual would invert the "
                                    "closure sign semantics and require re-deriving the bands.",
            "shape_semantics": "disagreement_closed_fraction = movement of the shape metric "
                               "toward unity; negative = pushed further away",
            "compatible": "closure in [0.5, 2]: prior-order effect of the residual's order",
            "tension": "closure in [0.1, 0.5) or (2, 10]: right direction, magnitude off",
            "excluded_leaning": "closure < 0.1 (negligible at prior level) or > 10 (dominant at "
                                "prior level); for shape, < 0.1 also covers worsening",
            "ranking_key": "brightness-axis minimum gap closure across the four events (weakest "
                           "link), descending; shape verdicts are recorded per event alongside "
                           "and enter the owner's reading, not this mechanical key",
            "gloss": {"compatible": "相容", "tension": "张力", "excluded_leaning": "排除倾向"},
        },
        "band_fractions_fixed_opac": band_fractions,
        # Keyed by candidate name so product citations resolve as
        # candidates[<name>]; ranking carries the mechanical order.
        "candidates": {item["candidate"]: item for item in candidates},
        "not_probed": [{
            "candidate": "free_troposphere_humidified",
            "reason": "The free troposphere is tabulated at 50% RH only and humidified scaling "
                      "needs the per-component OPAC data files, marked unavailable in "
                      "prior_parameters.json (row free_troposphere_reference/aod550_at_other_rh "
                      "and each boundary-layer type's aod550_at_rh_50_70_90 row carry the "
                      "unavailable status and the search trail); no sourced prior AOD exists, "
                      "so no probe runs and the candidate carries no verdict.",
            "direction_support": "The 677 nm same-event calibration direction supports 'more "
                                 "same-type extinction' (docs/ROUND2.md), but that calibrated "
                                 "value is not a prior and is excluded from this ranking.",
        }],
        "ranking": [item["candidate"] for item in ranked],
    }


def write_discrimination(out: Path, stats: dict) -> dict:
    discrimination = build_discrimination(stats)
    (out/"discrimination.json").write_text(json.dumps(discrimination, indent=2)+"\n")
    return discrimination


@lru_cache(maxsize=1)
def build_height_stats() -> dict:
    records = height_records()
    worst = max(record["relative_gap"] for record in records)
    if worst > MOMENT_TOLERANCE:
        raise RuntimeError(f"Zeroth-moment gate failed: worst relative gap {worst:.3e} "
                           f"exceeds {MOMENT_TOLERANCE:.0e}")
    return {
        "scope": "Contribution-path altitude statistics of the same weighted mixed kernel that "
                 "produced each aperture's committed prediction (kernel identity = plan design "
                 "point 1). P10/P50/P90 and band fractions describe the atmospheric branch only; "
                 "the vacuum branch (tangent above table top, transmission one) is recorded "
                 "separately and enters the zeroth-moment identity only. Single authoritative "
                 "implementation of four-aperture height statistics; round2 altitude_quantiles is "
                 "a historical product generator and its committed June values are the anchor.",
        "summary": {"record_count": len(records), "worst_relative_gap": worst},
        "kernel_identity": {
            "zeroth_moment_definition": "sum(weights) + vacuum branch == committed prediction",
            "relative_tolerance": MOMENT_TOLERANCE,
            "tolerance_semantics": "structural reuse makes equality hold by construction; the "
                                   "budget covers float summation-order differences only. Gate "
                                   "failure routing: aperture/kernel-identity errors fix the "
                                   "kernel; pure reproduction differences fix the tolerance "
                                   "definition with written justification, never silently.",
            "band_edges_m": [BAND_BELOW_TOP_M, BAND_OVER_TOP_M],
            "band_semantics": "below 2 km: no aerosol component in the model; 2-12 km: OPAC free "
                              "troposphere; over 12 km: no component (molecular above; vacuum "
                              "above the 120 km table top reported separately)",
        },
        "committed_sources": {str(path.relative_to(ROOT)): {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in [RESULTS/"round2/analysis.json", RESULTS/"hst2019/comparison.json",
                         RESULTS/"hst2019/comparison.csv", RESULTS/"dec2011/comparison.json",
                         RESULTS/"multiple_events/analysis.json",
                         OBSERVATIONS/"ugolnikov2011/geometry.json",
                         OBSERVATIONS/"hst2019/geometry.json"]},
        "production_module_sha256": {
            name: hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
            for name in ["aerosol.py", "atmosphere.py", "rays.py", "transport.py", "spectra.py"]},
        "records": records,
    }


def write_height_stats(out: Path) -> dict:
    stats = build_height_stats()
    out.mkdir(parents=True, exist_ok=True)
    (out/"height_stats.json").write_text(json.dumps(stats, indent=2)+"\n")
    fields = ["event", "readout", "variant", "wavelength_nm", "shadow_angle_deg",
              "p10_p50_p90_height_m", "band_fraction_below_2km", "band_fraction_2_to_12km",
              "band_fraction_over_12km", "vacuum_branch_fraction", "kernel_zeroth_moment",
              "committed_prediction", "relative_gap", "per_wavelength_max_relative_gap"]
    with (out/"height_stats.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(fields)
        for record in stats["records"]:
            writer.writerow([json.dumps(record[field]) if isinstance(record.get(field), list)
                             else record.get(field, "") for field in fields])
    return stats


def write_correlation(out: Path, stats: dict) -> dict:
    correlation = build_correlation(stats)
    (out/"correlation.json").write_text(json.dumps(correlation, indent=2)+"\n")
    return correlation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("outputs/extinction_gap_diagnosis"))
    args = parser.parse_args()
    stats = write_height_stats(args.out)
    correlation = write_correlation(args.out, stats)
    discrimination = write_discrimination(args.out, stats)
    print(json.dumps({"records": len(stats["records"]),
                      "max_relative_gap": max(record["relative_gap"] for record in stats["records"]),
                      "cross_event": {variant: {key: value for key, value in summary.items()
                                                if key != "pairwise_direction"}
                                      for variant, summary in correlation["cross_event"].items()},
                      "ranking": discrimination["ranking"]},
                     indent=2))


if __name__ == "__main__":
    main()
