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

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.spectra import ozone_cross_section, rayleigh_cross_section
from eclipse.transport import Geometry, circle_overlap, receiver_quadrature

from .compare_dec2011 import exposure_average
from .compare_hst2019 import G430L_MAX_NM, UMBRAL_QUANTITATIVE_MIN_NM
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
def hst_window_data() -> tuple:
    """Committed quantitative-window wavelengths and comparison.csv rows."""
    with (RESULTS/"hst2019/comparison.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    wavelengths = np.array([float(row["wavelength_nm"]) for row in rows])
    window = (wavelengths >= UMBRAL_QUANTITATIVE_MIN_NM) & (wavelengths <= G430L_MAX_NM)
    return wavelengths[window], [row for row, keep in zip(rows, window) if keep]


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("outputs/extinction_gap_diagnosis"))
    args = parser.parse_args()
    stats = write_height_stats(args.out)
    print(json.dumps({"records": len(stats["records"]),
                      "max_relative_gap": max(record["relative_gap"] for record in stats["records"])},
                     indent=2))


if __name__ == "__main__":
    main()
