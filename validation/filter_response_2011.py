"""Bandpass-integrated versus monochromatic comparison for the June 2011 readouts.

The script first rebuilds the monochromatic columns exactly as committed
(validation/results/round2/analysis.json and
validation/results/observations2011/comparison.json) and cross-checks them
point by point before any bandpass product is computed. Filter information
comes solely from validation/observations/ugolnikov2011/source.json
(filter_information, branch C: no width or response stated in the paper).
Rectangular and Gaussian passbands with FWHM 10/20/40 nm form a probe
envelope, not a parameter claim. No physical parameter is fitted.

Band observable (main definition): the published quantity is the in-eclipse
over out-of-eclipse flux ratio of the same lunar element through the filter,
so the out-of-eclipse reference carries the R(lambda)*S(lambda) weights and
    T_band(theta) = integral R*S*T(lambda,theta) dlambda / integral R*S dlambda,
with T the transmission_profile output (E_lambda/F_lambda) and S the AM0 solar
irradiance. Lunar albedo is treated as grey inside each narrow band — the same
assumption level as the existing monochromatic comparison. An R-only weighting
variant quantifies the sensitivity to this definition.

Deterministic outputs: no wall-clock or runtime fields; fixed JSON key order
and fixed CSV row/column order; PNG written without timestamp metadata.

Reproduce:
    uv run python -m validation.filter_response_2011 --out outputs/filter_response_2011
Snapshot gate: rerun to outputs/regen_check; analysis.json and comparison.csv
must equal the committed copies in validation/results/filter_response_2011/
byte for byte.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.spectra import solar_irradiance
from eclipse.transport import Geometry, transmission_profile
from .compare_observations import DATA, load_observations

BANDS_NM = (503.0, 677.0)
WIDTHS_NM = (10.0, 20.0, 40.0)
SHAPES = ("rect", "gauss")
NOMINAL_WIDTH_NM = 20.0
SIGMA_CUTOFF = 3.0
GAUSSIAN_TRUNCATED_TAIL_FRACTION = math.erfc(SIGMA_CUTOFF/math.sqrt(2.0))
BASE_STEP_NM = 1.0
REFINED_STEP_NM = 0.5
GRID_CONVERGENCE_THRESHOLD = 1e-3
OUTER_MIN_ANGLE_DEG = 0.400
# Classification thresholds are planning judgements pinned here so the
# generated conclusion fields stay deterministic and independently auditable.
MATERIAL_SHIFT_FRACTION = 0.5
OUTERMOST_REDUCTION_MATERIAL = 2.0/3.0
OUTERMOST_REDUCTION_NONE = 1.5
SHAPE_VS_WIDTH_THRESHOLD = 0.1
CROSS_REPRODUCE_RTOL = 1e-12
CROSS_ADJUDICATE_RTOL = 1e-9
ROUND2_REFERENCE = "validation/results/round2/analysis.json"
OBS2011_REFERENCE = "validation/results/observations2011/comparison.json"
VARIANT_TABLES = (("molecular", "Molecular v1 (no aerosols)"),
                  ("fixed_opac", "Fixed OPAC free-troposphere prior"))

CSV_FIELDS = ["wavelength_nm", "shadow_angle_deg", "receiver_radius_m",
              "observed_relative_brightness",
              "molecular_mono_prediction", "fixed_opac_mono_prediction",
              "molecular_mono_ratio", "fixed_opac_mono_ratio"]
for _model, _ in VARIANT_TABLES:
    for _shape in SHAPES:
        for _width in WIDTHS_NM:
            CSV_FIELDS.append(f"{_model}_band_{_shape}{_width:g}_ratio")
for _model, _ in VARIANT_TABLES:
    for _shape in SHAPES:
        for _width in WIDTHS_NM:
            CSV_FIELDS.append(f"{_model}_band_{_shape}{_width:g}_diff")
CSV_FIELDS += ["fixed_opac_band_gauss20_R_weight_only_ratio",
               "fixed_opac_band_gauss20_R_weight_only_diff"]


def gaussian_sigma_nm(fwhm_nm: float) -> float:
    return fwhm_nm/(2.0*math.sqrt(2.0*math.log(2.0)))


def response_values(wavelength_nm, center_nm: float, shape: str, width_nm: float):
    """Passband response on a wavelength grid; zero outside the support."""
    wl = np.asarray(wavelength_nm, dtype=float)
    if shape == "rect":
        return np.asarray(np.abs(wl-center_nm) <= width_nm/2.0, dtype=float)
    if shape == "gauss":
        sigma = gaussian_sigma_nm(width_nm)
        response = np.exp(-4.0*math.log(2.0)*((wl-center_nm)/width_nm)**2)
        response[np.abs(wl-center_nm) > SIGMA_CUTOFF*sigma] = 0.0
        return response
    raise ValueError(f"Unknown passband shape {shape!r}")


def support_grid(center_nm: float, half_extent_nm: float, step_nm: float):
    """Center-anchored wavelength grid covering the support at the given step."""
    points = int(math.ceil(half_extent_nm/step_nm - 1e-9))
    return center_nm + step_nm*np.arange(-points, points+1)


def band_mean(weights, values) -> float:
    """Normalized weighted integral: sum(w*v)/sum(w). Pure kernel, no I/O."""
    weights = np.asarray(weights, dtype=float)
    values = np.asarray(values, dtype=float)
    if weights.shape != values.shape or weights.ndim != 1 or weights.size < 2:
        raise ValueError("Kernel weights and values must be equal-length 1-D vectors")
    total = float(weights.sum())
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("Kernel weight sum must be positive and finite")
    return float(weights @ values)/total


def variant_grids(center_nm: float) -> dict:
    """Base and refined grids for every shape/width variant of one band."""
    grids = {}
    for shape in SHAPES:
        for width in WIDTHS_NM:
            half = width/2.0 if shape == "rect" else SIGMA_CUTOFF*gaussian_sigma_nm(width)
            grids[(shape, width, BASE_STEP_NM)] = support_grid(center_nm, half, BASE_STEP_NM)
            grids[(shape, width, REFINED_STEP_NM)] = support_grid(center_nm, half, REFINED_STEP_NM)
    return grids


def union_grid(grids: dict) -> np.ndarray:
    stacked = np.concatenate([np.asarray(grid, dtype=float) for grid in grids.values()])
    return np.unique(stacked)


def grid_indices(union: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Exact-match indices of a kernel grid inside the union grid."""
    indices = np.searchsorted(union, grid)
    if (np.any(indices >= union.size)
            or not np.all(union[np.minimum(indices, union.size-1)] == grid)):
        raise ValueError("Kernel grid points are not an exact subset of the union grid")
    return indices


def load_round2_predictions() -> dict:
    path = Path(__file__).resolve().parent.parent/ROUND2_REFERENCE
    report = json.loads(path.read_text())
    return {(float(r["wavelength_nm"]), float(r["shadow_angle_deg"])):
            {"molecular_prediction": float(r["molecular_prediction"]),
             "fixed_opac_prediction": float(r["fixed_opac_prediction"])}
            for r in report["records"]}


def load_obs2011_predictions() -> dict:
    path = Path(__file__).resolve().parent.parent/OBS2011_REFERENCE
    report = json.loads(path.read_text())
    predictions = {}
    for band in report["bands"]:
        for endpoint in ("innermost_readout", "outermost_readout"):
            row = band[endpoint]
            predictions[(float(row["wavelength_nm"]), float(row["shadow_angle_deg"]))] = {
                "model_event_300du": float(row["model_event_300du"])}
    return predictions


def environment_snapshot() -> dict:
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "platform": sys.platform}


def cross_check_gate(mono: dict, observations: list, adjudicated_rtol: float | None) -> dict:
    """Fail-closed reproduction gate for the monochromatic columns.

    rtol <= 1e-12 counts as reproduced; the middle band pauses and escalates to
    the owner unless --adjudicated-rtol covers the observed maximum (the
    adjudication must already be recorded in the task progress log); above 1e-9
    the run is blocked outright. No formal product is written in either
    stopping case because this gate runs before any output exists.
    """
    keys = [(float(row["wavelength_nm"]), float(row["shadow_angle_deg"]))
            for row in observations]
    molecular = {key: float(mono["molecular"][i]) for i, key in enumerate(keys)}
    opac = {key: float(mono["fixed_opac"][i]) for i, key in enumerate(keys)}
    round2 = load_round2_predictions()
    obs2011 = load_obs2011_predictions()
    references = {
        "round2_molecular": (
            ROUND2_REFERENCE, "records[].molecular_prediction",
            relative_differences(molecular, {k: v["molecular_prediction"] for k, v in round2.items()})),
        "round2_fixed_opac": (
            ROUND2_REFERENCE, "records[].fixed_opac_prediction",
            relative_differences(opac, {k: v["fixed_opac_prediction"] for k, v in round2.items()})),
        "observations2011_model_event_300du": (
            OBS2011_REFERENCE, "bands[].innermost_readout/outermost_readout.model_event_300du",
            relative_differences(molecular, {k: v["model_event_300du"] for k, v in obs2011.items()})),
    }
    summary, worst_value, worst_name, worst_key = {}, -1.0, None, None
    for name, (path, field, diffs) in references.items():
        maximum = max(diffs.values())
        summary[name] = {"reference": path, "field": field,
                         "point_count": len(diffs), "max_relative_difference": maximum}
        if maximum > worst_value:
            worst_value, worst_name = maximum, name
            worst_key = max(diffs, key=diffs.get)
    if worst_value <= CROSS_REPRODUCE_RTOL:
        return {"status": "reproduced", "threshold_reproduce": CROSS_REPRODUCE_RTOL,
                "threshold_adjudicate": CROSS_ADJUDICATE_RTOL, "references": summary,
                "adjudicated_rtol": None, "adjudication_reference": None}
    evidence = [
        "MONO CROSS-CHECK ADJUDICATION REQUIRED (no formal product written)",
        f"max rtol {worst_value:.6g} at {worst_key} via {worst_name}",
        f"environment: {environment_snapshot()}",
        "per-point rtol: "
        + str({str(k): f"{v:.6g}" for k, v in sorted(references[worst_name][2].items())}),
        "reproduction commands:",
        "  uv run python -m validation.filter_response_2011 --out outputs/filter_response_2011",
        "  uv run python -m validation.round2_analysis --out outputs/round2_check",
        f"references: {ROUND2_REFERENCE} (records[].molecular_prediction, "
        f"records[].fixed_opac_prediction); {OBS2011_REFERENCE} "
        "(bands[].innermost_readout/outermost_readout.model_event_300du)",
    ]
    if worst_value <= CROSS_ADJUDICATE_RTOL:
        if adjudicated_rtol is not None and adjudicated_rtol >= worst_value:
            return {"status": "adjudicated", "threshold_reproduce": CROSS_REPRODUCE_RTOL,
                    "threshold_adjudicate": CROSS_ADJUDICATE_RTOL, "references": summary,
                    "adjudicated_rtol": float(adjudicated_rtol),
                    "adjudication_reference": "owner adjudication recorded in the task progress log"}
        evidence.append("pause: owner adjudication required before any bandpass product is "
                        "produced; after recording the adjudication in the task progress log, "
                        "rerun with --adjudicated-rtol <value> covering the rtol above")
    else:
        evidence.append("block: rtol above the adjudication ceiling; investigate before rerunning")
    print("\n".join(evidence))
    raise RuntimeError("Monochromatic cross-check gate stopped the run; see evidence above")


def relative_differences(computed: dict, reference: dict) -> dict:
    """Relative differences at the readouts the reference contains (>=1)."""
    common = sorted(set(computed) & set(reference))
    if not common:
        raise ValueError("Cross-check reference shares no readout with the observations")
    return {key: abs(computed[key]-reference[key])/abs(reference[key]) for key in common}


def build_band_products(observations: list, geometry: Geometry, tables: dict,
                        adjudicated_rtol: float | None) -> dict:
    """Monochromatic gate first, then band-integrated columns.

    One transmission_profile call per (band, model variant) carries the union
    of all kernel grids (base and refined) at all readout radii.
    """
    angles = np.array([row["shadow_angle_deg"] for row in observations])
    wavelengths = np.array([row["wavelength_nm"] for row in observations])
    observed = np.array([row["observed_relative_brightness"] for row in observations])
    unique_wl, selection = np.unique(wavelengths, return_inverse=True)
    radius = geometry.moon_distance_m*np.tan(np.deg2rad(angles))
    diagonal = np.arange(len(angles))

    mono = {name: transmission_profile(table, radius, unique_wl, geometry=geometry)[diagonal, selection]
            for name, table in tables.items()}
    cross_check = cross_check_gate(mono, observations, adjudicated_rtol)

    band_values, sensitivity_values, grid_report = {}, {}, {}
    for band in BANDS_NM:
        grids = variant_grids(band)
        union = union_grid(grids)
        sampled_by_table = {}
        for name, table in tables.items():
            sampled = transmission_profile(table, radius, union, geometry=geometry)
            sampled_by_table[name] = sampled
            for shape in SHAPES:
                for width in WIDTHS_NM:
                    for step in (BASE_STEP_NM, REFINED_STEP_NM):
                        grid = grids[(shape, width, step)]
                        indices = grid_indices(union, grid)
                        weights = kernel_weights(grid, band, shape, width,
                                                 solar_irradiance(grid))
                        values = sampled[:, indices]
                        band_values[(name, band, shape, width, step)] = np.array(
                            [band_mean(weights, values[i]) for i in range(len(angles))])
        grid = grids[("gauss", NOMINAL_WIDTH_NM, BASE_STEP_NM)]
        indices = grid_indices(union, grid)
        weights = kernel_weights(grid, band, "gauss", NOMINAL_WIDTH_NM, None)
        values = sampled_by_table["fixed_opac"][:, indices]
        sensitivity_values[band] = np.array(
            [band_mean(weights, values[i]) for i in range(len(angles))])
        grid_report[band] = {"union_point_count": int(union.size),
                             "union_half_extent_nm": float(union[-1]-band)}
    return {"angles": angles, "wavelengths": wavelengths, "observed": observed,
            "radius": radius, "unique_wl": unique_wl, "mono": mono,
            "band_values": band_values, "sensitivity_values": sensitivity_values,
            "grid_report": grid_report, "cross_check": cross_check}


def kernel_weights(wavelength_nm, center_nm: float, shape: str, width_nm: float,
                   solar_weight) -> np.ndarray:
    weights = response_values(wavelength_nm, center_nm, shape, width_nm)
    if solar_weight is not None:
        weights = weights*np.asarray(solar_weight, dtype=float)
    return weights


def band_ratio(products, model: str, band: float, shape: str, width: float,
               step=BASE_STEP_NM) -> np.ndarray:
    values = products["band_values"][(model, band, shape, width, step)]
    return values/products["observed"]


def per_band_summary(products, band: float) -> dict:
    angles, observed = products["angles"], products["observed"]
    wavelengths = products["wavelengths"]
    mono_ratios = {model: products["mono"][model]/observed
                   for model, _ in VARIANT_TABLES}
    main = "fixed_opac"
    mask = wavelengths == band
    outer = mask & (angles >= OUTER_MIN_ANGLE_DEG)
    variants = {}
    for model, _ in VARIANT_TABLES:
        for shape in SHAPES:
            for width in WIDTHS_NM:
                difference = band_ratio(products, model, band, shape, width)-mono_ratios[model]
                variants[f"{model}_{shape}{width:g}"] = {
                    "max_abs_difference": float(np.max(np.abs(difference[mask]))),
                    "outer_max_abs_difference": float(np.max(np.abs(difference[outer])))}
    main_variants = {name: value for name, value in variants.items()
                     if name.startswith(main+"_")}
    outermost = mask & (angles == np.max(angles[mask]))
    mono_outermost = float(mono_ratios[main][outermost][0])
    best = min(float(np.abs(band_ratio(products, main, band, shape, width)[outermost][0]-1.0))
               for shape in SHAPES for width in WIDTHS_NM)/abs(mono_outermost-1.0)
    summary = {
        "readout_count": int(mask.sum()),
        "angle_range_deg": [float(np.min(angles[mask])), float(np.max(angles[mask]))],
        "mono_ratio_range_fixed_opac": [float(np.min(mono_ratios[main][mask])),
                                        float(np.max(mono_ratios[main][mask]))],
        "mono_ratio_range_molecular": [float(np.min(mono_ratios["molecular"][mask])),
                                       float(np.max(mono_ratios["molecular"][mask]))],
        "band_vs_mono_max_abs_difference": variants,
        "band_vs_mono_envelope_max": max(v["max_abs_difference"] for v in variants.values()),
        "main_variant_envelope_max": max(v["max_abs_difference"] for v in main_variants.values()),
        "outer_band_effect_max": max(v["outer_max_abs_difference"] for v in variants.values()),
        "main_variant_outer_effect_max": max(v["outer_max_abs_difference"]
                                             for v in main_variants.values()),
        "outer_mono_deviation_max": float(np.max(np.abs(mono_ratios[main][outer]-1.0))),
        "outermost_deviation_ratio_best_variant": best,
        "outermost_classification": classify_outermost(best),
    }
    inner = mask & (angles < OUTER_MIN_ANGLE_DEG)
    if np.any(inner):
        index = np.where(mask)[0][np.argmin(angles[mask])]
        effect = max(abs(float(band_ratio(products, main, band, shape, width)[index])
                         - float(mono_ratios[main][index]))
                     for shape in SHAPES for width in WIDTHS_NM)
        deviation = abs(float(mono_ratios[main][index])-1.0)
        summary["innermost_angle_deg"] = float(angles[index])
        summary["innermost_mono_deviation"] = deviation
        summary["innermost_band_effect_max"] = effect
        summary["innermost_classification"] = (
            "band definition cannot account for the innermost deep-point residual"
            if effect < MATERIAL_SHIFT_FRACTION*deviation else
            "band definition materially shifts the innermost deep-point ratio")
    return summary


def classify_outermost(best_ratio: float) -> str:
    if best_ratio < OUTERMOST_REDUCTION_MATERIAL:
        return ("band definition materially reduces the outermost model/observed deviation "
                f"(best variant deviation ratio {best_ratio:.3g} < "
                f"{OUTERMOST_REDUCTION_MATERIAL:.3g})")
    if best_ratio <= OUTERMOST_REDUCTION_NONE:
        return ("band definition changes the outermost model/observed deviation only marginally "
                f"(best variant deviation ratio {best_ratio:.3g})")
    return ("band definition does not reduce the outermost model/observed deviation "
            f"(best variant deviation ratio {best_ratio:.3g} > {OUTERMOST_REDUCTION_NONE:.3g})")


def shape_vs_width_summary(products) -> dict:
    shape_effect, width_effect = 0.0, 0.0
    for band in BANDS_NM:
        mask = products["wavelengths"] == band
        for width in WIDTHS_NM:
            difference = (band_ratio(products, "fixed_opac", band, "rect", width)
                          - band_ratio(products, "fixed_opac", band, "gauss", width))
            shape_effect = max(shape_effect, float(np.max(np.abs(difference)[mask])))
        for shape in SHAPES:
            difference = (band_ratio(products, "fixed_opac", band, shape, WIDTHS_NM[0])
                          - band_ratio(products, "fixed_opac", band, shape, WIDTHS_NM[-1]))
            width_effect = max(width_effect, float(np.max(np.abs(difference)[mask])))
    ratio = shape_effect/width_effect if width_effect > 0 else math.inf
    adequate = ratio < SHAPE_VS_WIDTH_THRESHOLD
    return {"max_shape_effect_same_width": shape_effect,
            "max_width_effect_10_to_40nm": width_effect,
            "shape_over_width_ratio": ratio,
            "threshold": SHAPE_VS_WIDTH_THRESHOLD,
            "shape_approximation_adequate": adequate,
            "statement": ("shape effect is below 10% of the width effect: rectangular/Gaussian "
                          "approximation is adequate at the envelope's present resolution"
                          if adequate else
                          "shape effect is not negligible against the width effect: the assumed "
                          "band shape is a material limitation of the envelope")}


def grid_convergence_report(products) -> dict:
    worst = 0.0
    for band in BANDS_NM:
        mask = products["wavelengths"] == band
        for model, _ in VARIANT_TABLES:
            for shape in SHAPES:
                for width in WIDTHS_NM:
                    base = band_ratio(products, model, band, shape, width, BASE_STEP_NM)
                    refined = band_ratio(products, model, band, shape, width, REFINED_STEP_NM)
                    worst = max(worst, float(np.max(np.abs(base/refined-1.0)[mask])))
    within = worst <= GRID_CONVERGENCE_THRESHOLD
    return {"base_step_nm": BASE_STEP_NM, "refined_step_nm": REFINED_STEP_NM,
            "max_relative_change_base_to_refined": worst,
            "threshold": GRID_CONVERGENCE_THRESHOLD, "within_threshold": within,
            "termination_reason": (
                "converged within threshold at the first refinement"
                if within else
                f"refined step {REFINED_STEP_NM:g} nm is below the AM0 table native resolution "
                "(1 nm); further halving cannot improve a linearly interpolated integrand, and "
                "the recorded value is the honest terminal value")}


def build_records(products) -> list:
    records = []
    for i, angle in enumerate(products["angles"]):
        band = float(products["wavelengths"][i])
        row = {"wavelength_nm": band,
               "shadow_angle_deg": float(angle),
               "receiver_radius_m": float(products["radius"][i]),
               "observed_relative_brightness": float(products["observed"][i]),
               "molecular_mono_prediction": float(products["mono"]["molecular"][i]),
               "fixed_opac_mono_prediction": float(products["mono"]["fixed_opac"][i]),
               "molecular_mono_ratio": float(products["mono"]["molecular"][i]/products["observed"][i]),
               "fixed_opac_mono_ratio": float(products["mono"]["fixed_opac"][i]/products["observed"][i])}
        for model, _ in VARIANT_TABLES:
            for shape in SHAPES:
                for width in WIDTHS_NM:
                    row[f"{model}_band_{shape}{width:g}_ratio"] = float(
                        band_ratio(products, model, band, shape, width)[i])
        for model, _ in VARIANT_TABLES:
            for shape in SHAPES:
                for width in WIDTHS_NM:
                    row[f"{model}_band_{shape}{width:g}_diff"] = (
                        row[f"{model}_band_{shape}{width:g}_ratio"]-row[f"{model}_mono_ratio"])
        sensitivity = float(products["sensitivity_values"][band][i]/products["observed"][i])
        row["fixed_opac_band_gauss20_R_weight_only_ratio"] = sensitivity
        row["fixed_opac_band_gauss20_R_weight_only_diff"] = (
            sensitivity-row["fixed_opac_mono_ratio"])
        records.append(row)
    return records


def compute_analysis(adjudicated_rtol: float | None = None):
    observations = load_observations()
    geometry_data = json.loads((DATA/"geometry.json").read_text())
    source = json.loads((DATA/"source.json").read_text())
    filter_info = source["filter_information"]
    geometry = Geometry(moon_distance_m=geometry_data["moon_distance_m"],
                        sun_distance_m=geometry_data["sun_distance_m"])
    tables = {"molecular": build_ray_table(Atmosphere()),
              "fixed_opac": build_ray_table(Atmosphere(),
                                            aerosol_profile=AerosolProfile.opac_free_troposphere())}
    products = build_band_products(observations, geometry, tables, adjudicated_rtol)
    summaries = {f"{band:g}": per_band_summary(products, band) for band in BANDS_NM}
    shape_width = shape_vs_width_summary(products)
    convergence = grid_convergence_report(products)
    records = build_records(products)
    sensitivity = max(abs(row["fixed_opac_band_gauss20_R_weight_only_diff"]) for row in records)
    per_band_conclusions = {}
    for band in BANDS_NM:
        summary = summaries[f"{band:g}"]
        per_band_conclusions[f"{band:g}"] = {
            "statement": (
                f"Band-integrated ratios differ from the main-variant monochromatic prediction "
                f"by at most {summary['main_variant_envelope_max']:.3g} in model/observed ratio "
                f"across the {{10,20,40}} nm rectangular/Gaussian envelope (both model variants "
                f"together: {summary['band_vs_mono_envelope_max']:.3g}); outer-region "
                f"(theta >= {OUTER_MIN_ANGLE_DEG:g} deg) main-variant monochromatic ratios deviate "
                f"from unity by up to {summary['outer_mono_deviation_max']:.3g}. "
                + summary["outermost_classification"])}
    report = {
        "task": "filter-response-2011",
        "observable_definition": {
            "published_quantity": "in-eclipse over out-of-eclipse brightness of the same lunar element through the filter",
            "band_observable": "T_band(theta) = integral R(lambda)*S(lambda)*T(lambda,theta) dlambda / integral R(lambda)*S(lambda) dlambda",
            "T_source": "eclipse.transport.transmission_profile (E_lambda/F_lambda)",
            "S_source": "eclipse.spectra.solar_irradiance (AM0); weights of the out-of-eclipse reference",
            "albedo_assumption": "lunar albedo treated as grey inside each narrow band; same assumption level as the committed monochromatic comparison",
            "sensitivity_variant": "R(lambda)-only weighting at the nominal width, main variant; quantifies the S(lambda)-weighting choice",
        },
        "provenance": {
            "observations_dir": "validation/observations/ugolnikov2011",
            "source_pdf_sha256": source["source_pdf_sha256"],
            "filter_branch": filter_info["branch"],
            "filter_information_authority": "validation/observations/ugolnikov2011/source.json filter_information",
            "filter_approximation_grade": filter_info["approximation_grade"],
            "geometry_source": "validation/observations/ugolnikov2011/geometry.json",
            "adopted_geometry": {"moon_distance_m": geometry_data["moon_distance_m"],
                                 "sun_distance_m": geometry_data["sun_distance_m"]},
            "production_module_sha256": {
                name: hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
                for name in ["aerosol.py", "atmosphere.py", "rays.py", "spectra.py", "transport.py"]},
            "fitted_parameters": [],
            "adjudicated_rtol": products["cross_check"]["adjudicated_rtol"],
            "adjudication_reference": products["cross_check"]["adjudication_reference"],
        },
        "model_variants": {name: label for name, label in VARIANT_TABLES},
        "band_variants": {
            "shapes": list(SHAPES), "widths_nm": list(WIDTHS_NM),
            "nominal_width_nm": NOMINAL_WIDTH_NM,
            "figure_draws": "nominal width only, both shapes; full envelope not drawn",
            "gaussian_truncation_sigma": SIGMA_CUTOFF,
            "gaussian_truncated_tail_fraction": GAUSSIAN_TRUNCATED_TAIL_FRACTION,
        },
        "grid": {f"{band:g}": products["grid_report"][band] for band in BANDS_NM},
        "grid_convergence": convergence,
        "definition_sensitivity_R_weight_only": {
            "variant": "fixed_opac gauss20 R-only weighting",
            "max_abs_ratio_difference_vs_main_definition": sensitivity},
        "monochromatic_cross_check": products["cross_check"],
        "per_band": summaries,
        "shape_vs_width": shape_width,
        "conclusion": {
            "outer_residuals": per_band_conclusions,
            "innermost_deep_point": summaries["677"].get(
                "innermost_classification", "no 677 nm readout inside the outer-region boundary"),
            "shape_vs_width": shape_width["statement"],
            "definition_sensitivity": (
                f"R-only versus R*S weighting changes the nominal band ratio by at most "
                f"{sensitivity:.3g}; the S(lambda)-weighting choice is "
                + ("negligible at the 1e-3 level" if sensitivity < 1e-3 else "quantified above")),
        },
        "limitations": [
            "Digitized smoothed published curves, not independent exposures; graph readout allowances are not observational errors",
            "Filter width and shape are assumed probes (branch C envelope), not measured properties",
            "Lunar albedo treated as grey within each passband",
            "Representative geocentric plane distance; no per-readout timing or topocentric geometry",
            "Two bands and one event do not validate full visible colour",
        ],
    }
    return report, records, products


def write_outputs(output_dir: Path, report: dict, records: list, products: dict):
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir/"comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    (output_dir/"analysis.json").write_text(json.dumps(report, indent=2)+"\n")
    style = {"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
             "axes.spines.right": False, "figure.facecolor": "white", "axes.facecolor": "white"}
    with plt.rc_context(style):
        fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), sharex=True,
                                 gridspec_kw={"wspace": .22})
        mono_ratio = products["mono"]["fixed_opac"]/products["observed"]
        for column, band in enumerate(BANDS_NM):
            ax = axes[column]
            mask = products["wavelengths"] == band
            angles = products["angles"][mask]
            ax.plot(angles, mono_ratio[mask], "o-", color="#ac4f31", ms=4, lw=1.6,
                    label="Monochromatic (nominal wavelength)")
            for shape, linestyle in (("rect", "--"), ("gauss", ":")):
                curve = band_ratio(products, "fixed_opac", band, shape, NOMINAL_WIDTH_NM)[mask]
                label = (f"Band-integrated, {NOMINAL_WIDTH_NM:g} nm "
                         f"{'rectangular' if shape == 'rect' else 'Gaussian'} FWHM")
                ax.plot(angles, curve, linestyle=linestyle, color="#227580", lw=1.6,
                        marker="s", ms=3, label=label)
            ax.axhline(1, color="#253d4a", lw=.8)
            ax.set_title(f"{band:g} nm / fixed OPAC prior", loc="left", fontsize=12, pad=10)
            ax.set_xlabel("Angular distance from shadow center (degrees)")
            ax.set_xlim(.12, .60)
            ax.grid(alpha=.18)
            ax.legend(frameon=False, fontsize=8, loc="upper left")
        axes[0].set_ylabel("Model / published observation")
        fig.suptitle("June 2011 eclipse: monochromatic versus band-integrated predictions",
                     fontsize=15, y=.98)
        fig.text(.085, .02,
                 "Assumed passbands (branch C envelope, nominal width drawn): grade and keyword search in "
                 "validation/observations/ugolnikov2011/source.json filter_information. No fitted parameters.",
                 fontsize=8.5, color="#535d67")
        fig.subplots_adjust(top=.86, bottom=.20, left=.07, right=.98)
        fig.savefig(output_dir/"comparison.png", dpi=180,
                    metadata={"Software": None, "CreationDate": None})
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("outputs/filter_response_2011"))
    parser.add_argument("--adjudicated-rtol", type=float, default=None,
                        help="owner-adjudicated cross-check rtol ceiling; applies only inside "
                             "the 1e-12 to 1e-9 adjudication band and must match the "
                             "adjudication recorded in the task progress log")
    args = parser.parse_args()
    report, records, products = compute_analysis(args.adjudicated_rtol)
    write_outputs(args.out, report, records, products)
    for band, summary in report["per_band"].items():
        print(f"{band} nm: band-vs-mono envelope max |difference| = "
              f"{summary['band_vs_mono_envelope_max']:.3g}; "
              f"{summary['outermost_classification']}")
    print(f"Monochromatic cross-check: {report['monochromatic_cross_check']['status']}")
    print(f"Results: {args.out.resolve()}")


if __name__ == "__main__":
    main()
