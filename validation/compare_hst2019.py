"""Fixed-model comparison of the 2019-01-21 HST/STIS eclipse spectra.

Extraction replicates the paper's own method (Youngblood et al. 2020 Sec.
II.2) on the raw archive: per-column median across the filled 52" slit of the
CR-rejected (or single-exposure) calibrated images, counts per second, no
background subtraction, no absolute flux calibration; in-eclipse spectra are
divided by the time-weighted average of the twelve full-moon G430L spectra.
The pipeline _sx1/_x1d products contribute only their wavelength solutions:
their point-source flux extraction places the background regions on the
filled-slit lunar spectrum and is never used for flux.

The paper divides each spectrum by a Lommel-Seeliger factor (F_LS in
0.4975-0.5021, incidence/emission angles closely coupled), which changes
ratios by <=1%, and divides only the penumbral spectra by a limb-darkening
factor (G430L I/I0 0.7-0.79); both omissions are recorded, so the umbral
main comparison aligns with the published convention while the auxiliary
penumbral branch carries a known ~1/0.75 scale difference against published
figures.

The model side is the fixed annular transport: transmission_profile at the
reconstructed shadow-axis distance of each sampled position, molecular-only
and fixed OPAC default variants, nothing fitted. Full-moon reference
transmission is unity by construction and is verified geometrically.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.transport import Geometry, geometric_transmission, transmission_profile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"validation/observations/hst2019"
CACHE = ROOT/"outputs/sources/hst2019"
MODEL_TABLE_MIN_NM = 300.0
G430L_MAX_NM = 570.0
# The paper reports no photons at <=4200 A during the umbral phase, so the
# umbral quantitative metrics use 430-570 nm; 300-420 nm is reported as a
# blue-end non-detection check instead.
UMBRAL_QUANTITATIVE_MIN_NM = 430.0
UMBRAL_BLUE_CHECK_MAX_NM = 420.0
BAND_HALFWIDTH_NM = 5.0
SLOPE_BANDS_NM = [(430.0, 570.0), (500.0, 570.0)]
# Physically invalid pixel categories only (lost data, bad detector, aperture
# feature, saturation). Hot pixels and cosmic-ray flags stay in: the median
# over ~1024 rows is robust to them, matching the paper's flag list.
DQ_EXCLUDE_MASK = 2 | 4 | 8 | 4096
MIN_VALID_FRACTION = 0.5


def column_medians(sci: np.ndarray, dq: np.ndarray, exptime_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Median counts/s per column over unflagged rows, with validity flags.

    The slit fills the detector height, so every unflagged row is part of the
    spectrum. Columns whose valid fraction drops below MIN_VALID_FRACTION are
    flagged unreliable (e.g. saturated penumbral bins).
    """
    sci = np.asarray(sci,dtype=float)
    dq = np.asarray(dq)
    if sci.shape!=dq.shape or sci.ndim!=2:
        raise ValueError("Science and DQ arrays must share one 2D shape")
    valid = (dq & DQ_EXCLUDE_MASK)==0
    counts = np.full(sci.shape[1],np.nan)
    valid_fraction = np.zeros(sci.shape[1])
    for column in range(sci.shape[1]):
        rows = sci[valid[:,column],column]
        valid_fraction[column] = len(rows)/sci.shape[0]
        if len(rows):
            counts[column] = np.median(rows)/exptime_s
    return counts,valid_fraction


def wavelength_solution(wavelength_table) -> np.ndarray:
    """Column wavelengths [nm] from the 1D product's WAVELENGTH column."""
    values = np.asarray(wavelength_table,dtype=float).ravel()
    if values.size!=1024 or np.any(np.diff(values)<=0):
        raise ValueError("Expected a strictly increasing 1024-sample wavelength solution")
    return values/10.0


def header_wavelength_solution(crval: float, cd1_1: float, crpix: float) -> np.ndarray:
    """Nominal linear wavelengths [nm] from the 2D image header, for audit."""
    columns_1based = np.arange(1,1025)
    return (crval+cd1_1*(columns_1based-crpix))/10.0


def window_mask(wavelength_nm: np.ndarray, low_nm: float, high_nm: float) -> np.ndarray:
    mask = (wavelength_nm>=low_nm)&(wavelength_nm<=high_nm)
    if not np.any(mask):
        raise ValueError(f"Empty wavelength window [{low_nm}, {high_nm}] nm")
    return mask


def time_weighted_mean(rates: list[np.ndarray], exptime_s: list[float]) -> np.ndarray:
    """Exposure-time-weighted mean rate; sums counts before dividing."""
    total = sum(exptime_s)
    if total<=0 or any(value<=0 for value in exptime_s):
        raise ValueError("Exposure times must be positive")
    return sum(rate*weight for rate,weight in zip(rates,exptime_s))/total


def band_median(values: np.ndarray, wavelength_nm: np.ndarray, center_nm: float,
                halfwidth_nm: float = BAND_HALFWIDTH_NM) -> float:
    mask = window_mask(wavelength_nm,center_nm-halfwidth_nm,center_nm+halfwidth_nm)
    return float(np.median(values[mask]))


def dataset_model_profile(table, member_records: list[dict], wavelengths_nm: np.ndarray,
                          radial_offset_m: float = 0.) -> np.ndarray:
    """Exposure-time-weighted transmission over the set's member positions."""
    weighted = np.zeros(len(wavelengths_nm))
    total = sum(member["texp_s"] for member in member_records)
    for member in member_records:
        geometry = Geometry(moon_distance_m=member["receiver_axial_distance_m"],
                            sun_distance_m=member["sun_site_distance_m"])
        radius = member["shadow_axis_distance_m"]+radial_offset_m
        weighted += member["texp_s"]*transmission_profile(table,[radius],wavelengths_nm,
                                                          geometry=geometry)[0]
    return weighted/total


def load_products(contract: dict, geometry: dict) -> dict:
    """Extract one spectrum per dataset and attach its geometry members."""
    from astropy.io import fits
    products = contract["mast_download"]["products"]
    membership = {**contract["cr_split_membership"],
                  **{dataset:[row["exposure_name"] for row in geometry["records"]
                              if row["dataset_id"]==dataset]
                     for dataset in contract["grouping"]["full_moon_out_of_eclipse"]}}
    records_by_name = {record["exposure_name"]:record for record in geometry["records"]}
    datasets = {}
    for dataset,exposures in membership.items():
        two_d = f"{dataset}_crj.fits" if f"{dataset}_crj.fits" in products else f"{dataset}_flt.fits"
        one_d = f"{dataset}_sx1.fits" if f"{dataset}_sx1.fits" in products else f"{dataset}_x1d.fits"
        with fits.open(CACHE/two_d) as hdus:
            counts,valid_fraction = column_medians(hdus["sci",1].data,hdus["dq",1].data,
                                                   float(hdus["sci",1].header["EXPTIME"]))
            header = hdus["sci",1].header
            nominal = header_wavelength_solution(header["CRVAL1"],header["CD1_1"],header["CRPIX1"])
        with fits.open(CACHE/one_d) as hdus:
            wavelength_nm = wavelength_solution(hdus[1].data["WAVELENGTH"][0])
        alignment = float(np.max(np.abs(wavelength_nm-nominal)))
        if alignment>20.0:
            raise RuntimeError(f"{dataset}: 1D wavelength solution disagrees with the image "
                               f"header solution by {alignment:.1f} nm")
        datasets[dataset] = {"exposures":exposures,"counts_per_s":counts,
                             "wavelength_nm":wavelength_nm,
                             "header_vs_1d_max_diff_nm":alignment,
                             "min_valid_fraction":float(valid_fraction.min()),
                             "members":[records_by_name[name] for name in exposures]}
    return datasets


def run_analysis():
    contract = json.loads((DATA/"source.json").read_text())
    geometry = json.loads((DATA/"geometry.json").read_text())
    datasets = load_products(contract,geometry)

    window = window_mask(datasets["ody904080"]["wavelength_nm"],MODEL_TABLE_MIN_NM,G430L_MAX_NM)
    for dataset,data in datasets.items():
        if not np.array_equal(np.where(window)[0],np.where(window_mask(
                data["wavelength_nm"],MODEL_TABLE_MIN_NM,G430L_MAX_NM))[0]):
            raise ValueError(f"{dataset}: wavelength grids disagree inside the window")
        if data["min_valid_fraction"]<MIN_VALID_FRACTION:
            raise RuntimeError(f"{dataset}: columns with <{MIN_VALID_FRACTION:.0%} valid pixels present")
    wavelengths = datasets["ody904080"]["wavelength_nm"][window]

    grouping = contract["grouping"]
    umbral_datasets = grouping["umbral_in_eclipse"]
    penumbral_datasets = grouping["penumbral_auxiliary"]
    in_eclipse_datasets = umbral_datasets+penumbral_datasets
    full_moon = grouping["full_moon_out_of_eclipse"]
    reference = time_weighted_mean([datasets[name]["counts_per_s"][window] for name in full_moon],
                                   [datasets[name]["members"][0]["texp_s"] for name in full_moon])
    for name in full_moon:
        radius = datasets[name]["members"][0]["shadow_axis_distance_m"]
        geometry_single = Geometry(moon_distance_m=datasets[name]["members"][0]["receiver_axial_distance_m"],
                                   sun_distance_m=datasets[name]["members"][0]["sun_site_distance_m"])
        if geometric_transmission([radius],geometry_single)[0] < 1-1e-12:
            raise RuntimeError(f"{name}: full-moon reference is not fully outside the geometric shadow")

    solar_scale = {}
    for dataset in in_eclipse_datasets:
        members = datasets[dataset]["members"]
        weight = sum(member["texp_s"] for member in members)
        sun_in = sum(member["texp_s"]/member["sun_site_distance_m"]**2 for member in members)/weight
        out_weight = sum(datasets[name]["members"][0]["texp_s"] for name in full_moon)
        sun_out = sum(datasets[name]["members"][0]["texp_s"]
                      /datasets[name]["members"][0]["sun_site_distance_m"]**2
                      for name in full_moon)/out_weight
        solar_scale[dataset] = float(sun_in/sun_out)
    for dataset,scale in solar_scale.items():
        if abs(scale-1)>0.001:
            raise RuntimeError(f"{dataset}: Sun-site range changed by more than 0.1% within the eclipse day")

    profiles = {"molecular":None,"fixed_opac":AerosolProfile.opac_free_troposphere()}
    observed = {}
    predictions = {}
    probes = {}
    for dataset in in_eclipse_datasets:
        observed[dataset] = datasets[dataset]["counts_per_s"][window]/reference
    for variant,profile in profiles.items():
        table = build_ray_table(Atmosphere(),aerosol_profile=profile)
        for dataset in in_eclipse_datasets:
            predictions[(variant,dataset)] = dataset_model_profile(
                table,datasets[dataset]["members"],wavelengths)*solar_scale[dataset]
        probes[variant] = {}
        for dataset in umbral_datasets:
            cases = {}
            for offset in geometry["spatial_sensitivity"]["sample_probe_offsets_m"]:
                cases[f"{int(offset/1000):+d}km"] = dataset_model_profile(
                    table,datasets[dataset]["members"],wavelengths,radial_offset_m=offset)
            envelope = geometry["spatial_sensitivity"]["pointing_track_error_km_max"]*1000.0
            cases["+track_envelope"] = dataset_model_profile(
                table,datasets[dataset]["members"],wavelengths,radial_offset_m=envelope)
            probes[variant][dataset] = cases

    metrics = {}
    for variant in profiles:
        metrics[variant] = {}
        for dataset in umbral_datasets:
            quantitative = window_mask(wavelengths,UMBRAL_QUANTITATIVE_MIN_NM,G430L_MAX_NM)
            ratio = predictions[(variant,dataset)][quantitative]/observed[dataset][quantitative]
            blue = window_mask(wavelengths,MODEL_TABLE_MIN_NM,UMBRAL_BLUE_CHECK_MAX_NM)
            metrics[variant][dataset] = {
                "model_to_observed_median":float(np.median(ratio)),
                "model_to_observed_min":float(ratio.min()),
                "model_to_observed_max":float(ratio.max()),
                "slope_model_vs_observed":{f"{num:g}_over_{den:g}":
                    (band_median(predictions[(variant,dataset)],wavelengths,num)
                     /band_median(predictions[(variant,dataset)],wavelengths,den))
                    /(band_median(observed[dataset],wavelengths,num)
                      /band_median(observed[dataset],wavelengths,den))
                    for num,den in SLOPE_BANDS_NM},
                "blue_end_observed_median":float(np.median(observed[dataset][blue])),
                "blue_end_observed_max":float(observed[dataset][blue].max()),
                "blue_end_model_median":float(np.median(predictions[(variant,dataset)][blue])),
                "blue_end_scope":"Scattered-light upper-limit region only: the observed blue end carries a "
                                 "long-wavelength scattered-light pedestal (observed ratios there exceed the "
                                 "numerically-zero model and peak at the 300 nm grid edge), so these numbers "
                                 "bound contamination, not umbral transmission; the paper likewise reports no "
                                 "photons <=4200 A in the umbral phase",
                "quantitative_scope":"430-570 nm window; the observed spectrum rises monotonically through "
                                     "this window (Rayleigh-slope region the paper detects), free of the "
                                     "blue-edge scattered-light bump"}
    metrics["penumbral_auxiliary_fixed_opac_only"] = {}
    for dataset in penumbral_datasets:
        quantitative = window_mask(wavelengths,MODEL_TABLE_MIN_NM,G430L_MAX_NM)
        ratio = predictions[("fixed_opac",dataset)][quantitative]/observed[dataset][quantitative]
        metrics["penumbral_auxiliary_fixed_opac_only"][dataset] = {
            "model_to_observed_median":float(np.median(ratio)),
            "model_to_observed_min":float(ratio.min()),
            "model_to_observed_max":float(ratio.max()),
            "scope":"Auxiliary only. Two recorded caveats: (1) published penumbral ratios additionally divide "
                    "by a limb-darkening factor (I/I0 0.7-0.79) not reproduced here, so published figures sit "
                    "about 1/0.75 above these raw ratios; (2) the model is numerically zero through the "
                    "Hartley-Huggins-saturated blue end while the observation carries scattered light, so "
                    "ratios below ~400 nm are meaningless and the min statistics reflect that region. "
                    "ody9a2020 flattens beyond ~530 nm, consistent with the saturation the paper reports for "
                    "its e=0.83 deg spectrum above 4500 A; treat its red end as contaminated."}

    report = {
        "event":"2019-01-21",
        "scope":"First spaceborne-receiver independent event; fixed atmosphere; conditional "
                "umbral/full-moon comparison at the published mean linear-track positions, not "
                "complete color or turquoise-band validation",
        "normalization":"In-eclipse counts/s divided by the time-weighted average of the twelve "
                        "full-moon G430L counts/s spectra; no background subtraction, no absolute flux "
                        "calibration, no re-exposure correction (single exposures and CR-rejected sets "
                        "are already per second); Sun-site range change within the day verified <0.1% "
                        "and applied",
        "wavelength_window_nm":[MODEL_TABLE_MIN_NM,float(wavelengths.max())],
        "umbral_quantitative_window_nm":[UMBRAL_QUANTITATIVE_MIN_NM,float(wavelengths.max())],
        "umbral_blue_check_window_nm":[MODEL_TABLE_MIN_NM,UMBRAL_BLUE_CHECK_MAX_NM],
        "extraction_method":"Per-column median across the filled 52"+"'"+'" slit of _crj (CR-SPLIT sets) '
                            'or _flt (single exposures) counts divided by exposure time; DQ exclusion '
                            'mask (lost data | bad detector | aperture feature | saturation) = '
                            f"{DQ_EXCLUDE_MASK}; wavelengths from the _sx1/_x1d solutions, audited "
                            "against the image headers",
        "literature_correction_convention":{
            "paper":"Youngblood et al. 2020 Sec. II.2: CalSTIS flt -> L.A.Cosmic -> x2d -> median across "
                    "the cross-dispersion direction excluding flagged pixels; in-eclipse spectrum divided "
                    "by an average of out-of-eclipse spectra; no absolute flux calibration for the "
                    "transmission spectra",
            "reproduced_here":"Same median extraction on _crj/_flt without x2d rectification (the filled "
                              "slit spans the full detector height, so rectification cannot change a "
                              "full-height median); CR handling via the pipeline _crj instead of L.A.Cosmic",
            "omitted_with_bounds":["Lommel-Seeliger per-exposure division: F_LS in 0.4975-0.5021 with "
                                   "closely coupled incidence/emission angles, <=1% effect on in/out ratios",
                                   "Penumbral-only limb-darkening division: G430L I/I0 0.7-0.79; published "
                                   "penumbral ratios are therefore about 1/0.75 above the raw ratios here"],
            "published_cross_check":"No machine-readable published spectrum exists; the paper's umbral "
                                    "statements used here are qualitative anchors: no photons <=4200 A, "
                                    "umbral G430L S/N 60-120, and Chappuis-band O3 detection >4000 A",
        },
        "profiles":{name:asdict(profile) if profile else None for name,profile in profiles.items()},
        "geometry_subset":{key:geometry[key] for key in
                           ["vector_construction","longitude_convention","epoch_contract",
                            "spatial_sensitivity","independent_checks","records"]},
        "solar_scale_in_over_out":solar_scale,
        "numerical_checks":{
            "full_moon_geometric_transmission":"verified >= 1-1e-12 for all twelve reference positions",
            "wavelength_alignment_max_nm":{name:datasets[name]["header_vs_1d_max_diff_nm"]
                                           for name in datasets},
            "min_valid_fraction":{name:datasets[name]["min_valid_fraction"] for name in datasets},
        },
        "metrics":metrics,
        "probe_cases":{variant:{dataset:{label:values.tolist() for label,values in cases.items()}
                                for dataset,cases in probes[variant].items()}
                       for variant in probes},
        "probe_scope":"+/-150 km sampled displacement mirrors the December 2011 radial probe; the track-error "
                      "envelope case (up to ~1042 km, Youngblood et al. pointing budget) brackets the "
                      "same-region premise and is recorded, not folded into the headline metrics",
        "limitations":[
            "Pointing wanders within the paper's documented budget (up to 10 arcmin track approximation "
            "error), so in/out sampled lunar regions can differ; lunar albedo cancellation is conditional "
            "on the published mean linear-track coordinates and degrades within the recorded envelope",
            "The model treats the solar disk as uniform; the paper's solar center-to-limb residual applies "
            "to the umbral branch uncorrected (paper: possible factor <=2 on umbral transmission spectra)",
            "No background/scattered-light subtraction; STIS CCD long-wavelength scattered light adds a "
            "broad pedestal that only partially cancels in the in/out ratio, as in the paper",
            "G430L ends at 570 nm and the model tables start at 300 nm: no red-end Chappuis core and no "
            "289-300 nm edge; turquoise-band color claims are outside this source's reach",
            "Auxiliary penumbral branch omits the published limb-darkening division; see its scope field",
        ],
        "production_module_sha256":{name:hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
                                    for name in ["aerosol.py","atmosphere.py","rays.py","transport.py","spectra.py"]},
    }
    records = []
    for index,wavelength in enumerate(wavelengths):
        record = {"wavelength_nm":float(wavelength),
                  "fullmoon_rate_avg_cps":float(reference[index])}
        for dataset in in_eclipse_datasets:
            record[f"{dataset}_rate_cps"] = float(datasets[dataset]["counts_per_s"][window][index])
            record[f"{dataset}_obs_ratio"] = float(observed[dataset][index])
            if dataset in umbral_datasets:
                for variant in profiles:
                    record[f"{dataset}_{variant}_prediction"] = float(predictions[(variant,dataset)][index])
                    record[f"{dataset}_{variant}_model_to_observed"] = float(
                        predictions[(variant,dataset)][index]/observed[dataset][index])
        records.append(record)
    report["records"] = records
    return report,datasets,window,observed,predictions,probes


def write_analysis(out: Path):
    report,datasets,window,observed,predictions,probes = run_analysis()
    out.mkdir(parents=True,exist_ok=True)
    (out/"comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    with (out/"comparison.csv").open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(report["records"][0]),lineterminator="\n")
        writer.writeheader();writer.writerows(report["records"])
    with (DATA/"spectra.csv").open("w",newline="") as stream:
        fields = ["wavelength_nm"]+[f"{name}_rate_cps" for name in datasets]
        rows = []
        for index in np.where(window)[0]:
            row = {"wavelength_nm":float(datasets["ody904080"]["wavelength_nm"][index])}
            for name,data in datasets.items():
                row[f"{name}_rate_cps"] = float(data["counts_per_s"][index])
            rows.append(row)
        writer = csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    wavelengths = np.array([record["wavelength_nm"] for record in report["records"]])
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False}):
        fig,axes = plt.subplots(2,2,figsize=(11,7.3),gridspec_kw={"hspace":.35})
        for dataset,color in [("ody902030","#a66b44"),("ody902040","#4a6c96")]:
            axes[0,0].plot(wavelengths,observed[dataset],"-",color=color,alpha=.75,
                           label=f"observed umbral/full-moon {dataset}")
            for variant,style in [("molecular","--"),("fixed_opac","-")]:
                axes[0,0].plot(wavelengths,predictions[(variant,dataset)],style,color=color,lw=1.2,
                               label=f"{variant} model {dataset}" if variant=="fixed_opac" else None)
        axes[0,0].set_yscale("log");axes[0,0].set_ylim(1e-9,1)
        axes[0,0].set_ylabel("Umbral / full-moon spectrum");axes[0,0].legend(frameon=False,fontsize=7)
        axes[0,0].set_title("Umbral comparison, 300-570 nm window",loc="left")
        for variant,color in [("molecular","#a66b44"),("fixed_opac","#227580")]:
            for dataset in ("ody902030","ody902040"):
                ratio = np.array([record[f"{dataset}_{variant}_model_to_observed"]
                                  for record in report["records"]])
                axes[1,0].plot(wavelengths,ratio,"-",color=color,alpha=.8,
                               label=f"{variant} {dataset}" if dataset=="ody902030" else None)
        axes[1,0].axhline(1,color="#253d4a",lw=.8)
        axes[1,0].set_ylabel("Model / observed");axes[1,0].set_xlabel("Wavelength (nm)")
        axes[1,0].legend(frameon=False,fontsize=7);axes[1,0].set_ylim(0.05,20)
        for dataset,color in [("ody902030","#a66b44"),("ody9a2010","#ad8540")]:
            axes[0,1].plot(wavelengths,np.array([record[f"{dataset}_rate_cps"]
                          for record in report["records"]]),color=color,alpha=.8,label=dataset)
        axes[0,1].plot(wavelengths,np.array([record["fullmoon_rate_avg_cps"]
                       for record in report["records"]]),color="#253d4a",label="full-moon average")
        axes[0,1].set_yscale("log");axes[0,1].set_ylabel("Counts / s")
        axes[0,1].legend(frameon=False,fontsize=7);axes[0,1].set_title("Extracted slit spectra",loc="left")
        for dataset,color in [("ody9a2010","#ad8540"),("ody9a2020","#a64d4d")]:
            axes[1,1].plot(wavelengths,observed[dataset],color=color,alpha=.8,
                           label=f"penumbral (aux) {dataset}")
            axes[1,1].plot(wavelengths,predictions[("fixed_opac",dataset)],"--",color=color,lw=1.1)
        axes[1,1].set_yscale("log");axes[1,1].set_ylim(1e-5,1)
        axes[1,1].set_ylabel("Penumbral / full-moon (raw)");axes[1,1].set_xlabel("Wavelength (nm)")
        axes[1,1].legend(frameon=False,fontsize=7)
        axes[1,1].set_title("Auxiliary penumbral branch (no CLV division)",loc="left")
        for ax in axes.flat: ax.grid(alpha=.18)
        fig.suptitle("2019-01-21 HST/STIS G430L / fixed-atmosphere comparison",fontsize=12)
        fig.savefig(out/"comparison.png",dpi=160,bbox_inches="tight");plt.close(fig)
    for dataset in ("ody902030","ody902040"):
        for variant in ("molecular","fixed_opac"):
            metric = report["metrics"][variant][dataset]
            print(f"{dataset} {variant}: model/observed median "
                  f"{metric['model_to_observed_median']:.3g} "
                  f"over 430-570 nm; blue-end observed max {metric['blue_end_observed_max']:.2e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("outputs/hst2019"))
    write_analysis(parser.parse_args().out)
