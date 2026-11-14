"""Time-average the fixed model over Yan et al.'s full U1 exposure.

Near-Tycho nominal ephemeris geometry is independently checked against Fig. 1.
The +/-150km radial probe describes sensitivity, not a pointing confidence band.
Neither atmospheric parameters nor brightness normalization are fitted.
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
from scipy.integrate import simpson

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.report import BASE_RC
from eclipse.transport import Geometry, geometric_transmission, transmission_profile

DATA = Path(__file__).resolve().parent/"observations/yan2011dec"


def exposure_average(values, seconds):
    """Average linear spectral irradiance, before ratios, logarithms or color."""
    t = np.asarray(seconds,dtype=float)
    values = np.asarray(values,dtype=float)
    if t.ndim!=1 or len(t)<3 or np.any(np.diff(t)<=0) or values.shape[0]!=len(t):
        raise ValueError("Ordered exposure samples and matching spectral rows required")
    return simpson(values,x=t,axis=0)/(t[-1]-t[0])


def predict_samples(table, rows, wavelengths, radial_offset_m=0., axial_offset_m=0.):
    values = []
    for row in rows:
        geometry = Geometry(moon_distance_m=row["receiver_axial_distance_m"]+axial_offset_m,
                            sun_distance_m=row["sun_earth_distance_m"])
        radius = row["shadow_radius_m"]+radial_offset_m
        # Same-region U1/B2 cancels F_lambda and lunar reflectance only under
        # the paper's assumptions; retain the small Sun-site range difference.
        values.append(transmission_profile(table,[radius],wavelengths,geometry=geometry)[0]
                      * row["solar_flux_relative_to_b2"])
    return np.array(values)


def run_analysis():
    geometry_data = json.loads((DATA/"geometry.json").read_text())
    with (DATA/"points.csv").open() as stream:
        observations = [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]
    wavelengths = np.array([row["wavelength_nm"] for row in observations])
    observed = np.array([row["observed_relative_brightness"] for row in observations])
    rows = geometry_data["records"]
    seconds = np.array([row["seconds_since_start"] for row in rows])
    if seconds[0]!=0 or seconds[-1]!=1800:
        raise ValueError("Comparison must cover the entire U1 exposure")
    offsets = geometry_data["spatial_sensitivity_radial_offsets_m"]
    profiles = {"molecular":None,"fixed_opac":AerosolProfile.opac_free_troposphere()}
    averages,instantaneous,checks,descriptive = {},{},{},{}
    b2 = geometry_data["bright_reference"]
    bright_geometry = Geometry(moon_distance_m=b2["receiver_axial_distance_m"],sun_distance_m=b2["sun_earth_distance_m"])
    if geometric_transmission([b2["shadow_radius_m"]],bright_geometry)[0] < 1-1e-12:
        raise RuntimeError("B2 reference is not fully outside the geometric shadow")
    for variant,profile in profiles.items():
        table = build_ray_table(Atmosphere(),aerosol_profile=profile)
        fine_table = build_ray_table(Atmosphere(),aerosol_profile=profile,height_step_m=50)
        cases = {}
        for offset in offsets:
            spectra = predict_samples(table,rows,wavelengths,offset)
            cases[str(int(offset))] = exposure_average(spectra,seconds)
            if offset==0:
                instantaneous[variant] = spectra
        nominal = cases["0"]
        averages[variant] = cases
        coarse_time = exposure_average(instantaneous[variant][::2],seconds[::2])
        fine_height = exposure_average(predict_samples(fine_table,rows,wavelengths),seconds)
        midpoint = instantaneous[variant][len(rows)//2]
        simultaneous_rows = [{**row,"shadow_radius_m":row["simultaneous_geometry_radius_m"]} for row in rows]
        # This diagnostic changes only rho: it isolates the radius sensitivity
        # to dropping light-time, rather than claiming a complete second geometry.
        simultaneous = exposure_average(predict_samples(table,simultaneous_rows,wavelengths),seconds)
        depth_probe = [exposure_average(predict_samples(table,rows,wavelengths,axial_offset_m=delta),seconds)
                       for delta in [-150_000.,150_000.]]
        checks[variant] = {
            "max_relative_change_time120s_to60s":float(np.max(np.abs(coarse_time/nominal-1))),
            "max_relative_change_height100m_to50m":float(np.max(np.abs(fine_height/nominal-1))),
            "midpoint_to_exposure_average_ratio":(midpoint/nominal).tolist(),
            "simultaneous_radius_only_to_retarded_ratio":(simultaneous/nominal).tolist(),
            "max_relative_change_axial_distance_plus_minus150km":float(np.max(np.abs(np.array(depth_probe)/nominal-1))),
            "bright_reference_geometric_transmission":1.,
        }
        if checks[variant]["max_relative_change_time120s_to60s"]>.001:
            raise RuntimeError("Time averaging is not converged within0.1%")
        if checks[variant]["max_relative_change_height100m_to50m"]>.001:
            raise RuntimeError("Ray-height resolution is not converged within0.1%")
        ratio = nominal/observed
        descriptive[variant] = {
            "model_to_observed_min_max":[float(ratio.min()),float(ratio.max())],
            "readout_rms_ln_model_to_observed":float(np.sqrt(np.mean(np.log(ratio)**2))),
            "450_to750_transmission_ratio_model_over_observed":float((nominal[0]/nominal[-1])/(observed[0]/observed[-1])),
            "500_to750_transmission_ratio_model_over_observed":float((nominal[1]/nominal[-1])/(observed[1]/observed[-1])),
            "450_to750_shape_ratio_in_same_offset_probes":[float((case[0]/case[-1])/(observed[0]/observed[-1])) for case in cases.values()],
            "smallest_model_to_observed_across_radial_probes":float(np.min(np.array(list(cases.values()))/observed)),
        }
    records = []
    for i,obs in enumerate(observations):
        record = {"wavelength_nm":wavelengths[i],"observed_curve_ratio":observed[i],
                  "graph_readout_log10_halfwidth":obs["graph_readout_log10_halfwidth"]}
        for variant,cases in averages.items():
            nominal = cases["0"][i]
            record[variant+"_prediction"] = float(nominal)
            record[variant+"_model_to_observed"] = float(nominal/observed[i])
            record[variant+"_radial_probe_min"] = float(min(case[i] for case in cases.values()))
            record[variant+"_radial_probe_max"] = float(max(case[i] for case in cases.values()))
        records.append(record)
    report = {"event":"2011-12-10","scope":"Independent event; fixed atmosphere; conditional near-Tycho continuum comparison, not exact fibre reconstruction or complete color validation",
              "normalization":"Fig15 already(U1/B2)/300; no repeated exposure correction, fitted scale, geometry shift or atmospheric fit",
              "profiles":{name:asdict(profile) if profile else None for name,profile in profiles.items()},
              "geometry":geometry_data,"temporal_method":"Simpson integration of linear spectral irradiance on31samples across1800s; solar flux scaled to B2 source-site range",
              "numerical_checks":checks,"descriptive_diagnostics":descriptive,
              "radial_probe_cases":{variant:{offset:values.tolist() for offset,values in cases.items()}
                                    for variant,cases in averages.items()},
              "diagnostic_scope":"Correlated graph readouts, not independent exposures; RMSlog ratio is descriptive, not chi-square. Pairwise ratios diagnose continuum shape, not full observed color.",
              "sensitivity_scope":"+/-150km radial displacement is a chosen sensitivity tied to the schematic reading scale, not a confidence interval or complete pointing/systematic uncertainty",
              "limitations":["Exact fibre coordinates and tracking are not published; nominal Tycho center independently agrees with Fig1",
                             "No exposure-level observational uncertainties or covariance",
                             "Same-region lunar albedo and spectrograph response cancellation follow the paper; lunar BRDF and solar limb variation not modeled",
                             "Local airmass mismatch0.04 and time-varying atmosphere assumed small, not independently removed here",
                             "Current model omits NO2, O2/CIA and H2O lines; selected continuum points avoid dominant O2 bands but residual weak features remain",
                             "Seven sparse continuum readouts do not define a reliable measured color or turquoise-fringe width"],
              "production_module_sha256":{name:hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
                                          for name in ["aerosol.py","atmosphere.py","rays.py","transport.py","spectra.py"]},
              "records":records}
    return report,seconds,instantaneous


def write_analysis(out):
    report,seconds,instantaneous = run_analysis()
    out.mkdir(parents=True,exist_ok=True)
    (out/"comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    records = report["records"]
    with (out/"comparison.csv").open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(records[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(records)
    with (out/"time_samples.csv").open("w",newline="") as stream:
        fields = ["observed_utc","seconds_since_start","shadow_radius_m","wavelength_nm","molecular","fixed_opac"]
        writer = csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");writer.writeheader()
        for j,row in enumerate(report["geometry"]["records"]):
            for i,point in enumerate(records):
                writer.writerow({**{key:row[key] for key in fields[:3]},"wavelength_nm":point["wavelength_nm"],
                                 **{key:values[j,i] for key,values in instantaneous.items()}})
    wl = np.array([r["wavelength_nm"] for r in records])
    obs = np.array([r["observed_curve_ratio"] for r in records])
    with plt.rc_context(BASE_RC):
        fig,axes = plt.subplots(2,2,figsize=(11,7.3),gridspec_kw={"height_ratios":[2,1.3],"hspace":.35})
        allowances = np.array([r["graph_readout_log10_halfwidth"] for r in records])
        axes[0,0].errorbar(wl,obs,yerr=[obs-obs*10**-allowances,obs*10**allowances-obs],fmt="o",color="#253d4a",ms=4,label="Observed Fig.15 / graph allowance only")
        for variant,color,label in [("molecular","#a66b44","Molecular v1"),("fixed_opac","#227580","Fixed OPAC prior")]:
            nominal = np.array([r[variant+"_prediction"] for r in records])
            lo = np.array([r[variant+"_radial_probe_min"] for r in records])
            hi = np.array([r[variant+"_radial_probe_max"] for r in records])
            axes[0,0].plot(wl,nominal,".-",color=color,label=label)
            axes[1,0].plot(wl,nominal/obs,".-",color=color)
            if variant=="fixed_opac":
                axes[0,0].fill_between(wl,lo,hi,color=color,alpha=.16,label="+/-150km radial probe (not CI)")
                axes[1,0].fill_between(wl,lo/obs,hi/obs,color=color,alpha=.16)
        axes[0,0].set_yscale("log");axes[0,0].set_ylabel("Eclipsed / bright continuum ratio")
        axes[0,0].legend(frameon=False,fontsize=7.6);axes[0,0].set_title("Full 1800 s exposure average",loc="left")
        axes[1,0].axhline(1,color="#253d4a",lw=.8);axes[1,0].set_ylabel("Model / observed");axes[1,0].set_xlabel("Wavelength (nm)")
        rows = report["geometry"]["records"]
        axes[0,1].plot(seconds/60,np.array([r["shadow_radius_m"] for r in rows])/1000,color="#253d4a",label="Retarded JPL near-Tycho geometry")
        check = report["geometry"]["independent_figure_check"]
        axes[0,1].plot([0,30],np.array(check["plot_start_end_radius_m"])/1000,"o",mfc="none",color="#a66b44",label="Independent Fig.1 schematic endpoints")
        axes[0,1].set_ylabel("Target distance from shadow axis (km)");axes[0,1].set_xlabel("Minutes from U1 start");axes[0,1].legend(frameon=False,fontsize=7.6)
        for index,color in [(0,"#4a6c96"),(3,"#ad8540"),(6,"#a64d4d")]:
            series = instantaneous["fixed_opac"][:,index]
            avg = next(r["fixed_opac_prediction"] for r in records if r["wavelength_nm"]==wl[index])
            axes[1,1].plot(seconds/60,series/avg,color=color,label=f"{wl[index]:g}nm")
        axes[1,1].axhline(1,color="#253d4a",lw=.8);axes[1,1].set_ylabel("Instantaneous / exposure average");axes[1,1].set_xlabel("Minutes from U1 start");axes[1,1].legend(frameon=False,fontsize=7.6,ncol=3)
        for ax in axes.flat: ax.grid(alpha=.18)
        fig.suptitle("2011-12-10 / independent spectral comparison with fixed atmosphere",fontsize=12)
        fig.savefig(out/"comparison.png",dpi=160,bbox_inches="tight");plt.close(fig)
    for row in records:
        print(f"{row['wavelength_nm']:g}nm: molecular/observed={row['molecular_model_to_observed']:.3f}, fixed/observed={row['fixed_opac_model_to_observed']:.3f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("outputs/dec2011"))
    write_analysis(parser.parse_args().out)
