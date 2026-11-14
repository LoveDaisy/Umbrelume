"""Compare the unchanged first model with digitized June 2011 observations.

No scale factor, spatial offset or physical parameter is fitted. Missing filter
responses make effective-wavelength predictions provisional, not a formal test
with observational error bars. Numerical completion is independent of agreement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import eclipse
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.report import BASE_RC
from eclipse.transport import Geometry, transmission_profile

DATA = Path(__file__).resolve().parent/"observations"/"ugolnikov2011"


def load_observations(data_dir=DATA):
    with (data_dir/"points.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    return [{key: float(value) for key,value in row.items()} for row in rows]


def compute_comparison(data_dir=DATA):
    observations = load_observations(data_dir)
    provenance = json.loads((data_dir/"source.json").read_text())
    distances = json.loads((data_dir/"geometry.json").read_text())
    geometry = Geometry(moon_distance_m=distances["moon_distance_m"],
                        sun_distance_m=distances["sun_distance_m"])
    angles = np.array([row["shadow_angle_deg"] for row in observations])
    wavelength = np.array([row["wavelength_nm"] for row in observations])
    unique_wavelengths, select = np.unique(wavelength, return_inverse=True)
    table = build_ray_table(Atmosphere())

    def prediction(ray_table, geom):
        radius = geom.moon_distance_m*np.tan(np.deg2rad(angles))
        values = transmission_profile(ray_table, radius, unique_wavelengths, geometry=geom)
        return values[np.arange(len(angles)), select]

    event = prediction(table, geometry)
    mean_geometry = prediction(table, Geometry())
    ozone250 = prediction(build_ray_table(Atmosphere(ozone_du=250)), geometry)
    refined = prediction(build_ray_table(Atmosphere(), height_step_m=50), geometry)
    event_start = Geometry(moon_distance_m=distances["moon"]["ranges"][0]["distance_m"],
                           sun_distance_m=distances["sun"]["ranges"][0]["distance_m"])
    event_stop = Geometry(moon_distance_m=distances["moon"]["ranges"][-1]["distance_m"],
                          sun_distance_m=distances["sun"]["ranges"][-1]["distance_m"])
    start_prediction, stop_prediction = prediction(table,event_start), prediction(table,event_stop)
    records = []
    for i, row in enumerate(observations):
        obs = row["observed_relative_brightness"]
        records.append({**row,
                        "receiver_radius_m": geometry.moon_distance_m*np.tan(np.deg2rad(angles[i])),
                        "model_event_300du": event[i], "model_mean_geometry_300du": mean_geometry[i],
                        "model_event_250du": ozone250[i], "model_event_300du_50m": refined[i],
                        "model_to_observed_ratio": event[i]/obs,
                        "model_log10_residual": np.log10(event[i]/obs)})
    bands = []
    for band in unique_wavelengths:
        subset = [row for row in records if row["wavelength_nm"] == band]
        ratio = np.array([row["model_to_observed_ratio"] for row in subset])
        bands.append({"effective_wavelength_nm": band,
                      "angle_range_deg": [subset[0]["shadow_angle_deg"], subset[-1]["shadow_angle_deg"]],
                      "number_of_curve_readouts": len(subset),
                      "median_model_to_observed_ratio": float(np.median(ratio)),
                      "range_model_to_observed_ratio": [float(ratio.min()),float(ratio.max())],
                      "innermost_readout": subset[0], "outermost_readout": subset[-1]})
    paired = []
    for angle in sorted(set(angles[select==0]) & set(angles[select==1])):
        blue = next(row for row in records if row["shadow_angle_deg"] == angle and row["wavelength_nm"] == 503)
        red = next(row for row in records if row["shadow_angle_deg"] == angle and row["wavelength_nm"] == 677)
        observed_ratio = blue["observed_relative_brightness"]/red["observed_relative_brightness"]
        model_ratio = blue["model_event_300du"]/red["model_event_300du"]
        paired.append({"shadow_angle_deg": angle, "observed_T503_over_T677": observed_ratio,
                       "model_T503_over_T677": model_ratio,
                       "model_to_observed_color_ratio": model_ratio/observed_ratio})
    dense_angle = np.linspace(0.12,0.60,301)
    dense_radius = geometry.moon_distance_m*np.tan(np.deg2rad(dense_angle))
    dense = transmission_profile(table,dense_radius,unique_wavelengths,geometry=geometry)
    dense250 = transmission_profile(build_ray_table(Atmosphere(ozone_du=250)),dense_radius,
                                    unique_wavelengths,geometry=geometry)
    report = {
        "comparison_status": "Provisional observed-curve comparison; substantial inner-umbra overprediction",
        "source": provenance, "geometry": distances,
        "production_module_sha256": {
            name: hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
            for name in ["atmosphere.py", "spectra.py", "rays.py", "transport.py"]
        },
        "predicted_observable": "E_lambda/F_lambda at stated effective wavelength, monochromatic proxy for unknown bandpass",
        "fixed_baseline": {"ozone_du":300,"density_scale":1,"ray_height_step_m":100,
                           "refraction_reference_nm":550,"aerosols":False,"in_scattered_light":False,
                           "solar_disk":"uniform", "atmosphere":"US1976 + Gaussian ozone"},
        "fitted_parameters": [],
        "bands": bands, "paired_transmission_colors": paired,
        "numerical_and_geometry_sensitivity": {
            "max_relative_prediction_change_height100_to50m":float(np.max(np.abs(event/refined-1))),
            "max_relative_prediction_change_1920_to1940":float(np.max(np.abs(start_prediction/stop_prediction-1))),
            "max_relative_prediction_change_300du_to250du":float(np.max(np.abs(ozone250/event-1))),
            "max_relative_prediction_change_mean_to_event_geometry":float(np.max(np.abs(mean_geometry/event-1))),
        },
        "interpretation": "The fixed clean-atmosphere model is close in scale at the outermost sampled "
                          "western points but increasingly overbright inward. This motivates investigating "
                          "additional low-altitude extinction; it does not uniquely identify aerosols.",
        "limitations": [
            "Readouts from smoothed published curves, not independent original measurements",
            "No instrumental response curves or original pointwise uncertainties were released with this paper",
            "Monochromatic effective-wavelength comparison, not a filter-integrated prediction",
            "Representative geocentric plane distance; no topocentric or lunar-surface reconstruction",
            "Paper's atmospheric comparison used an equatorial temperature profile, not US1976",
            "Two spectral bands do not validate full visible color or the turquoise fringe",
            "Spherical atmosphere cannot predict east-west weather-dependent asymmetry",
        ],
    }
    return report, records, dense_angle, unique_wavelengths, dense, dense250


def write_comparison(output_dir: Path, data_dir=DATA):
    report, records, angles, wavelengths, model, model250 = compute_comparison(data_dir)
    output_dir.mkdir(parents=True,exist_ok=True)
    with (output_dir/"comparison.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=records[0].keys(),lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    (output_dir/"comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    with plt.rc_context(BASE_RC):
        fig,axes=plt.subplots(2,2,figsize=(11.8,8.0),sharex="col",
                              gridspec_kw={"height_ratios":[2.2,1.2],"hspace":.12,"wspace":.23})
        for j,band in enumerate(wavelengths):
            rows=[row for row in records if row["wavelength_nm"]==band]
            x=np.array([row["shadow_angle_deg"] for row in rows])
            observed=np.array([row["observed_relative_brightness"] for row in rows])
            ratio=np.array([row["model_to_observed_ratio"] for row in rows])
            ax=axes[0,j]
            ax.semilogy(angles,model[:,j],color="#ae4d31",lw=2,label="First model / 300 DU ozone")
            ax.semilogy(angles,model250[:,j],color="#ae4d31",ls="--",lw=1.4,
                        label="250 DU sensitivity (not fitted)")
            ax.semilogy(x,observed,"o-",color="#283e4a",ms=4,lw=1.3,
                        label="Published observed curve (digitized)")
            ax.set_title(f"{band:g} nm effective wavelength / western umbra",loc="left",fontsize=12,pad=10)
            ax.grid(alpha=.18); ax.legend(frameon=False,fontsize=8,loc="upper left")
            axes[1,j].plot(x,ratio,"o-",color="#ae4d31",ms=4,lw=1.6)
            axes[1,j].plot(x,[row["model_event_250du"]/row["observed_relative_brightness"] for row in rows],
                           color="#ae4d31",ls="--",lw=1.2)
            axes[1,j].axhline(1,color="#283e4a",lw=1)
            axes[1,j].grid(alpha=.18)
            axes[1,j].set_xlabel("Angular distance from shadow center (degrees)")
            axes[1,j].set_xlim(.12,.60)
        axes[0,0].set_ylabel("Eclipsed / uneclipsed surface brightness")
        axes[1,0].set_ylabel("Model / published observation")
        fig.suptitle("June 15, 2011 lunar eclipse: first observational comparison",fontsize=17,y=.98)
        fig.text(.085,.045,"No scale, position or physical parameters fitted. Source: Ugolnikov et al. (2011), Fig. 2.",fontsize=9)
        fig.text(.085,.023,"Unknown filter responses and measurement errors: effective-wavelength, graph-readout comparison only.",
                 fontsize=9,color="#535d67")
        fig.subplots_adjust(top=.88,bottom=.13,left=.085,right=.98)
        fig.savefig(output_dir/"comparison.png",dpi=180)
        plt.close(fig)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("outputs/observations2011"))
    args=parser.parse_args()
    report=write_comparison(args.out)
    for band in report["bands"]:
        low,high=band["range_model_to_observed_ratio"]
        print(f"{band['effective_wavelength_nm']:g} nm: model / observed = {low:.2f}–{high:.2f}")
    print("No fitted parameters; digitized curves, unknown observational uncertainties.")
    print(f"Results: {args.out.resolve()}")


if __name__=="__main__":
    main()
