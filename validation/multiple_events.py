"""Fixed-model photopic comparisons; one measured event and one flagged context.

This is a neutral-Moon approximation to reflected photometry, not a surface/BRDF
retrieval. Never pool the 2007 extrapolation with independent observed minima.
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
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere, LAYER_BOUNDARIES_M
from eclipse.rays import build_ray_table
from eclipse.spectra import cie_xyz, solar_irradiance, rayleigh_cross_section
from eclipse.transport import Geometry, MOON_RADIUS_M, ring_fraction, transmission_profile

DATA = Path(__file__).resolve().parent/"observations/hernitschek2008"


def aperture_quadrature(center_m, aperture_radius_m, order=128):
    """Area-average weights of axisymmetric illumination over an offset circle.

    For shadow-centered radius r, covered azimuth is 2*pi*f(r); divide the
    annular area by pi*R**2. No post-hoc normalization is applied. Splitting at
    the containment boundary handles circles that contain the shadow axis.
    """
    if not np.isfinite(center_m) or center_m < 0:
        raise ValueError("Center must be finite and nonnegative")
    if not np.isfinite(aperture_radius_m) or aperture_radius_m <= 0 or order < 2:
        raise ValueError("Positive aperture radius and quadrature order >=2 required")
    lower = max(0., center_m-aperture_radius_m)
    upper = center_m+aperture_radius_m
    edges = sorted(set([lower, upper, max(lower, aperture_radius_m-center_m)]))
    nodes, weights = leggauss(order)
    radii, areas = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        # A sin² change of variable removes the square-root endpoint behavior
        # of partial-circle azimuth coverage without renormalizing the area.
        theta = (nodes+1)*np.pi/4
        r = lo+(hi-lo)*np.sin(theta)**2
        dr = 2*(hi-lo)*np.sin(theta)*np.cos(theta)*np.pi/4*weights
        radii.extend(r)
        areas.extend(2*r*dr/aperture_radius_m**2*ring_fraction(center_m,r,aperture_radius_m))
    return np.array(radii), np.array(areas)


def photopic_response(wavelength_nm, observer_rayleigh=False):
    """Solar-weighted nominal V(lambda); optional *one* observer airmass.

    Rayleigh-only observer weighting is a sensitivity probe, not a full model
    of the published instruments, ground atmosphere, or lunar reflectance.
    """
    response = solar_irradiance(wavelength_nm)*cie_xyz(wavelength_nm)[:,1]
    if observer_rayleigh:
        atmosphere = Atmosphere()
        edges = np.unique(np.r_[0., LAYER_BOUNDARIES_M, atmosphere.top_m])
        edges = edges[edges <= atmosphere.top_m]
        column = sum(quad(atmosphere.number_density,lo,hi,epsabs=0,epsrel=1e-9)[0]
                     for lo,hi in zip(edges[:-1],edges[1:]))
        response *= np.exp(-column*rayleigh_cross_section(wavelength_nm))
    return response


def predict(table, geometry, center, radius, order, wavelength_step):
    wavelengths = np.arange(380.,780.+wavelength_step/2,wavelength_step)
    samples, area = aperture_quadrature(center,radius,order)
    transmission = transmission_profile(table,samples,wavelengths,geometry=geometry)
    mean_spectrum = area @ transmission
    values = {}
    for name, ground in [("nominal_v",False),("one_airmass_rayleigh_proxy",True)]:
        response = photopic_response(wavelengths,ground)
        values[name] = float(np.trapezoid(mean_spectrum*response,wavelengths)
                             / np.trapezoid(response,wavelengths))
    return values, float(area.sum())


def run_analysis():
    with (DATA/"points.csv").open() as stream:
        observations = list(csv.DictReader(stream))
    geometry_data = json.loads((DATA/"geometry.json").read_text())
    geometries = {row["event_date"]:row for row in geometry_data["events"]}
    profiles = {"molecular":None,"fixed_opac":AerosolProfile.opac_free_troposphere()}
    tables = {name:build_ray_table(Atmosphere(),aerosol_profile=profile)
              for name,profile in profiles.items()}
    records, checks = [], []
    for obs in observations:
        event = obs["event_date"]
        distances = geometries[event]
        geometry = Geometry(moon_distance_m=distances["moon_distance_m"],
                            sun_distance_m=distances["sun_distance_m"])
        center = geometry.moon_distance_m*np.tan(np.deg2rad(float(obs["geocentric_minimum_phase_angle_arcmin"])/60))
        # Published field is a diameter. Whole-disk measurement was corrected
        # for sky/aperture dilution, so use the physical lunar disk for 1990.
        radius = (MOON_RADIUS_M if obs["aperture_scope"]=="whole_lunar_disk" else
                  geometry.moon_distance_m*np.tan(np.deg2rad(float(obs["field_diameter_arcmin"])/120)))
        full, full_err = float(obs["full_luminance_cd_m2"]), float(obs["full_luminance_reported_error_cd_m2"])
        low, low_err = float(obs["minimum_luminance_cd_m2"]), float(obs["minimum_luminance_reported_error_cd_m2"])
        observed = low/full
        for variant, table in tables.items():
            coarse, area = predict(table,geometry,center,radius,128,5.)
            refined, refined_area = predict(table,geometry,center,radius,256,2.5)
            checks.append({"event_date":event,"variant":variant,
                           "area_weight_sum_order128":area,"area_weight_sum_order256":refined_area,
                           "max_relative_change_order128_5nm_to_order256_2p5nm":
                               max(abs(coarse[key]/refined[key]-1) for key in coarse)})
            for response, prediction in refined.items():
                records.append({"event_date":event,"comparison_role":obs["comparison_role"],
                                "minimum_status":obs["minimum_status"],"aperture_scope":obs["aperture_scope"],
                                "aperture_radius_m":radius,"shadow_center_offset_m":center,
                                "observed_minimum_to_reference":observed,
                                "reported_value_corner_lower":(low-low_err)/(full+full_err),
                                "reported_value_corner_upper":(low+low_err)/(full-full_err),
                                "variant":variant,"response":response,"prediction":prediction,
                                "model_to_observed":prediction/observed})
    if max(abs(check["area_weight_sum_order256"]-1) for check in checks)>1e-6:
        raise RuntimeError("Aperture area integration did not close")
    if max(check["max_relative_change_order128_5nm_to_order256_2p5nm"] for check in checks)>.01:
        raise RuntimeError("Photometry refinement changed predictions by >1%")
    report = {"scope":"Fixed model, no fitted atmospheric parameter or brightness scaling; "
                      "neutral uniform lunar screen and approximate minimum epoch.",
              "independent_new_measured_events":["1990-02-09"],
              "context_only_extrapolated_minima":["2007-03-03"],
              "profiles":{name:asdict(profile) if profile else None for name,profile in profiles.items()},
              "geometry":geometry_data,"numerical_checks":checks,
              "response_contract":"380-780nm solar-weighted CIE y; second response retains one "
                                  "vertical Rayleigh airmass as sensitivity only. Neither includes "
                                  "actual instrument response, lunar spectral albedo/BRDF or full ground extinction.",
              "reported_bounds_contract":"Extremal ratios of printed plus/minus values, not a confidence "
                                         "interval; no covariance or complete systematics available.",
              "production_module_sha256":{name:hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
                                          for name in ["aerosol.py","atmosphere.py","rays.py","transport.py","spectra.py"]},
              "records":records}
    return report


def write_analysis(out):
    report = run_analysis()
    out.mkdir(parents=True,exist_ok=True)
    (out/"analysis.json").write_text(json.dumps(report,indent=2)+"\n")
    records = report["records"]
    with (out/"comparison.csv").open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(records[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(records)
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False}):
        fig,axes = plt.subplots(1,2,figsize=(10.5,5),sharey=True)
        for ax,event in zip(axes,["1990-02-09","2007-03-03"]):
            rows = [r for r in records if r["event_date"]==event]
            first = rows[0]
            observed = first["observed_minimum_to_reference"]
            ax.axhline(observed,color="#253d4a",label="Published minimum/reference")
            ax.axhspan(first["reported_value_corner_lower"],first["reported_value_corner_upper"],
                       color="#253d4a",alpha=.12,label="Printed-value bounds (not CI)")
            for i,variant in enumerate(["molecular","fixed_opac"]):
                subset = [r for r in rows if r["variant"]==variant]
                for j,row in enumerate(subset):
                    ax.plot(i+(j-.5)*.13,row["prediction"],"o" if j==0 else "s",
                            color="#a66b44" if i==0 else "#227580",
                            ms=6,mfc="none" if event=="2007-03-03" else None,
                            label=("Nominal V(lambda)" if j==0 else "+ one Rayleigh airmass (proxy)") if i==0 else None)
            ax.set_xticks([0,1],["Molecular v1","Fixed OPAC prior"])
            ax.set_yscale("log");ax.grid(axis="y",alpha=.18)
            ax.set_title(event+"\n"+("Observed whole disk" if event=="1990-02-09" else "Extrapolated central aperture / context only"),fontsize=10)
        axes[0].set_ylabel("Minimum / uneclipsed photopic brightness")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,frameon=False,fontsize=8,loc="lower center",ncol=2)
        fig.suptitle("Fixed atmosphere across additional eclipse events; neutral-Moon approximation",fontsize=11)
        fig.tight_layout(rect=(0,.11,1,1));fig.savefig(out/"comparison.png",dpi=160);plt.close(fig)
    for row in records:
        print(row["event_date"],row["variant"],row["response"],f"model/observed={row['model_to_observed']:.3f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("outputs/multiple-events"))
    write_analysis(parser.parse_args().out)
