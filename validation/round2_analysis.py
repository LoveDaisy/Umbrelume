"""Discriminate missing physics and compare one aerosol-extinction increment.

The OPAC prediction is fixed independently of the eclipse data. A separate
single-parameter AOD calibration uses only 677 nm, with 503 nm withheld. Neither
curve readouts nor their graph readout allowances are independent error bars.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import roots_jacobi
from scipy.integrate import quad

import eclipse
from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.spectra import rayleigh_cross_section, ozone_cross_section
from eclipse.report import BASE_RC
from eclipse.transport import Geometry, receiver_quadrature, transmission_profile, ring_fraction
from .compare_observations import DATA, load_observations


def fit_red_aod(table, profile, radius_red_m, observed_red, geometry):
    """Calibrate exactly one column parameter, receiving no blue observations."""
    unit_tau = table.aerosol_layer_taus[0]/profile.aod550

    def scaled(aod):
        return replace(table, aerosol_layer_taus=(unit_tau*aod,),
                       aerosol_profiles=(replace(profile,aod550=float(aod)),))

    def objective(aod):
        prediction = transmission_profile(scaled(aod),radius_red_m,[677.],geometry=geometry)[:,0]
        return float(np.mean(np.log(prediction/observed_red)**2))

    optimum = minimize_scalar(objective,bounds=(0.,.1),method="bounded",
                              options={"xatol":1e-9})
    if not optimum.success:
        raise RuntimeError("Red-band AOD calibration failed")
    return scaled(optimum.x), {"aod550":float(optimum.x),"bounds":[0.,.1],
                               "objective_mean_squared_ln_ratio":float(optimum.fun)}


def limb_darkened_probe(table, radius, wavelengths, geometry, order):
    """Empirical solar source diagnostic using a positive uniform-disk mixture.

    Hestroffer & Magnan1998 Table2 (NL). I(mu)=mu**alpha; whole-disk flux
    normalized to the same AM0. The mixture variable t=mu^2 has Beta(alpha/2,2)
    density, with disk radius R_s*sqrt(1-t). No production source change is made.
    """
    alpha = np.array([np.interp(w,[492.905,519.930,669.400,700.875],
                                [.570,.538,.407,.386]) for w in wavelengths])
    output = np.zeros((len(radius),len(wavelengths)))
    for j, power in enumerate(alpha):
        nodes,weights = roots_jacobi(order,1,power/2-1)
        weights = weights/weights.sum()
        for node,weight in zip(nodes,weights):
            source = replace(geometry,sun_radius_m=geometry.sun_radius_m*np.sqrt((1-node)/2))
            output[:,j] += weight*transmission_profile(table,radius,[wavelengths[j]],geometry=source)[:,0]
    return output


def altitude_quantiles(table, radius, wavelength, geometry):
    b,area = receiver_quadrature(table,float(radius),geometry)
    tau = (np.interp(b,table.impact_m,table.air_column_m2)*rayleigh_cross_section(wavelength)
           + np.interp(b,table.impact_m,table.ozone_column_m2)*ozone_cross_section(wavelength))
    flux = area*np.exp(-tau)
    heights = np.interp(b,table.impact_m,table.altitude_m)
    cdf = np.cumsum(flux)/flux.sum()
    return np.interp([.1,.5,.9],cdf,heights).tolist()


def source_kernel_check():
    """Verify the source mixture against direct angular integration on a ring."""
    records=[]
    for power in [.5580466234967622,.4019293089753773]:
        exponent=power/2
        for r,s in [(0.,.8),(.2,.3),(.6,.6),(.9,.6),(1.4,.8)]:
            if r*s==0:
                direct=(exponent+1)*max(1-(r-s)**2,0)**exponent
            else:
                limit=(np.pi if r+s<=1 else
                       np.arccos(np.clip((r*r+s*s-1)/(2*r*s),-1,1)))
                direct=(exponent+1)/np.pi*quad(
                    lambda phi:max(1-r*r-s*s+2*r*s*np.cos(phi),0)**exponent,
                    0,limit,epsabs=1e-11,limit=200)[0]
            points=[(1-q*q)**exponent for q in [abs(r-s),r+s] if 0<q<1]
            mixture=(exponent+1)*quad(
                lambda z:float(ring_fraction(r,s,np.sqrt(max(1-z**(1/exponent),1e-30)))),
                0,1,points=points,epsabs=1e-11,limit=200)[0]
            records.append({"alpha":power,"receiver_radius":r,"ring_radius":s,
                            "direct_angular":direct,"positive_disk_mixture":mixture,
                            "absolute_difference":abs(direct-mixture)})
    error=max(row["absolute_difference"] for row in records)
    if error>1e-7:
        raise RuntimeError("Source-mixture diagnostic failed independent angular check")
    return {"pass":True,"max_absolute_difference":error,"cases":records}


def run_analysis():
    observations = load_observations()
    geometry_data = json.loads((DATA/"geometry.json").read_text())
    geometry = Geometry(moon_distance_m=geometry_data["moon_distance_m"],
                        sun_distance_m=geometry_data["sun_distance_m"])
    angles = np.array([row["shadow_angle_deg"] for row in observations])
    wl = np.array([row["wavelength_nm"] for row in observations])
    obs = np.array([row["observed_relative_brightness"] for row in observations])
    wavelengths,selection = np.unique(wl,return_inverse=True)
    radius = geometry.moon_distance_m*np.tan(np.deg2rad(angles))
    molecular = build_ray_table(Atmosphere())
    prior = AerosolProfile.opac_free_troposphere()
    aerosol_table = build_ray_table(Atmosphere(),aerosol_profile=prior)
    no_aerosol = transmission_profile(molecular,radius,wavelengths,geometry=geometry)
    with_aerosol = transmission_profile(aerosol_table,radius,wavelengths,geometry=geometry)
    red = wl==677
    calibrated_table,fit = fit_red_aod(aerosol_table,prior,radius[red],obs[red],geometry)
    calibrated = transmission_profile(calibrated_table,radius,wavelengths,geometry=geometry)
    ld64 = limb_darkened_probe(molecular,radius,wavelengths,geometry,64)
    ld128 = limb_darkened_probe(molecular,radius,wavelengths,geometry,128)
    density_proxy = transmission_profile(build_ray_table(Atmosphere(density_scale=.9605)),
                                         radius,wavelengths,geometry=geometry)
    refined_table = build_ray_table(Atmosphere(),aerosol_profile=prior,height_step_m=50)
    refined = transmission_profile(refined_table,radius,wavelengths,geometry=geometry)
    records=[]
    variants={"molecular":no_aerosol,"fixed_opac":with_aerosol,"red_calibrated":calibrated,
              "limb_darkened_diagnostic":ld128,"density_proxy":density_proxy}
    for i,row in enumerate(observations):
        record={"wavelength_nm":wl[i],"shadow_angle_deg":angles[i],
                "observed_curve_brightness":obs[i]}
        for name,value in variants.items():
            record[name+"_prediction"]=float(value[i,selection[i]])
            record[name+"_model_to_observed"]=float(value[i,selection[i]]/obs[i])
        record["molecular_tangent_height_p10_p50_p90_m"] = altitude_quantiles(
            molecular,radius[i],wl[i],geometry)
        records.append(record)
    scores={}
    for name,values in variants.items():
        diagonal=values[np.arange(len(obs)),selection]
        scores[name]={"all_readouts_rms_ln_ratio":float(np.sqrt(np.mean(np.log(diagonal/obs)**2)))}
        for band in wavelengths:
            mask=wl==band
            scores[name][f"{band:g}nm_rms_ln_ratio"] = float(np.sqrt(np.mean(np.log(diagonal[mask]/obs[mask])**2)))
    report={
        "selected_increment":"Free-tropospheric aerosol extinction on the existing curved paths",
        "identification_scope":"Largest tested omitted effect with the required altitude dependence; "
                               "not a unique retrieval of the real2011aerosol/cloud profile",
        "fixed_prior":asdict(prior),
        "prior_source":"Hess,Koepke&Schult1998 OPAC Tables3,5; independent climatological/model prior, not measured2011profile",
        "source_url":"https://lweb.cfa.harvard.edu/HITRAN/HITRAN2012/Aerosols/papers/hess_opac_1998.pdf",
        "calibration":{**fit,"fitted_parameters":["aod550"],"fit_band_nm":677,
                       "fit_readout_count":int(red.sum()),"withheld_band_nm":503,
                       "withheld_readout_count":int((~red).sum()),"default_uses_calibration":False,
                       "scope":"Same-event cross-band check, not independent observational validation. "
                               "Effective AOD can absorb unmodeled source/profile/cloud effects."},
        "geometry":geometry_data,
        "descriptive_scores":scores,
        "score_scope":"RMS natural-log ratio of correlated published-curve readouts; not chi-square or statistical errors",
        "limb_darkening_diagnostic":{
            "source":"Hestroffer&Magnan1998 A&A333338 Table2 Neckel-Labs fits",
            "url":"https://www.physics.hmc.edu/faculty/esin/a101/labs/limbdark.pdf",
            "tabulated_rows_nm_alpha":[[492.905,.570],[519.930,.538],[669.400,.407],[700.875,.386]],
            "normalization":"I(mu)=mu**alpha; full-disk average2/(alpha+2); equal integrated flux",
            "method":"Normalized positive mixture of uniform disks, t=mu²~Beta(alpha/2,2)",
            "independent_ring_kernel_check":source_kernel_check(),
            "max_relative_difference_order64_to128":float(np.max(np.abs(ld64/ld128-1))),
            "limitation":"Fit calibrated mainly within solar r<=0.9; extreme limb extrapolated. Diagnostic only."},
        "density_proxy_scope":"Uniform density_scale=.9605 is a sensitivity proxy, not a tropical temperature profile",
        "fixed_prior_height100_to50m_max_relative_change":float(np.max(np.abs(with_aerosol/refined-1))),
        "production_module_sha256":{name:hashlib.sha256((Path(eclipse.__file__).parent/name).read_bytes()).hexdigest()
                                    for name in ["aerosol.py","atmosphere.py","rays.py","transport.py","spectra.py"]},
        "observational_limitations":["Unknown filter responses and measurement uncertainties",
                                     "Digitized smoothed curves, not independent exposures",
                                     "Only two visible bands, no direct turquoise-fringe validation"],
        "records":records,
    }
    dense_angles=np.linspace(.12,.60,301)
    dense_radius=geometry.moon_distance_m*np.tan(np.deg2rad(dense_angles))
    dense={"molecular":transmission_profile(molecular,dense_radius,wavelengths,geometry=geometry),
           "fixed_opac":transmission_profile(aerosol_table,dense_radius,wavelengths,geometry=geometry),
           "red_calibrated":transmission_profile(calibrated_table,dense_radius,wavelengths,geometry=geometry)}
    return report,records,dense_angles,wavelengths,dense


def write_analysis(out):
    report,records,angles,wavelengths,dense=run_analysis()
    out.mkdir(parents=True,exist_ok=True)
    (out/"analysis.json").write_text(json.dumps(report,indent=2)+"\n")
    fields=[key for key in records[0] if not isinstance(records[0][key],list)]
    with (out/"comparison.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n",extrasaction="ignore")
        writer.writeheader();writer.writerows(records)
    styles={"molecular":("#a7aeb4","--","Molecular v1"),
            "fixed_opac":("#ac4f31","-","Fixed OPAC prior / AOD550=0.013"),
            "red_calibrated":("#227580",":",f"AOD calibrated on677nm only / {report['calibration']['aod550']:.4f}")}
    with plt.rc_context(BASE_RC):
        fig,axes=plt.subplots(2,2,figsize=(11.8,8),sharex="col",gridspec_kw={"height_ratios":[2.2,1.2],"hspace":.12})
        for j,band in enumerate(wavelengths):
            subset=[row for row in records if row["wavelength_nm"]==band]
            x=[row["shadow_angle_deg"] for row in subset]
            y=[row["observed_curve_brightness"] for row in subset]
            axes[0,j].semilogy(x,y,"o-",color="#253d4a",ms=4,label="Published observed curve (digitized)")
            for name,(color,style,label) in styles.items():
                axes[0,j].semilogy(angles,dense[name][:,j],color=color,ls=style,lw=1.8,label=label)
                axes[1,j].plot(x,[row[name+"_model_to_observed"] for row in subset],color=color,ls=style,lw=1.8)
            axes[0,j].set_title(f"{band:g}nm / "+("withheld from calibration" if band==503 else "calibration band"),loc="left")
            axes[0,j].legend(frameon=False,fontsize=7.5,loc="upper left")
            axes[1,j].axhline(1,color="#253d4a",lw=.8)
            axes[1,j].set_xlim(.12,.60)
            axes[1,j].set_xlabel("Angular distance from shadow center (degrees)")
            for ax in axes[:,j]:ax.grid(alpha=.18)
        axes[0,0].set_ylabel("Eclipsed / uneclipsed surface brightness")
        axes[1,0].set_ylabel("Model / observed curve")
        fig.suptitle("One increment: free-tropospheric aerosol extinction",fontsize=17,y=.98)
        fig.text(.085,.045,"Solid: independent fixed prior. Dotted: same-event red-band AOD calibration; never used as default.",fontsize=9)
        fig.text(.085,.022,"Unknown filter responses and measurement errors. Two-band readouts do not validate full visible color.",fontsize=9,color="#535d67")
        fig.subplots_adjust(top=.88,bottom=.13,left=.085,right=.98,wspace=.23)
        fig.savefig(out/"comparison.png",dpi=180);plt.close(fig)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("outputs/round2"))
    args=parser.parse_args()
    report=write_analysis(args.out)
    print(json.dumps({"fixed_prior":report["fixed_prior"],"calibration":report["calibration"],
                      "scores":report["descriptive_scores"]},indent=2))


if __name__=="__main__":main()
