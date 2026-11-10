"""Independent aerosol-path validation by Bessel integrals, Cartesian ODE and flux.

Run ``uv run python -m validation.check_aerosol --out outputs/validation/aerosol.json``.
These checks concern the declared aerosol transport, not the adequacy of its
vertical profile or the fit to any lunar eclipse.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.special import k1e

from eclipse.aerosol import AerosolColumn, AerosolProfile
from eclipse.atmosphere import (
    Atmosphere, EARTH_RADIUS_M, LAYER_BOUNDARIES_M, REFERENCE_NUMBER_DENSITY,
)
from eclipse.rays import build_ray_table
from eclipse.spectra import ozone_cross_section, rayleigh_cross_section, refractive_index_air
from eclipse.transport import Geometry, signed_landing_radius, transmission_profile


def straight_exponential_checks():
    """Check spherical straight rays against an exact modified-Bessel integral.

    For beta(z)=beta0 exp(-z/H), the full infinite straight-ray optical depth
    is 2 beta0 r_t exp(-h/H) k1e(r_t/H). A 120 km upper bound omits less than
    ~1e-9 of the tested columns, avoiding an enormous exp(R/H) intermediate.
    """
    atmosphere = Atmosphere(density_scale=0, ozone_du=0)
    heights = np.array([0., 5000., 11000., 25000.])
    results = []
    for scale in [2000., 5000.]:
        profile = AerosolProfile(aod550=.05, bottom_m=0., top_m=120000.,
                                 scale_height_m=scale)
        table = build_ray_table(atmosphere, aerosol_profile=profile, height_step_m=500,
                                quadrature_order=64)
        beta0 = profile.aod550 / (scale * -np.expm1(-profile.top_m/scale))
        tangent_radius = EARTH_RADIUS_M + heights
        exact = 2*beta0*tangent_radius*np.exp(-heights/scale)*k1e(tangent_radius/scale)
        actual = np.interp(heights, table.altitude_m, table.aerosol_layer_taus[0])
        difference = np.abs(actual/exact - 1)
        results.append({"scale_height_m":scale, "tangent_altitude_m":heights.tolist(),
                        "bessel_tau550":exact.tolist(), "production_tau550":actual.tolist(),
                        "relative_difference":difference.tolist(),
                        "pass":bool(np.all(difference < 3e-7))})
    return {"pass":all(row["pass"] for row in results), "profiles":results}


def cartesian_aerosol_reference(atmosphere, profiles, tangent_altitude_m):
    """Integrate one outward half-ray; split every aerosol support discontinuity.

    The Cartesian ray equation is independent of the production turning-point
    quadrature. Each layer's exponential extinction is evaluated directly from
    its AOD normalization and profile parameters, independently of
    extinction550(), and one ODE state accumulates per layer. Each radial
    segment carries a fixed per-layer active/inactive flag, so an ODE step
    across a cutoff cannot integrate the wrong side of a jump.
    """
    profiles = tuple(profiles)
    if not 0 <= tangent_altitude_m < atmosphere.top_m:
        raise ValueError("Tangent must lie inside atmosphere")
    coefficient = (float(refractive_index_air(550.))-1)/REFERENCE_NUMBER_DENSITY
    beta_bottoms = [p.aod550/(p.scale_height_m * -np.expm1(
        -(p.top_m-p.bottom_m)/p.scale_height_m)) for p in profiles]
    boundaries = np.r_[LAYER_BOUNDARIES_M,
                       [edge for p in profiles for edge in (p.bottom_m, p.top_m)],
                       atmosphere.top_m]
    boundaries = np.unique(boundaries[(boundaries > tangent_altitude_m)
                                     & (boundaries <= atmosphere.top_m)])
    state = np.array([EARTH_RADIUS_M+tangent_altitude_m, 0., 0., 1.,
                      *[0.]*len(profiles)])
    first_impact = ((1+coefficient*float(atmosphere.number_density(tangent_altitude_m)))
                    * state[0])
    invariant_error = 0.
    direction_error = 0.
    steps = 0
    lower = tangent_altitude_m
    for upper in boundaries:
        active = [upper > p.bottom_m and lower < p.top_m for p in profiles]

        def derivative(_, value):
            position, direction = value[:2], value[2:4]
            radius = np.linalg.norm(position)
            altitude = max(float(radius-EARTH_RADIUS_M),0.)
            density = float(atmosphere.number_density(altitude))
            n = 1+coefficient*density
            gradient = (coefficient*float(atmosphere.density_gradient(altitude))/n
                        * position/radius)
            acceleration = gradient-direction*np.dot(direction,gradient)
            extinctions = [beta_bottoms[i]*np.exp(-(altitude-profiles[i].bottom_m)
                                                  /profiles[i].scale_height_m)
                           if active[i] else 0. for i in range(len(profiles))]
            return np.r_[direction,acceleration,extinctions]

        def exit_shell(_, value):
            return np.linalg.norm(value[:2])-(EARTH_RADIUS_M+upper)

        exit_shell.terminal = True
        exit_shell.direction = 1
        solution = solve_ivp(derivative, (0.,4*(EARTH_RADIUS_M+atmosphere.top_m)),state,
                             events=exit_shell,method="DOP853",rtol=1e-10,
                             atol=[1e-6,1e-6,1e-13,1e-13,*([1e-12]*len(profiles))],
                             max_step=3000.)
        if not solution.success or not len(solution.t_events[0]):
            raise RuntimeError("Cartesian aerosol reference failed to exit radial shell")
        states = solution.y
        radii = np.hypot(states[0],states[1])
        n = 1+coefficient*atmosphere.number_density(np.maximum(radii-EARTH_RADIUS_M,0.))
        impact = n*(states[0]*states[3]-states[1]*states[2])
        invariant_error = max(invariant_error,float(np.max(np.abs(impact/first_impact-1))))
        direction_error = max(direction_error,float(np.max(np.abs(np.hypot(states[2],states[3])-1))))
        state = states[:,-1]
        steps += len(solution.t)
        lower = float(upper)
    return {"tau550":float(2*np.sum(state[4:])),
            "bending_rad":float(-2*np.arctan2(state[2],state[3])),
            "max_invariant_relative_error":invariant_error,
            "max_direction_norm_error":direction_error,"accepted_steps":steps}


def curved_aerosol_checks():
    atmosphere = Atmosphere()
    profiles = (AerosolProfile(), AerosolProfile.opac_stratospheric_background())
    column = AerosolColumn(profiles)
    table = build_ray_table(atmosphere,aerosol_profile=column)
    results = []
    for altitude in [1000.,2000.,5000.,11000.,12000.,13000.,25000.,35000.]:
        reference = cartesian_aerosol_reference(atmosphere,profiles,altitude)
        actual = float(np.interp(altitude,table.altitude_m,
                                 np.sum(table.aerosol_layer_taus,axis=0)))
        difference = abs(actual/reference["tau550"]-1) if reference["tau550"] else abs(actual)
        results.append({"tangent_altitude_m":altitude,"production_tau550":actual,
                        "cartesian_reference":reference,"relative_or_zero_absolute_difference":difference,
                        "pass":bool(difference < 2e-6 and reference["max_invariant_relative_error"] < 2e-8)})
    return {"pass":all(row["pass"] for row in results),"profile":asdict(column),"rays":results}


def aerosol_energy_check(radius_samples=2049):
    """Compare receiver-plane missing flux with independently sampled pupil flux."""
    profile = AerosolProfile()
    geometry = Geometry()
    table = build_ray_table(Atmosphere(),aerosol_profile=profile)
    wavelengths = np.array([450.,550.,677.])
    landing = signed_landing_radius(table,geometry)
    endpoint = max(float(np.abs(landing).max()),geometry.magnification*table.impact_m[-1])
    radius = np.linspace(0,endpoint+1.02*geometry.solar_blur_m,radius_samples)
    transmission = transmission_profile(table,radius,wavelengths,geometry=geometry)
    receiver_deficit = 2*np.pi*simpson(radius[:,None]*(1-transmission),x=radius,axis=0)
    subsamples = np.linspace(0,1,33)
    b = np.unique((table.impact_m[:-1,None]*(1-subsamples)
                   +table.impact_m[1:,None]*subsamples).ravel())
    tau = np.interp(b,table.impact_m,table.air_column_m2)[:,None]*rayleigh_cross_section(wavelengths)
    tau += np.interp(b,table.impact_m,table.ozone_column_m2)[:,None]*ozone_cross_section(wavelengths)
    for layer_tau, layer_profile in zip(table.aerosol_layer_taus, table.aerosol_profiles):
        tau += (np.interp(b,table.impact_m,layer_tau)[:,None]
                *layer_profile.spectral_factor(wavelengths))
    m = geometry.magnification
    transmitted_area = 2*np.pi*m*m*np.trapezoid(b[:,None]*np.exp(-tau),x=b,axis=0)
    pupil_deficit = np.pi*(m*table.impact_m[-1])**2-transmitted_area
    difference = np.abs(receiver_deficit/pupil_deficit-1)
    return {"pass":bool(np.all(difference < 2e-5)),"wavelength_nm":wavelengths.tolist(),
            "receiver_deficit_m2":receiver_deficit.tolist(),"pupil_deficit_m2":pupil_deficit.tolist(),
            "relative_difference":difference.tolist(),"radius_samples":radius_samples,
            "profile":asdict(profile)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=Path("outputs/validation/aerosol.json"))
    args = parser.parse_args()
    checks = {"straight_bessel":straight_exponential_checks(),
              "curved_cartesian":curved_aerosol_checks(),"global_energy":aerosol_energy_check()}
    report = {"pass":all(value["pass"] for value in checks.values()),"checks":checks,
              "scope":"Numerical transport closure; not aerosol-profile or observational validation"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report["pass"] else 1)


if __name__=="__main__":
    main()
