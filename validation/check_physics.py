"""Independent Cartesian-ray and global-energy checks for the eclipse model.

Run ``uv run python -m validation.check_physics --out results/physics.json``.
The Cartesian reference integrates the ray equation in path length rather than
using the production turning-point quadrature. The global check compares a
receiver-plane integral with an incident-pupil energy budget, independent of
where any individual refracted ray lands.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp

from eclipse.atmosphere import (
    Atmosphere, EARTH_RADIUS_M, LAYER_BOUNDARIES_M, REFERENCE_NUMBER_DENSITY,
)
from eclipse.rays import build_ray_table
from eclipse.spectra import ozone_cross_section, rayleigh_cross_section, refractive_index_air
from eclipse.transport import Geometry, signed_landing_radius, transmission_profile


def cartesian_reference(
    atmosphere: Atmosphere,
    tangent_altitude_m: float,
    reference_wavelength_nm: float = 550.0,
) -> dict:
    """Integrate one outward half-ray; symmetry supplies the inward half.

    Position is Cartesian in the ray plane. Direction evolves according to
    du/ds = grad(log n) - u [u dot grad(log n)]. Columns are accumulated in
    scaled units to keep the ODE error norm meaningful for every component.
    """
    if not 0 <= tangent_altitude_m < atmosphere.top_m:
        raise ValueError("Reference ray tangent must lie inside the atmosphere.")
    coefficient = (float(refractive_index_air(reference_wavelength_nm)) - 1) / REFERENCE_NUMBER_DENSITY
    ozone_scale = 1e18
    tangent_radius = EARTH_RADIUS_M + tangent_altitude_m
    top_radius = EARTH_RADIUS_M + atmosphere.top_m

    def derivative(_, state):
        position, direction = state[:2], state[2:4]
        radial_distance = np.linalg.norm(position)
        altitude = max(float(radial_distance - EARTH_RADIUS_M), 0.0)
        density = float(atmosphere.number_density(altitude))
        refractive_index = 1 + coefficient * density
        radial_gradient = coefficient * float(atmosphere.density_gradient(altitude)) / refractive_index
        gradient = radial_gradient * position / radial_distance
        direction_derivative = gradient - direction * np.dot(direction, gradient)
        return np.array([
            *direction, *direction_derivative,
            density / REFERENCE_NUMBER_DENSITY,
            float(atmosphere.ozone_density(altitude)) / ozone_scale,
        ])

    def exit_atmosphere(_, state):
        return np.linalg.norm(state[:2]) - top_radius

    exit_atmosphere.terminal = True
    exit_atmosphere.direction = 1
    solution = solve_ivp(
        derivative, (0, 4 * top_radius),
        [tangent_radius, 0, 0, 1, 0, 0],
        events=exit_atmosphere, method="DOP853", rtol=2e-10,
        atol=[1e-6, 1e-6, 1e-13, 1e-13, 1e-5, 1e-5], max_step=3000,
    )
    if not solution.success or not len(solution.t_events[0]):
        raise RuntimeError("Independent Cartesian reference did not reach the atmosphere top.")
    states = solution.y
    radii = np.hypot(states[0], states[1])
    n = 1 + coefficient * atmosphere.number_density(np.maximum(radii - EARTH_RADIUS_M, 0))
    invariant = n * (states[0] * states[3] - states[1] * states[2])
    final = states[:, -1]
    return {
        "altitude_m": float(tangent_altitude_m),
        "bending_rad": float(-2 * np.arctan2(final[2], final[3])),
        "air_column_m2": float(2 * final[4] * REFERENCE_NUMBER_DENSITY),
        "ozone_column_m2": float(2 * final[5] * ozone_scale),
        "max_invariant_relative_error": float(np.max(np.abs(invariant / invariant[0] - 1))),
        "max_direction_norm_error": float(np.max(np.abs(np.hypot(states[2], states[3]) - 1))),
        "accepted_steps": int(len(solution.t)),
    }


def ray_checks(*, height_step_m: float = 100) -> dict:
    """Compare radial production quadrature with independent Cartesian ODEs."""
    atmosphere = Atmosphere()
    table = build_ray_table(atmosphere, height_step_m=height_step_m, quadrature_order=48)
    results = []
    tropopause = float(LAYER_BOUNDARIES_M[1])
    geometry = Geometry()
    coefficient = (float(refractive_index_air(550)) - 1) / REFERENCE_NUMBER_DENSITY
    for altitude in [5000.0, 11000.0, tropopause, 15000.0, 25000.0, 40000.0]:
        reference = cartesian_reference(atmosphere, altitude)
        comparison = {}
        for name in ["bending_rad", "air_column_m2", "ozone_column_m2"]:
            actual = float(np.interp(altitude, table.altitude_m, getattr(table, name)))
            comparison[name] = {
                "radial_quadrature": actual,
                "cartesian_ode": reference[name],
                "relative_difference": abs(actual / reference[name] - 1),
            }
        passed = all(item["relative_difference"] < 2e-5 for item in comparison.values())
        passed &= reference["max_invariant_relative_error"] < 2e-8
        passed &= reference["max_direction_norm_error"] < 2e-9
        impact = (EARTH_RADIUS_M + altitude) * (1 + coefficient * float(atmosphere.number_density(altitude)))
        landing = geometry.magnification * impact - geometry.moon_distance_m * reference["bending_rad"]
        results.append({**reference, "signed_landing_m": float(landing),
                        "comparisons": comparison, "pass": bool(passed)})
    # The derivative discontinuity in this idealised atmosphere produces a
    # small real fold of the landing map immediately below the tropopause.
    # An integrator must retain every branch rather than silently sorting it.
    landing_drop = results[2]["signed_landing_m"] - results[1]["signed_landing_m"]
    return {"pass": all(item["pass"] for item in results) and landing_drop < 0,
            "tropopause_landing_drop_m": landing_drop, "rays": results}


def energy_check(*, extinction: bool = True, radius_samples: int = 2049) -> dict:
    """Close whole-plane missing flux against occulted plus extinguished flux.

    Each source-disk convolution has unit integral. Hence for each wavelength,
    integral 2*pi*r*(1-E/F) dr equals pi*(m*b_top)^2 minus the transmitted
    atmospheric annular input area. This identity is independent of bending.
    """
    atmosphere = Atmosphere()
    geometry = Geometry()
    table = build_ray_table(atmosphere, height_step_m=100, quadrature_order=48)
    wavelengths = np.array([450.0, 550.0, 650.0])
    landing = signed_landing_radius(table, geometry)
    support_radius = max(float(np.max(np.abs(landing))), geometry.magnification * table.impact_m[-1])
    radius = np.linspace(0, support_radius + 1.02 * geometry.solar_blur_m, radius_samples)
    transmission = transmission_profile(
        table, radius, wavelengths, geometry=geometry, extinction=extinction,
        integration_order=8,
    )
    receiver_deficit = 2 * np.pi * simpson(radius[:, None] * (1 - transmission), x=radius, axis=0)
    m = geometry.magnification
    if extinction:
        # Refine only the optical-column interpolation, independently of the
        # receiver integrator's angular geometry and Gauss-Legendre weights.
        nodes = np.linspace(0, 1, 17)
        b = (table.impact_m[:-1, None] * (1 - nodes)
             + table.impact_m[1:, None] * nodes).ravel()
        b = np.unique(b)
        air = np.interp(b, table.impact_m, table.air_column_m2)
        ozone = np.interp(b, table.impact_m, table.ozone_column_m2)
        tau = air[:, None] * rayleigh_cross_section(wavelengths)
        tau += ozone[:, None] * ozone_cross_section(wavelengths)
        transmitted_area = 2 * np.pi * m**2 * np.trapezoid(b[:, None] * np.exp(-tau), x=b, axis=0)
        pupil_deficit = np.pi * (m * table.impact_m[-1])**2 - transmitted_area
    else:
        pupil_deficit = np.full(len(wavelengths), np.pi * (m * table.impact_m[0])**2)
    difference = np.abs(receiver_deficit / pupil_deficit - 1)
    return {
        "pass": bool(np.all(difference < 2e-5)),
        "extinction": extinction,
        "wavelength_nm": wavelengths.tolist(),
        "receiver_deficit_m2": receiver_deficit.tolist(),
        "pupil_deficit_m2": pupil_deficit.tolist(),
        "relative_difference": difference.tolist(),
        "radius_samples": radius_samples,
        "support_endpoint_m": float(radius[-1]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/physics.json"))
    args = parser.parse_args()
    checks = {"cartesian_rays": ray_checks(), "energy_extinction": energy_check(),
              "energy_transparent": energy_check(extinction=False)}
    passed = all(item["pass"] for item in checks.values())
    report = {"pass": passed, "checks": checks,
              "scope": "Numerical closure of the declared model; not observational validation."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
