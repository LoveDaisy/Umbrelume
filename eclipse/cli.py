"""Command-line entry point for reproducible eclipse simulations."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import numpy as np

from . import __version__
from .atmosphere import Atmosphere, EARTH_RADIUS_M
from .color import spectrum_to_xyz
from .report import write_report
from .transport import Geometry, simulate, signed_landing_radius


def _positive(value):
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be positive and finite")
    return number


def _nonnegative(value):
    number = float(value)
    if not np.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative and finite")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("simulate", help="Compute spectra, ray tables and eclipse figures")
    run.add_argument("--out", type=Path, default=Path("outputs/baseline"))
    run.add_argument("--ozone-du", type=_nonnegative, default=300.0)
    run.add_argument("--density-scale", type=_nonnegative, default=1.0)
    run.add_argument("--moon-distance-km", type=_positive, default=384400.0)
    run.add_argument("--moon-offset-km", type=_nonnegative, default=2860.0,
                     help="Moon center's distance from the shadow axis")
    run.add_argument("--exposure-stops", type=float, default=8.0)
    run.add_argument("--height-step-m", type=_positive, default=100.0)
    run.add_argument("--radial-step-km", type=_positive, default=20.0)
    run.add_argument("--radius-max-km", type=_positive, default=9000.0)
    run.add_argument("--wavelength-step-nm", type=_positive, default=5.0)
    run.add_argument("--no-plots", action="store_true", help="Write numerical results only")
    args = parser.parse_args(argv)
    if not np.isfinite(args.exposure_stops) or abs(args.exposure_stops) > 100:
        parser.error("--exposure-stops must be finite and between -100 and 100")
    if args.wavelength_step_nm > 100:
        parser.error("--wavelength-step-nm must be at most 100")
    if not args.no_plots and args.radius_max_km < args.moon_offset_km + 1737.4:
        parser.error("--radius-max-km must cover the entire lunar disk")
    radius = np.linspace(0, args.radius_max_km*1000,
                         int(np.ceil(args.radius_max_km/args.radial_step_km))+1)
    wl = np.linspace(380, 780, int(np.ceil(400/args.wavelength_step_nm))+1)
    atmosphere = Atmosphere(ozone_du=args.ozone_du, density_scale=args.density_scale)
    geometry = Geometry(moon_distance_m=args.moon_distance_km*1000)
    start = time.perf_counter()
    result = simulate(radius, wl, atmosphere=atmosphere, geometry=geometry,
                      height_step_m=args.height_step_m)
    elapsed = time.perf_counter()-start
    relative_y = spectrum_to_xyz(result.irradiance, wl)[:, 1] / spectrum_to_xyz(
        result.reference_irradiance, wl)[1]
    table = result.ray_table
    umbra = geometry.magnification*EARTH_RADIUS_M - geometry.solar_blur_m
    penumbra = geometry.magnification*EARTH_RADIUS_M + geometry.solar_blur_m
    metadata = {
        "version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model": "Finite-distance paraxial annular transport with spherical refracted ray tables",
        "atmosphere": asdict(atmosphere),
        "geometry": asdict(geometry),
        "numerics": {"height_step_m": args.height_step_m, "path_quadrature_order": 48,
                     "receiver_integration_order": 8, "radius_count": len(radius),
                     "wavelength_count": len(wl), "ray_count": len(table.altitude_m),
                     "wavelength_range_nm": [float(wl[0]), float(wl[-1])],
                     "wavelength_step_nm": float(wl[1]-wl[0]),
                     "radial_step_m": float(radius[1]-radius[0])},
        "geometric_umbra_radius_m": umbra,
        "geometric_penumbra_radius_m": penumbra,
        "moon_offset_m": args.moon_offset_km*1000,
        "axis_relative_photopic_irradiance": float(relative_y[0]),
        "axis_dimming_stops": float(-np.log2(relative_y[0])) if relative_y[0] > 0 else None,
        "simulation_seconds": elapsed,
        "approximations": [
            "Spherical standard dry atmosphere; Gaussian ozone height profile",
            "Achromatic refracted paths at 550 nm; spectral Rayleigh and ozone extinction",
            "Uniform solar disk, finite-distance first-order paraxial propagation",
            "No aerosol/cloud extinction and no scattered light entering receiver",
            "Neutral shadow-plane lunar screen; no lunar reflectance spectrum or BRDF",
            "Photopic colorimetry with global solar white; no observer atmosphere or dark adaptation",
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out/"spectra.npz", radius_m=radius, wavelength_nm=wl,
                        irradiance_W_m2_nm=result.irradiance,
                        irradiance_no_ozone_W_m2_nm=result.irradiance_no_ozone,
                        reference_irradiance_W_m2_nm=result.reference_irradiance)
    np.savez_compressed(args.out/"rays.npz", tangent_altitude_m=table.altitude_m,
                        vacuum_impact_parameter_m=table.impact_m,
                        bending_rad=table.bending_rad,
                        air_column_m2=table.air_column_m2,
                        ozone_column_m2=table.ozone_column_m2,
                        signed_landing_radius_m=signed_landing_radius(table, geometry))
    (args.out/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    csv = np.column_stack((radius/1000, relative_y, result.irradiance / result.reference_irradiance))
    header = "radius_km,relative_photopic_irradiance," + ",".join(f"T_{w:g}nm" for w in wl)
    np.savetxt(args.out/"profile.csv", csv, delimiter=",", header=header, comments="")
    if not args.no_plots:
        write_report(radius, wl, result.irradiance, result.reference_irradiance,
                     result.irradiance_no_ozone, args.out,
                     moon_offset_m=args.moon_offset_km*1000,
                     exposure_stops=args.exposure_stops, metadata=metadata)
    print(f"Computed {len(radius)} radii x {len(wl)} wavelengths from {len(table.altitude_m)} forward rays "
          f"in {elapsed:.2f} s")
    print(f"Shadow-axis relative photopic irradiance: {relative_y[0]:.6g} "
          f"({metadata['axis_dimming_stops']:.2f} stops below uneclipsed sunlight)"
          if relative_y[0] > 0 else "Shadow axis: no transmitted light")
    print(f"Results: {args.out.resolve()}")


if __name__ == "__main__":
    main()
