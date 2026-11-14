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
from .aerosol import (
    AerosolColumn, AerosolProfile, FREE_TROPOSPHERE_SOURCE, OPAC_DEFAULT_LAYER_SOURCES,
    OPAC_PDF_SHA256, PRIOR_PARAMETERS_JSON, RH_MAX, RH_MIN, ft_growth_factor,
    ft_humidity_provenance,
)
from .color import EXPOSURE_STOP_LIMIT, relative_photopic
from .rays import QUADRATURE_ORDER
from .report import write_report
from .transport import MOON_RADIUS_M, RECEIVER_INTEGRATION_ORDER, Geometry, simulate, signed_landing_radius


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


def _aerosol_layer_metadata(aerosol_prior, layer_sources, rh=0.5,
                            ft_aod550_base=0.013):
    """Per-layer metadata records with their pinned provenance.

    ``layer_sources`` is produced by the same construction branch that built
    ``aerosol_prior``, so every layer pairs with its own provenance. OPAC
    layers carry their prior_parameters.json rows and the OPAC PDF sha256; a
    calibration-lane layer carries the registry's own provenance verbatim
    (``source_class == "calibrated_registry"``). The free-troposphere layer's
    source record additionally carries the humidity-axis provenance at the
    effective ``rh`` (the 50% reference default included, so the schema is
    uniform; old consumers ignore the new keys).
    """
    layers = (aerosol_prior.layers if isinstance(aerosol_prior, AerosolColumn)
              else (aerosol_prior,))
    records = []
    for layer, source in zip(layers, layer_sources, strict=True):
        if source.get("source_class") == "calibrated_registry":
            records.append({**asdict(layer), "source": source})
        else:
            record_source = {**source,
                             "prior_parameters_json": PRIOR_PARAMETERS_JSON,
                             "opac1998_sha256": OPAC_PDF_SHA256}
            if source.get("candidate") == FREE_TROPOSPHERE_SOURCE["candidate"]:
                record_source.update(ft_humidity_provenance(rh, ft_aod550_base))
            records.append({**asdict(layer), "source": record_source})
    return records


CALIBRATED_PRESETS_JSON = Path(__file__).parent/"data"/"calibrated_presets.json"


def load_calibrated_preset(name: str) -> tuple[AerosolColumn, dict]:
    """Lazy loader for one calibration-lane preset.

    Only the ``calibrated:<name>`` branch of ``--aerosol`` ever touches the
    registry; the default path performs no registry IO. Provenance is
    mandatory at field level: an entry without the ``kind='calibrated'`` and
    ``calibrated_not_physical_prior=true`` markers, or without ``sha256`` /
    ``inversion_product``, fails loud here (the sha-vs-file equality itself is
    anchored pytest-side). The name-to-event mapping lives solely in the
    entry's ``event`` field. Any malformed structure surfaces as a ValueError
    carrying the preset name, so the CLI fails through its one parser.error
    path.
    """
    try:
        registry = json.loads(CALIBRATED_PRESETS_JSON.read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f"calibrated preset registry is unreadable ({error})") from None
    try:
        presets = (registry or {}).get("presets") or {}
        entry = presets.get(name)
    except AttributeError as error:
        raise ValueError(f"calibrated preset registry is malformed ({error})") from None
    if entry is None:
        raise KeyError(f"calibrated preset {name!r} is not in "
                       f"eclipse/data/calibrated_presets.json")
    try:
        provenance = entry.get("provenance") or {}
        if (provenance.get("kind") != "calibrated"
                or provenance.get("calibrated_not_physical_prior") is not True):
            raise ValueError(f"calibrated preset {name!r}: provenance must carry "
                             f"kind='calibrated' and "
                             f"calibrated_not_physical_prior=true")
        for field in ("sha256", "inversion_product"):
            if not provenance.get(field):
                raise ValueError(f"calibrated preset {name!r}: provenance field "
                                 f"{field!r} is missing (mandatory provenance, "
                                 f"lane discipline)")
        if not isinstance(entry.get("event"), str) or not entry["event"]:
            raise ValueError(f"calibrated preset {name!r}: 'event' must be a "
                             f"non-empty string (the sole name-to-event mapping)")
        layer = AerosolProfile(**entry["increment_layer"])
    except (TypeError, AttributeError, KeyError) as error:
        raise ValueError(f"calibrated preset {name!r}: malformed registry entry "
                         f"({error})") from None
    column = AerosolColumn((*AerosolColumn.opac_default().layers, layer))
    source = {"source_class": "calibrated_registry",
              "registry_path": "eclipse/data/calibrated_presets.json",
              "preset": name, "event": entry.get("event"), **provenance}
    return column, source


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, prog="lunar-eclipse")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("simulate", help="Compute spectra, ray tables and eclipse figures")
    run.add_argument("--out", type=Path, default=Path("outputs/current"))
    run.add_argument("--ozone-du", type=_nonnegative, default=300.0)
    run.add_argument("--density-scale", type=_nonnegative, default=1.0)
    run.add_argument("--aerosol", default="opac",
                     help="Aerosol preset: two-layer OPAC default (free troposphere "
                          "2-12 km + stratospheric background 12-35 km), "
                          "free-troposphere-only (historical preset), "
                          "molecular-only counterfactual, or 'calibrated:<name>' "
                          "to load a calibration-lane preset from "
                          "eclipse/data/calibrated_presets.json (calibrated, "
                          "not a physical prior; never a default)")
    run.add_argument("--aod550", type=_nonnegative,
                     help="Override the free-troposphere layer's 550 nm vertical AOD "
                          "(default 0.013; the stratospheric layer stays at 0.005); "
                          "0 disables the whole aerosol component")
    run.add_argument("--rh", type=float, default=None,
                     help="Free-troposphere relative humidity for the pinned OPAC "
                          "humidity growth (finite, within [0.5, 0.99]; default "
                          "0.5 = the published 50%% reference state, bit-exact "
                          "unchanged). Legal only with --aerosol opac or opac-ft; "
                          "it acts on the free-troposphere layer only")
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
    if not np.isfinite(args.exposure_stops) or abs(args.exposure_stops) > EXPOSURE_STOP_LIMIT:
        parser.error(f"--exposure-stops must be finite and between "
                     f"{-EXPOSURE_STOP_LIMIT:g} and {EXPOSURE_STOP_LIMIT:g}")
    if args.wavelength_step_nm > 100:
        parser.error("--wavelength-step-nm must be at most 100")
    if not args.no_plots and args.radius_max_km < args.moon_offset_km + MOON_RADIUS_M / 1000:
        parser.error("--radius-max-km must cover the entire lunar disk")
    calibrated_name = None
    if args.aerosol.startswith("calibrated:"):
        calibrated_name = args.aerosol[len("calibrated:"):]
        if not calibrated_name:
            parser.error("--aerosol calibrated:<name> requires a preset name")
    elif args.aerosol not in ("opac", "opac-ft", "none"):
        parser.error("--aerosol must be opac, opac-ft, none, or calibrated:<name>")
    if args.aerosol == "none" and args.aod550 is not None:
        parser.error("--aod550 requires --aerosol opac or opac-ft")
    if calibrated_name is not None and args.aod550 is not None:
        parser.error("--aod550 cannot be combined with --aerosol "
                     "calibrated:<name>: the preset fully determines the column")
    if args.rh is not None:
        if not np.isfinite(args.rh) or not RH_MIN <= args.rh <= RH_MAX:
            parser.error(f"--rh must be finite within [{RH_MIN}, {RH_MAX}] "
                         f"(the pinned humidity domain)")
        if calibrated_name is not None:
            parser.error("--rh cannot be combined with --aerosol calibrated:<name>: "
                         "the preset fully determines the column")
        if args.aerosol == "none" or args.aod550 == 0:
            parser.error("--rh requires a nonzero free-troposphere aerosol "
                         "component (--aerosol opac or opac-ft)")
    rh = 0.5 if args.rh is None else args.rh
    ft_aod550_base = 0.013 if args.aod550 is None else args.aod550
    radius = np.linspace(0, args.radius_max_km*1000,
                         int(np.ceil(args.radius_max_km/args.radial_step_km))+1)
    wl = np.linspace(380, 780, int(np.ceil(400/args.wavelength_step_nm))+1)
    atmosphere = Atmosphere(ozone_du=args.ozone_du, density_scale=args.density_scale)
    if calibrated_name is not None:
        try:
            aerosol_prior, preset_source = load_calibrated_preset(calibrated_name)
        except (KeyError, ValueError) as error:
            parser.error(str(error))
        layer_sources = (*OPAC_DEFAULT_LAYER_SOURCES, preset_source)
    elif args.aerosol == "none" or args.aod550 == 0:
        # AOD=0 is the molecular-only counterfactual: the whole aerosol
        # component is disabled and reproduces --aerosol none exactly.
        aerosol_prior, layer_sources = None, ()
    elif args.aerosol == "opac":
        aerosol_prior = AerosolColumn.opac_default(ft_aod550=ft_aod550_base, rh=rh)
        layer_sources = OPAC_DEFAULT_LAYER_SOURCES
    else:
        aerosol_prior = AerosolProfile.opac_free_troposphere(
            aod550=ft_aod550_base, rh=rh)
        layer_sources = (FREE_TROPOSPHERE_SOURCE,)
    if aerosol_prior is None:
        aerosol_text = "No aerosol extinction"
    elif calibrated_name is not None:
        aerosol_text = (f"Calibrated preset '{calibrated_name}': the production "
                        f"two-layer OPAC default plus an inversion increment fitted "
                        f"to the {preset_source['event']} brightness residual; "
                        f"calibrated, not a physical prior, and never a default")
    elif args.aerosol == "opac":
        aerosol_text = ("OPAC free-troposphere (2-12 km) and stratospheric-background "
                        "(12-35 km) aerosol extinction only; no boundary-layer component")
    else:
        aerosol_text = ("OPAC free-troposphere aerosol extinction only; no boundary "
                        "layer or stratospheric aerosol preset")
    if rh != 0.5:
        aerosol_text += (f"; free-troposphere humidity growth at rh={rh:g} "
                         f"(pinned OPAC WASO growth, factor "
                         f"{ft_growth_factor(rh):.6g} on the 50% reference state)")
    geometry = Geometry(moon_distance_m=args.moon_distance_km*1000)
    start = time.perf_counter()
    result = simulate(radius, wl, atmosphere=atmosphere, geometry=geometry,
                      height_step_m=args.height_step_m, aerosol_profile=aerosol_prior)
    elapsed = time.perf_counter()-start
    relative_y = relative_photopic(result.irradiance, wl, result.reference_irradiance)
    table = result.ray_table
    umbra = geometry.magnification*EARTH_RADIUS_M - geometry.solar_blur_m
    penumbra = geometry.magnification*EARTH_RADIUS_M + geometry.solar_blur_m
    metadata = {
        "version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model": "Finite-distance paraxial annular transport with spherical refracted ray tables",
        "atmosphere": asdict(atmosphere),
        "aerosol": ({"layers": _aerosol_layer_metadata(aerosol_prior, layer_sources,
                                                       rh=rh,
                                                       ft_aod550_base=ft_aod550_base)}
                    if aerosol_prior is not None else None),
        "aerosol_preset": args.aerosol,
        "aerosol_aod_source": ("Disabled" if aerosol_prior is None else
                                "Calibrated inversion increment; not an independent "
                                "prior" if calibrated_name is not None else
                                "OPAC literature prior" if args.aod550 is None else
                                "User override; not an independent prior"),
        "geometry": asdict(geometry),
        "numerics": {"height_step_m": args.height_step_m, "path_quadrature_order": QUADRATURE_ORDER,
                     "receiver_integration_order": RECEIVER_INTEGRATION_ORDER, "radius_count": len(radius),
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
            aerosol_text,
            "No cloud extinction or scattered light entering receiver",
            "Neutral shadow-plane lunar screen; no lunar reflectance spectrum or BRDF",
            "Photopic colorimetry with global solar white; no observer atmosphere or dark adaptation",
        ],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out/"spectra.npz", radius_m=radius, wavelength_nm=wl,
                        irradiance_W_m2_nm=result.irradiance,
                        irradiance_no_ozone_W_m2_nm=result.irradiance_no_ozone,
                        reference_irradiance_W_m2_nm=result.reference_irradiance)
    ray_records = {"tangent_altitude_m": table.altitude_m,
                   "vacuum_impact_parameter_m": table.impact_m,
                   "bending_rad": table.bending_rad,
                   "air_column_m2": table.air_column_m2,
                   "ozone_column_m2": table.ozone_column_m2}
    if table.aerosol_layer_taus is not None:
        ray_records["aerosol_layer_taus"] = np.stack(table.aerosol_layer_taus)
    else:
        # Molecular-only tables keep the historical zero-array key untouched.
        ray_records["aerosol_tau550"] = np.zeros_like(table.altitude_m)
    ray_records["signed_landing_radius_m"] = signed_landing_radius(table, geometry)
    np.savez_compressed(args.out/"rays.npz", **ray_records)
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
