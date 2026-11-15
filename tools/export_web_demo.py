"""Export a deterministic, production-derived browser interpretation outlet.

The core remains authoritative: this tool only serializes production ray columns,
receiver quadrature and colorimetry. No display exposure is baked into the assets.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

# Direct script invocation must never pick an editable install in another worktree.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from eclipse import __version__
from eclipse.aerosol import AerosolColumn, AerosolProfile, active_aerosol_layers, ft_growth_factor
from eclipse.atmosphere import Atmosphere, EARTH_RADIUS_M
from eclipse.cli import CALIBRATED_PRESETS_JSON, load_calibrated_preset
from eclipse.color import BRADFORD, D65, XYZ_TO_SRGB, spectrum_to_xyz
from eclipse.rays import QUADRATURE_ORDER, build_ray_table
from eclipse.spectra import cie_xyz, ozone_cross_section, rayleigh_cross_section, solar_irradiance
from eclipse.transport import (
    AU_M, Geometry, MOON_RADIUS_M, RECEIVER_INTEGRATION_ORDER,
    circle_overlap, receiver_quadrature,
)

SCHEMA_VERSION = 1
MAX_DATA_BYTES = 32 * 2**20


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def radius_grid() -> np.ndarray:
    """20 km baseline, with only cells failing the fixed spatial gate refined."""
    return np.unique(np.concatenate((
        np.arange(0, 9_000_001, 20_000, dtype=float),
        [2_010_000., 2_030_000.],
        np.arange(4_570_000., 4_720_000., 20_000.),
        np.arange(4_615_000., 4_670_000., 10_000.),
    )))


def preset_columns() -> dict:
    presets = {
        'none': (None, {'kind': 'counterfactual', 'label': 'Molecular only'}),
        'opac-ft': (AerosolProfile.opac_free_troposphere(),
                    {'kind': 'literature_prior', 'label': 'OPAC free troposphere'}),
        'opac': (AerosolColumn.opac_default(),
                 {'kind': 'literature_prior', 'label': 'OPAC two-layer prior'}),
    }
    registry = json.loads(CALIBRATED_PRESETS_JSON.read_text())
    for name in registry['presets']:
        column, source = load_calibrated_preset(name)
        presets[f'calibrated:{name}'] = (column, {
            'kind': 'calibrated', 'label': f'{source["event"]} — fitted, not a prior',
            'provenance': source,
        })
    return presets


def input_fingerprints() -> dict:
    paths = sorted({
        *ROOT.glob('eclipse/*.py'), *ROOT.glob('eclipse/data/*'),
        *ROOT.glob('web/*'), ROOT/'tools/export_web_demo.py',
        ROOT/'pyproject.toml', ROOT/'uv.lock',
    })
    return {str(p.relative_to(ROOT)): sha256(p.read_bytes()) for p in paths if p.is_file()}


def pack_arrays(arrays: dict) -> tuple[bytes, dict]:
    parts, fields, offset = [], {}, 0
    for name, values in arrays.items():
        dtype = '<u4' if np.asarray(values).dtype.kind in 'iu' else '<f8'
        a = np.ascontiguousarray(values, dtype=dtype)
        padding = (-offset) % 8
        parts.append(b'\0' * padding)
        offset += padding
        fields[name] = {'shape': list(a.shape), 'dtype': dtype,
                        'offset': offset, 'byteLength': a.nbytes}
        parts.append(a.tobytes())
        offset += a.nbytes
    return b''.join(parts), fields


def common_arrays(radius: np.ndarray, wave: np.ndarray, geometry: Geometry) -> dict:
    reference = solar_irradiance(wave) * (
        AU_M/(geometry.sun_distance_m + geometry.moon_distance_m))**2
    ref_xyz = spectrum_to_xyz(reference, wave)
    adaptation = np.linalg.solve(BRADFORD, np.diag(
        (BRADFORD @ D65)/(BRADFORD @ (ref_xyz/ref_xyz[1]))) @ BRADFORD)
    trap = np.concatenate(([.5], np.ones(len(wave)-2), [.5])) * 5
    return {
        'radius': radius, 'wave': wave, 'reference': reference,
        'cie': cie_xyz(wave)*trap[:, None],
        'rgbMatrix': XYZ_TO_SRGB @ adaptation / ref_xyz[1],
        'referenceXYZ': ref_xyz,
    }


def family_arrays(column, radius: np.ndarray, wave: np.ndarray, geometry: Geometry) -> dict:
    # Unit layers are built through the production integrator. Their support and
    # spectrum are untouched; only their amplitude is deferred to the browser.
    layers = active_aerosol_layers(column)
    unit = AerosolColumn(tuple(replace(layer, aod550=1.) for layer in layers)) if layers else None
    table = build_ray_table(Atmosphere(), aerosol_profile=unit)
    samples, weights, offsets = [], [], [0]
    for r in radius:
        bb, area = receiver_quadrature(table, float(r), geometry)
        samples.append(bb)
        weights.append(area)
        offsets.append(offsets[-1]+len(bb))
    nodes, index = np.unique(np.concatenate(samples), return_inverse=True)
    components = [np.interp(nodes, table.impact_m, table.air_column_m2),
                  np.interp(nodes, table.impact_m, table.ozone_column_m2)/300.]
    factors = [rayleigh_cross_section(wave), ozone_cross_section(wave)]
    for tau, profile in zip(table.aerosol_layer_taus or (), table.aerosol_profiles, strict=True):
        components.append(np.interp(nodes, table.impact_m, tau))
        factors.append(profile.spectral_factor(wave))
    a, m = geometry.solar_blur_m, geometry.magnification
    return {
        'offsets': np.array(offsets), 'index': index,
        'weight': np.concatenate(weights), 'b': nodes,
        'components': np.array(components).T, 'factors': np.array(factors),
        'outer': 1-circle_overlap(radius, a, m*table.impact_m[-1])/(np.pi*a*a),
        'altitude': table.altitude_m, 'impact': table.impact_m,
        'bending': table.bending_rad,
    }


def export_site(out: Path) -> dict:
    out = out.absolute()
    if out.is_symlink() or (out.exists() and (not out.is_dir() or any(out.iterdir()))):
        raise ValueError('--out must be an absent or empty, non-symlink directory')
    # Everything is staged next to the destination; failures never publish a
    # half-site. No downloads, generated scientific snapshots or input rewrites.
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.web-export-', dir=out.parent) as temp:
        stage = Path(temp)/'site'
        shutil.copytree(ROOT/'web', stage)
        (stage/'data').mkdir()
        radius, wave, geometry = radius_grid(), np.arange(380., 781., 5.), Geometry()
        manifest = {
            'schema_version': SCHEMA_VERSION, 'software_version': __version__,
            'generator': 'uv run python tools/export_web_demo.py --out <empty-dir>',
            'inputs_sha256': input_fingerprints(),
            'geometry': asdict(geometry), 'moon_radius_m': MOON_RADIUS_M,
            'boundaries_m': [geometry.magnification*EARTH_RADIUS_M+s*geometry.solar_blur_m
                             for s in (-1, 1)],
            'numerics': {'height_step_m': 100., 'ray_quadrature_order': QUADRATURE_ORDER,
                         'receiver_integration_order': RECEIVER_INTEGRATION_ORDER},
            'units': {'radius': 'm', 'b': 'm', 'altitude': 'm', 'impact': 'm',
                      'bending': 'rad', 'wave': 'nm', 'reference': 'W m^-2 nm^-1',
                      'components': ['air molecules m^-2', 'ozone molecules m^-2 DU^-1',
                                     'per-layer tau550 per vertical AOD550'],
                      'weight': 'conserved annular area / solar disk area'},
            'axes': {'components': ['unique_b', 'component'], 'factors': ['component', 'wavelength'],
                     'cie': ['wavelength', 'XYZ'], 'rgbMatrix': ['RGB', 'XYZ']},
            'domain': {'radius_m': [0, 9_000_000], 'du': [0, 600], 'aod': [0, .2],
                       'positive_aod_min': .001, 'rh': [.5, .99], 'ev': [-4, 20]},
            'defaults': {'preset': 'opac', 'du': 300, 'aod': .013, 'rh': .5,
                         'ev': 8, 'x_m': 2_860_000, 'y_m': 0},
            'rh_knots': [[r, ft_growth_factor(r)] for r in (.5, .7, .99)],
            'presets': {}, 'families': {},
            'source_notice': 'Source-code MIT license does not license all source data. '
                             'See SOURCES.md for per-source terms and verification status.',
        }

        def save(name, arrays):
            raw, fields = pack_arrays(arrays)
            encoded = gzip.compress(raw, mtime=0)
            path = f'data/{name}.bin.gz'
            (stage/path).write_bytes(encoded)
            return {'file': path, 'sha256': sha256(encoded), 'byteLength': len(encoded),
                    'decodedByteLength': len(raw), 'encoding': 'gzip', 'fields': fields}

        manifest['common'] = save('common', common_arrays(radius, wave, geometry))
        seen = {}
        for name, (column, meta) in preset_columns().items():
            arrays = family_arrays(column, radius, wave, geometry)
            # Exact byte identity, not approximate height/node merging. All
            # calibrated presets share a geometry but retain independent amplitudes.
            raw, _ = pack_arrays(arrays)
            identity = sha256(raw)
            if identity not in seen:
                family = f'family-{len(seen)}'
                seen[identity] = family
                manifest['families'][family] = save(family, arrays)
            layers = active_aerosol_layers(column)
            manifest['presets'][name] = {**meta, 'family': seen[identity],
                                        'amplitudes': [1., 300.] + [p.aod550 for p in layers],
                                        'layers': [asdict(p) for p in layers]}
        total = manifest['common']['byteLength'] + sum(
            f['byteLength'] for f in manifest['families'].values())
        if total > MAX_DATA_BYTES:
            raise ValueError(f'data budget exceeded: {total} > {MAX_DATA_BYTES}')
        manifest['data_bytes'] = total
        (stage/'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
        shutil.copyfile(ROOT/'eclipse/data/SOURCES.md', stage/'SOURCES.md')
        shutil.copyfile(ROOT/'LICENSE', stage/'LICENSE.txt')
        if out.exists():
            out.rmdir()  # Checked empty above; a concurrent writer makes this fail closed.
        stage.rename(out)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        manifest = export_site(args.out)
    except (ValueError, OSError) as error:
        parser.exit(1, f'export failed: {error}\n')
    print(json.dumps({'site': str(args.out), 'families': len(manifest['families']),
                      'data_bytes': manifest['data_bytes']}))


if __name__ == '__main__':
    main()
