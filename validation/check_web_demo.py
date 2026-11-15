"""Independent browser/core cross-check; failures block deployment.

Run separately from pytest. References are freshly integrated through production
APIs, never reconstructed from the exported unit-column/kernel tables.
"""
from __future__ import annotations

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import platform
import threading

import numpy as np

from eclipse.aerosol import AerosolColumn, AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.cli import load_calibrated_preset
from eclipse.color import linear_display_rgb, relative_photopic, spectrum_to_xyz
from eclipse.rays import build_ray_table
from eclipse.spectra import solar_irradiance
from eclipse.transport import AU_M, Geometry, transmission_profile
from tools.export_web_demo import radius_grid

ROOT = Path(__file__).resolve().parents[1]


def production(parameters, radius):
    name, aod, rh = parameters['preset'], parameters['aod'], parameters['rh']
    if name == 'none' or (name in ('opac', 'opac-ft') and aod == 0):
        column = None
    elif name == 'opac':
        column = AerosolColumn.opac_default(aod, rh)
    elif name == 'opac-ft':
        column = AerosolProfile.opac_free_troposphere(aod, rh)
    else:
        column = load_calibrated_preset(name.removeprefix('calibrated:'))[0]
    g, wave = Geometry(), np.arange(380., 781., 5.)
    table = build_ray_table(Atmosphere(ozone_du=parameters['du']), aerosol_profile=column)
    reference = solar_irradiance(wave)*(AU_M/(g.sun_distance_m+g.moon_distance_m))**2
    t = transmission_profile(table, radius, wave)
    spectra = t*reference
    return {
        'transmission': t, 'xyz': spectrum_to_xyz(spectra, wave)/spectrum_to_xyz(reference, wave)[1],
        'rgb': linear_display_rgb(spectra, wave, reference),
        'relativeY': relative_photopic(spectra, wave, reference),
        'reference': reference,
    }


def compare(actual, target, atol, rtol):
    a = np.asarray(actual, dtype=float).reshape(target.shape)
    if not np.all(np.isfinite(a)):
        raise ValueError('nonfinite browser result')
    difference = np.abs(a-target)
    score = difference/(atol+rtol*np.abs(target))
    index = np.unravel_index(int(score.argmax()), score.shape)
    return {'pass': bool(np.all(score <= 1)), 'samples': a.size,
            'max_abs': float(difference.max()), 'worst_tolerance_ratio': float(score[index]),
            'worst_index': [int(i) for i in index]}


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def check_site(site: Path) -> dict:
    from playwright.sync_api import sync_playwright
    site = site.resolve()
    manifest = json.loads((site/'manifest.json').read_text())
    report = {'pass': False, 'mobile_real_device': 'NOT VERIFIED', 'checks': [],
              'host': platform.platform(), 'site': str(site), 'schema_version': manifest['schema_version']}
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(QuietHandler, directory=str(site.parent)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            report['browser_version'] = browser.version
            page = browser.new_page(viewport={'width': 1280, 'height': 900})
            errors = []
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(f'http://127.0.0.1:{server.server_port}/{site.name}/')
            page.evaluate('''async () => {
                window.kernel=await import('./kernel.js');
                window.data=await kernel.loadSite(new URL('./',location.href));
                window.map=kernel.pixelMap(data.common.radius,data.manifest.moon_radius_m,2860000,0);
            }''')
            cases = [dict(preset=n,du=du,aod=aod,rh=rh) for n,du,aod,rh in [
                ('opac',300,.013,.5), ('opac',0,.2,.99), ('opac',600,.001,.7),
                ('opac-ft',175,.03,.85), ('none',300,0,.5), ('opac',300,0,.99),
                *[(n,300,.013,.5) for n in manifest['presets'] if n.startswith('calibrated:')],
            ]]
            for parameters in cases:
                reference = production(parameters, radius_grid())
                actual = page.evaluate('''async parameters => {
                    const result=await kernel.reconstruct(data,parameters);
                    return Object.fromEntries(Object.entries(result).map(([k,v])=>[k,Array.from(v)]));
                }''', parameters)
                quantities = {key: compare(actual[key], reference[key], 2e-11 if key == 'rgb' else 2e-12, 2e-8)
                              for key in ('transmission', 'xyz', 'rgb', 'relativeY')}
                report['checks'].append({'parameters': parameters, 'quantities': quantities})
            report['timing'] = page.evaluate('''async () => {
                const ctx=document.querySelector('canvas').getContext('2d'),samples=[];
                for(let n=0;n<230;n++) {
                    const start=performance.now();
                    const field=await kernel.reconstruct(data,{preset:'opac',du:n%3*300,aod:.001+n%20*.01,rh:[.5,.7,.99][n%3]});
                    const computed=performance.now();
                    kernel.draw(ctx,field.rgb,map,8);
                    const end=performance.now();
                    if(n>=30)samples.push({compute:computed-start,draw:end-computed,total:end-start,communication:0});
                }
                return samples;
            }''')
            values = np.array([s['total'] for s in report['timing']])
            report['timing_summary'] = {'p95_ms': float(np.percentile(values, 95)), 'max_ms': float(values.max()),
                                        'pass': bool(values.max() <= 100), 'samples': len(values)}
            report['console_errors'] = errors
            report['pass'] = bool(report['checks'] and not errors and report['timing_summary']['pass'] and all(
                q['pass'] for row in report['checks'] for q in row['quantities'].values()))
            browser.close()
    finally:
        server.shutdown()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = check_site(args.site)
    except Exception as error:
        report = {'pass': False, 'error': f'{type(error).__name__}: {error}'}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('checks', 'timing')}))
    raise SystemExit(0 if report['pass'] else 1)


if __name__ == '__main__':
    main()
