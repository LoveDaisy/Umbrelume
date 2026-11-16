"""Independent browser/core gate. Nonzero exit prevents publishing the site.

References are fresh production integrations and real CLI products, never
reconstructions from the exported unit-column/W tables. This is deliberately
outside pytest: installing a browser is not a prerequisite for the Python gate.
"""
from __future__ import annotations

import argparse
import base64
from functools import partial
import gzip
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import threading

import numpy as np

from eclipse.aerosol import AerosolColumn, AerosolProfile
from eclipse.atmosphere import Atmosphere, EARTH_RADIUS_M
from eclipse.cli import load_calibrated_preset
from eclipse.color import finish_display, linear_display_rgb, relative_photopic, spectrum_to_xyz
from eclipse.rays import build_ray_table
from eclipse.spectra import solar_irradiance
from eclipse.transport import AU_M, Geometry, MOON_RADIUS_M, transmission_profile
from tools.export_web_demo import input_fingerprints, radius_grid

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
        'relativeY': relative_photopic(spectra, wave, reference), 'reference': reference,
    }


def compare(actual, target, atol, rtol):
    a = np.asarray(actual, dtype=float)
    if a.size != target.size or not target.size:
        raise ValueError(f'empty/mismatched comparison: {a.size} != {target.size}')
    a = a.reshape(target.shape)
    if not np.all(np.isfinite(a)):
        raise ValueError('nonfinite browser result')
    difference = np.abs(a-target)
    score = difference/(atol+rtol*np.abs(target))
    index = np.unravel_index(int(score.argmax()), score.shape)
    return {'pass': bool(np.all(score <= 1)), 'samples': a.size,
            'max_abs': float(difference.max()), 'worst_tolerance_ratio': float(score[index]),
            'worst_index': [int(i) for i in index]}


def parameter_cases(presets):
    rows = []
    for name in ('opac', 'opac-ft'):
        for du in (0, 300, 600):
            for aod in (0, .001, .013, .2):
                for rh in ((.5,) if aod == 0 else (.5, .7, .99)):
                    rows.append(dict(preset=name, du=du, aod=aod, rh=rh))
        rows.extend(dict(preset=name, du=175, aod=.03, rh=rh) for rh in (.6, .85))
    for name in presets:
        if name in ('opac', 'opac-ft'):
            continue
        rows.extend(dict(preset=name, du=du, aod=.013, rh=.5) for du in (0, 300, 600))
    rng = np.random.default_rng(260_116)
    for _ in range(8):
        rows.append(dict(preset=str(rng.choice(['opac', 'opac-ft'])), du=float(rng.uniform(0, 600)),
                         aod=float(rng.uniform(.001, .2)), rh=float(rng.uniform(.5, .99))))
    return rows


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def evaluate_async(page, expression, argument=None, timeout_ms=30_000):
    """Bound application promises too: Playwright's action timeout excludes evaluate."""
    return page.evaluate('''async argument => {
        let timer;
        try {
            return await Promise.race([
                (EXPRESSION)(argument.value),
                new Promise((_,reject)=>{timer=setTimeout(()=>reject(new Error('Browser evaluation deadline exceeded')),argument.timeout);})
            ]);
        } finally { clearTimeout(timer); }
    }'''.replace('EXPRESSION', expression), {'value': argument, 'timeout': timeout_ms})


def await_ready(page):
    return evaluate_async(page, 'async()=>await demo.ready')


def browser_field(page, parameters):
    return evaluate_async(page, '''async parameters => {
        await demo.update(parameters);
        const result=Object.fromEntries(Object.entries(demo.field).map(([k,v])=>[k,Array.from(v)]));
        const {finish}=await import('./kernel.js');
        result.displays=Object.fromEntries([-4,0,8,12,20].map(ev=>[ev,Array.from({length:demo.common.radius.length},(_,r)=>finish(Array.from(demo.field.rgb.slice(r*3,r*3+3)),ev))]));
        return result;
    }''', parameters)


def quantities(actual, reference):
    result = {key: compare(actual[key], reference[key], 2e-11 if key == 'rgb' else 2e-12, 2e-8)
              for key in ('transmission', 'xyz', 'rgb', 'relativeY')}
    if any(np.any(np.asarray(actual[key]) < 0) for key in ('transmission', 'xyz', 'relativeY')):
        raise ValueError('negative irradiance/transmission/XYZ')
    for ev in (-4, 0, 8, 12, 20):
        result[f'display_EV{ev}'] = compare(
            actual['displays'][str(ev)], finish_display(reference['rgb'], ev), 1/255, 0)
    return result


def passed_quantities(rows):
    return bool(rows) and all(q['pass'] for row in rows for q in row['quantities'].values())


def cli_checks(page, presets, radius):
    rows = []
    cases = [dict(preset=name, du=300, aod=.013, rh=.5) for name in presets]
    cases.append(dict(preset='opac', du=175, aod=.03, rh=.85))
    for parameters in cases:
        with tempfile.TemporaryDirectory(prefix='web-cli-') as tmp:
            command = [sys.executable, '-m', 'eclipse', 'simulate', '--out', tmp, '--no-plots',
                       '--aerosol', parameters['preset'], '--ozone-du', str(parameters['du'])]
            if parameters['preset'] in ('opac', 'opac-ft') and parameters != dict(preset=parameters['preset'],du=300,aod=.013,rh=.5):
                command += ['--aod550', str(parameters['aod']), '--rh', str(parameters['rh'])]
            subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
            csv = np.loadtxt(Path(tmp)/'profile.csv', delimiter=',', skiprows=1)
            with np.load(Path(tmp)/'spectra.npz') as archive:
                t = archive['irradiance_W_m2_nm']/archive['reference_irradiance_W_m2_nm']
                ref = archive['reference_irradiance_W_m2_nm'].copy()
                np.testing.assert_array_equal(archive['radius_m'], csv[:, 0]*1000)
            reference = production(parameters, csv[:, 0]*1000)
            # This also verifies fresh API vs real CSV/NPZ, not just browser vs API.
            np.testing.assert_allclose(t, csv[:, 2:], atol=2e-12, rtol=2e-8)
            np.testing.assert_allclose(t, reference['transmission'], atol=2e-12, rtol=2e-8)
            np.testing.assert_array_equal(csv[:, 1], reference['relativeY'])
            actual = browser_field(page, parameters)
            indices = np.searchsorted(radius, csv[:, 0]*1000)
            np.testing.assert_array_equal(radius[indices], csv[:, 0]*1000)
            web_t = np.array(actual['transmission']).reshape(len(radius), 81)[indices]
            web_y = np.array(actual['relativeY'])[indices]
            common_ref = page.evaluate('Array.from(demo.common.reference)')
            rows.append({'parameters': parameters, 'quantities': {
                'profile_T': compare(web_t, csv[:, 2:], 2e-12, 2e-8),
                'profile_Y': compare(web_y, csv[:, 1], 2e-12, 2e-8),
                'SI_reference': compare(common_ref, ref, 1e-14, 2e-12),
            }})
    return rows


def spatial_checks(page, presets, radius):
    g = Geometry()
    boundaries = g.magnification*EARTH_RADIUS_M + np.array([-1, 1])*g.solar_blur_m
    probes = np.unique(np.concatenate(((radius[1:]+radius[:-1])/2,
        np.arange(4_400_000., 4_800_001., 2000.),
        *(b+np.array([-1000., 0., 1000.]) for b in boundaries))))
    cases = [dict(preset='opac', du=du, aod=aod, rh=rh)
             for du, aod, rh in [(300, .013, .5), (0, .013, .5), (600, .013, .5),
                                (300, 0, .5), (300, .2, .5), (300, .013, .99)]]
    cases += [dict(preset=n, du=300, aod=.013, rh=.5) for n in presets if n.startswith('calibrated:')]
    rows = []
    for parameters in cases:
        reference = production(parameters, probes)
        browser_field(page, parameters)
        actual = page.evaluate('''async probes => {
            const {sample,finish}=await import('./kernel.js');
            return {y:probes.map(r=>sample(demo.field.relativeY,demo.common.radius,r,1)[0]),
                displays:Object.fromEntries([-4,0,8,12,20].map(ev=>[ev,probes.map(r=>finish(sample(demo.field.rgb,demo.common.radius,r),ev))]))};
        }''', probes.tolist())
        checks = {'relativeY': compare(actual['y'], reference['relativeY'], 5e-9, .005)}
        for ev in (-4, 0, 8, 12, 20):
            checks[f'RGB_EV{ev}'] = compare(actual['displays'][str(ev)], finish_display(reference['rgb'], ev), 3/255, 0)
        row = {'parameters': parameters, 'quantities': checks}
        for q in checks.values():
            q['worst_radius_m'] = float(probes[q['worst_index'][0]])
        rows.append(row)
    return rows


def pixel_checks(page, radius):
    parameters = dict(preset='opac', du=300, aod=.013, rh=.5)
    reference = production(parameters, radius)
    browser_field(page, parameters)
    g, size = Geometry(), 560
    positions = [0., 2_860_000., g.magnification*EARTH_RADIUS_M-g.solar_blur_m,
                 g.magnification*EARTH_RADIUS_M+g.solar_blur_m,
                 9_000_000.-MOON_RADIUS_M]
    # The geometric penumbra center lies outside the legal whole-moon domain;
    # touch its boundary with the outer moon edge instead, without clamping samples.
    positions[3] -= MOON_RADIUS_M
    yy, xx = np.mgrid[:size, :size]
    dx, dy = ((xx+.5)/size*2-1)*MOON_RADIUS_M, ((yy+.5)/size*2-1)*MOON_RADIUS_M
    mask = dx*dx+dy*dy <= MOON_RADIUS_M**2
    rows = []
    for center in positions:
        pixel_radius = np.hypot(dx+center, dy)
        linear = np.stack([np.interp(pixel_radius, radius, reference['rgb'][:, c]) for c in range(3)], axis=-1)
        for ev in (-4, 0, 8, 12, 20):
            encoded = page.evaluate('''async values => {
                demo.display(values);
                const bytes=document.querySelector('#moon').getContext('2d').getImageData(0,0,560,560).data;
                return await new Promise(resolve=>{const reader=new FileReader();reader.onload=()=>resolve(reader.result.split(',')[1]);reader.readAsDataURL(new Blob([bytes]));});
            }''', {'x_m': center, 'y_m': 0, 'ev': ev})
            actual = np.frombuffer(base64.b64decode(encoded), dtype=np.uint8).reshape(size, size, 4)
            expected = np.rint(finish_display(linear, ev)*255).astype(np.uint8)
            error = np.abs(actual[:, :, :3].astype(int)-expected.astype(int))[mask]
            mask_ok = np.array_equal(actual[:, :, 3] == 255, mask)
            rows.append({'center_m': float(center), 'ev': ev, 'max_code_error': int(error.max()),
                         'mask_equal': bool(mask_ok), 'samples': int(mask.sum()*3),
                         'pass': bool(mask_ok and error.max() <= 1)})
    page.evaluate('demo.display({x_m:2860000,y_m:0,ev:8})')
    # Explicit zero and vacuum/empty-support display samples.
    synthetic = page.evaluate('''async()=>{const {finish}=await import('./kernel.js');return [finish([0,0,0],8),finish([1,1,1],20)];}''')
    np.testing.assert_allclose(synthetic, [[0, 0, 0], [1, 1, 1]], atol=1/255, rtol=0)
    return rows


def interaction_checks(page, screenshot_dir):
    checks = {}
    page.locator('#reset').click()
    evaluate_async(page, '()=>demo.update({preset:"opac",du:300,aod:.013,rh:.5})')
    before = page.evaluate('Array.from(demo.field.transmission)')
    page.locator('#ev').fill('12')
    page.locator('#ev').dispatch_event('input')
    page.locator('#moon').focus()
    page.keyboard.press('ArrowRight')
    checks['display_raw_unchanged'] = before == page.evaluate('Array.from(demo.field.transmission)')
    checks['keyboard_moved'] = page.evaluate('demo.state.x_m===2910000')
    # Actual Pointer Events, not merely changing DOM attributes.
    rect = page.locator('#moon').bounding_box()
    page.mouse.move(rect['x']+rect['width']/2, rect['y']+rect['height']/2)
    page.mouse.down()
    page.mouse.move(rect['x']+rect['width']/2+10, rect['y']+rect['height']/2+5)
    page.mouse.up()
    checks['pointer_moved'] = page.evaluate('demo.state.x_m>2910000 && demo.state.y_m>0')
    page.locator('#aod').fill('0')
    page.locator('#aod').dispatch_event('input')
    evaluate_async(page, '()=>demo.update({aod:0})')
    checks['zero_disables_rh'] = page.locator('#rh').is_disabled()
    zero = page.evaluate('Array.from(demo.field.transmission)')
    evaluate_async(page, '()=>demo.update({preset:"none"})')
    checks['zero_restores_none'] = zero == page.evaluate('Array.from(demo.field.transmission)')
    checks['none_disables_aod'] = page.locator('#aod').is_disabled()
    calibrated = page.evaluate('Object.keys(demo.manifest.presets).find(k=>k.startsWith("calibrated:"))')
    page.locator('#preset').select_option(calibrated)
    evaluate_async(page, '(p)=>demo.update({preset:p})', calibrated)
    checks['calibration_identity'] = page.locator('#aod').is_disabled() and page.locator('#rh').is_disabled() and '已拟合' in page.locator('#identity').inner_text()
    evaluate_async(page, '()=>demo.update({preset:"opac",aod:.013,rh:.7,du:300})')
    checks['nonzero_restores_rh'] = not page.locator('#rh').is_disabled()
    evaluate_async(page, '()=>demo.update({preset:"opac-ft"})')
    checks['ft_identity_excludes_stratosphere'] = '不含平流层' in page.locator('#identity').inner_text()
    try:
        evaluate_async(page, '()=>new Promise(()=>{})', timeout_ms=50)
    except Exception as error:
        checks['unsettled_promise_times_out'] = 'Browser evaluation deadline exceeded' in str(error)
    else:
        checks['unsettled_promise_times_out'] = False
    flood = evaluate_async(page, '''async()=>{
        const promises=[];
        for(let i=0;i<20;i++)promises.push(demo.update({du:i*30}));
        const q=demo.queue,results=await Promise.all(promises);
        return {bounded:q.active!==null&&q.pending!==null,final:demo.state.du===570,
            stale:results.filter(r=>r.superseded).length===19};
    }''')
    checks.update({f'queue_{k}': v for k, v in flood.items()})
    checks['invalid_rejected'] = evaluate_async(page, '''async()=>{
        for(const p of [{du:NaN},{du:601},{aod:-1},{aod:.0001},{rh:.49},{preset:'unknown'}]) {
            try{await demo.update(p);return false;}catch{}
        }
        try{demo.display({x_m:9000001});return false;}catch{}
        return true;
    }''')
    page.locator('#reset').click()
    evaluate_async(page, '()=>demo.update({preset:"opac",du:300,aod:.013,rh:.5})')
    page.screenshot(path=str(screenshot_dir/'desktop.png'), full_page=True)
    page.set_viewport_size({'width': 390, 'height': 844})
    checks['no_mobile_width_overflow'] = page.evaluate('document.documentElement.scrollWidth<=innerWidth')
    page.screenshot(path=str(screenshot_dir/'mobile-layout.png'), full_page=True)
    page.set_viewport_size({'width': 1280, 'height': 900})
    return checks


def fault_checks(browser, url, site, manifest):
    results = []
    default = dict(preset='opac', du=300, aod=.013, rh=.5)
    reference = production(default, radius_grid())
    cases = ['schema', 'missing', 'checksum', 'shape', 'W', 'aerosol', 'color']
    for fault in cases:
        m = json.loads(json.dumps(manifest))
        body = None
        descriptor = m['common'] if fault == 'color' else m['families'][m['presets']['opac']['family']]
        if fault == 'schema':
            m['schema_version'] = -1
        elif fault == 'shape':
            descriptor['fields']['components']['shape'][0] += 1
        elif fault in ('W', 'aerosol', 'color'):
            raw = bytearray(gzip.decompress((site/descriptor['file']).read_bytes()))
            field = descriptor['fields'][{'W': 'weight', 'aerosol': 'factors', 'color': 'rgbMatrix'}[fault]]
            values = np.frombuffer(raw, dtype=field['dtype'], offset=field['offset'], count=field['byteLength']//8)
            if fault == 'W':
                values[int(values.argmax())] *= 1.2
            elif fault == 'aerosol':
                values[2*81:3*81] *= 1.2
            else:
                values[0] *= 1.2
            body = gzip.compress(raw, mtime=0)
            descriptor['sha256'] = hashlib.sha256(body).hexdigest()
            descriptor['byteLength'] = len(body)
        context = browser.new_context()
        page = context.new_page()
        page.route('**/manifest.json', lambda route, request, m=m: route.fulfill(json=m))
        if fault in ('missing', 'checksum'):
            page.route('**/'+descriptor['file'], lambda route, request, f=fault: route.fulfill(status=404 if f == 'missing' else 200, body=b'broken asset'))
        elif body is not None:
            page.route('**/'+descriptor['file'], lambda route, request, body=body: route.fulfill(body=body, content_type='application/gzip'))
        page.goto(url)
        if fault in ('schema', 'missing', 'checksum', 'shape'):
            rejected = evaluate_async(page, 'async()=>{try{await demo.ready;return false;}catch{return true;}}')
            disabled = page.locator('#du').is_disabled()
            results.append({'fault': fault, 'ready_rejected': bool(rejected),
                            'controls_disabled': disabled, 'pass': bool(rejected and disabled)})
        else:
            await_ready(page)
            actual = browser_field(page, default)
            checks = quantities(actual, reference)
            results.append({'fault': fault, 'pass': not all(q['pass'] for q in checks.values()), 'detected_by': checks})
        context.close()
    return results


def buffer_inventory(manifest):
    decoded = sum(x['decodedByteLength'] for x in manifest['families'].values())
    largest = max(x['decodedByteLength'] for x in manifest['families'].values())
    largest_packed = max(x['byteLength'] for x in manifest['families'].values())
    nn = max(x['fields']['b']['shape'][0] for x in manifest['families'].values())
    nr = len(radius_grid())
    rows = {
        'all_decoded_families': decoded, 'worker_and_main_common': 2*manifest['common']['decodedByteLength'],
        'largest_exponential_workspace': nn*81*8, 'two_fields_during_transfer': 2*nr*(81+3+3+1)*8,
        'pixel_map': 560*560*12, 'three_RGBA_buffers': 3*560*560*4,
        # Conservatively include transient Blob/Response copies during preload.
        'preload_transient_decoded_copies': 2*largest, 'preload_transient_gzip_copies': 3*largest_packed,
    }
    total = sum(rows.values())
    return {'buffers_bytes': rows, 'conservative_peak_bytes': total, 'limit_bytes': 192*2**20,
            'pass': total <= 192*2**20, 'scope': 'Explicit numerical/image/transfer buffers, not whole browser RSS or JS/DOM heap'}


def check_site(site: Path, out: Path) -> dict:
    from playwright.sync_api import sync_playwright
    site = site.resolve()
    manifest = json.loads((site/'manifest.json').read_text())
    if manifest['inputs_sha256'] != input_fingerprints():
        raise ValueError('Stale site: source/data/lock fingerprint mismatch; regenerate before validation')
    radius = radius_grid()
    report = {'pass': False, 'mobile_real_device': 'NOT VERIFIED', 'checks': [],
              'host': platform.platform(), 'schema_version': manifest['schema_version'],
              'manifest_sha256': hashlib.sha256((site/'manifest.json').read_bytes()).hexdigest(),
              'inputs_sha256': manifest['inputs_sha256'], 'buffers': buffer_inventory(manifest)}
    def checkpoint(stage):
        report['stage'] = stage
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2)+'\n')
        print(f'WEB_GATE_STAGE: {stage}', flush=True)

    checkpoint('initialization')
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(QuietHandler, directory=str(site.parent)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    screenshot_dir = out.parent/(out.stem+'-screens')
    screenshot_dir.mkdir(parents=True, exist_ok=True)
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            report['browser_version'] = browser.version
            page = browser.new_page(viewport={'width': 1280, 'height': 900})
            errors = []
            page.on('pageerror', lambda error: errors.append(str(error)))
            url = f'http://127.0.0.1:{server.server_port}/{site.name}/'
            page.goto(url)
            await_ready(page)
            actual_radius = page.evaluate('Array.from(demo.common.radius)')
            np.testing.assert_array_equal(actual_radius, radius)
            np.testing.assert_array_equal(page.evaluate('Array.from(demo.common.wave)'), np.arange(380, 781, 5))
            for parameters in parameter_cases(manifest['presets']):
                reference = production(parameters, radius)
                actual = browser_field(page, parameters)
                checks = quantities(actual, reference)
                for name, quantity in checks.items():
                    quantity['worst_radius_m'] = float(radius[quantity['worst_index'][0]])
                    if name == 'transmission':
                        quantity['worst_wavelength_nm'] = float(380+5*quantity['worst_index'][1])
                report['checks'].append({'parameters': parameters, 'quantities': checks})
            checkpoint('same-node matrix complete')
            report['cli'] = cli_checks(page, manifest['presets'], radius)
            checkpoint('CLI complete')
            report['spatial'] = spatial_checks(page, manifest['presets'], radius)
            checkpoint('spatial complete')
            report['pixels'] = pixel_checks(page, radius)
            checkpoint('canvas pixels complete')
            report['interactions'] = interaction_checks(page, screenshot_dir)
            checkpoint('interactions complete')
            report['fault_injection'] = fault_checks(browser, url, site, manifest)
            checkpoint('fault injection complete')
            # This measures the real app/Worker/Canvas path, including cross-family
            # preset changes, not a standalone hot helper or low-resolution preview.
            report['timing'] = evaluate_async(page, '''async()=>{
                const samples=[],presets=Object.keys(demo.manifest.presets);
                for(let n=0;n<230;n++) {
                    const result=await demo.update({preset:presets[n%presets.length],du:n%3*300,
                        aod:n%7===0?0:.001+n%20*.01,rh:[.5,.7,.99][n%3]});
                    if(n>=30)samples.push(result.timing);
                }
                return samples;
            }''')
            values = np.array([s['total'] for s in report['timing']])
            report['timing_summary'] = {'p95_ms': float(np.percentile(values, 95)), 'max_ms': float(values.max()),
                                        'pass': bool(values.size == 200 and values.max() <= 100), 'samples': len(values)}
            report['display_timing'] = page.evaluate('''()=>{
                const samples=[];
                for(let n=0;n<230;n++) {
                    const start=performance.now();demo.display({ev:n%25-4,x_m:2860000+n*1000,y_m:0});
                    if(n>=30)samples.push(performance.now()-start);
                }
                return samples;
            }''')
            report['display_timing_p95_ms'] = float(np.percentile(report['display_timing'], 95))
            report['console_errors'] = errors
            report['pass'] = bool(not errors and report['buffers']['pass'] and report['timing_summary']['pass'] and
                report['display_timing_p95_ms'] <= 33 and passed_quantities(report['checks']) and
                passed_quantities(report['cli']) and passed_quantities(report['spatial']) and
                all(row['pass'] for row in report['pixels']) and all(report['interactions'].values()) and
                all(row['pass'] for row in report['fault_injection']))
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = check_site(args.site, args.out)
    except Exception as error:
        report = json.loads(args.out.read_text()) if args.out.exists() else {}
        report.update({'pass': False, 'error': f'{type(error).__name__}: {error}'})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k in ('pass', 'error', 'browser_version', 'timing_summary', 'display_timing_p95_ms')}))
    raise SystemExit(0 if report['pass'] else 1)


if __name__ == '__main__':
    main()
