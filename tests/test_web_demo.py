"""Lightweight exporter/contract tests; the real browser gate runs separately."""
import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from eclipse.aerosol import AerosolColumn, ft_growth_factor
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.transport import Geometry, receiver_quadrature
from tools.export_web_demo import (
    MAX_DATA_BYTES, export_site, pack_arrays, preset_columns, radius_grid, sha256,
)
from validation.check_web_demo import buffer_inventory, compare, parameter_cases


@pytest.fixture(scope='module')
def exported(tmp_path_factory):
    out = tmp_path_factory.mktemp('web')/'site'
    manifest = export_site(out)
    return out, manifest


def unpack(site, descriptor):
    encoded = (site/descriptor['file']).read_bytes()
    assert sha256(encoded) == descriptor['sha256']
    assert len(encoded) == descriptor['byteLength']
    raw = gzip.decompress(encoded)
    assert len(raw) == descriptor['decodedByteLength']
    fields = {}
    for name, field in descriptor['fields'].items():
        dtype = np.dtype(field['dtype'])
        count = int(np.prod(field['shape']))
        assert count*dtype.itemsize == field['byteLength']
        fields[name] = np.frombuffer(raw, dtype=dtype, offset=field['offset'], count=count).reshape(field['shape'])
    return fields


def test_export_contract_and_exact_sharing(exported):
    site, m = exported
    assert m['schema_version'] == 1
    assert m['data_bytes'] <= MAX_DATA_BYTES
    assert buffer_inventory(m)['pass']
    common = unpack(site, m['common'])
    np.testing.assert_array_equal(common['radius'], radius_grid())
    np.testing.assert_array_equal(common['wave'], np.arange(380, 781, 5))
    assert common['reference'].shape == (81,)
    assert m['defaults']['preset'] == 'opac'
    assert set(m['presets']) == set(preset_columns())
    assert len(m['families']) == 4
    assert len({p['family'] for n, p in m['presets'].items() if n.startswith('calibrated:')}) == 1
    for rh, value in m['rh_knots']:
        assert value == ft_growth_factor(rh)
    for name, p in m['presets'].items():
        family = unpack(site, m['families'][p['family']])
        assert family['components'].shape[1] == len(p['amplitudes']) == len(p['layers'])+2
        assert family['index'].shape == family['weight'].shape
        assert family['offsets'][-1] == len(family['index'])
        assert np.all(family['weight'] >= 0)
        assert family['outer'][-1] == 1
        assert family['offsets'][-1] == family['offsets'][-2]  # Vacuum: no supported atmospheric ring.
        if name.startswith('calibrated:'):
            assert p['kind'] == 'calibrated'
            assert p['provenance']['calibrated_not_physical_prior']
            assert p['provenance']['sha256']


def test_export_byte_determinism_and_only_out(exported, tmp_path):
    site, _ = exported
    second = tmp_path/'second'
    export_site(second)
    files = sorted(p.relative_to(site) for p in site.rglob('*') if p.is_file())
    assert files == sorted(p.relative_to(second) for p in second.rglob('*') if p.is_file())
    for relative in files:
        assert (site/relative).read_bytes() == (second/relative).read_bytes(), relative
    assert sorted(p.name for p in tmp_path.iterdir()) == ['second']


def test_export_rejects_existing_output_before_writing(tmp_path):
    out = tmp_path/'not-empty'
    out.mkdir()
    marker = out/'keep.txt'
    marker.write_text('do not overwrite')
    with pytest.raises(ValueError, match='empty'):
        export_site(out)
    assert marker.read_text() == 'do not overwrite'
    link = tmp_path/'link'
    link.symlink_to(out, target_is_directory=True)
    with pytest.raises(ValueError, match='symlink'):
        export_site(link)
    marker2 = tmp_path/'file'
    marker2.write_text('preserve')
    with pytest.raises(ValueError):
        export_site(marker2)
    assert marker2.read_text() == 'preserve'


def test_array_endian_offsets_and_finite_shape():
    raw, fields = pack_arrays({'ids': np.array([1, 2, 3]), 'weights': np.array([[.125, .5]])})
    assert fields['ids']['dtype'] == '<u4'
    assert fields['weights']['offset'] == 16
    assert raw[12:16] == b'\0'*4
    np.testing.assert_array_equal(np.frombuffer(raw, dtype='<f8', offset=16), [.125, .5])


def test_geometry_invariance_is_scoped_to_active_grid_family():
    base = build_ray_table(Atmosphere(), aerosol_profile=AerosolColumn.opac_default())
    changed = build_ray_table(Atmosphere(ozone_du=600), aerosol_profile=AerosolColumn.opac_default(.2, .99))
    for field in ('altitude_m', 'impact_m', 'bending_rad'):
        np.testing.assert_array_equal(getattr(base, field), getattr(changed, field))
    for r in (0, 2e6, 4.65e6, 8e6, 9e6):
        old = receiver_quadrature(base, r, Geometry())
        new = receiver_quadrature(changed, r, Geometry())
        for a, b in zip(old, new, strict=True):
            np.testing.assert_array_equal(a, b)
    none = build_ray_table(Atmosphere())
    assert len(none.altitude_m) != len(base.altitude_m)
    _, a, b = np.intersect1d(base.altitude_m, none.altitude_m, return_indices=True)
    np.testing.assert_array_equal(base.impact_m[a], none.impact_m[b])
    np.testing.assert_array_equal(base.bending_rad[a], none.bending_rad[b])


def test_cross_gate_never_accepts_zero_or_wrong_samples():
    with pytest.raises(ValueError, match='empty'):
        compare([], np.array([]), 2e-12, 2e-8)
    with pytest.raises(ValueError, match='mismatched'):
        compare([1], np.array([1, 2]), 2e-12, 2e-8)
    with pytest.raises(ValueError, match='nonfinite'):
        compare([np.nan], np.array([1]), 2e-12, 2e-8)
    assert not compare([1.01], np.array([1.]), 2e-12, 2e-8)['pass']
    cases = parameter_cases(preset_columns())
    assert len(cases) == 87
    assert {row['preset'] for row in cases} == set(preset_columns())
    assert json.loads(json.dumps(cases)) == cases


def test_workflow_failure_blocks_deploy_and_all_core_changes_are_covered():
    import yaml
    root = Path(__file__).resolve().parents[1]
    workflow = yaml.load((root/'.github/workflows/web-demo.yml').read_text(), Loader=yaml.BaseLoader)
    triggers = workflow['on']
    assert 'pull_request' in triggers and triggers['push']['branches'] == ['main']
    assert 'paths' not in triggers['pull_request'] if isinstance(triggers['pull_request'], dict) else True
    assert 'paths' not in triggers['push']
    validate, deploy = workflow['jobs']['validate'], workflow['jobs']['deploy']
    assert deploy['needs'] == 'validate'
    assert 'always()' not in deploy['if']
    assert "github.ref == 'refs/heads/main'" in deploy['if']
    assert deploy['permissions'] == {'contents': 'read', 'pages': 'write', 'id-token': 'write'}
    assert workflow['permissions'] == {'contents': 'read'}
    assert deploy['environment']['name'] == 'github-pages'
    assert deploy['concurrency']['cancel-in-progress'] == 'false'
    gate = [step for step in validate['steps'] if 'validation.check_web_demo' in step.get('run', '')]
    assert len(gate) == 1 and 'continue-on-error' not in gate[0]
    # With needs' default success guard, a failed gate cannot make a deploy job
    # eligible. No always()/continue-on-error path bypasses this dependency.
    assert all(step.get('continue-on-error', 'false') == 'false' for step in validate['steps'])
    commands = '\n'.join(step.get('run', '') for step in deploy['steps'])
    assert 'export_web_demo.py' not in commands
    assert 'source_sha' in commands and 'manifest_sha256' in commands and "report['pass']" in commands
    for job in (validate, deploy):
        for step in job['steps']:
            if 'uses' in step:
                action, revision = step['uses'].split('@')
                assert len(revision) == 40 and all(c in '0123456789abcdef' for c in revision), action
