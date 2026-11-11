"""Lightweight figure gates. Full CLI double-rendering is an on-demand check."""
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
from PIL import Image, PngImagePlugin
import pytest

from tools import render_mainline as render

ROOT = Path(__file__).resolve().parents[1]


def test_cli_state_configuration():
    for state in render.STATES:
        args = render.cli_args(state, "<raw>/state")
        assert args[0] == "simulate"
        assert args[args.index("--aerosol") + 1] == state
        assert ("--rh" in args) == (state != "none")
        if state != "none":
            assert args[args.index("--rh") + 1] == "0.5"
        assert "--aod550" not in args
        assert "--no-plots" not in args
        assert args[args.index("--exposure-stops") + 1] == "12"


def test_clean_png_preserves_pixels_and_removes_actual_chunks(tmp_path):
    original, clean = tmp_path / "original.png", tmp_path / "clean.png"
    pixels = np.arange(60, dtype=np.uint8).reshape(3, 5, 4)
    text = PngImagePlugin.PngInfo()
    text.add_text("Software", "volatile renderer")
    text.add_itxt("Date", "volatile wall clock")
    Image.fromarray(pixels).save(original, pnginfo=text, dpi=(170, 170))
    assert "tEXt" in render.png_chunks(original)
    render.clean_png(original, clean)
    assert set(render.png_chunks(clean)) == {"IHDR", "IDAT", "IEND"}
    with Image.open(clean) as result:
        assert result.mode == "RGBA"
        np.testing.assert_array_equal(result, pixels)
    again = tmp_path / "again.png"
    render.clean_png(original, again)
    assert clean.read_bytes() == again.read_bytes()


def test_recursive_volatile_removal():
    assert render.stable_value({"created_utc": "now", "keep": [
        {"simulation_seconds": 3, "nested": {"created_utc": "now", "value": 4}}
    ]}) == {"keep": [{"nested": {"value": 4}}]}


def test_spectrum_table_is_direct_division_at_exact_cli_radii(tmp_path):
    radius = np.array([0, 1000, 4e6, 4.6e6])
    wl = np.array([380, 500, 780])
    reference = np.array([2., 3., 4.])
    irradiance = np.arange(12).reshape(4, 3) * 0.001
    path = tmp_path / "spectra.npz"
    np.savez(path, radius_m=radius, wavelength_nm=wl,
             reference_irradiance_W_m2_nm=reference, irradiance_W_m2_nm=irradiance)
    table = render.read_spectra(path)
    np.testing.assert_array_equal(table["transmission"], irradiance[[0, 2, 3]] / reference)
    np.savez(path, radius_m=radius + 1, wavelength_nm=wl,
             reference_irradiance_W_m2_nm=reference, irradiance_W_m2_nm=irradiance)
    with pytest.raises(ValueError, match="does not contain radius"):
        render.read_spectra(path)


def test_failed_cli_never_publishes_or_changes_existing_output(tmp_path, monkeypatch):
    out = tmp_path / "bundle"
    commands = []

    def fail_cli(command, **kwargs):
        assert kwargs["cwd"] == ROOT
        commands.append(command)
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, str(ROOT / "eclipse/__init__.py"))
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(render.subprocess, "run", fail_cli)
    with pytest.raises(subprocess.CalledProcessError):
        render.generate(out)
    assert not out.exists()
    assert not list(tmp_path.glob(".mainline-*"))
    assert commands[-1][1:3] == ["-m", "eclipse"]
    out.mkdir()
    (out / "user.txt").write_text("keep")
    with pytest.raises(ValueError, match="empty or nonexistent"):
        render.generate(out)
    assert (out / "user.txt").read_text() == "keep"


def test_wrong_worktree_import_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(render.subprocess, "run", lambda *args, **kw:
                        subprocess.CompletedProcess(args, 0, "/elsewhere/eclipse/__init__.py"))
    with pytest.raises(RuntimeError, match="another worktree"):
        render.generate(tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()


def test_committed_bundle_provenance_and_inputs():
    directory = ROOT / "figures/mainline"
    meta = render.validate_bundle(directory)
    assert {p.name for p in directory.iterdir()} == {*render.FIGURES, "metadata.json", "README.md"}
    assert meta["inputs_sha256"] == render.input_fingerprints()
    # Resolve fingerprints against the source tree, not a self-reported HEAD.
    for path, fingerprint in meta["inputs_sha256"].items():
        assert not Path(path).is_absolute()
        assert render.sha256(ROOT / path) == fingerprint
    serialized = json.dumps(meta)
    assert "created_utc" not in serialized
    assert "simulation_seconds" not in serialized
    assert str(ROOT) not in serialized
    assert "/private/" not in serialized
    assert "/tmp/" not in serialized
    for state, record in meta["states"].items():
        model, report = record["model"], record["report"]
        assert model["moon_offset_m"] == report["moon_offset_m"] == 2860000
        assert model["geometry"] == meta["states"]["none"]["model"]["geometry"]
        assert model["atmosphere"] == meta["states"]["none"]["model"]["atmosphere"]
        assert model["numerics"]["radial_step_m"] == 20000
        assert model["numerics"]["radius_count"] == 451
        assert model["numerics"]["wavelength_count"] == 81
        assert model["numerics"]["height_step_m"] == 100
        assert model["numerics"]["wavelength_step_nm"] == 5
        table = record["spectral_table"]
        np.testing.assert_array_equal(table["wavelength_nm"], np.arange(380, 781, 5))
        assert table["radii_m"][-1] > model["geometric_umbra_radius_m"]
        values = np.asarray(table["transmission"])
        assert values.min() >= meta["plot"]["ylim"][0]
        assert values.max() <= meta["plot"]["ylim"][1]
        if state != "none":
            assert len(model["aerosol"]["layers"]) == (1 if state == "opac-ft" else 2)
    with Image.open(directory / render.STATES["none"]["figure"]) as im:
        size = im.size
    for spec in render.STATES.values():
        with Image.open(directory / spec["figure"]) as im:
            assert im.size == size


@pytest.fixture
def bundle_copy(tmp_path):
    out = tmp_path / "bundle"
    shutil.copytree(ROOT / "figures/mainline", out)
    return out


@pytest.mark.parametrize("defect", [
    "missing_figure", "missing_model", "missing_white", "missing_command",
    "wrong_cli_state", "wrong_model_state", "wrong_figure_state", "unequal_exposure",
    "missing_aerosol_source", "volatile_nested", "png_metadata",
])
def test_bundle_gate_rejects_corruptions(bundle_copy, defect):
    path = bundle_copy / "metadata.json"
    meta = json.loads(path.read_text())
    state = meta["states"]["opac"]
    if defect == "missing_figure":
        (bundle_copy / "l1_molecular.png").unlink()
    elif defect == "missing_model":
        del state["model"]
    elif defect == "missing_white":
        del state["report"]["display"]["white_reference"]
    elif defect == "missing_command":
        del state["command"]
    elif defect == "wrong_cli_state":
        state["argv"][state["argv"].index("--aerosol") + 1] = "none"
    elif defect == "wrong_model_state":
        state["model"]["aerosol_preset"] = "opac-ft"
    elif defect == "wrong_figure_state":
        meta["figures"]["l25_stratosphere.png"]["states"] = ["opac-ft"]
    elif defect == "unequal_exposure":
        state["report"]["display"]["exposure_stops"] = 8
    elif defect == "missing_aerosol_source":
        del state["model"]["aerosol"]["layers"][0]["source"]
    elif defect == "volatile_nested":
        state["model"]["created_utc"] = "not deterministic"
    elif defect == "png_metadata":
        image_path = bundle_copy / "l1_molecular.png"
        info = PngImagePlugin.PngInfo()
        info.add_text("Software", "not stripped")
        with Image.open(image_path) as im:
            im.save(image_path, pnginfo=info)
        # Defeat the hash gate deliberately to test the independent chunk gate.
        meta["figures"][image_path.name]["sha256"] = render.sha256(image_path)
    path.write_text(json.dumps(meta))
    with pytest.raises((ValueError, KeyError, FileNotFoundError)):
        render.validate_bundle(bundle_copy)
