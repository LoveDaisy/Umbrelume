"""Lightweight figure gates. Full CLI double-rendering is an on-demand check."""
import json
from pathlib import Path
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
