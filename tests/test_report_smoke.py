"""End-to-end smoke test of the plotting CLI path.

Every other CLI test runs with ``--no-plots``; this is the only mechanical
check that report.py's matplotlib rendering path executes and produces the
full published file set with a sane report.json.
"""
import json
import subprocess
import sys

EXPECTED_PLOTS = {"profile.png", "spectra.png", "moon.png", "moon_chromaticity.png",
                  "ozone_comparison.png", "moon_exposures.png"}


def test_cli_full_run_writes_plots_and_report(tmp_path):
    out = tmp_path / "run"
    # Coarse grid keeps the full plotting run within the pytest time budget;
    # the default 9000 km radius covers the offset lunar disk.
    result = subprocess.run(
        [sys.executable, "-m", "eclipse", "simulate", "--out", str(out),
         "--radial-step-km", "200", "--wavelength-step-nm", "40",
         "--height-step-m", "500"],
        capture_output=True, text=True, timeout=60, check=True)
    assert "Results:" in result.stdout
    report = json.loads((out / "report.json").read_text())
    assert EXPECTED_PLOTS == set(report["files"]) - {"report.json"}
    for name in EXPECTED_PLOTS:
        assert (out / name).stat().st_size > 0
    meta = report["model_metadata"]
    assert isinstance(meta["version"], str) and meta["version"]
    for field in ("moon_offset_m", "moon_radius_m",
                  "minimum_relative_photopic_irradiance"):
        assert isinstance(report[field], float)
    assert set(report["display"]) >= {"white_reference", "exposure_stops",
                                      "luminance_reference"}
