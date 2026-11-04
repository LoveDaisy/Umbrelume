"""Separate ray-table and spectral-resolution checks of the baseline model."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from eclipse.color import linear_display_rgb, spectrum_to_xyz
from eclipse.transport import simulate


def compare(coarse, fine):
    def quantities(result):
        y = spectrum_to_xyz(result.irradiance, result.wavelength_nm)[:, 1]
        y /= spectrum_to_xyz(result.reference_irradiance, result.wavelength_nm)[1]
        rgb = linear_display_rgb(result.irradiance, result.wavelength_nm,
                                 result.reference_irradiance, normalize_luminance=True)
        return y, rgb
    y_coarse, rgb_coarse = quantities(coarse)
    y_fine, rgb_fine = quantities(fine)
    relative = np.abs(y_coarse/y_fine - 1)
    color = np.max(np.abs(rgb_coarse-rgb_fine), axis=1)
    return {
        "pass": bool(relative.max() < 0.01 and color.max() < 0.02),
        "max_relative_photopic_irradiance_difference": float(relative.max()),
        "photopic_worst_radius_m": float(coarse.radius_m[np.argmax(relative)]),
        "max_equal_Y_linear_sRGB_channel_difference": float(color.max()),
        "color_worst_radius_m": float(coarse.radius_m[np.argmax(color)]),
        "photopic_tolerance": 0.01,
        "color_tolerance": 0.02,
    }


def run_checks():
    # Dense sampling around the faint cold edge, plus the entire shadow.
    radius = np.unique(np.r_[np.arange(0, 9_000_001, 20_000),
                             np.arange(4_400_000, 4_800_001, 2_000)])
    wl = np.arange(380, 781, 5)
    base = simulate(radius, wl, height_step_m=100)
    fine = simulate(radius, wl, height_step_m=50)
    finer = simulate(radius, wl, height_step_m=25)
    spectral = simulate(radius, np.arange(380, 780.1, 2.5), height_step_m=50)
    checks = {"height_100m_vs_50m": compare(base, fine),
              "height_50m_vs_25m": compare(fine, finer),
              "spectrum_5nm_vs_2p5nm": compare(fine, spectral)}
    return {"pass": all(value["pass"] for value in checks.values()),
            "checks": checks,
            "radius_count": len(radius),
            "height_spacing_note": "Above 60 km the table uses 500 m spacing in all runs.",
            "scope": "Numerical resolution at sampled receiver radii; excludes physical-model error, "
                     "finite-angle error and screen-pixel interpolation error."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("outputs/validation/convergence.json"))
    args = parser.parse_args()
    result = run_checks()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["pass"] else 1)


if __name__ == "__main__":
    main()
