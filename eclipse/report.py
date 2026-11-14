"""Reproducible scientific plots and neutral-disk eclipse visualizations.

The disk is a uniform, spectrally neutral screen sampling the shadow plane.
It deliberately omits lunar BRDF, incidence-angle shading, topography and the
observer's atmosphere. Photopic colorimetry does not predict dark-adapted vision.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .color import finish_display, linear_display_rgb, relative_photopic, spectrum_to_xyz
from .style import BASE_RC, HOUSE_STYLE
from .transport import MOON_RADIUS_M


_ON = "#ad502e"
_OFF = "#2d797a"


def _validate_inputs(radius_m, wavelength_nm, irradiance, reference_irradiance,
                     irradiance_no_ozone):
    radius = np.asarray(radius_m, dtype=float)
    wl = np.asarray(wavelength_nm, dtype=float)
    spectra = np.asarray(irradiance, dtype=float)
    reference = np.asarray(reference_irradiance, dtype=float)
    no_ozone = np.asarray(irradiance_no_ozone, dtype=float)
    if radius.ndim != 1 or radius.size < 2 or not np.all(np.isfinite(radius)):
        raise ValueError("radius_m must be a finite one-dimensional grid with at least two points")
    if radius[0] < 0 or not np.all(np.diff(radius) > 0):
        raise ValueError("radius_m must be nonnegative and strictly increasing")
    if spectra.shape != (radius.size, wl.size) or no_ozone.shape != spectra.shape:
        raise ValueError("irradiance arrays must have shape (n_radius, n_wavelength)")
    if reference.shape != (wl.size,):
        raise ValueError("reference_irradiance must have shape (n_wavelength,)")
    # Colorimetry checks wavelength order and finite, nonnegative spectra.
    spectrum_to_xyz(spectra, wl)
    spectrum_to_xyz(no_ozone, wl)
    if spectrum_to_xyz(reference, wl)[1] <= 0:
        raise ValueError("reference_irradiance must have positive CIE Y")
    return radius, wl, spectra, reference, no_ozone


def _disk_image(radius, linear_rgb, moon_offset_m, moon_radius_m, exposure_stops,
                size=560):
    """Interpolate linear light onto the screen before display encoding."""
    axis = np.linspace(-1.08, 1.08, size) * moon_radius_m
    x, y = np.meshgrid(axis, axis)
    distance = np.hypot(x + moon_offset_m, y)
    mask = x*x + y*y <= moon_radius_m*moon_radius_m
    if distance[mask].min() < radius[0] or distance[mask].max() > radius[-1]:
        raise ValueError("radius_m does not cover the requested lunar disk")
    mapped = np.stack([np.interp(distance, radius, linear_rgb[:, j]) for j in range(3)], axis=-1)
    mapped[~mask] = 0.0
    return finish_display(mapped, exposure_stops)


def _save(fig, path):
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return path


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def write_report(radius_m, wavelength_nm, irradiance, reference_irradiance,
                 irradiance_no_ozone, output_dir, *, moon_offset_m=2860000.0,
                 moon_radius_m=MOON_RADIUS_M, exposure_stops=8.0, metadata=None) -> list[Path]:
    """Write plots, common-exposure disk images, diagnostics and report.json.

    Input irradiance is a normal-to-axis, shadow-plane spectral irradiance in
    W m^-2 nm^-1. A single scene-wide solar white is used for every panel. The
    ozone-off comparison must otherwise use the identical atmosphere/geometry.
    """
    radius, wl, spectra, reference, no_ozone = _validate_inputs(
        radius_m, wavelength_nm, irradiance, reference_irradiance, irradiance_no_ozone
    )
    if not np.isfinite(moon_offset_m) or not np.isfinite(moon_radius_m) or moon_radius_m <= 0:
        raise ValueError("moon offset must be finite and moon radius must be positive")
    # Validate before creating output files; derived exposure panels remain
    # within the same supported range at the extreme user setting.
    finish_display(np.zeros(3), exposure_stops)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    relative_y = relative_photopic(spectra, wl, reference)
    relative_y_off = relative_photopic(no_ozone, wl, reference)
    rgb_linear = linear_display_rgb(spectra, wl, reference)
    off_linear = linear_display_rgb(no_ozone, wl, reference)
    diagnostic_linear = linear_display_rgb(spectra, wl, reference, normalize_luminance=True)
    moon = _disk_image(radius, rgb_linear, moon_offset_m, moon_radius_m, exposure_stops)
    moon_diagnostic = _disk_image(radius, diagnostic_linear, moon_offset_m, moon_radius_m, 0.0)
    moon_off = _disk_image(radius, off_linear, moon_offset_m, moon_radius_m, exposure_stops)
    files = []

    with plt.rc_context(HOUSE_STYLE):
        fig, axes = plt.subplots(4, 1, figsize=(10.6, 7.3), sharex=True,
                                 gridspec_kw={"height_ratios": [4.5, 0.5, 0.5, 0.5], "hspace": 0.36})
        ax = axes[0]
        ax.semilogy(radius / 1e3, np.maximum(relative_y, 1e-15), color=_ON, lw=2,
                    label="Standard ozone")
        ax.semilogy(radius / 1e3, np.maximum(relative_y_off, 1e-15), color=_OFF,
                    lw=1.6, ls="--", label="Ozone removed")
        ax.axhline(1, color="#aab1b8", lw=0.8)
        ax.set_ylabel("Photopic irradiance / uneclipsed sunlight")
        ax.set_ylim(max(1e-9, min(relative_y.min(), relative_y_off.min()) * 0.6), 1.5)
        ax.grid(axis="y", alpha=0.2)
        ax.legend(frameon=False, loc="lower right")
        ax.set_title("Across Earth's shadow at the Moon's distance", loc="left", fontsize=16, pad=14)
        extent = [radius[0] / 1e3, radius[-1] / 1e3, 0, 1]
        strip_radius = np.linspace(radius[0], radius[-1], max(1200, len(radius)))
        def strip_colors(linear, exposure=0):
            uniform = np.stack([np.interp(strip_radius, radius, linear[:, j])
                                for j in range(3)], axis=-1)
            return finish_display(uniform, exposure)
        strips = [
            (strip_colors(rgb_linear, exposure_stops), f"Ozone / {exposure_stops:+g} EV"),
            (strip_colors(off_linear, exposure_stops), f"No ozone / {exposure_stops:+g} EV"),
            (strip_colors(diagnostic_linear), "Ozone / equal Y*"),
        ]
        for ax, (strip, label) in zip(axes[1:], strips):
            ax.imshow(strip[None, :, :], extent=extent, aspect="auto", origin="lower")
            ax.set_yticks([])
            ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9, labelpad=14)
            ax.spines[["left", "bottom"]].set_visible(False)
        axes[-1].set_xlabel("Distance from shadow axis (km)")
        fig.text(0.125, 0.015, "*Equal-Y colors are a chromaticity diagnostic; brightness has been removed.",
                 fontsize=8.5, color="#68727b")
        files.append(_save(fig, out / "profile.png"))

        candidates = np.array([0.0, 2e6, 4e6, 4.6e6, 5.2e6, 6.5e6])
        indices = np.unique([np.abs(radius-r).argmin() for r in candidates if radius[0] <= r <= radius[-1]])
        if indices.size < 2:
            indices = np.unique(np.linspace(0, radius.size-1, 4).astype(int))
        colors = plt.colormaps["viridis"](np.linspace(0.1, 0.9, len(indices)))
        fig, axes = plt.subplots(2, 1, figsize=(10.6, 8.2), sharex=True)
        for idx, color in zip(indices, colors):
            label = f"r = {radius[idx]/1e3:,.0f} km"
            axes[0].semilogy(wl, np.maximum(spectra[idx], 1e-18), color=color, label=label, lw=1.6)
            transmission = np.divide(spectra[idx], reference, out=np.zeros_like(reference), where=reference > 0)
            axes[1].semilogy(wl, np.maximum(transmission, 1e-15), color=color, lw=1.6)
        axes[0].set_title("Spectra across the eclipse", loc="left", fontsize=16, pad=14)
        axes[0].set_ylabel(r"Spectral irradiance (W m$^{-2}$ nm$^{-1}$)")
        axes[0].legend(frameon=False, ncol=3, fontsize=9)
        axes[1].set_ylabel("Spectral irradiance / solar reference")
        axes[1].set_xlabel("Wavelength (nm)")
        for ax in axes:
            ax.grid(alpha=0.17)
        fig.tight_layout(h_pad=2.0)
        files.append(_save(fig, out / "spectra.png"))

        for name, image, title, subtitle in [
            ("moon.png", moon, f"Common exposure: {exposure_stops:+g} EV",
             "Photopic colors; neutral screen sampling the shadow plane"),
            ("moon_chromaticity.png", moon_diagnostic, "Equal-Y chromaticity diagnostic",
             "Brightness removed at every radius; not a photometric lunar image"),
        ]:
            fig, ax = plt.subplots(figsize=(6.8, 7.25))
            ax.imshow(image, origin="lower")
            ax.set_axis_off()
            ax.set_title(title, loc="left", fontsize=15, pad=13)
            fig.text(0.13, 0.07, subtitle, fontsize=8.5, color="#68727b")
            fig.text(0.13, 0.045, f"Moon center {moon_offset_m/1e3:,.0f} km from shadow axis; white = solar reference",
                     fontsize=8.5, color="#68727b")
            files.append(_save(fig, out / name))

        fig, axes = plt.subplots(1, 2, figsize=(11, 6))
        for ax, image, title in zip(axes, [moon, moon_off], ["Standard ozone", "Ozone removed"]):
            ax.imshow(image, origin="lower")
            ax.set_title(title, loc="left", fontsize=14)
            ax.set_axis_off()
        fig.suptitle(f"Ozone counterfactual / identical {exposure_stops:+g} EV exposure", fontsize=16, y=0.96)
        fig.text(0.13, 0.06, "Same density, refraction, solar source and display white; uniform neutral screen.",
                 fontsize=9, color="#68727b")
        files.append(_save(fig, out / "ozone_comparison.png"))

        exposures = [0.0, exposure_stops, min(exposure_stops + 4.0, 100.0)]
        fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.1))
        for ax, exposure in zip(axes, exposures):
            image = _disk_image(radius, rgb_linear, moon_offset_m, moon_radius_m, exposure)
            ax.imshow(image, origin="lower")
            ax.set_title(f"{exposure:+g} EV", loc="left", fontsize=14)
            ax.set_axis_off()
        fig.suptitle("One physical spectrum field, three global exposures", fontsize=16, y=0.99)
        fig.text(0.13, 0.04, "Above display range: common RGB-channel compression; photopic, solar-white-adapted visualization.",
                 fontsize=9, color="#68727b")
        files.append(_save(fig, out / "moon_exposures.png"))

    summary = {
        "model_metadata": metadata or {},
        "display": {
            "white_reference": "Uneclipsed solar spectrum, Bradford-adapted to D65",
            "exposure_stops": float(exposure_stops),
            "luminance_reference": "CIE Y of uneclipsed solar irradiance = 1",
            "highlights": "Common max-channel compression above linear RGB = 1",
            "gamut": "Negative linear sRGB channels clipped to zero",
            "vision": "Photopic CIE 1931; no mesopic or dark-adaptation model",
            "moon_geometry": "Uniform neutral screen at lunar distance; no lunar BRDF or incidence shading",
            "diagnostic": "moon_chromaticity.png normalizes CIE Y at each radius; not photometric",
        },
        "moon_offset_m": float(moon_offset_m),
        "moon_radius_m": float(moon_radius_m),
        "minimum_relative_photopic_irradiance": float(relative_y.min()),
        "axis_relative_photopic_irradiance": float(relative_y[0]) if radius[0] == 0 else None,
        "sampled_spectra": [
            {"radius_m": float(radius[i]), "relative_photopic_irradiance": float(relative_y[i]),
             "relative_photopic_irradiance_no_ozone": float(relative_y_off[i])}
            for i in indices
        ],
        "files": [p.name for p in files],
    }
    summary_path = out / "report.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=_json_value) + "\n")
    files.append(summary_path)
    return files
