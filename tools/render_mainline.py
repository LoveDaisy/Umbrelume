"""Offline mainline figures, using only the current worktree's production CLI.

Run from the repository root: uv run python tools/render_mainline.py --out <new-dir>
Only five PNGs and one deterministic provenance document are published. No
observational pipeline or solver is imported here; staging owns all CLI writes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import tempfile
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, ft2font
import numpy as np
from PIL import Image

from eclipse.report import HOUSE_STYLE

ROOT = Path(__file__).resolve().parents[1]
# One configuration expands into execution, provenance and test expectations.
COMMON_ARGS = (
    "--ozone-du", "300", "--density-scale", "1",
    "--moon-distance-km", "384400", "--moon-offset-km", "2860",
    "--exposure-stops", "12", "--height-step-m", "100",
    "--radial-step-km", "20", "--radius-max-km", "9000",
    "--wavelength-step-nm", "5",
)
STATES = {
    "none": {"figure": "l1_molecular.png", "label": "L1 / molecular",
             "color": "#2a78d6", "linestyle": "-"},
    "opac-ft": {"figure": "l2_free_troposphere.png", "label": "L2 / + free troposphere",
                "color": "#eb6834", "linestyle": "--"},
    "opac": {"figure": "l25_stratosphere.png", "label": "L2.5 / + stratosphere",
             "color": "#1baf7a", "linestyle": ":"},
}
SAMPLE_RADII_M = (0, 4_000_000, 4_600_000)
PLOT = {"radii_m": list(SAMPLE_RADII_M), "quantity": "E_lambda / F_lambda",
        "x_unit": "nm", "y_unit": "dimensionless", "yscale": "log",
        "ylim": [1e-12, 1.5], "uneclipsed_reference": 1,
        "normalization": "none; no peak or equal-Y normalization",
        "figsize_inches": [12.4, 5.4], "dpi": 170,
        "font": "DejaVu Sans", "series": STATES}
FIGURES = {spec["figure"]: {"states": [state], "source": "moon.png"}
           for state, spec in STATES.items()}
FIGURES.update({
    "shadow_profile.png": {"states": ["opac"], "source": "profile.png"},
    "transmission_ladder.png": {"states": list(STATES), "source": "spectra.npz"},
})
VOLATILE = {"created_utc", "simulation_seconds"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stable_value(value):
    """Remove volatile keys recursively, including report.model_metadata copies."""
    if isinstance(value, dict):
        return {k: stable_value(v) for k, v in value.items() if k not in VOLATILE}
    if isinstance(value, list):
        return [stable_value(v) for v in value]
    return value


def cli_args(state, out):
    args = ["simulate", *COMMON_ARGS, "--aerosol", state]
    if state != "none":
        args += ["--rh", "0.5"]
    return [*args, "--out", str(out)]


def clean_png(source, target):
    """Re-encode pixels, not display transforms; discard all ancillary metadata."""
    with Image.open(source) as original:
        if original.mode not in ("RGB", "RGBA"):
            raise ValueError(f"expected RGB/RGBA PNG: {source}")
        clean = Image.frombytes(original.mode, original.size, original.tobytes())
        clean.save(target, format="PNG", optimize=False, compress_level=9)
        with Image.open(target) as result:
            if result.size != original.size or result.tobytes() != original.tobytes():
                raise ValueError("PNG cleaning changed pixels")


def png_chunks(path):
    """Read actual chunk types, not only the Pillow-exposed metadata subset."""
    data = Path(path).read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"not a PNG: {path}")
    chunks, offset = [], 8
    while offset < len(data):
        length = int.from_bytes(data[offset:offset+4], "big")
        kind = data[offset+4:offset+8].decode("ascii")
        chunks.append(kind)
        offset += 12 + length
    if offset != len(data) or not chunks or chunks[-1] != "IEND":
        raise ValueError(f"invalid PNG chunks: {path}")
    return chunks


def read_spectra(path):
    """Expose the exact plotted values as a numeric table in metadata as well."""
    with np.load(path, allow_pickle=False) as data:
        radius = data["radius_m"]
        wl = data["wavelength_nm"]
        reference = data["reference_irradiance_W_m2_nm"]
        irradiance = data["irradiance_W_m2_nm"]
        if (irradiance.shape != (len(radius), len(wl)) or reference.shape != wl.shape
                or not np.all(np.isfinite(irradiance)) or np.any(irradiance < 0)
                or not np.all(np.isfinite(reference)) or np.any(reference <= 0)
                or not np.all(np.isfinite(wl)) or not np.all(np.diff(wl) > 0)):
            raise ValueError("invalid CLI spectral arrays")
        indices = []
        for sample in SAMPLE_RADII_M:
            matches = np.flatnonzero(radius == sample)
            if len(matches) != 1:
                raise ValueError(f"CLI grid does not contain radius {sample}")
            indices.append(int(matches[0]))
        return {"wavelength_nm": wl.tolist(), "radii_m": list(SAMPLE_RADII_M),
                "transmission": (irradiance[indices] / reference).tolist()}


def plot_transmission(tables, umbra_m, target):
    with plt.rc_context(HOUSE_STYLE):
        fig, axes = plt.subplots(1, 3, figsize=PLOT["figsize_inches"], sharey=True)
        for panel, (ax, radius) in enumerate(zip(axes, SAMPLE_RADII_M, strict=True)):
            for state, spec in STATES.items():
                table = tables[state]
                ax.semilogy(table["wavelength_nm"], table["transmission"][panel],
                            color=spec["color"], ls=spec["linestyle"], lw=2,
                            label=spec["label"])
            ax.axhline(1, color="#898781", lw=0.8)
            ax.text(0.03, 0.96, "Uneclipsed sunlight: T = 1", transform=ax.transAxes,
                    va="top", fontsize=8, color="#52514e")
            region = "umbra" if radius < umbra_m else "penumbra"
            ax.set_title(f"r = {radius/1000:,.0f} km / {region}", loc="left", fontsize=11)
            ax.set_ylim(*PLOT["ylim"])
            ax.set_xlim(380, 780)
            ax.set_xlabel("Wavelength (nm)")
            ax.grid(axis="y", color="#e1e0d9", lw=0.5)
        axes[0].set_ylabel("Spectral irradiance / uneclipsed sunlight (E / F)")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.91),
                   frameon=False, ncol=3)
        fig.suptitle("Three physical states, one absolute transmission scale",
                     x=0.065, ha="left", fontsize=16)
        fig.text(0.065, 0.035,
                 "Fixed radial samples, not event sites. No per-curve normalization. "
                 "Line styles also identify states; numeric table in metadata.json.",
                 fontsize=8, color="#52514e")
        fig.subplots_adjust(left=0.075, right=0.985, bottom=0.17, top=0.77, wspace=0.12)
        fig.savefig(target, dpi=PLOT["dpi"])
        plt.close(fig)


def input_fingerprints():
    # All production source files, actual runtime tables, and the generation tool.
    paths = [*sorted((ROOT / "eclipse").glob("*.py")),
             *sorted((ROOT / "eclipse/data").glob("*.csv")),
             ROOT / "eclipse/data/SOURCES.md", ROOT / "uv.lock", Path(__file__).resolve()]
    return {str(path.relative_to(ROOT)): sha256(path) for path in paths}


def software_versions():
    return {"python": platform.python_version(), "platform": platform.system(),
            "machine": platform.machine(), "zlib": zlib.ZLIB_RUNTIME_VERSION,
            "freetype": ft2font.__freetype_version__,
            "font_sha256": sha256(font_manager.findfont(PLOT["font"])),
            **{name: importlib.metadata.version(name)
               for name in ("numpy", "scipy", "matplotlib", "pillow")}}


def validate_bundle(directory):
    """Shared publication/readback gate; never runs a solver."""
    directory = Path(directory)
    meta = json.loads((directory / "metadata.json").read_text())
    if (meta["schema"] != "mainline-figures-v1" or meta["plot"] != PLOT
            or set(meta["states"]) != set(STATES)
            or set(meta["figures"]) != set(FIGURES)
            or not meta["reproduce_command"] or not meta["inputs_sha256"]
            or not meta["software"]):
        raise ValueError("incomplete bundle provenance")
    common_display = meta["states"]["none"]["report"]["display"]
    common_geometry = meta["states"]["none"]["model"]["geometry"]
    for state, record in meta["states"].items():
        if record["argv"] != cli_args(state, f"<raw>/{state}"):
            raise ValueError(f"CLI parameter mismatch: {state}")
        if record["command"] != shlex.join(["uv", "run", "lunar-eclipse", *record["argv"]]):
            raise ValueError(f"CLI command mismatch: {state}")
        model, report = record["model"], record["report"]
        if (model["aerosol_preset"] != state or report["display"] != common_display
                or model["geometry"] != common_geometry
                or report["display"]["exposure_stops"] != 12.0
                or not all(report["display"].get(k) for k in (
                    "white_reference", "luminance_reference", "highlights", "gamut",
                    "vision", "moon_geometry", "diagnostic"))
                or not model["atmosphere"] or not model["numerics"]
                or not model["approximations"]
                or record["overrides"] != {"aod550": None,
                                          "rh": None if state == "none" else 0.5}):
            raise ValueError(f"model/display provenance mismatch: {state}")
        if state == "none":
            if model["aerosol"] is not None or record["rh_semantics"] != "not applicable":
                raise ValueError("molecular state cannot carry aerosol or RH")
        elif (not model["aerosol"]["layers"]
              or model["aerosol_aod_source"] != "OPAC literature prior"
              or not all(layer["source"] for layer in model["aerosol"]["layers"])):
            raise ValueError(f"aerosol provenance missing: {state}")
        table = record["spectral_table"]
        if table["radii_m"] != list(SAMPLE_RADII_M):
            raise ValueError("spectral radius mismatch")
        values = np.asarray(table["transmission"])
        if (values.shape != (3, len(table["wavelength_nm"]))
                or not np.all(np.isfinite(values)) or np.any(values < 0)):
            raise ValueError("invalid transmission table")
    for name, provenance in FIGURES.items():
        record = meta["figures"][name]
        if any(record[k] != value for k, value in provenance.items()):
            raise ValueError(f"figure state mismatch: {name}")
        if record["sha256"] != sha256(directory / name):
            raise ValueError(f"figure hash mismatch: {name}")
        if set(png_chunks(directory / name)) != {"IHDR", "IDAT", "IEND"}:
            raise ValueError(f"PNG metadata not stripped: {name}")
    if stable_value(meta) != meta:
        raise ValueError("volatile metadata in bundle")
    return meta


def generate(out):
    out = Path(out).absolute()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"output must be empty or nonexistent: {out}")
    # Check the actual interpreter import, not an assumed editable installation.
    probe = subprocess.run([sys.executable, "-c", "import eclipse; print(eclipse.__file__)"],
                           cwd=ROOT, check=True, capture_output=True, text=True)
    if Path(probe.stdout.strip()).resolve() != ROOT / "eclipse/__init__.py":
        raise RuntimeError("interpreter imports eclipse from another worktree")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".mainline-", dir=out.parent) as tmp:
        staging = Path(tmp)
        publish = staging / "publish"
        publish.mkdir()
        records, tables = {}, {}
        # Keep matplotlib cache/config writes inside this invocation's staging too.
        env = {**os.environ, "MPLCONFIGDIR": str(staging / "mpl-cache"),
               "PYTHONHASHSEED": "0"}
        with matplotlib.rc_context():
            for state, spec in STATES.items():
                raw = staging / state
                subprocess.run([sys.executable, "-m", "eclipse", *cli_args(state, raw)],
                               cwd=ROOT, env=env, check=True)
                model = stable_value(json.loads((raw / "metadata.json").read_text()))
                report = stable_value(json.loads((raw / "report.json").read_text()))
                if report.pop("model_metadata") != model:
                    raise ValueError("CLI report and metadata disagree")
                report.pop("files")  # Unpublished helper PNGs are not bundle references.
                tables[state] = read_spectra(raw / "spectra.npz")
                argv = cli_args(state, f"<raw>/{state}")
                records[state] = {"argv": argv,
                                  "command": shlex.join(["uv", "run", "lunar-eclipse", *argv]),
                                  "overrides": {"aod550": None,
                                                "rh": None if state == "none" else 0.5},
                                  "rh_semantics": "not applicable" if state == "none" else
                                                  "explicit published 50% reference; FT only",
                                  "model": model, "report": report,
                                  "spectral_table": tables[state]}
                clean_png(raw / "moon.png", publish / spec["figure"])
                if state == "opac":
                    clean_png(raw / "profile.png", publish / "shadow_profile.png")
            if any(table["wavelength_nm"] != tables["none"]["wavelength_nm"]
                   for table in tables.values()):
                raise ValueError("CLI states have different wavelength grids")
            plot_transmission(tables, records["opac"]["model"]["geometric_umbra_radius_m"],
                              staging / "transmission.png")
            clean_png(staging / "transmission.png", publish / "transmission_ladder.png")
        meta = {"schema": "mainline-figures-v1",
                "reproduce_command": "uv run python tools/render_mainline.py --out <new-dir>",
                "states": records, "plot": PLOT, "inputs_sha256": input_fingerprints(),
                "software": software_versions(),
                "figures": {name: {**source, "sha256": sha256(publish / name)}
                            for name, source in FIGURES.items()}}
        (publish / "metadata.json").write_text(
            json.dumps(meta, indent=2, sort_keys=True, allow_nan=False) + "\n")
        validate_bundle(publish)
        if out.exists():
            out.rmdir()  # Only the previously checked empty directory; never recursive.
        publish.rename(out)  # Same filesystem: no half-published bundle on CLI failure.
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    print(f"Published: {generate(args.out)}")


if __name__ == "__main__":
    main()
