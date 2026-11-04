"""Read seven observed continuum points from Yan et al. Fig. 15.

One-time provenance tool (PyMuPDF is not a simulation dependency):
uv run --with pymupdf==2.1 python tools/extract_yan2011.py

The fixed PDF embeds a raster figure. Only neutral black pixels of the
OBSERVED curve are used; the red/yellow/blue/green model curves are excluded.
These are graph readouts, not released numerical spectra or original exposures.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pymupdf


URL = "https://arxiv.org/pdf/1405.4780v1"
PDF_SHA256 = "bf6771ce3b4a01ecd80e195d2f1ca82897f10bf78ee51a8690295249024b2105"
# Coordinates refer to the native 4096 x 2037 RGB image, origin at upper left.
IMAGE_XREF = 273
X_400_NM, X_1000_NM = 359.5, 3736.5
Y_1E_MINUS2, PIXELS_PER_DECADE = 365.5, 375.0
# Narrow windows, independently checked against the figure, exclude text,
# ticks, the green residual, and the blue skylight model.
READOUTS = {450: 1101, 500: 989, 550: 994, 600: 1010,
            650: 816, 700: 694, 750: 644}


def extract(pdf: Path, out: Path, qa: Path | None = None):
    digest = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if digest != PDF_SHA256:
        raise ValueError(f"PDF differs from audited source: {digest}")
    document = pymupdf.open(pdf)
    if len(document) != 12 or IMAGE_XREF not in {x[0] for x in document[10].get_images()}:
        raise ValueError("Unexpected Fig. 15 source structure")
    pixmap = pymupdf.Pixmap(document, IMAGE_XREF)
    if (pixmap.width, pixmap.height, pixmap.n) != (4096, 2037, 3):
        raise ValueError("Unexpected native figure dimensions/color space")
    rgb = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(2037, 4096, 3)
    black = rgb.max(axis=2) < 80
    records = []
    for wavelength, y_window_center in READOUTS.items():
        x = round(X_400_NM + (wavelength-400)/(1000-400)*(X_1000_NM-X_400_NM))
        y_start, y_stop = y_window_center-30, y_window_center+31
        y_hits = np.where(black[y_start:y_stop, x-3:x+4])[0] + y_start
        if not len(y_hits) or y_hits.min() == y_start or y_hits.max() == y_stop-1:
            raise ValueError(f"Observed curve absent or clipped at {wavelength} nm")
        y = float(np.median(y_hits))
        # Include the full observed stroke/window spread plus two raster pixels.
        halfwidth_y = max(y-y_hits.min(), y_hits.max()-y) + 2.0
        records.append({
            "wavelength_nm": wavelength,
            "observed_relative_brightness": 10**(-2-(y-Y_1E_MINUS2)/PIXELS_PER_DECADE),
            "graph_readout_wavelength_halfwidth_nm": 3.5*600/(X_1000_NM-X_400_NM),
            "graph_readout_log10_halfwidth": halfwidth_y/PIXELS_PER_DECADE,
            "image_x_pixel": x, "image_y_pixel": y,
            "stroke_y_min_pixel": int(y_hits.min()), "stroke_y_max_pixel": int(y_hits.max()),
        })
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    if qa:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.imshow(rgb)
        for row in records:
            x, y = row["image_x_pixel"], row["image_y_pixel"]
            ax.plot(x, y, marker="o", ms=10, mfc="none", mec="magenta")
            ax.annotate(str(row["wavelength_nm"]), (x, y), xytext=(0, -22),
                        textcoords="offset points", color="magenta", ha="center")
        ax.set_axis_off()
        qa.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(qa, dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(f"Extracted {len(records)} observed graph readouts to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=Path("outputs/sources/yan2011dec.pdf"))
    parser.add_argument("--out", type=Path,
                        default=Path("validation/observations/yan2011dec/points.csv"))
    parser.add_argument("--qa", type=Path)
    args = parser.parse_args()
    if not args.pdf.exists():
        args.pdf.parent.mkdir(parents=True, exist_ok=True)
        with urlopen(URL, timeout=60) as response:
            args.pdf.write_bytes(response.read())
    extract(args.pdf, args.out, args.qa)


if __name__ == "__main__":
    main()
