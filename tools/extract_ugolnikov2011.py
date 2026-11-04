"""Read observed curves from Fig. 2 of the fixed arXiv 1106.6178v2 PDF.

One-time provenance tool, not a runtime dependency:
uv run --with pymupdf==2.1 python tools/extract_ugolnikov2011.py
Only published OBSERVED solid curves are extracted. The PDF's dashed model
curves are not reference measurements. Original exposures are not available here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import pymupdf


URL = "https://arxiv.org/pdf/1106.6178v2"
SHA256 = "366b6bd302937530f52d8fc3768be26ebb75657dba5072a3bd07d0733b67c043"
AXES = {"x_zero_pt": 135.0, "x_0p8_deg_pt": 506.4,
        "y_1e_minus5_pt": 289.32, "y_1e_minus2_pt": 61.08}
CURVES = {503: (228, [0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575]),
          677: (229, [0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.575])}


def vertical_intersections(drawing, x):
    """Intersect a vertical readout line with the filled polyline's edges."""
    intersections = []
    for kind, start, stop in drawing["items"]:
        if kind != "l":
            raise ValueError("Unexpected drawing command in fixed-source curve")
        if start.x == stop.x:
            continue
        fraction = (x-start.x)/(stop.x-start.x)
        if 0 <= fraction <= 1:
            intersections.append(start.y + fraction*(stop.y-start.y))
    if len(intersections) < 2:
        raise ValueError("Readout does not cross both edges of the observed curve")
    return min(intersections), max(intersections)


def extract(pdf: Path, out: Path):
    raw = pdf.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != SHA256:
        raise ValueError(f"PDF SHA256 differs from audited source: {digest}")
    document = pymupdf.open(pdf)
    page = document[3]
    if len(document) != 5 or len(page.get_drawings()) != 234:
        raise ValueError("Unexpected source page structure; re-audit curve selection")
    drawings = page.get_drawings()
    x_scale = (AXES["x_0p8_deg_pt"]-AXES["x_zero_pt"])/0.8
    y_scale = (AXES["y_1e_minus5_pt"]-AXES["y_1e_minus2_pt"])/3
    records = []
    for wavelength, (index, angles) in CURVES.items():
        curve = drawings[index]
        if curve["type"] != "f" or curve["fill"] != (0.0, 0.0, 0.0):
            raise ValueError("Expected the audited solid black filled curve")
        for angle in angles:
            x = AXES["x_zero_pt"] + angle*x_scale
            y_min, y_max = vertical_intersections(curve, x)
            y = (y_min+y_max)/2
            log_transmission = -5 + (AXES["y_1e_minus5_pt"]-y)/y_scale
            records.append({
                "wavelength_nm": wavelength, "shadow_angle_deg": angle,
                "observed_relative_brightness": 10**log_transmission,
                "graph_readout_angle_halfwidth_deg": 1/x_scale,
                "graph_readout_log10_halfwidth": ((y_max-y_min)/2+1)/y_scale,
                "pdf_center_x_pt": x, "pdf_center_y_pt": y,
                "pdf_stroke_y_min_pt": y_min, "pdf_stroke_y_max_pt": y_max,
            })
    out.mkdir(parents=True, exist_ok=True)
    with (out/"points.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    manifest = {
        "paper": "Ugolnikov, Maslov & Korotkiy (2011), arXiv:1106.6178v2",
        "url": URL, "source_pdf_sha256": digest, "source_pdf_bytes": len(raw),
        "figure": 2, "page_one_based": 4, "pymupdf_version": pymupdf.VersionBind,
        "axes": AXES,
        "observed_solid_path_indices_zero_based": {str(k): v[0] for k,v in CURVES.items()},
        "source_status": "Digitized published observed curves, NOT original exposures or a released numerical table",
        "measurement": "Same lunar surface element's in-eclipse brightness divided by fully out-of-eclipse brightness",
        "direction": "Western umbra, 30 degrees north of equator, as labeled in Fig. 2",
        "readout_definition": "Half PDF stroke thickness plus 1 PDF point vertically; 1 point horizontally. "
                              "Graph readout allowance only, NOT observational uncertainty or a confidence interval.",
        "missing_information": ["Filter response curves and FWHM", "Pointwise measurement uncertainties",
                                "Original CCD data and reduction products", "Exact time attached to each curve coordinate"],
        "excluded": "867 nm is outside the present 380–780 nm baseline range",
        "number_of_curve_readouts": len(records),
        "independence": "Readouts sample smoothed published curves and are not statistically independent exposures",
    }
    (out/"source.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"Extracted {len(records)} curve readouts from audited PDF to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=Path("outputs/sources/ugolnikov2011.pdf"))
    parser.add_argument("--out", type=Path, default=Path("validation/observations/ugolnikov2011"))
    args = parser.parse_args()
    if not args.pdf.exists():
        args.pdf.parent.mkdir(parents=True, exist_ok=True)
        with urlopen(URL, timeout=60) as response:
            args.pdf.write_bytes(response.read())
    extract(args.pdf, args.out)


if __name__ == "__main__":
    main()
