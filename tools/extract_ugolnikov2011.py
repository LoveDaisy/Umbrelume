"""Read observed curves from Fig. 2 of the fixed arXiv 1106.6178v2 PDF.

One-time provenance tool, not a runtime dependency:
uv run --with pymupdf==2.1 python tools/extract_ugolnikov2011.py --branch C
Only published OBSERVED solid curves are extracted. The PDF's dashed model
curves are not reference measurements. Original exposures are not available here.

The ``--branch`` argument carries the H1 adjudication made from a human reading
of this fixed PDF (recorded in the task progress log): branch A = measured
filter response curves stated, branch B = FWHM/bandwidth stated without curves,
branch C = no filter width or response information stated. The guard re-checks
that claim against the extracted text so the adjudication is mechanically
reproducible; contradicting hits abort instead of silently switching branches.
Branch A/B currently abort even when their anchors are found — producing A/B
output first requires an owner adjudication of the carrier form and sampling
method. Making ``--branch`` required intentionally retires the previous
no-argument invocation; the curve-extraction path itself is unchanged.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from urllib.request import urlopen

import pymupdf


URL = "https://arxiv.org/pdf/1106.6178v2"
SHA256 = "366b6bd302937530f52d8fc3768be26ebb75657dba5072a3bd07d0733b67c043"
AXES = {"x_zero_pt": 135.0, "x_0p8_deg_pt": 506.4,
        "y_1e_minus5_pt": 289.32, "y_1e_minus2_pt": 61.08}
CURVES = {503: (228, [0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575]),
          677: (229, [0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.575])}

# Branch C trigger keywords: terms whose presence would contradict the
# adjudication that the paper states no filter width or response information.
# The paper is in English; the list is literal substrings of the normalized
# (lowercase, dehyphenated, whitespace-collapsed) full text.
FILTER_KEYWORDS = ["fwhm", "bandwidth", "band width", "bandpass", "band-pass",
                   "passband", "pass band", "half-width", "halfwidth", "half width",
                   "central wavelength", "center wavelength", "centre wavelength",
                   "wavelength shift", "wavelength offset", "interference filter",
                   "response curve", "filter response", "spectral response",
                   "transmission curve", "transmission profile", "spectral transmission",
                   "filter transmission", "spectral bandpass"]
BRANCH_A_ANCHORS = ["response curve", "filter response", "spectral response",
                    "transmission curve", "spectral transmission", "filter transmission"]
BRANCH_B_ANCHORS = ["fwhm", "bandwidth", "band width", "bandpass", "band-pass",
                    "passband", "pass band", "half-width", "halfwidth", "half width"]
# Human-adjudicated on 2026-11-10 from the fixed PDF (progress log DECISION
# entry); fixed constant so regeneration stays deterministic.
ADJUDICATION_DATE = "2026-11-10"
FILTER_STATEMENTS = [
    {"page_one_based": 1, "section": "abstract",
     "quote": "three narrow spectral bands with effective wavelengths equal to 503, 677, and 867 nm"},
    {"page_one_based": 2, "section": "2. Observations",
     "quote": "SAO observations were conducted by two CCD-cameras with narrow filters with the "
              "effective wavelengths 503 and 677 nm."},
    {"page_one_based": 2, "section": "2. Observations",
     "quote": "Data obtained at 677 and 867 nm were also used to correct the effects of weak "
              "“red tail” of 503 nm filter."},
]


def normalized_text(document) -> str:
    """Lowercased, dehyphenated, whitespace-collapsed text of every page."""
    raw = "\n".join(page.get_text("text") for page in document)
    text = raw.lower().replace("-\n", "").replace("- ", "-")
    return re.sub(r"\s+", " ", text)


def keyword_hits(text: str, keywords: list[str]) -> dict[str, int]:
    """Hyphen-insensitive literal counts.

    Each keyword is counted on the normalized text and on a hyphen-stripped
    copy, so a natively hyphenated spelling in the body ("band-width",
    "pass-band") is matched by the space-form keyword; the recorded count is
    the larger of the two forms.
    """
    dehyphenated = text.replace("-", "")
    return {keyword: max(text.count(keyword),
                         dehyphenated.count(keyword.replace("-", "")))
            for keyword in keywords}


def filter_information(text: str, branch: str, digest: str) -> dict:
    """Branch-aware filter-information record; contradicted claims abort."""
    if branch == "C":
        hits = keyword_hits(text, FILTER_KEYWORDS)
        found = {keyword: count for keyword, count in hits.items() if count}
        if found:
            raise ValueError(
                "Branch C contradicted: filter-parameter keywords present in the "
                f"fixed PDF: {found}; re-adjudicate H1 before regenerating")
        return {
            "branch": "C",
            "status": "not stated in paper",
            "adjudicated_date": ADJUDICATION_DATE,
            "method": "pymupdf==2.1 full-text extraction of the fixed-SHA PDF, case-insensitive "
                      "literal keyword search over body text, figure captions and table text; "
                      "figure pages 3-5 additionally inspected visually",
            "search": {
                "keywords": FILTER_KEYWORDS,
                "scope": "full extracted text of all 5 pages (body, figure captions, table text); "
                         "figure raster content covered by manual page inspection recorded in the "
                         "task progress log",
                "normalization": "lowercase; hyphenated line breaks joined; whitespace collapsed; "
                                 "keywords additionally matched against a hyphen-stripped copy so "
                                 "natively hyphenated spelling variants (e.g. band-width) are "
                                 "covered; recorded count is the max of both forms",
                "hit_counts": hits,
            },
            "filter_statements_present": FILTER_STATEMENTS,
            "approximation_grade": "Assumed rectangular and Gaussian passbands, FWHM 10/20/40 nm "
                                   "probe envelope; filter width and shape not measured and not "
                                   "stated in the paper",
            "source_pdf_sha256": digest,
            "upgrade_note": "Requesting the measured filter response curves from the authors is "
                            "recorded as an upgrade path only; no external contact was made",
        }
    if branch == "A":
        found = {keyword: count for keyword, count
                 in keyword_hits(text, BRANCH_A_ANCHORS).items() if count}
        if not found:
            raise ValueError("Branch A claimed but no filter response-curve statement found in "
                             "the fixed PDF; re-adjudicate H1 before regenerating")
        raise ValueError("Branch A response-curve statement present; data sampling requires the "
                         "owner-adjudicated carrier form and sampling method (plan Step 1/2), "
                         f"anchors found: {found}")
    if branch == "B":
        found = {keyword: count for keyword, count
                 in keyword_hits(text, BRANCH_B_ANCHORS).items() if count}
        if not found:
            raise ValueError("Branch B claimed but no FWHM/bandwidth statement found in the fixed "
                             "PDF; re-adjudicate H1 before regenerating")
        raise ValueError("Branch B FWHM/bandwidth statement present; extraction of stated values "
                         "requires owner adjudication (plan Step 2 branch B), anchors found: "
                         f"{found}")
    raise ValueError(f"Unknown branch {branch!r}; expected one of A, B, C")


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


def extract(pdf: Path, out: Path, branch: str):
    raw = pdf.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != SHA256:
        raise ValueError(f"PDF SHA256 differs from audited source: {digest}")
    document = pymupdf.open(pdf)
    page = document[3]
    if len(document) != 5 or len(page.get_drawings()) != 234:
        raise ValueError("Unexpected source page structure; re-audit curve selection")
    info = filter_information(normalized_text(document), branch, digest)
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
        "filter_information": info,
        "missing_information": list(MISSING_INFORMATION[branch]),
        "excluded": "867 nm is outside the present 380–780 nm baseline range",
        "number_of_curve_readouts": len(records),
        "independence": "Readouts sample smoothed published curves and are not statistically independent exposures",
    }
    (out/"source.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"Extracted {len(records)} curve readouts from audited PDF to {out} (branch {branch})")


MISSING_INFORMATION = {
    "C": ["Filter response curves and FWHM", "Pointwise measurement uncertainties",
          "Original CCD data and reduction products", "Exact time attached to each curve coordinate"],
    "B": ["Filter response curves (FWHM is recorded in filter_information)",
          "Pointwise measurement uncertainties",
          "Original CCD data and reduction products", "Exact time attached to each curve coordinate"],
    "A": ["Pointwise measurement uncertainties",
          "Original CCD data and reduction products", "Exact time attached to each curve coordinate"],
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=Path("outputs/sources/ugolnikov2011.pdf"))
    parser.add_argument("--out", type=Path, default=Path("validation/observations/ugolnikov2011"))
    parser.add_argument("--branch", required=True, choices=["A", "B", "C"],
                        help="H1 adjudication from the human reading of the fixed PDF")
    args = parser.parse_args()
    if not args.pdf.exists():
        args.pdf.parent.mkdir(parents=True, exist_ok=True)
        with urlopen(URL, timeout=60) as response:
            args.pdf.write_bytes(response.read())
    extract(args.pdf, args.out, args.branch)


if __name__ == "__main__":
    main()
