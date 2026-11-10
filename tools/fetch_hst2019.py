"""Download the HST/STIS 2019-01-21 eclipse products used by validation/compare_hst2019.py.

Every product is fetched once into a cache directory and pinned by SHA256; a
re-run verifies the cached bytes against the recorded contract instead of
re-downloading. No MAST metadata API is used: both the v0.1 and v0.3 invoke
endpoints returned HTTP 404 at retrieval time, so provenance is anchored to
the product-byte SHA256 plus the FITS headers' observation record. Product
filenames containing hexadecimal dataset suffixes are served by MAST fully
lowercased (verified by probe, see source.json).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"validation/observations/hst2019"
DOWNLOAD_BASE = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/{name}"
PAPER_URL = "https://arxiv.org/pdf/2008.01837v1"
PAPER_NAME = "youngblood2020_arxiv2008.01837v1.pdf"

UMBRAL_DATASETS = ["ody902030", "ody902040"]
PENUMBRAL_DATASETS = ["ody9a2010", "ody9a2020"]
FULLMOON_DATASETS = ["ody904080", "ody904090", "ody9040a0", "ody9040b0", "ody9040c0",
                     "ody9040d0", "ody9040e0", "ody9040f0", "ody9040g0", "ody9040h0",
                     "ody9040i0", "ody9040j0"]
# Table 1 (Youngblood et al. 2020) CR-SPLIT membership per observation set.
CR_SPLIT_MEMBERS = {"ody902030": ["ody902alq", "ody902amq", "ody902aoq"],
                    "ody902040": ["ody902apq", "ody902aqq"],
                    "ody9a2010": ["ody9a2atq", "ody9a2avq", "ody9a2axq"],
                    "ody9a2020": ["ody9a2ayq", "ody9a2azq", "ody9a2b0q"]}

# Products consumed by the extraction pipeline: CR-rejected combined 2D
# spectra + wavelength arrays for the CR-SPLIT sets, single-exposure flt +
# x1d for the full-moon G430L reference.
MANIFEST = []
for dataset in UMBRAL_DATASETS + PENUMBRAL_DATASETS:
    MANIFEST += [(f"{dataset}_crj.fits", "2D science spectrum, CR-rejected combination"),
                 (f"{dataset}_sx1.fits", "1D product, used only for its wavelength solution")]
for dataset in FULLMOON_DATASETS:
    MANIFEST += [(f"{dataset}_flt.fits", "2D science spectrum, single exposure"),
                 (f"{dataset}_x1d.fits", "1D product, used only for its wavelength solution")]


def fetch(name: str, cache: Path) -> bytes:
    path = cache/name
    if path.exists():
        return path.read_bytes()
    url = DOWNLOAD_BASE.format(name=name) if name.endswith(".fits") else PAPER_URL
    with urlopen(url, timeout=180) as response:
        raw = response.read()
    if len(raw) < 1000:
        raise RuntimeError(f"Download of {name} returned only {len(raw)} bytes")
    # Validate before caching (same discipline as the Horizons fetch paths):
    # an HTTP-200 error page must not poison the cache with non-FITS bytes.
    if name.endswith(".fits") and not raw[:80].lstrip().startswith(b"SIMPLE  ="):
        raise RuntimeError(f"Download of {name} is not a FITS file (missing SIMPLE header)")
    cache.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    return raw


def header_record(path: Path) -> dict:
    with fits.open(path) as hdus:
        primary = hdus[0].header
        sci = hdus["sci", 1].header
        record = {"rootname": primary.get("ROOTNAME"), "proposal_id": primary.get("PROPOSID"),
                  "target": primary.get("TARGNAME"), "instrument": primary.get("INSTRUME"),
                  "grating": primary.get("OPT_ELEM"), "aperture": primary.get("APERTURE"),
                  "cr_split": primary.get("CR-SPLIT"),
                  "exptime_s": sci.get("EXPTIME"), "date_obs": sci.get("DATE-OBS"),
                  "time_obs": sci.get("TIME-OBS"), "ra_aper_deg": sci.get("RA_APER"),
                  "dec_aper_deg": sci.get("DEC_APER"), "pa_aper_deg": sci.get("PA_APER"),
                  "orientat_deg": sci.get("ORIENTAT"), "bunit": sci.get("BUNIT"),
                  "crval1_angstrom": sci.get("CRVAL1"), "cd1_1_angstrom_per_pix": sci.get("CD1_1")}
        if path.name.endswith(("_sx1.fits", "_x1d.fits")):
            wavelength = hdus[1].data["WAVELENGTH"][0]
            record["wavelength_min_angstrom"] = float(wavelength.min())
            record["wavelength_max_angstrom"] = float(wavelength.max())
            record["wavelength_samples"] = int(wavelength.size)
        return record


def check_header(name: str, record: dict) -> None:
    if not name.endswith(".fits"):
        return
    if record["proposal_id"] != 15674 or record["instrument"] != "STIS":
        raise RuntimeError(f"{name}: unexpected proposal or instrument")
    if record["grating"] != "G430L" or record["aperture"] != "52X2":
        raise RuntimeError(f"{name}: unexpected grating or aperture {record['grating']}/{record['aperture']}")
    if record["exptime_s"] is None or record["exptime_s"] <= 0:
        raise RuntimeError(f"{name}: missing exposure time")


def build_contract(cache: Path, existing: dict | None) -> dict:
    products = {}
    for name, role in MANIFEST:
        raw = fetch(name, cache)
        path = cache/name
        digest = hashlib.sha256(raw).hexdigest()
        record = header_record(path)
        check_header(name, record)
        if existing is not None and name in existing["mast_download"]["products"]:
            pinned = existing["mast_download"]["products"][name]
            if pinned["sha256"] != digest or pinned["bytes"] != len(raw):
                raise RuntimeError(f"{name}: cached bytes differ from the pinned contract")
            products[name] = pinned
            continue
        products[name] = {"role": role, "sha256": digest, "bytes": len(raw), "header": record}
    paper_raw = fetch(PAPER_NAME, cache)
    paper_sha = hashlib.sha256(paper_raw).hexdigest()
    if existing is not None and existing["arxiv_pdf_sha256"] != paper_sha:
        raise RuntimeError("Paper PDF bytes differ from the pinned contract")
    contract = {
        "source_status": "Imported STIS archive spectra with reproducible extraction; fixed-model comparison produced by validation/compare_hst2019.py",
        "retrieved_date": existing["retrieved_date"] if existing else "2026-11-10",
        "citation": "Youngblood, A., Arney, G. N., Garcia Munoz, A., Stocke, J. T., & France, K. (2020), "
                    "'The Hubble Space Telescope's near-UV and optical transmission spectrum of Earth as an exoplanet', arXiv:2008.01837v1",
        "arxiv_pdf_url": PAPER_URL,
        "arxiv_pdf_sha256": paper_sha,
        "observatory": "HST/STIS CCD, G430L (2889-5696 A on chip), 52X2 slit, HST-GO/DD-15674, above-atmosphere receiver (no telluric airmass)",
        "mast_download": {
            "base_url": DOWNLOAD_BASE,
            "note": "MAST metadata API (v0.1 and v0.3 invoke) returned HTTP 404 at retrieval time; provenance is anchored to "
                    "product-byte SHA256 and the FITS header observation record. Product filenames containing hexadecimal "
                    "dataset suffixes are served fully lowercased (e.g. ody9040a0_flt.fits).",
            "products": products,
        },
        "cr_split_membership": CR_SPLIT_MEMBERS,
        "grouping": {
            "umbral_in_eclipse": UMBRAL_DATASETS,
            "penumbral_auxiliary": PENUMBRAL_DATASETS,
            "full_moon_out_of_eclipse": FULLMOON_DATASETS,
            "ruling": "Plan default A2 confirmed against the archive: umbral (5 G430L exposures in 2 CR-SPLIT sets) and "
                      "full-moon (12 single G430L exposures) both exist, so the main ratio is umbral/full-moon; penumbral "
                      "sets are auxiliary and are not part of the headline comparison.",
        },
        "product_availability_survey": {
            "recorded_utc_date": "2026-11-10",
            "confirmed_entry_points": ["ody902030_flt.fits", "ody9a2010_flt.fits", "ody904080_flt.fits"],
            "pipeline_1d_products": "_sx1 present for the four CR-SPLIT sets; _x1d present for the twelve single-exposure full-moon sets; _x2d present only for ody904080",
            "literature_extracted_values": "No HLSP or machine-readable transmission-spectra table was found for this program; "
                                           "the published spectra exist only as paper figures, so the extraction replicates the paper's own method on the raw archive",
            "point_source_products_unused": "The pipeline _sx1/_x1d flux columns use point-source extraction with background regions placed on the "
                                            "filled-slit lunar spectrum; they are used only for wavelength calibration, never for flux.",
        },
        "what_this_source_cannot_support": [
            "No red-end coverage beyond 570 nm: G430L ends at 5700 A, so the ~600 nm Chappuis-band center and any turquoise-band color cannot be validated",
            "Model tables start at 300 nm; the observed 289-300 nm edge is excluded and never extrapolated",
            "Full visible color: the G230LB spectra (1685-3060 A) are outside the model wavelength range and are not extracted",
            "Pointing drift within and between exposures (linear tracks, gyro guiding) makes the in/out sampled lunar regions differ; "
            "lunar albedo cancellation is conditional on the nominal pointing record, see geometry.json",
        ],
        "extraction_products": {
            "per_exposure_spectra": "spectra.csv",
            "comparison": "validation/results/hst2019/comparison.json",
            "geometry": "geometry.json",
            "documentation": "docs/MULTI_EVENT_OBSERVATIONS.md",
        },
    }
    return contract


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=ROOT/"outputs/sources/hst2019")
    parser.add_argument("--out", type=Path, default=DATA/"source.json")
    args = parser.parse_args()
    existing = json.loads(args.out.read_text()) if args.out.exists() else None
    contract = build_contract(args.cache, existing)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(contract, indent=2)+"\n")
    total = sum(p["bytes"] for p in contract["mast_download"]["products"].values())
    print(f"{len(contract['mast_download']['products'])} FITS products + paper pinned; {total/1e6:.1f} MB cached under {args.cache}")
    print(f"contract written to {args.out}")
    for name, product in contract["mast_download"]["products"].items():
        if name.endswith(("_crj.fits", "_flt.fits")):
            header = product["header"]
            print(f"{name}: exptime={header['exptime_s']}s {header['date_obs']}T{header['time_obs']} "
                  f"ra={header['ra_aper_deg']:.5f} dec={header['dec_aper_deg']:.5f}")


if __name__ == "__main__":
    main()
