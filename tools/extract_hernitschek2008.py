"""Check published numerical table transcription against a pinned HTML snapshot."""
from __future__ import annotations

import argparse
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "validation/observations/hernitschek2008"


class Tables(HTMLParser):
    def __init__(self):
        super().__init__()
        self.depth = 0
        self.rows = []
        self.row = None
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.depth += 1
        if tag == "tr" and self.depth:
            self.row = []
        if tag in ("td", "th") and self.row is not None:
            self.cell = []

    def handle_data(self, value):
        if self.cell is not None:
            self.cell.append(value)

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None
        if tag == "tr" and self.row is not None:
            self.rows.append(self.row)
            self.row = None
        if tag == "table":
            self.depth -= 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, default=ROOT / "outputs/sources/hernitschek2008.html")
    args = parser.parse_args()
    metadata = json.loads((DATA / "source.json").read_text())
    if not args.html.exists():
        args.html.parent.mkdir(parents=True, exist_ok=True)
        with urlopen(metadata["url"], timeout=60) as response:
            args.html.write_bytes(response.read())
    raw = args.html.read_bytes()
    if hashlib.sha256(raw).hexdigest() != metadata["source_html_sha256"]:
        raise SystemExit("HTML snapshot hash differs; review the new publisher page before updating provenance")
    tables = Tables()
    tables.feed(raw.decode("utf-8"))
    values = {}
    for row in tables.rows:
        if len(row) in (4, 5):
            pair = row[-2:]
            if row[0] in values and values[row[0]] != pair:
                raise SystemExit(f"Conflicting duplicated table row: {row[0]}")
            values[row[0]] = pair
    columns = {
        "Full Moon luminance": ("full_luminance_cd_m2", "full_luminance_reported_error_cd_m2"),
        "Minimum luminance": ("minimum_luminance_cd_m2", "minimum_luminance_reported_error_cd_m2"),
        "Linear dynamic range": ("published_dynamic_range",),
        "Full Moon brightness": ("full_brightness_mag_vis", "full_brightness_reported_error_mag"),
        "Minimum brightness": ("minimum_brightness_mag_vis", "minimum_brightness_reported_error_mag"),
        "Magnitude drop": ("magnitude_drop", "magnitude_drop_reported_error"),
        "Earth’s penumbral radius": ("geocentric_penumbra_radius_arcmin",),
        "Earth’s umbral radius": ("geocentric_umbra_radius_arcmin",),
        "Lunar radius at mid eclipse": ("geocentric_lunar_radius_arcmin",),
        "Smallest phase angle": ("geocentric_minimum_phase_angle_arcmin",),
        "Greatest umbral magnitude": ("greatest_umbral_magnitude",),
    }
    with (DATA / "points.csv").open(newline="") as stream:
        points = list(csv.DictReader(stream))
    if len(points) != 2:
        raise SystemExit("Expected exactly two event records")
    checks = 0
    for label, keys in columns.items():
        for event_index, point in enumerate(points):
            numbers = re.findall(r"[-−]?\d+(?:\.\d+)?", values[label][event_index].replace("− ", "−"))
            if len(numbers) != len(keys):
                raise SystemExit(f"Unexpected source format for {label}")
            for key, number in zip(keys, numbers):
                if float(point[key]) != float(number.replace("−", "-")):
                    raise SystemExit(f"Transcription mismatch: {point['event_date']} {key}")
                checks += 1
    print(f"Verified {checks} photometric/geometry numbers against SHA-256-pinned Tables 1–2")


if __name__ == "__main__":
    main()
