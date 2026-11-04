"""Fetch geocentric Sun/Moon ranges for the visible June 2011 eclipse data.

Source: JPL Horizons observer quantity 20. No credentials or extra dependencies.
The returned apparent ranges approximate the model's axis-normal plane distance.
The observation curves combine a short time window rather than one exposure.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

AU_M = 149_597_870_700.0


def fetch(body, cache):
    parameters = {"format": "json", "COMMAND": f"'{body}'", "OBJ_DATA": "'NO'",
                  "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'500@399'",
                  "START_TIME": "'2011-06-15 19:20'", "STOP_TIME": "'2011-06-15 19:40'",
                  "STEP_SIZE": "'10m'", "QUANTITIES": "'20'", "CSV_FORMAT": "'YES'"}
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(parameters)
    path = cache/f"horizons_{body}.json"
    if path.exists():
        raw = path.read_bytes()
    else:
        with urlopen(url, timeout=60) as response:
            raw = response.read()
        cache.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    response = json.loads(raw)
    if "error" in response:
        raise ValueError(response["error"])
    text = response["result"]
    if "GEOCENTRIC" not in text or "DE441" not in text:
        raise ValueError("Unexpected Horizons center or ephemeris")
    table = text.split("$$SOE")[1].split("$$EOE")[0]
    rows = []
    for row in csv.reader(table.strip().splitlines()):
        rows.append({"utc": row[0].strip(), "distance_au": float(row[3]),
                     "distance_m": float(row[3])*AU_M})
    return {"target": body, "url": url, "parameters": parameters,
            "response_sha256": hashlib.sha256(raw).hexdigest(),
            "api_signature": response["signature"], "ranges": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("outputs/sources"))
    parser.add_argument("--out", type=Path, default=Path("validation/observations/ugolnikov2011/geometry.json"))
    args = parser.parse_args()
    moon = fetch("301", args.cache)
    sun = fetch("10", args.cache)
    time = "2011-Jun-15 19:30"
    result = {
        "source": "JPL Horizons DE441, geocentric observer, quantity20 apparent range",
        "api_docs": "https://ssd-api.jpl.nasa.gov/doc/horizons.html",
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "adopted_time_utc": time,
        "moon_distance_m": next(r["distance_m"] for r in moon["ranges"] if r["utc"] == time),
        "sun_distance_m": next(r["distance_m"] for r in sun["ranges"] if r["utc"] == time),
        "moon": moon, "sun": sun,
        "approximation": "One representative geocentric range for the first15minutes of totality; "
                         "no per-readout timing, topocentric or lunar-surface depth reconstruction",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(f"Saved distances at {time}: Moon {result['moon_distance_m']/1000:.3f} km, "
          f"Sun {result['sun_distance_m']/AU_M:.9f} AU")


if __name__ == "__main__":
    main()
