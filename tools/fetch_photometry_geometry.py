"""Fetch fixed geocentric ranges at published greatest-eclipse epochs.

These epochs approximate photometric minima; the paper gives no exposure times.
No atmospheric or photometric parameters are inferred from the observations.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

AU_M = 149_597_870_700.0
DATA = Path(__file__).resolve().parents[1]/"validation/observations/hernitschek2008"


def fetch(body, epoch, cache):
    start = datetime.fromisoformat(epoch.replace("Z", "+00:00"))
    parameters = {"format": "json", "COMMAND": f"'{body}'", "OBJ_DATA": "'NO'",
                  "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'500@399'",
                  "START_TIME": "'"+start.strftime("%Y-%m-%d %H:%M:%S")+"'",
                  "STOP_TIME": "'"+(start+timedelta(minutes=1)).strftime("%Y-%m-%d %H:%M:%S")+"'",
                  "STEP_SIZE": "'1m'", "TIME_DIGITS": "'SECONDS'",
                  "QUANTITIES": "'20'", "CSV_FORMAT": "'YES'"}
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(parameters)
    path = cache/f"horizons_{start.date()}_{body}.json"
    if path.exists():
        raw = path.read_bytes()
    else:
        with urlopen(url, timeout=45) as response:
            raw = response.read()
        cache.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    response = json.loads(raw)
    if "error" in response:
        raise ValueError(response["error"])
    result = response["result"]
    if "GEOCENTRIC" not in result or "DE441" not in result:
        raise ValueError("Unexpected center or ephemeris")
    rows = list(csv.reader(result.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()))
    expected = start.strftime("%Y-%b-%d %H:%M:%S")
    if rows[0][0].strip() != expected:
        raise ValueError("Unexpected epoch")
    return {"target": body, "url": url, "parameters": parameters,
            "response_sha256": hashlib.sha256(raw).hexdigest(),
            "api_signature": response["signature"], "utc": expected,
            "distance_m": float(rows[0][3])*AU_M}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("outputs/sources"))
    parser.add_argument("--out", type=Path, default=DATA/"geometry.json")
    args = parser.parse_args()
    events = []
    with (DATA/"points.csv").open() as stream:
        for observation in csv.DictReader(stream):
            epoch = observation["greatest_eclipse_utc"]
            moon = fetch("301", epoch, args.cache)
            sun = fetch("10", epoch, args.cache)
            events.append({"event_date": observation["event_date"], "adopted_time_utc": epoch,
                           "moon_distance_m": moon["distance_m"], "sun_distance_m": sun["distance_m"],
                           "moon": moon, "sun": sun})
    report = {"source": "JPL Horizons DE441, geocentric observer, quantity20 apparent range",
              "api_docs": "https://ssd-api.jpl.nasa.gov/doc/horizons.html",
              "recorded_utc": datetime.now(timezone.utc).isoformat(),
              "approximation": "Greatest-eclipse epoch substitutes for unpublished photometric-minimum epoch; "
                               "geocentric ranges approximate the axis-normal receiver plane; "
                               "no topocentric surface reconstruction.",
              "events": events}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2)+"\n")
    for event in events:
        print(event["event_date"], event["moon_distance_m"]/1000, "km", event["sun_distance_m"]/AU_M, "AU")


if __name__ == "__main__":
    main()
