"""Recover near-Tycho eclipse geometry without fitting atmospheric parameters.

Horizons handles the DE441 orbit and high-precision MEAN_ME lunar orientation.
Light-time corrected Earth/Sun vectors are evaluated from a nominal Tycho site
at the lunar reflection epoch. Fig. 1 independently checks this nominal target;
the exact fibre coordinates remain unknown and require spatial sensitivity.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"validation/observations/yan2011dec"
OBS_START = datetime.fromisoformat("2011-12-10T14:08:25+00:00")
OBS_END = OBS_START+timedelta(seconds=1800)
BRIGHT_MID = datetime.fromisoformat("2011-12-10T18:10:27+00:00")
DOWNLEG_SECONDS = 1.32
SITE = {"longitude_deg_east":-11.22,"latitude_deg_north":-43.30,"altitude_km":0.}


def shadow_coordinates(earth_from_site, sun_from_site):
    """SI receiver radius, axial distance and Earth-Sun distance from vectors.

    Both input vectors share the lunar site origin and ICRF axes, in metres.
    With LT corrections their endpoints are the retarded Earth/Sun positions.
    r=|e x u| and D=e.u for sunward unit vector u=(s-e)/|s-e|.
    """
    e,s = np.asarray(earth_from_site),np.asarray(sun_from_site)
    source_distance = np.linalg.norm(s-e,axis=-1)
    sunward = (s-e)/source_distance[...,None]
    radius = np.linalg.norm(np.cross(e,sunward),axis=-1)
    axial_distance = np.sum(e*sunward,axis=-1)
    return radius,axial_distance,source_distance


def fetch(request, cache):
    label,body,correction,start,end = request
    fmt = lambda t:t.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
    coords = ",".join(f"{SITE[key]:g}" for key in ["longitude_deg_east","latitude_deg_north","altitude_km"])
    parameters = {"format":"json","COMMAND":f"'{body}'","OBJ_DATA":"'NO'",
                  "MAKE_EPHEM":"'YES'","EPHEM_TYPE":"'VECTORS'","CENTER":"'coord@301'",
                  "COORD_TYPE":"'GEODETIC'","SITE_COORD":f"'{coords}'",
                  "START_TIME":"'"+fmt(start)+"'","STOP_TIME":"'"+fmt(end)+"'",
                  "STEP_SIZE":"'1m'","TIME_TYPE":"'UT'","REF_PLANE":"'FRAME'",
                  "REF_SYSTEM":"'ICRF'","VEC_CORR":f"'{correction}'","VEC_TABLE":"'2'",
                  "OUT_UNITS":"'KM-S'","CSV_FORMAT":"'YES'"}
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(parameters)
    path = cache/f"dec2011_{label}.json"
    if path.exists():
        raw = path.read_bytes()
    else:
        with urlopen(url,timeout=45) as response:
            raw = response.read()
        cache.mkdir(parents=True,exist_ok=True)
        path.write_bytes(raw)
    response = json.loads(raw)
    if "error" in response:
        raise ValueError(response["error"])
    result = response["result"]
    if not all(text in result for text in ["DE441","MEAN_ME (high precision)","ICRF"]):
        raise ValueError("Unexpected ephemeris, lunar orientation, frame or site")
    site_match = re.search(r"Center geodetic\s*:\s*([\d.-]+),\s*([\d.-]+)",result)
    if not site_match or not np.allclose([float(x) for x in site_match.groups()],
                                        [SITE["longitude_deg_east"]%360,SITE["latitude_deg_north"]],rtol=0,atol=1e-6):
        raise ValueError("Unexpected lunar site coordinates")
    rows = list(csv.reader(result.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()))
    vectors = np.array([[float(x)*1000 for x in row[2:5]] for row in rows])
    expected_count = int((end-start).total_seconds()/60)+1
    if len(rows)!=expected_count:
        raise ValueError("Unexpected time sampling")
    first_epoch = datetime.strptime(rows[0][1].strip().removeprefix("A.D. "),"%Y-%b-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
    if abs((first_epoch-start).total_seconds())>.001:
        raise ValueError("Cached epoch differs from requested epoch")
    jd = np.array([float(row[0]) for row in rows])
    if not np.allclose(np.diff(jd)*86400,60.,rtol=0,atol=.0001):
        raise ValueError("Unexpected time step")
    meta = {"url":url,"parameters":parameters,"response_sha256":hashlib.sha256(raw).hexdigest(),
            "api_signature":response["signature"],"first_calendar_utc":rows[0][1].strip(),
            "last_calendar_utc":rows[-1][1].strip(),"row_count":len(rows)}
    return label,vectors,meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=ROOT/"outputs/sources")
    parser.add_argument("--out",type=Path,default=DATA/"geometry.json")
    args = parser.parse_args()
    shift = timedelta(seconds=DOWNLEG_SECONDS)
    start,end = OBS_START-shift,OBS_END-shift
    bright = BRIGHT_MID-shift
    requests = [(f"u1_{name}_{correction.lower()}",body,correction,start,end)
                for correction in ["LT","NONE"] for name,body in [("earth","399"),("sun","10")]]
    requests += [(f"b2_{name}_lt",body,"LT",bright,bright+timedelta(minutes=1))
                 for name,body in [("earth","399"),("sun","10")]]
    with ThreadPoolExecutor(max_workers=3) as pool:
        replies = list(pool.map(lambda request:fetch(request,args.cache),requests))
    vectors = {label:vec for label,vec,meta in replies}
    metadata = {label:meta for label,vec,meta in replies}
    r,d,s = shadow_coordinates(vectors["u1_earth_lt"],vectors["u1_sun_lt"])
    r_geometric,_,_ = shadow_coordinates(vectors["u1_earth_none"],vectors["u1_sun_none"])
    bright_r,bright_d,bright_s = shadow_coordinates(vectors["b2_earth_lt"][0],vectors["b2_sun_lt"][0])
    bright_sun_range = float(np.linalg.norm(vectors["b2_sun_lt"][0]))
    records = []
    for i in range(len(r)):
        records.append({"observed_utc":(OBS_START+timedelta(minutes=i)).isoformat(),
                        "reflection_utc":(start+timedelta(minutes=i)).isoformat(),
                        "seconds_since_start":60*i,"shadow_radius_m":float(r[i]),
                        "receiver_axial_distance_m":float(d[i]),"sun_earth_distance_m":float(s[i]),
                        "sun_site_distance_m":float(np.linalg.norm(vectors["u1_sun_lt"][i])),
                        "solar_flux_relative_to_b2":float((bright_sun_range/np.linalg.norm(vectors["u1_sun_lt"][i]))**2),
                        "simultaneous_geometry_radius_m":float(r_geometric[i])})
    figure = json.loads((DATA/"figure1_geometry.json").read_text())
    center = np.array(figure["shadow_center_px"])
    endpoints = np.array([figure["trace_start_px"],figure["trace_end_px"]])
    plot_radii = np.linalg.norm(endpoints-center,axis=1)/figure["adopted_moon_radius_px"]*1_737_400
    tolerance = 150_000.
    endpoint_error = np.abs(np.array([r[0],r[-1]])-plot_radii)
    if np.max(endpoint_error)>tolerance:
        raise RuntimeError("Ephemeris target disagrees with the published fibre-trace schematic")
    report = {"source":"JPL Horizons DE441; ICRF vectors; lunar MEAN_ME high-precision site orientation",
              "api_docs":"https://ssd-api.jpl.nasa.gov/doc/horizons.html",
              "manual":"https://ssd.jpl.nasa.gov/horizons/manual.html",
              "recorded_utc":datetime.now(timezone.utc).isoformat(),
              "nominal_site":{**SITE,"source":"https://planetarynames.wr.usgs.gov/Feature/6163",
                              "scope":"Tycho crater center as nominal near-Tycho fibre location, not a measured exact pointing"},
              "epoch_contract":{"observed_start_utc":OBS_START.isoformat(),"observed_end_utc":OBS_END.isoformat(),
                                "exposure_s":1800,"reflection_epoch_shift_s":-DOWNLEG_SECONDS,
                                "scope":"Approximate Moon-to-ground downleg; 396000km/c about1.32s. Site-dependent correction <=0.022s. Earth/Sun vectors use LT (not LT+S) to retain physical ray coordinates."},
              "vector_construction":"At lunar reflection epoch t, e=Earth(t-LT_E)-site(t), s=Sun(t-LT_S)-site(t). Sunward u=(s-e)/|s-e|; rho=|e cross u|, D=e dot u, S=|s-e|. No apparent-sky stellar aberration rotation applied.",
              "independent_figure_check":{"source":"Yan Fig.1, native raster on PDF page2",
                                          "plot_start_end_radius_m":plot_radii.tolist(),
                                          "ephemeris_start_end_radius_m":[float(r[0]),float(r[-1])],
                                          "absolute_endpoint_difference_m":endpoint_error.tolist(),
                                          "schematic_reading_allowance_m":tolerance,
                                          "scope":"Cross-check, not an observational pointing confidence interval; no atmosphere-derived tangent heights used"},
              "spatial_sensitivity_radial_offsets_m":[-tolerance,0.,tolerance],
              "bright_reference":{"mid_observed_utc":BRIGHT_MID.isoformat(),"reflection_utc":bright.isoformat(),
                                  "shadow_radius_m":float(bright_r),"receiver_axial_distance_m":float(bright_d),
                                  "sun_earth_distance_m":float(bright_s),"sun_site_distance_m":bright_sun_range},
              "requests":metadata,"records":records}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+"\n")
    print(f"31 samples; shadow radius {r.min()/1000:.3f}..{r.max()/1000:.3f} km")
    print("Fig.1 endpoint differences (km):",endpoint_error/1000)
    print("Simultaneous vs retarded radii max difference (km):",np.max(np.abs(r-r_geometric))/1000)


if __name__ == "__main__":
    main()
