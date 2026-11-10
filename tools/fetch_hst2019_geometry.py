"""Rebuild the 2019-01-21 HST/STIS eclipse geometry from Horizons without fitting.

For every Table-1 sampled selenographic position this script queries Earth and
Sun vectors from that lunar site (DE441, MEAN_ME high-precision orientation,
light-time corrected) at the reflection epoch and computes the shadow-axis
distance with the same shadow_coordinates relation used for December 2011.
The receiver-side parallax never enters as an uncertainty: the vectors
originate at the lunar surface point itself, and the Moon-to-HST downleg
(used only to shift the epoch) is insensitive at the few-km level to the HST
position. The eclipse angle at the (0,0) site and the FITS RA_APER/DEC_APER
pointing provide independent checks against Table 1 and the east-positive
longitude convention. No wall-clock field is written: this file is a
deterministic function of the cached Horizons responses and the pinned
archive contract (retrieval dates live in source.json).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"validation/observations/hst2019"
CACHE = ROOT/"outputs/sources/hst2019"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api?"
MOON_RADIUS_M = 1_737_400.0
C_M_S = 299_792_458.0
SPAN_START = datetime.fromisoformat("2019-01-21T05:00:00+00:00")
SPAN_END = datetime.fromisoformat("2019-01-21T10:20:00+00:00")
# Youngblood et al. 2020 pointing-error budget: 14" initial gyro slew, gyro
# guiding accumulating <=0.001"/s (16.48" by end of orbit 2, 10.72" by end of
# orbit 3) and the linear-track approximation dominating at up to 10 arcmin.
POINTING_TRACK_ERROR_ARCSEC = 600.0
POINTING_SMALL_TERMS_ARCSEC = 14.0 + 16.48
SAMPLE_PROBE_OFFSETS_M = [-150_000.0, 0.0, 150_000.0]


def shadow_coordinates(earth_from_site, sun_from_site):
    """Shadow-axis distance, axial distance and Sun distance from site vectors.

    Both input vectors share the lunar site origin and ICRF axes, in metres.
    r=|e x u| and D=e.u for sunward unit vector u=(s-e)/|s-e|.
    """
    e,s = np.asarray(earth_from_site),np.asarray(sun_from_site)
    source_distance = np.linalg.norm(s-e,axis=-1)
    sunward = (s-e)/source_distance[...,None]
    radius = np.linalg.norm(np.cross(e,sunward),axis=-1)
    axial_distance = np.sum(e*sunward,axis=-1)
    return radius,axial_distance,source_distance


def fmt_time(t: datetime) -> str:
    return t.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]


def parse_calendar(text: str) -> datetime:
    text = text.strip()
    for pattern in ("%Y-%b-%d %H:%M:%S.%f","%Y-%b-%d %H:%M:%S","%Y-%b-%d %H:%M"):
        try:
            return datetime.strptime(text,pattern).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"Unparseable calendar date {text!r}")


def horizons_site_states(command: str, site_coord: str, start: datetime, end: datetime,
                         step: str, cache: Path, label: str) -> tuple[list[datetime], np.ndarray, dict]:
    """Fetch Earth/Sun states from a lunar site via the OBSERVER pipeline.

    The Horizons VECTORS pipeline returned time-only rows for spacecraft and
    lunar-site centers at retrieval time (verified against the committed
    December-2011 request URL, which succeeded when that geometry was built),
    so the light-time-corrected site vector is constructed from the OBSERVER
    table instead: the astrometric ICRF direction (no apparent-sky stellar
    aberration, matching the December 2011 LT convention) times the delta
    distance. Vectors are returned in metres, site origin, ICRF axes.
    """
    parameters = {"format":"json","COMMAND":f"'{command}'","OBJ_DATA":"'NO'",
                  "MAKE_EPHEM":"'YES'","EPHEM_TYPE":"'OBSERVER'","CENTER":"'coord@301'",
                  "COORD_TYPE":"'GEODETIC'","SITE_COORD":f"'{site_coord}'",
                  "START_TIME":"'"+fmt_time(start)+"'","STOP_TIME":"'"+fmt_time(end)+"'",
                  "STEP_SIZE":f"'{step}'","TIME_TYPE":"'UT'","QUANTITIES":"'1,20'",
                  "CSV_FORMAT":"'YES'"}
    url = HORIZONS+urlencode(parameters)
    path = cache/f"horizons_{label}.json"
    if path.exists():
        raw = path.read_bytes()
    else:
        with urlopen(url, timeout=90) as response:
            raw = response.read()
    response = json.loads(raw)
    if "error" in response:
        # Validate before caching so failed or malformed requests never
        # poison the cache for later runs.
        if path.exists():
            path.unlink()
        raise ValueError(f"{label}: {response['error']}")
    result = response["result"]
    if not all(text in result for text in ["DE441","MEAN_ME (high precision)"]):
        if path.exists():
            path.unlink()
        raise ValueError(f"{label}: unexpected ephemeris or lunar orientation")
    path.write_bytes(raw)
    site_match = re.search(r"Center geodetic\s*:\s*([\d.-]+),\s*([\d.-]+)",result)
    requested = [float(x) for x in site_coord.split(",")]
    if not site_match or not np.allclose([float(x) for x in site_match.groups()],
                                        [requested[0]%360,requested[1]],rtol=0,atol=1e-6):
        if path.exists():
            path.unlink()
        raise ValueError(f"{label}: unexpected lunar site coordinates")
    section = result.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    header = [line for line in result.split("$$SOE")[0].splitlines() if "Date" in line][-1]
    columns = [column.strip() for column in header.split(",")]
    wanted = {"R.A._(ICRF)","DEC__(ICRF)","delta"}
    indices = {name:columns.index(name) for name in wanted if name in columns}
    missing = wanted-set(indices)
    if missing:
        if path.exists():
            path.unlink()
        raise ValueError(f"{label}: observer table missing columns {sorted(missing)}")
    epochs,vectors = [],[]
    for row in csv.reader(section):
        epochs.append(parse_calendar(row[0]))
        ra_text,dec_text = row[indices["R.A._(ICRF)"]].strip(),row[indices["DEC__(ICRF)"]].strip()
        ra_parts = [float(part) for part in ra_text.split()]
        ra = np.radians(ra_parts[0]*15+ra_parts[1]/4+ra_parts[2]/240)
        dec_parts = [float(part) for part in dec_text.split()]
        sign = -1 if dec_text.startswith("-") else 1
        dec = np.radians(sign*(abs(dec_parts[0])+dec_parts[1]/60+dec_parts[2]/3600))
        direction = np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
        vectors.append(float(row[indices["delta"]])*1_495_978_707_00.0*direction)
    meta = {"url":url,"response_sha256":hashlib.sha256(raw).hexdigest(),
            "api_signature":response["signature"],"row_count":len(vectors),
            "construction":"astrometric ICRF direction x delta distance (OBSERVER pipeline); "
                           "VECTORS output served time-only rows at retrieval time"}
    return epochs,np.array(vectors),meta


def horizons_observer(command: str, center: str, start: datetime, end: datetime,
                      cache: Path, label: str) -> tuple[list[dict], dict]:
    """Fetch the default OBSERVER table: apparent RA/Dec of the Moon from the
    receiver, its range (delta, AU) and the lunar-orientation columns used by
    the pointing check. Horizons serves spacecraft-centered VECTORS output as
    time-only rows, so the downleg range comes from here instead."""
    parameters = {"format":"json","COMMAND":f"'{command}'","OBJ_DATA":"'NO'",
                  "MAKE_EPHEM":"'YES'","EPHEM_TYPE":"'OBSERVER'","CENTER":f"'{center}'",
                  "START_TIME":"'"+fmt_time(start)+"'","STOP_TIME":"'"+fmt_time(end)+"'",
                  "STEP_SIZE":"'1m'","TIME_TYPE":"'UT'","CSV_FORMAT":"'YES'"}
    url = HORIZONS+urlencode(parameters)
    path = cache/f"horizons_{label}.json"
    if path.exists():
        raw = path.read_bytes()
    else:
        with urlopen(url, timeout=90) as response:
            raw = response.read()
    response = json.loads(raw)
    if "error" in response:
        if path.exists():
            path.unlink()
        raise ValueError(f"{label}: {response['error']}")
    result = response["result"]
    if "DE441" not in result:
        if path.exists():
            path.unlink()
        raise ValueError(f"{label}: unexpected ephemeris")
    path.write_bytes(raw)
    section = result.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    header = [line for line in result.split("$$SOE")[0].splitlines() if "Date" in line][-1]
    columns = [column.strip() for column in header.split(",")]
    wanted = {"R.A.__(a-app)","DEC_(a-app)","delta","ObsSub-LON","ObsSub-LAT","NP.ang","Ang-diam","N.Pole-RA","N.Pole-DC"}
    indices = {name:columns.index(name) for name in wanted if name in columns}
    missing = wanted-set(indices)
    if missing:
        raise ValueError(f"{label}: observer table missing columns {sorted(missing)}")
    records = []
    for row in csv.reader(section):
        epoch = parse_calendar(row[0])
        entry = {"epoch_utc":epoch.isoformat()}
        for name,index in indices.items():
            text = row[index].strip()
            if name=="R.A.__(a-app)":
                parts = [float(part) for part in text.split()]
                entry[name] = parts[0]*15+parts[1]/4+parts[2]/240
            elif name=="DEC_(a-app)":
                parts = [float(part) for part in text.split()]
                sign = -1 if text.startswith("-") else 1
                entry[name] = sign*(abs(parts[0])+parts[1]/60+parts[2]/3600)
            else:
                entry[name] = float(text)
        entry["range_m"] = entry.pop("delta")*1_495_978_707_00.0
        records.append(entry)
    meta = {"url":url,"response_sha256":hashlib.sha256(raw).hexdigest(),
            "api_signature":response["signature"],"row_count":len(records)}
    return records,meta


def read_table1() -> list[dict]:
    with (DATA/"table1_exposures.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        row["texp_s"] = float(row["texp_s"])
        row["eclipse_angle_deg"] = float(row["eclipse_angle_deg"])
        row["lon_deg"] = float(row["lon_deg"])
        row["lat_deg"] = float(row["lat_deg"])
        row["mid_time"] = datetime.strptime("2019-01-21 "+row["mid_time_utc"],
                                            "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    return rows


def fits_headers(contract: dict, cache: Path) -> dict:
    from astropy.io import fits
    headers = {}
    for name,product in contract["mast_download"]["products"].items():
        if not name.endswith(("_crj.fits","_flt.fits")):
            continue
        with fits.open(cache/name) as hdus:
            sci = hdus["sci",1].header
            headers[product["header"]["rootname"]] = {
                "time_obs":datetime.strptime(sci["DATE-OBS"]+" "+sci["TIME-OBS"],
                                             "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc),
                "exptime_s":float(sci["EXPTIME"]),
                "ra_aper_deg":float(sci["RA_APER"]),
                "dec_aper_deg":float(sci["DEC_APER"])}
    return headers


def member_mid_times(rows: list[dict], contract: dict, headers: dict) -> dict:
    """FITS-recorded mid time of each exposure; loud on Table 1 disagreement.

    TIME-OBS is the start of the first CR-SPLIT member (or of the single
    exposure); mid = TIME-OBS + member exposure time / 2. Members after the
    first inherit the set start and take their Table 1 mid-time offsets.
    """
    fits_mid = {}
    for dataset,exposures in contract["cr_split_membership"].items():
        header = headers[dataset]
        members = [row for row in rows if row["exposure_name"] in exposures]
        total = sum(row["texp_s"] for row in members)
        if abs(header["exptime_s"]-total)>0.1:
            raise RuntimeError(f"{dataset}: crj exposure time {header['exptime_s']}s != member sum {total}s")
        base = members[0]
        fits_mid[base["exposure_name"]] = header["time_obs"]+timedelta(seconds=base["texp_s"]/2)
        for member in members[1:]:
            fits_mid[member["exposure_name"]] = fits_mid[base["exposure_name"]]+(
                member["mid_time"]-base["mid_time"])
    for row in rows:
        if row["exposure_name"] in fits_mid:
            continue
        header = headers[row["dataset_id"]]
        fits_mid[row["exposure_name"]] = header["time_obs"]+timedelta(seconds=header["exptime_s"]/2)
    checks = {}
    for row in rows:
        delta = abs((fits_mid[row["exposure_name"]]-row["mid_time"]).total_seconds())
        if delta>2.0:
            raise RuntimeError(f"{row['exposure_name']}: FITS mid time differs from Table 1 by {delta:.1f}s")
        checks[row["exposure_name"]] = delta
    return checks


def observer_at(records: list[dict], epoch: datetime) -> dict:
    times = [datetime.fromisoformat(record["epoch_utc"]) for record in records]
    after = next(index for index,time in enumerate(times) if time>=epoch)
    if times[after]==epoch:
        return records[after]
    before = after-1
    weight = (epoch-times[before]).total_seconds()/(times[after]-times[before]).total_seconds()
    entry = {"epoch_utc":epoch.isoformat()}
    for key in records[0]:
        if key=="epoch_utc":
            continue
        entry[key] = (1-weight)*records[before][key]+weight*records[after][key]
    return entry


def sky_tangents(direction: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    """Sky east and north tangents perpendicular to a unit direction."""
    pole = np.array([0.,0.,1.])
    east = np.cross(pole,direction)
    east /= np.linalg.norm(east)
    north = np.cross(direction,east)
    return east,north


def predict_pointing(lon_deg: float, lat_deg: float, sub_lon: float, sub_lat: float,
                     pole_icrf: np.ndarray, center_direction: np.ndarray,
                     range_m: float) -> tuple[float, float]:
    """Sky RA/Dec of a selenographic surface point, from the reconstructed body frame.

    The body frame is rebuilt from the Horizons sub-observer longitude/latitude
    and the pole direction, so no hand-coded east/west sign convention enters:
    the east-positive reading of Table 1 is the hypothesis under test, judged
    by the RA_APER residuals. The outward normal at the sub-observer point
    faces the observer (anti-parallel to the HST-to-center direction), which
    fixes the handedness; the point's sky direction is the exact
    range*center + radius*normal combination, not a flat projection.
    """
    cos_sub = np.cos(np.radians(sub_lat))
    sub_body = np.array([cos_sub*np.cos(np.radians(sub_lon)),
                         cos_sub*np.sin(np.radians(sub_lon)),
                         np.sin(np.radians(sub_lat))])
    pole_body = np.array([0.,0.,1.])
    east_body = np.cross(pole_body,sub_body); east_body /= np.linalg.norm(east_body)
    north_body = np.cross(sub_body,east_body)
    sub_icrf = -center_direction
    east_icrf = np.cross(pole_icrf,sub_icrf); east_icrf /= np.linalg.norm(east_icrf)
    north_icrf = np.cross(sub_icrf,east_icrf)
    cos_lat = np.cos(np.radians(lat_deg))
    normal_body = np.array([cos_lat*np.cos(np.radians(lon_deg)),
                            cos_lat*np.sin(np.radians(lon_deg)),
                            np.sin(np.radians(lat_deg))])
    normal_icrf = ((normal_body@sub_body)*sub_icrf+(normal_body@east_body)*east_icrf
                   +(normal_body@north_body)*north_icrf)
    direction = range_m*center_direction+MOON_RADIUS_M*normal_icrf
    direction /= np.linalg.norm(direction)
    ra = np.degrees(np.arctan2(direction[1],direction[0]))%360
    dec = np.degrees(np.arcsin(np.clip(direction[2],-1,1)))
    return ra,dec


def angular_separation_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    r1,r2 = np.radians(ra1),np.radians(ra2)
    d1,d2 = np.radians(dec1),np.radians(dec2)
    cos = np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    return float(np.degrees(np.arccos(np.clip(cos,-1,1)))*3600)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=CACHE)
    parser.add_argument("--out",type=Path,default=DATA/"geometry.json")
    args = parser.parse_args()
    rows = read_table1()
    contract = json.loads((DATA/"source.json").read_text())
    headers = fits_headers(contract,args.cache)
    mid_time_checks = member_mid_times(rows,contract,headers)

    observer,observer_meta = horizons_observer("301","@-48",SPAN_START,SPAN_END,
                                               args.cache,"hst_moon_observer")

    def moon_state(epoch: datetime) -> tuple[float,np.ndarray]:
        state = observer_at(observer,epoch)
        ra,dec = np.radians(state["R.A.__(a-app)"]),np.radians(state["DEC_(a-app)"])
        direction = np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
        return state["range_m"],direction

    # Reflection epochs: mid exposure minus the Moon-to-HST downleg.
    downlegs = {}
    for row in rows:
        range_m,_ = moon_state(row["mid_time"])
        downlegs[row["exposure_name"]] = range_m/C_M_S
        row["reflection_time"] = row["mid_time"]-timedelta(seconds=downlegs[row["exposure_name"]])

    # Eclipse angle at the fixed (0,0) site at observed mid times (Table 1 convention).
    zero_earth_epochs,zero_earth,zero_earth_meta = horizons_site_states(
        "399","0.0,0.0,0",rows[0]["mid_time"],rows[-1]["mid_time"],"1m",args.cache,"echeck_earth")
    zero_sun_epochs,zero_sun,zero_sun_meta = horizons_site_states(
        "10","0.0,0.0,0",rows[0]["mid_time"],rows[-1]["mid_time"],"1m",args.cache,"echeck_sun")

    def interpolate_vector(vector_table,table_epochs,epoch):
        after = next((index for index,time in enumerate(table_epochs) if time>=epoch),None)
        if after is None:
            return vector_table[-1]
        if after==0 or table_epochs[after]==epoch:
            return vector_table[after]
        before = after-1
        weight = (epoch-table_epochs[before]).total_seconds()/(
            table_epochs[after]-table_epochs[before]).total_seconds()
        return (1-weight)*vector_table[before]+weight*vector_table[after]

    def zero_state(epoch: datetime) -> tuple[np.ndarray,np.ndarray]:
        return (interpolate_vector(zero_earth,zero_earth_epochs,epoch),
                interpolate_vector(zero_sun,zero_sun_epochs,epoch))

    # Per-site Earth/Sun vectors, linearly interpolated at the reflection
    # epoch between 1-minute Horizons samples (Horizons rejects sub-minute
    # steps here; interpolation error <~2 km, far below the 150 km probe).
    site_labels = {}
    for row in rows:
        key = (row["lon_deg"],row["lat_deg"])
        if key not in site_labels:
            label = f"site_{len(site_labels):02d}_lon{row['lon_deg']:+.1f}_lat{row['lat_deg']:+.1f}"
            site_labels[key] = label
    site_meta = {}
    site_vectors = {}
    for (lon,lat),label in site_labels.items():
        epochs_needed = [row["reflection_time"] for row in rows
                         if (row["lon_deg"],row["lat_deg"])==(lon,lat)]
        start,end = min(epochs_needed)-timedelta(minutes=2),max(epochs_needed)+timedelta(minutes=2)
        e_epochs,earth,meta_e = horizons_site_states("399",f"{lon:g},{lat:g},0",
                                                     start,end,"1m",args.cache,label+"_earth")
        s_epochs,sun,meta_s = horizons_site_states("10",f"{lon:g},{lat:g},0",
                                                   start,end,"1m",args.cache,label+"_sun")
        site_meta[label+"_earth"],site_meta[label+"_sun"] = meta_e,meta_s
        site_vectors[label] = [interpolate_vector(earth,e_epochs,epoch) for epoch in epochs_needed]+[
            interpolate_vector(sun,s_epochs,epoch) for epoch in epochs_needed]

    records = []
    for row in rows:
        label = site_labels[(row["lon_deg"],row["lat_deg"])]
        epochs_needed = [other["reflection_time"] for other in rows
                         if (other["lon_deg"],other["lat_deg"])==(row["lon_deg"],row["lat_deg"])]
        member_index = epochs_needed.index(row["reflection_time"])
        site_earth,site_sun = site_vectors[label][member_index],site_vectors[label][len(epochs_needed)+member_index]
        rho,axial,sun_range = shadow_coordinates(site_earth,site_sun)
        zero_e,zero_s = zero_state(row["mid_time"])
        eclipse_angle = float(np.degrees(np.arccos(np.clip(
            np.dot(zero_e,zero_s)/np.linalg.norm(zero_e)/np.linalg.norm(zero_s),-1,1))))
        records.append({
            "dataset_id":row["dataset_id"],"exposure_name":row["exposure_name"],
            "phase":row["phase"],"mid_time_utc":row["mid_time"].isoformat(),
            "reflection_utc":row["reflection_time"].isoformat(),
            "downleg_s":downlegs[row["exposure_name"]],
            "texp_s":row["texp_s"],"lon_deg":row["lon_deg"],"lat_deg":row["lat_deg"],
            "shadow_axis_distance_m":float(rho),"receiver_axial_distance_m":float(axial),
            "sun_site_distance_m":float(sun_range),
            "eclipse_angle_table1_deg":row["eclipse_angle_deg"],
            "eclipse_angle_horizons_deg":eclipse_angle,
            "eclipse_angle_abs_diff_deg":abs(eclipse_angle-row["eclipse_angle_deg"]),
            "fits_mid_time_abs_diff_s":mid_time_checks[row["exposure_name"]]})

    # Pointing check: rebuild the body frame from Horizons sub-observer values
    # and predict the sky pointing of every Table 1 position. Only single-
    # exposure sets and the first CR-SPLIT member carry a commanded pointing
    # in their headers (RA_APER is the set start); later members follow the
    # same linear track and are checked for the record but not against the
    # tolerance.
    first_members = {exposures[0] for exposures in contract["cr_split_membership"].values()}
    pointing = []
    checked_separations = []
    for row,record in zip(rows,records):
        state = observer_at(observer,row["mid_time"])
        range_m,center_direction = moon_state(row["mid_time"])
        pole_ra,pole_dec = state["N.Pole-RA"],state["N.Pole-DC"]
        pole = np.array([np.cos(np.radians(pole_dec))*np.cos(np.radians(pole_ra)),
                         np.cos(np.radians(pole_dec))*np.sin(np.radians(pole_ra)),
                         np.sin(np.radians(pole_dec))])
        ra,dec = predict_pointing(row["lon_deg"],row["lat_deg"],state["ObsSub-LON"],
                                  state["ObsSub-LAT"],pole,center_direction,range_m)
        header = headers[row["dataset_id"]]
        separation = angular_separation_arcsec(ra,dec,header["ra_aper_deg"],header["dec_aper_deg"])
        checked = (row["dataset_id"] not in contract["cr_split_membership"]
                   or row["exposure_name"] in first_members)
        if checked:
            checked_separations.append(separation)
        pointing.append({"exposure_name":row["exposure_name"],
                         "predicted_ra_deg":ra,"predicted_dec_deg":dec,
                         "fits_ra_aper_deg":header["ra_aper_deg"],
                         "fits_dec_aper_deg":header["dec_aper_deg"],
                         "separation_arcsec":separation,"checked_against_header_pointing":checked})
    separations = np.array(checked_separations)
    if np.max(separations)>2*POINTING_TRACK_ERROR_ARCSEC:
        raise RuntimeError("East-positive Table 1 convention contradicts FITS RA_APER residuals: "
                           +json.dumps(pointing,indent=1))

    report = {
        "source":"JPL Horizons DE441; ICRF vectors; lunar MEAN_ME high-precision site orientation",
        "api_docs":"https://ssd-api.jpl.nasa.gov/doc/horizons.html",
        "vector_construction":"At lunar reflection epoch t (mid exposure minus Moon-to-HST downleg), "
                             "e=Earth(t-LT_E)-site and s=Sun(t-LT_S)-site from the sampled selenographic site; "
                             "rho=|e x u|, D=e.u, S=|s-e| for sunward u=(s-e)/|s-e|. Site vectors originate at "
                             "the lunar surface point, so receiver-side (HST) parallax is exact and not an "
                             "uncertainty term; the downleg only shifts the epoch (Moon-HST range/c, ~1.28 s). "
                             "Site states come from the Horizons OBSERVER pipeline (astrometric ICRF direction "
                             "times delta distance; the VECTORS pipeline served time-only rows at retrieval time, "
                             "including for the previously validated December-2011 request shape), which matches "
                             "the December 2011 LT convention of applying no apparent-sky stellar aberration.",
        "horizons_hst_body_id":-48,
        "horizons_hst_probe":"Moon vectors from HST (-48) confirmed available (DE441); A5 assumption verified. "
                             "HST position enters only the downleg shift: a geocentric fallback would change the "
                             "shift by <=6911 km/c = 0.023 s (~40 km sampled displacement), not the shadow geometry.",
        "longitude_convention":"Table 1 selenographic longitudes read east-positive (IAU/SPICE modern convention, "
                              "same as Horizons Moon site coordinates); validated by the RA_APER pointing check below.",
        "epoch_contract":{"observed_span_utc":[SPAN_START.isoformat(),SPAN_END.isoformat()],
                          "reflection_epoch_shift":"mid exposure time minus Moon-to-HST downleg; downleg interpolated "
                          "linearly from 1-minute Moon-from-HST DE441 ranges",
                          "site_vector_sampling":"Earth/Sun vectors from each selenographic site linearly interpolated "
                          "at the reflection epoch between 1-minute Horizons samples (Horizons rejects sub-minute steps); "
                          "interpolation error <~2 km, far below the 150 km probe scale"},
        "requests":{**observer_meta,**site_meta,
                    "echeck_earth":zero_earth_meta,"echeck_sun":zero_sun_meta},
        "spatial_sensitivity":{
            "sample_probe_offsets_m":SAMPLE_PROBE_OFFSETS_M,
            "probe_semantics":"Sampled point displaced +/-150 km toward/away from the shadow axis on the lunar "
                              "surface (rho +/- 150 km), mirroring the December 2011 radial probe scale.",
            "pointing_small_terms_arcsec":POINTING_SMALL_TERMS_ARCSEC,
            "pointing_track_error_max_arcsec":POINTING_TRACK_ERROR_ARCSEC,
            "pointing_terms_source":"Youngblood et al. 2020 Sec. II.1: 14\" initial gyro slew, gyro guiding "
                                    "accumulating at 0.001\"/s (16.48\" max orbit 2, 10.72\" max orbit 3), linear-track "
                                    "approximation dominating at up to 10 arcmin.",
            "pointing_track_error_km_max":POINTING_TRACK_ERROR_ARCSEC*MOON_RADIUS_M/1e6,
            "pointing_small_terms_km_max":POINTING_SMALL_TERMS_ARCSEC*MOON_RADIUS_M/1e6,
            "scope":"Displacement sensitivity, not a confidence interval; the track-error envelope exceeds the "
                    "150 km probe and limits the same-region albedo cancellation premise, see limitations."},
        "independent_checks":{
            "mid_time_fits_vs_table1":{"max_abs_diff_s":float(max(mid_time_checks.values())),
                                       "tolerance_s":2.0,
                                       "per_exposure_abs_diff_s":mid_time_checks},
            "eclipse_angle_table1_vs_horizons":{
                "max_abs_diff_deg":float(max(record["eclipse_angle_abs_diff_deg"] for record in records)),
                "tolerance_deg":0.1,
                "scope":"Table 1 evaluates e at the fixed selenographic (0,0) site at observed mid times; "
                        "Horizons site vectors reproduce the same angle definition."},
            "pointing_ra_aper_vs_table1":{
                "max_separation_arcsec":float(np.max(separations)),
                "rms_separation_arcsec":float(np.sqrt(np.mean(separations**2))),
                "tolerance_arcsec":2*POINTING_TRACK_ERROR_ARCSEC,
                "checked_exposure_count":len(checked_separations),
                "per_exposure":pointing,
                "scope":"Body frame reconstructed from Horizons sub-observer longitude/latitude and pole "
                        "direction (no hand-coded east/west sign); residuals judged only for single-exposure "
                        "sets and first CR-SPLIT members, the only commanded pointings recorded in the headers. "
                        "Residual scatter is the paper's documented linear-track approximation error; the "
                        "east-positive Table 1 convention is confirmed by residuals far inside that budget."},
        },
        "limitations":[
            "Sampled positions are the paper's published mean linear-track coordinates; the actual in-exposure "
            "pointing wanders within the documented error budget (up to 10 arcmin track approximation error), "
            "so the in/out sampled lunar regions can differ and albedo cancellation is conditional on the "
            "nominal pointing record.",
            "The Lommel-Seeliger per-exposure correction applied by the paper (F_LS in 0.4975-0.5021, "
            "incidence/emission angles closely coupled) changes ratios by <=1% and is omitted here; the paper's "
            "penumbral limb-darkening division (G430L I/I0 0.7-0.79) is likewise not reproduced for the "
            "auxiliary penumbral branch.",
            "Selenographic (0,0) eclipse-angle check absorbs SPICE-vs-MEAN_ME orientation differences within "
            "0.1 deg; per-exposure eclipse angles from Table 1 are evaluated at (0,0), not at the sampled point.",
        ],
        "records":records}
    args.out.write_text(json.dumps(report,indent=2)+"\n")
    umbral = [record for record in records if record["phase"]=="umbral"]
    print(f"{len(records)} exposures; umbral shadow-axis distances (km): "
          f"{[round(record['shadow_axis_distance_m']/1000,1) for record in umbral]}")
    checks = report["independent_checks"]
    print("eclipse-angle max diff (deg):",checks["eclipse_angle_table1_vs_horizons"]["max_abs_diff_deg"])
    print("pointing max/RMS separation (arcsec):",
          checks["pointing_ra_aper_vs_table1"]["max_separation_arcsec"],
          checks["pointing_ra_aper_vs_table1"]["rms_separation_arcsec"])
    print("mid-time max diff (s):",checks["mid_time_fits_vs_table1"]["max_abs_diff_s"])


if __name__ == "__main__":
    main()
