"""Synthetic-FITS extraction and pointing/geometry checks for the HST 2019 comparison."""
import csv
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from tools.fetch_hst2019_geometry import member_mid_times,predict_pointing,shadow_coordinates
from validation.compare_hst2019 import (band_median,column_medians,header_wavelength_solution,
                                        time_weighted_mean,window_mask,wavelength_solution)

DATA = Path(__file__).resolve().parents[1]/"validation/observations/hst2019"


@pytest.mark.parametrize("rotation",[np.eye(3),Rotation.from_rotvec([.4,-.7,1.2]).as_matrix()])
def test_shadow_axis_distance_for_site_on_the_sun_earth_axis(rotation):
    rho,d,s = 3.7e6,3.96e8,1.47e11
    earth = rotation @ np.array([-rho,0.,d])
    sun = rotation @ np.array([-rho,0.,d+s])
    assert shadow_coordinates(earth,sun)==pytest.approx((rho,d,s),rel=1e-12)


def test_predicted_pointing_east_limb_lands_on_the_sky_west_side():
    # Mare Crisium (selenographic east) appears on the sky-west side of the
    # disk; this pins the handedness of the reconstructed body frame.
    range_m,radius = 3.84e8,1.7374e6
    ra,dec = predict_pointing(90.,0.,0.,0.,np.array([0.,0.,1.]),np.array([1.,0.,0.]),range_m)
    angular_offset = np.degrees(radius/range_m)
    assert ra==pytest.approx(360-angular_offset,abs=1e-3)
    assert dec==pytest.approx(0.,abs=1e-6)


def test_predicted_pointing_north_point_stays_on_the_sky_north_side():
    range_m,radius = 3.84e8,1.7374e6
    ra,dec = predict_pointing(0.,45.,0.,0.,np.array([0.,0.,1.]),np.array([1.,0.,0.]),range_m)
    assert dec>0.1 and dec<0.2
    assert abs(ra-360)<1e-6 or abs(ra)<1e-6


@pytest.mark.parametrize("rotation",[Rotation.from_rotvec([.9,-.4,2.1]).as_matrix()])
def test_predicted_pointing_separation_is_rotation_invariant(rotation):
    range_m,radius = 3.84e8,1.7374e6
    center = rotation @ np.array([1.,0.,0.])
    pole = rotation @ np.array([0.,0.,1.])
    for lon,lat in [(30.,0.),(0.,-20.)]:
        ra,dec = predict_pointing(lon,lat,0.,0.,pole,center,range_m)
        point = np.array([np.cos(np.radians(dec))*np.cos(np.radians(ra)),
                          np.cos(np.radians(dec))*np.sin(np.radians(ra)),np.sin(np.radians(dec))])
        separation = np.degrees(np.arccos(np.clip(point@center,-1,1)))
        surface_angle = np.radians(np.hypot(lon,lat))
        expected = np.degrees(np.arctan2(radius*np.sin(surface_angle),
                                         range_m-radius*np.cos(surface_angle)))
        assert separation==pytest.approx(expected,abs=1e-6)


def test_column_medians_median_robust_and_counts_per_second():
    sci = np.tile(np.arange(1,5)*10.0,(16,1))  # 16 rows, 4 columns
    dq = np.zeros_like(sci,dtype=np.int16)
    dq[3,1] = 16  # hot pixel with an extreme value must not move the median
    sci[3,1] = 9e4
    counts,valid_fraction = column_medians(sci,dq,2.0)
    assert counts==pytest.approx([5.,10.,15.,20.])
    assert np.all(valid_fraction==1.0)


def test_column_medians_flags_low_validity_and_dead_columns():
    sci = np.full((16,3),100.0)
    dq = np.zeros((16,3),dtype=np.int16)
    dq[:12,1] = 4096  # saturated: only 4/16 rows remain
    dq[:,2] = 4  # bad detector column: no rows remain
    counts,valid_fraction = column_medians(sci,dq,1.0)
    assert counts[0]==pytest.approx(100.)
    assert counts[1]==pytest.approx(100.)
    assert np.isnan(counts[2])
    assert valid_fraction[1]==pytest.approx(4/16)
    assert valid_fraction[2]==0.


def test_wavelength_solution_converts_and_validates():
    values = 2889.0+2.746*np.arange(1024)
    solution = wavelength_solution(values)
    assert solution[0]==pytest.approx(288.9)
    assert solution[-1]==pytest.approx(569.8158)
    with pytest.raises(ValueError):
        wavelength_solution(values[::-1])
    with pytest.raises(ValueError):
        wavelength_solution(values[:1023])


def test_header_wavelength_solution_matches_crval_at_crpix():
    solution = header_wavelength_solution(4300.0,2.746,519.263)
    columns_1based = np.arange(1,1025)
    assert np.interp(519.263,columns_1based,solution)==pytest.approx(430.0)


def test_window_mask_clamps_to_positive_finite_intersection():
    wavelengths = np.array([250.,300.,400.,569.9,600.])
    mask = window_mask(wavelengths,300.,570.)
    assert mask.tolist()==[False,True,True,True,False]
    with pytest.raises(ValueError):
        window_mask(np.array([100.,200.]),300.,570.)


def test_time_weighted_mean_sums_counts_before_dividing():
    assert time_weighted_mean([np.array([10.,20.]),np.array([30.,40.])],[1.,3.]) \
        ==pytest.approx([25.,35.])


def test_band_median_uses_only_the_band_window():
    wavelengths = np.linspace(400.,600.,201)
    values = wavelengths**2
    assert band_median(values,wavelengths,500.,2.) \
        ==pytest.approx(np.median(values[(wavelengths>=498)&(wavelengths<=502)]))


def _rows():
    base = datetime(2019,1,21,4,57,tzinfo=timezone.utc)
    rows = []
    for index,offset in enumerate([75.,375.]):
        rows.append({"dataset_id":"d1","exposure_name":f"e{index+1}","texp_s":150.,
                     "mid_time":base+timedelta(seconds=offset)})
    rows.append({"dataset_id":"d2","exposure_name":"f1","texp_s":2.,
                 "mid_time":datetime(2019,1,21,10,0,tzinfo=timezone.utc)})
    return rows


def test_member_mid_times_matches_table_and_rejects_disagreement():
    contract = {"cr_split_membership":{"d1":["e1","e2"]}}
    headers = {"d1":{"time_obs":datetime(2019,1,21,4,57,tzinfo=timezone.utc),"exptime_s":300.},
               "d2":{"time_obs":datetime(2019,1,21,9,59,59,tzinfo=timezone.utc),"exptime_s":2.}}
    checks = member_mid_times(_rows(),contract,headers)
    assert max(checks.values())<=2.0
    headers_bad = {"d1":{"time_obs":datetime(2019,1,21,4,57,tzinfo=timezone.utc),"exptime_s":300.},
                   "d2":{"time_obs":datetime(2019,1,21,10,0,10,tzinfo=timezone.utc),"exptime_s":2.}}
    with pytest.raises(RuntimeError):
        member_mid_times(_rows(),contract,headers_bad)


def test_committed_contract_internally_consistent():
    contract = json.loads((DATA/"source.json").read_text())
    products = contract["mast_download"]["products"]
    grouping = contract["grouping"]
    with (DATA/"table1_exposures.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    texp = {row["exposure_name"]:float(row["texp_s"]) for row in rows}
    for dataset in grouping["umbral_in_eclipse"]+grouping["penumbral_auxiliary"]:
        assert dataset in contract["cr_split_membership"]
        assert f"{dataset}_crj.fits" in products and f"{dataset}_sx1.fits" in products
        total = sum(texp[name] for name in contract["cr_split_membership"][dataset])
        header_texp = products[f"{dataset}_crj.fits"]["header"]["exptime_s"]
        assert abs(total-header_texp)<=0.1
    for dataset in grouping["full_moon_out_of_eclipse"]:
        assert dataset not in contract["cr_split_membership"]
        assert f"{dataset}_flt.fits" in products and f"{dataset}_x1d.fits" in products
        row = next(row for row in rows if row["dataset_id"]==dataset)
        header_texp = products[f"{dataset}_flt.fits"]["header"]["exptime_s"]
        assert abs(float(row["texp_s"])-header_texp)<=0.1
    for name,product in products.items():
        if name.endswith(("_crj.fits","_flt.fits","_sx1.fits","_x1d.fits")):
            assert product["header"]["grating"]=="G430L"
            assert product["header"]["proposal_id"]==15674
