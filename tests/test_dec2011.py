"""Geometry and exposure integration against independently known analytic cases."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from tools.fetch_dec2011_geometry import shadow_coordinates
from validation.compare_dec2011 import exposure_average


@pytest.mark.parametrize("rotation",[np.eye(3),Rotation.from_rotvec([.4,-.7,1.2]).as_matrix()])
def test_site_distance_to_known_earth_sun_axis(rotation):
    # Earth at0, Sun at(0,0,S), lunar site at(rho,0,-D), with arbitrary ICRF axes.
    rho,d,s = 3.7e6,3.96e8,1.47e11
    earth = rotation @ np.array([-rho,0.,d])
    sun = rotation @ np.array([-rho,0.,d+s])
    recovered = shadow_coordinates(earth,sun)
    assert recovered==pytest.approx((rho,d,s),rel=1e-12)


def test_exposure_average_is_linear_irradiance_integral_not_midpoint_or_log_mean():
    seconds = np.linspace(0,1800,31)
    x = seconds/1800
    # Known continuous averages, deliberately varying chromaticity with time.
    values = np.column_stack([1+x*x,2+4*x**3,3-2*x+x*x])
    assert exposure_average(values,seconds)==pytest.approx([4/3,3,7/3],rel=1e-13)
    assert not np.allclose(exposure_average(values,seconds),values[15])
