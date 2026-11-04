"""Independent area moments protect aperture geometry in photometric comparisons."""
import numpy as np
import pytest

from validation.multiple_events import aperture_quadrature


@pytest.mark.parametrize("center,radius",[(0.,2.),(.5,2.),(2.,2.),(4.,2.)])
def test_offset_aperture_matches_cartesian_disk_area_and_second_moment(center,radius):
    r,area = aperture_quadrature(center,radius,256)
    # On a uniform Cartesian disk centered (c,0), <x²+y²>=c²+R²/2.
    assert area.sum()==pytest.approx(1.,abs=5e-8)
    assert area @ r**2==pytest.approx(center**2+radius**2/2,rel=5e-8)
    assert np.all(area>0)
