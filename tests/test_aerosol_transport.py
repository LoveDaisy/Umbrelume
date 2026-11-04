"""Independent aerosol transport references and physical invariants."""
from dataclasses import replace

import numpy as np

from eclipse.aerosol import AerosolProfile
from eclipse.atmosphere import Atmosphere
from eclipse.rays import build_ray_table
from eclipse.transport import transmission_profile
from validation.check_aerosol import (
    aerosol_energy_check, curved_aerosol_checks, straight_exponential_checks,
)


def test_straight_exponential_slant_tau_matches_exact_bessel_integral():
    report=straight_exponential_checks()
    assert report["pass"],report


def test_curved_slant_tau_matches_independent_cartesian_ode():
    report=curved_aerosol_checks()
    assert report["pass"],report


def test_aerosol_receiver_flux_matches_incident_pupil_budget():
    report=aerosol_energy_check(radius_samples=1025)
    assert report["pass"],report


def test_zero_aerosol_exactly_restores_molecular_tables_and_irradiance():
    atmosphere=Atmosphere()
    baseline=build_ray_table(atmosphere)
    zero=build_ray_table(atmosphere,aerosol_profile=replace(AerosolProfile(),aod550=0.))
    for field in ["altitude_m","impact_m","bending_rad","air_column_m2","ozone_column_m2"]:
        np.testing.assert_array_equal(getattr(zero,field),getattr(baseline,field))
    assert zero.aerosol_tau550 is None or np.all(zero.aerosol_tau550==0.)
    radius=np.array([0.,1e6,2.5e6,4e6,4.6e6,5e6,9e6])
    wavelengths=np.array([450.,550.,677.])
    np.testing.assert_array_equal(transmission_profile(zero,radius,wavelengths),
                                  transmission_profile(baseline,radius,wavelengths))


def test_aerosol_changes_no_molecular_column_or_ray_geometry():
    atmosphere=Atmosphere()
    baseline=build_ray_table(atmosphere)
    aerosol=build_ray_table(atmosphere,aerosol_profile=AerosolProfile(bottom_m=2134.,top_m=12345.))
    common,i,j=np.intersect1d(baseline.altitude_m,aerosol.altitude_m,return_indices=True)
    assert len(common)==len(baseline.altitude_m)
    for field in ["impact_m","bending_rad","air_column_m2","ozone_column_m2"]:
        np.testing.assert_allclose(getattr(baseline,field)[i],getattr(aerosol,field)[j],rtol=2e-12,atol=0.)
    assert 2134. in aerosol.altitude_m and 12345. in aerosol.altitude_m


def test_fixed_550nm_anchor_is_independent_of_requested_wavelength_grid():
    table=build_ray_table(Atmosphere(),aerosol_profile=AerosolProfile())
    radius=np.array([0.,1e6,3e6,4.5e6,5e6,9e6])
    whole=transmission_profile(table,radius,[450.,550.,677.])
    sliced=transmission_profile(table,radius,[550.])
    np.testing.assert_array_equal(table.aerosol_profile.spectral_factor([550.]),
                                  table.aerosol_profile.spectral_factor([450.,550.,677.])[1:2])
    # A one-column BLAS reduction can round differently from a matrix product.
    np.testing.assert_allclose(sliced[:,0],whole[:,1],rtol=2e-14,atol=0.)


def test_extinction_increases_only_attenuation_and_can_be_disabled():
    aerosol=build_ray_table(Atmosphere(),aerosol_profile=AerosolProfile())
    radius=np.linspace(0,9e6,61)
    wavelengths=np.array([450.,550.,677.])
    # Compare on the same height/impact mesh: nonzero aerosol adds cutoff
    # refinement nodes, which slightly change interpolation even if disabled.
    plain=transmission_profile(aerosol,radius,wavelengths,aerosol=False)
    attenuated=transmission_profile(aerosol,radius,wavelengths)
    assert np.all(attenuated<=plain+1e-14)
    assert np.any(attenuated<.9*plain)
    np.testing.assert_allclose(transmission_profile(aerosol,radius,wavelengths,extinction=False),
                               transmission_profile(aerosol,radius,wavelengths,extinction=False,aerosol=False),
                               rtol=0.,atol=1e-14)
