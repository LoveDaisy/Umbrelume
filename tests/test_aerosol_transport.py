"""Independent aerosol transport references and physical invariants."""
from dataclasses import replace

import numpy as np

from eclipse.aerosol import AerosolColumn, AerosolProfile
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
    assert zero.aerosol_layer_taus is None
    assert zero.aerosol_profiles == ()
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
    np.testing.assert_array_equal(table.aerosol_profiles[0].spectral_factor([550.]),
                                  table.aerosol_profiles[0].spectral_factor([450.,550.,677.])[1:2])
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


def test_single_layer_column_builds_the_same_table_as_the_bare_profile():
    atmosphere=Atmosphere()
    profile=AerosolProfile.opac_free_troposphere()
    bare=build_ray_table(atmosphere,aerosol_profile=profile)
    wrapped=build_ray_table(atmosphere,aerosol_profile=AerosolColumn((profile,)))
    for field in ["altitude_m","impact_m","bending_rad","air_column_m2","ozone_column_m2"]:
        np.testing.assert_array_equal(getattr(bare,field),getattr(wrapped,field))
    np.testing.assert_array_equal(bare.aerosol_layer_taus[0],wrapped.aerosol_layer_taus[0])
    assert bare.aerosol_profiles == wrapped.aerosol_profiles == (profile,)


def test_two_layer_table_is_the_layerwise_sum_of_single_layer_tables():
    """Decomposition identity at shared tangent nodes, bit-exact.

    The per-tangent quadrature depends only on the layer's own segment
    structure, never on the rest of the height grid, so the two-layer table's
    per-layer tau values equal the single-layer tables' values exactly at
    every shared tangent height.
    """
    atmosphere=Atmosphere()
    ft=AerosolProfile.opac_free_troposphere()
    strat=AerosolProfile.opac_stratospheric_background()
    two=build_ray_table(atmosphere,aerosol_profile=AerosolColumn((ft,strat)))
    ft_only=build_ray_table(atmosphere,aerosol_profile=ft)
    strat_only=build_ray_table(atmosphere,aerosol_profile=strat)
    # Both shared layer boundaries sit on the two-layer grid as exact nodes.
    for edge in (12_000.,35_000.):
        assert edge in two.altitude_m
    common=ft_only.altitude_m[np.isin(ft_only.altitude_m,strat_only.altitude_m)]
    mask_two=np.isin(two.altitude_m,common)
    total=two.aerosol_layer_taus[0][mask_two]+two.aerosol_layer_taus[1][mask_two]
    ft_at=ft_only.aerosol_layer_taus[0][np.isin(ft_only.altitude_m,common)]
    strat_at=strat_only.aerosol_layer_taus[0][np.isin(strat_only.altitude_m,common)]
    np.testing.assert_array_equal(total,ft_at+strat_at)
    # Support discipline through the ray geometry: a ray tangent at or above
    # a layer's top never enters it, so the FT tau vanishes for tangents
    # >= 12 km and the strat tau vanishes for tangents >= 35 km; tangents
    # below a layer's bottom still cross it on the way up and carry tau.
    heights=common
    ft_tau_at_height=ft_at
    assert np.all(ft_tau_at_height[heights>=12000.]==0.)
    assert np.all(ft_tau_at_height[heights<12000.]>0.)
    assert np.all(strat_at[heights>=35000.]==0.)
    assert np.all(strat_at[heights<35000.]>0.)


def test_two_layer_transmission_accumulates_both_layers():
    radius=np.array([0.,1e6,2.5e6,4e6,4.6e6,5e6,9e6])
    wavelengths=np.array([450.,550.,677.])
    atmosphere=Atmosphere()
    molecular=build_ray_table(atmosphere)
    two=build_ray_table(atmosphere,aerosol_profile=AerosolColumn.opac_default())
    ft_only=build_ray_table(atmosphere,aerosol_profile=AerosolProfile.opac_free_troposphere())
    molecular_t=transmission_profile(molecular,radius,wavelengths)
    ft_t=transmission_profile(ft_only,radius,wavelengths)
    two_t=transmission_profile(two,radius,wavelengths)
    assert np.all(two_t<=ft_t+1e-15)
    assert np.all(ft_t<=molecular_t+1e-15)
    # Both layers' spectral factors exceed one below the 500 nm break and
    # fall below one above it, so the dimming the stratospheric layer adds
    # is stronger in the blue channel than in the red one. The ordering is
    # asserted where the strat dimming is substantial: near the penumbra
    # edge the sampled rays tangent above 35 km and the different table
    # grids' interpolation noise swamp the spectral ordering.
    dimmed=ft_t[:,0]<1.
    blue_dimming=ft_t[dimmed,0]/two_t[dimmed,0]
    red_dimming=ft_t[dimmed,2]/two_t[dimmed,2]
    assert np.all(blue_dimming>1.) and np.all(red_dimming>1.)
    strong=(blue_dimming>1.05)&(red_dimming>1.05)
    assert np.any(strong)
    assert np.all(blue_dimming[strong]>red_dimming[strong])
