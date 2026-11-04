"""Physical checks using independent propagation and energy accounting."""
from validation.check_physics import energy_check, ray_checks


def test_radial_paths_match_independent_cartesian_ode():
    report = ray_checks()
    assert report["pass"], report


def test_global_receiver_flux_matches_incident_pupil_budget():
    report = energy_check()
    assert report["pass"], report


def test_transparent_refraction_preserves_flux():
    report = energy_check(extinction=False)
    assert report["pass"], report
