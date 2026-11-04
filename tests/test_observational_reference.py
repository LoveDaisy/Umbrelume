"""Independent rough visual anchors and contracts for the published curve readout."""
import json

import pytest

from validation.compare_observations import DATA, load_observations


@pytest.mark.parametrize("wavelength,angle,visual_value", [
    (503,.400,2.2e-5), (503,.500,1.2e-4), (503,.575,3.0e-4),
    (677,.150,1.6e-4), (677,.400,5.0e-4), (677,.575,1.5e-3),
])
def test_curve_values_match_independently_read_figure_anchors(wavelength,angle,visual_value):
    # Rounded readings of the rendered source figure, deliberately independent
    # of PDF path coordinates. This catches selecting the dashed model instead.
    row=next(row for row in load_observations()
             if row["wavelength_nm"]==wavelength and row["shadow_angle_deg"]==angle)
    assert row["observed_relative_brightness"]==pytest.approx(visual_value,rel=.06)


def test_reference_keeps_observations_separate_from_model_and_measurement_errors():
    source=json.loads((DATA/"source.json").read_text())
    assert "NOT original" in source["source_status"]
    assert "NOT observational uncertainty" in source["readout_definition"]
    assert len(load_observations())==18
    assert set(row["wavelength_nm"] for row in load_observations())=={503,677}


def test_event_geometry_is_independently_consistent_with_rounded_nasa_ephemeris():
    geometry=json.loads((DATA/"geometry.json").read_text())
    # NASA's event sheet at greatest eclipse gives horizontal parallax58'33.0''
    # (~374509km); observations occur~43min earlier, giving a slightly smaller D.
    assert 374_200_000<geometry["moon_distance_m"]<374_600_000
    assert 1.015<geometry["sun_distance_m"]/149_597_870_700<1.017
    assert geometry["moon"]["ranges"][0]["distance_m"]<geometry["moon_distance_m"]
    assert geometry["moon"]["ranges"][-1]["distance_m"]>geometry["moon_distance_m"]
