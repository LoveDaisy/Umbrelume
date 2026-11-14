"""Synthetic bandpass-kernel checks and observation data contracts.

Pure numpy tests only: no ray tables, no model calls, no network. The kernel
functions are exact analytic targets on constant and linear transmissions.
"""
import json
import math

import numpy as np
import pytest

from validation import filter_response_2011 as filter_response
from validation.filter_response_2011 import (
    BANDS_NM, GAUSSIAN_TRUNCATED_TAIL_FRACTION, NOMINAL_WIDTH_NM, SIGMA_CUTOFF,
    SHAPES, WIDTHS_NM, band_mean, classify_outermost, cross_check_gate,
    gaussian_sigma_nm, response_values,
    support_grid, union_grid, grid_indices, variant_grids,
    load_obs2011_predictions, load_round2_predictions)
from validation.compare_observations import DATA


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("width", WIDTHS_NM)
def test_kernel_leaves_constant_transmission_invariant(shape, width):
    grid = support_grid(677.0, SIGMA_CUTOFF*gaussian_sigma_nm(width) if shape == "gauss"
                        else width/2.0, 1.0)
    weights = response_values(grid, 677.0, shape, width)
    constant = np.full(grid.shape, 0.0125)
    assert band_mean(weights, constant) == pytest.approx(0.0125, rel=1e-12)


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("width", WIDTHS_NM)
def test_symmetric_kernel_reduces_linear_transmission_to_center_value(shape, width):
    half = width/2.0 if shape == "rect" else SIGMA_CUTOFF*gaussian_sigma_nm(width)
    grid = support_grid(503.0, half, 1.0)
    weights = response_values(grid, 503.0, shape, width)
    linear = 0.4 + 3.1e-4*(grid-503.0)
    assert band_mean(weights, linear) == pytest.approx(0.4, rel=1e-12, abs=1e-15)


def test_band_mean_is_scale_invariant_and_rejects_degenerate_weights():
    grid = support_grid(677.0, 10.0, 1.0)
    weights = response_values(grid, 677.0, "rect", 20.0)
    values = np.linspace(0.9, 1.1, grid.size)
    assert band_mean(2.0*weights, values) == pytest.approx(band_mean(weights, values), rel=1e-14)
    with pytest.raises(ValueError):
        band_mean(np.zeros_like(weights), values)


def test_rectangular_response_support_is_inclusive_and_hard_edged():
    grid = np.array([482.9, 483.0, 503.0, 513.1, 523.0, 523.1])
    response = response_values(grid, 503.0, "rect", 40.0)
    assert response.tolist() == [0.0, 1.0, 1.0, 1.0, 1.0, 0.0]


def test_gaussian_response_peak_cutoff_and_truncation_tail():
    assert response_values(np.array([677.0]), 677.0, "gauss", 20.0)[0] == 1.0
    sigma = gaussian_sigma_nm(20.0)
    beyond = response_values(np.array([677.0+3.1*sigma]), 677.0, "gauss", 20.0)
    assert beyond[0] == 0.0
    # Two-sided mass beyond the +/-3 sigma support, analytic and bounded.
    assert GAUSSIAN_TRUNCATED_TAIL_FRACTION == pytest.approx(math.erfc(3.0/math.sqrt(2.0)),
                                                             rel=1e-15)
    assert GAUSSIAN_TRUNCATED_TAIL_FRACTION < 0.003


def test_variant_grids_nest_exactly_inside_the_union_grid():
    for band in BANDS_NM:
        grids = variant_grids(band)
        union = union_grid(grids)
        for grid in grids.values():
            indices = grid_indices(union, grid)
            assert np.all(union[indices] == grid)
        # The refined widest Gaussian support bounds every other kernel grid.
        widest = grids[("gauss", max(WIDTHS_NM), 0.5)]
        assert union.size == widest.size


def test_committed_filter_information_contract():
    source = json.loads((DATA/"source.json").read_text())
    info = source["filter_information"]
    assert info["branch"] == "C"
    assert info["status"] == "not stated in paper"
    assert info["approximation_grade"].startswith("Assumed rectangular and Gaussian")
    assert info["source_pdf_sha256"] == source["source_pdf_sha256"]
    assert len(info["search"]["keywords"]) >= 20
    assert all(count == 0 for count in info["search"]["hit_counts"].values())
    assert info["search"]["hit_counts"].keys() == set(info["search"]["keywords"])
    assert len(info["filter_statements_present"]) >= 2
    assert source["missing_information"][0] == "Filter response curves and FWHM"


def test_committed_cross_check_references_contain_expected_fields():
    round2 = load_round2_predictions()
    obs2011 = load_obs2011_predictions()
    assert len(round2) == 18
    assert all("molecular_prediction" in value and "fixed_opac_prediction" in value
               for value in round2.values())
    assert len(obs2011) == 4
    assert all("model_event_300du" in value for value in obs2011.values())


def test_classify_outermost_thresholds_are_pinned():
    # Audit anchor for the thresholds quoted in docs/records/OBSERVATIONAL_VALIDATION.md.
    assert "materially reduces" in classify_outermost(0.5)
    assert "marginally" in classify_outermost(2.0/3.0)
    assert "marginally" in classify_outermost(1.5)
    assert "does not reduce" in classify_outermost(1.5 + 1e-12)


def _gate_with_predictions(monkeypatch, molecular, opac, obs2011_value):
    key = (503.0, 0.400)
    monkeypatch.setattr(filter_response, "load_round2_predictions",
                        lambda: {key: {"molecular_prediction": molecular,
                                       "fixed_opac_prediction": opac}})
    monkeypatch.setattr(filter_response, "load_obs2011_predictions",
                        lambda: {key: {"model_event_300du": obs2011_value}})
    observations = [{"wavelength_nm": 503.0, "shadow_angle_deg": 0.400}]
    mono = {"molecular": np.array([molecular]), "fixed_opac": np.array([opac])}
    return mono, observations


def test_cross_check_gate_reproduced_on_exact_match(monkeypatch):
    mono, observations = _gate_with_predictions(monkeypatch, 1.0, 1.0, 1.0)
    result = cross_check_gate(mono, observations, None)
    assert result["status"] == "reproduced"
    assert result["adjudicated_rtol"] is None


def test_cross_check_gate_middle_band_pauses_without_adjudication(monkeypatch, capsys):
    mono, observations = _gate_with_predictions(monkeypatch, 1.0*(1+5e-10), 1.0, 1.0)
    with pytest.raises(RuntimeError, match="gate stopped the run"):
        cross_check_gate(mono, observations, None)
    assert "pause: owner adjudication required" in capsys.readouterr().out


def test_cross_check_gate_middle_band_accepts_covering_adjudication(monkeypatch):
    mono, observations = _gate_with_predictions(monkeypatch, 1.0*(1+5e-10), 1.0, 1.0)
    result = cross_check_gate(mono, observations, 1e-9)
    assert result["status"] == "adjudicated"
    assert result["adjudicated_rtol"] == 1e-9


def test_cross_check_gate_middle_band_rejects_non_covering_adjudication(monkeypatch, capsys):
    mono, observations = _gate_with_predictions(monkeypatch, 1.0*(1+5e-10), 1.0, 1.0)
    with pytest.raises(RuntimeError, match="gate stopped the run"):
        cross_check_gate(mono, observations, 1e-11)
    assert "pause: owner adjudication required" in capsys.readouterr().out


def test_cross_check_gate_blocks_above_ceiling(monkeypatch, capsys):
    mono, observations = _gate_with_predictions(monkeypatch, 1.0*(1+1e-6), 1.0, 1.0)
    with pytest.raises(RuntimeError, match="gate stopped the run"):
        cross_check_gate(mono, observations, 1e-3)
    assert "block: rtol above the adjudication ceiling" in capsys.readouterr().out
