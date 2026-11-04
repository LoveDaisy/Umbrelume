"""Tests of colorimetric units, shared white and global exposure."""
import numpy as np
import pytest

from eclipse.color import (
    display_rgb, encode_srgb, finish_display, linear_display_rgb,
    spectrum_to_xyz, xyz_to_linear_srgb,
)
from eclipse.spectra import cie_xyz, solar_irradiance


@pytest.fixture
def source():
    wl = np.arange(380.0, 781.0, 5.0)
    return wl, solar_irradiance(wl)


def test_colorimetric_integration_and_linearity(source):
    wl, solar = source
    expected = np.array([np.trapezoid(solar * cie_xyz(wl)[:, j], x=wl) for j in range(3)])
    actual = spectrum_to_xyz(np.stack([solar, solar * 0.2]), wl)
    np.testing.assert_allclose(actual[0], expected, rtol=1e-14)
    np.testing.assert_allclose(actual[1], expected * 0.2, rtol=1e-14)
    assert actual.shape == (2, 3)


def test_d65_srgb_white():
    white = xyz_to_linear_srgb([0.9504559270516716, 1.0, 1.0890577507598784])
    np.testing.assert_allclose(white, [1.0, 1.0, 1.0], atol=1e-12)


def test_reference_white_and_black(source):
    wl, solar = source
    np.testing.assert_allclose(display_rgb(solar, wl, solar), np.ones(3), atol=1e-12)
    np.testing.assert_array_equal(display_rgb(solar * 0.0, wl, solar), np.zeros(3))


def test_exposure_preserves_relative_light_before_display_limit(source):
    wl, solar = source
    values = np.stack([solar * 0.01, solar * 0.02])
    linear = linear_display_rgb(values, wl, solar)
    np.testing.assert_allclose(linear[1], 2.0 * linear[0], rtol=1e-13)
    np.testing.assert_allclose(display_rgb(values, wl, solar, exposure_stops=1),
                               encode_srgb(linear * 2), atol=1e-13)
    assert np.all(display_rgb(values, wl, solar)[1] > display_rgb(values, wl, solar)[0])


def test_different_spectrum_is_not_individually_whitened(source):
    wl, solar = source
    red_filter = np.exp(-((wl - 650) / 50) ** 2)
    rgb = display_rgb(solar * red_filter * 0.1, wl, solar)
    assert rgb[0] > 3 * rgb[2]
    np.testing.assert_allclose(display_rgb(solar * red_filter * 0.01, wl, solar,
                                          normalize_luminance=True),
                               display_rgb(solar * red_filter * 0.1, wl, solar,
                                           normalize_luminance=True), atol=1e-13)


def test_highlight_compression_preserves_channel_ratios():
    expected = encode_srgb([1.0, 0.5, 0.25])
    np.testing.assert_allclose(finish_display([8, 4, 2]), expected)


@pytest.mark.parametrize("wl, spectra", [
    ([500, 400], [1, 1]), ([400, 400], [1, 1]), ([400], [1]),
    ([400, np.nan], [1, 1]), ([400, 500], [-1, 1]),
    ([400, 500], [np.inf, 1]), ([400, 500], [1, 2, 3]),
])
def test_invalid_spectral_input(wl, spectra):
    with pytest.raises(ValueError):
        spectrum_to_xyz(spectra, wl)


def test_invalid_reference_and_exposure(source):
    wl, solar = source
    with pytest.raises(ValueError, match="positive"):
        display_rgb(solar, wl, solar * 0)
    with pytest.raises(ValueError, match="one-dimensional"):
        display_rgb(solar, wl, np.stack([solar, solar]))
    with pytest.raises(ValueError, match="exposure"):
        display_rgb(solar, wl, solar, exposure_stops=np.nan)
