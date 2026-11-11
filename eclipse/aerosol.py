"""Fixed, independent OPAC free-tropospheric and stratospheric-background priors.

The default describes a climatological/model component, not a measured profile
for the June 2011 eclipse. Extinction includes scattering out of the direct
beam and absorption. No scattering-into-beam or aerosol refractivity is added.
The production default composes two layers: the OPAC free troposphere (2-12 km)
and the OPAC stratospheric background (12-35 km).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Provenance of the pinned prior rows. The numeric values themselves are
# anchored to these rows by pytest (json-anchor test); the SHA identifies the
# OPAC source paper recorded in sources.json alongside the parameter rows.
PRIOR_PARAMETERS_JSON = (
    "validation/results/extinction_gap_diagnosis/prior_parameters.json")
OPAC_PDF_SHA256 = ("68a20ba9488db1e0f65a9976d4ee46e6c879f89054fcff583ba67b"
                   "b5d68d425c")
FREE_TROPOSPHERE_SOURCE = {
    "candidate": "free_troposphere_reference",
    "rows": ("aod550_2_to_12km_reference", "angstrom_0.35_to_0.5um",
             "angstrom_0.5_to_0.8um", "scale_height_m"),
}
STRATOSPHERE_BACKGROUND_SOURCE = {
    "candidate": "stratosphere_background",
    "rows": ("aod550", "angstrom_0.35_to_0.5um", "angstrom_0.5_to_0.8um",
             "layer_bottom_m", "layer_top_m", "scale_height_m"),
}
# Layer order of ``AerosolColumn.opac_default()`` paired with each layer's
# pinned prior-row provenance. CLI metadata and validation products consume
# this single pairing, so the layer order and its source rows are written in
# exactly one place and cannot drift apart when the default column evolves.
OPAC_DEFAULT_LAYER_SOURCES = (FREE_TROPOSPHERE_SOURCE, STRATOSPHERE_BACKGROUND_SOURCE)


# --- Humidity axis (pinned OPAC 1998 component growth; task-humidity-axis) ------
#
# The free-troposphere layer carries an optional humidity growth: the OPAC 1998
# WASO component's 550 nm extinction per particle is pinned at three
# relative-humidity knots (recovered component files, sha256-pinned in
# ``humidity_growth.json``), composed into a column-AOD growth factor with the
# paper's FT composition (Sec 5a, printed p. 842: WASO number mixing ratio 0.6
# of a 730 particles/cm3 sea-level density over the 2-12 km layer with 8 km
# scale height) and the published 50%-RH reference AOD550 (Table 3 FT row
# 0.013). The non-growing soot+insoluble remainder is pinned by subtraction
# (``Q_REST``): soot and insoluble particles are hydrophobic (the OPAC dataset
# pins them only as 00 files), so their extinction does not grow with RH.
# At the 50% reference state the factor is exactly 1.0 and the default output
# is bit-exact unchanged. The growth rescales only the column AOD (J1): the
# layer's vertical shape and its broken-power-law spectral factor never change,
# and the ratio is invariant to the vertical profile reference because
# numerator and denominator share the same geometry.
#
# Transcription product with the pinned quotes and preregistered semantics:
# ``validation/results/prior_source_survey/humidity_axis_priors.json``.
# The humidity section is a removable component: deleting it (and restoring
# the factories' ``rh`` parameters) restores the previous state exactly.

HUMIDITY_GROWTH_JSON = (
    "validation/results/prior_source_survey/humidity_growth.json")
HUMIDITY_AXIS_PRIORS_JSON = (
    "validation/results/prior_source_survey/humidity_axis_priors.json")

FT_AOD550_PUBLISHED = 0.013
RH_MIN = 0.5
RH_MAX = 0.99
WASO_NUMBER_MIXING_RATIO = 0.6
FT_SEA_LEVEL_NUMBER_DENSITY_PER_CM3 = 730.0
# Paper FT vertical geometry (Table 5 FT row H = 8 km; Table 3 FT row layer
# bounds 2-12 km; Sec 5b: "the upper boundary of this layer is always situated
# at 12 km"). The repository FT layer uses the same geometry (mechanically
# compared in tests); the growth factor's geometry uses the paper's
# sea-level-referenced profile N(z) = N(0)*exp(-z/H).
FT_BOTTOM_KM = 2.0
FT_TOP_KM = 12.0
FT_SCALE_HEIGHT_KM = 8.0
# Vertical column integral of exp(-z/H) over the FT layer, in km:
# int_2km^12km exp(-z/8km) dz = 8*(exp(-2/8) - exp(-12/8)).
FT_GEOMETRY_FACTOR_KM = FT_SCALE_HEIGHT_KM*(
    math.exp(-FT_BOTTOM_KM/FT_SCALE_HEIGHT_KM)
    - math.exp(-FT_TOP_KM/FT_SCALE_HEIGHT_KM))

# Pinned WASO ext550 per particle [km^-1 per (particle cm^-3)] at the growth
# knots, transcribed from ``humidity_growth.json`` (WASO component; values and
# sha256 are anchored to that product by pytest). The dry 0% knot lies outside
# the pinned RH domain and is not consumed.
WASO_EXT550_PER_PARTICLE_KM = {
    0.5: 6.366e-06,
    0.7: 7.913e-06,
    0.99: 4.740000000000001e-05,
}
WASO_EXT550_SHA256 = {
    0.5: "eab5e8e27a470b3ee6c350698a1d6e4b34f407d3b3ae28f1fda20a237039579e",
    0.7: "33a13edffc37e6d8348ae5521fd03f8d36311878f2436649c7f3c9f02a2bbb80",
    0.99: "bd8209c423030b555910dbb4b1beb54f3d70301c60ca434db0f484573ef430fc",
}

_GROWTH_KNOTS = tuple(sorted(WASO_EXT550_PER_PARTICLE_KM))


def _linear_interp(x: float, knots: tuple, values: tuple) -> float:
    """Piecewise-linear interpolation in the lerp form ``v0 + w*(v1 - v0)``.

    ``x`` must lie within ``[knots[0], knots[-1]]``; knot queries return the
    knot value exactly (weight 0 or the exact Sterbenz subtraction property).
    """
    for lo, hi, lo_value, hi_value in zip(knots, knots[1:], values, values[1:]):
        if x <= hi:
            weight = (x-lo)/(hi-lo)
            return lo_value + weight*(hi_value-lo_value)
    return values[-1]


def _waso_column_aod550(rh: float) -> float:
    """E_waso(rh): the WASO share of the FT column AOD550 at RH = rh.

    Linear interpolation of the pinned ext550 knots in RH (the preregistered
    semantics), times the paper's composition and vertical geometry. Units
    cancel to a dimensionless column contribution: (particle cm^-3) times
    (km^-1 per particle cm^-3) times km.
    """
    if not _GROWTH_KNOTS[0] <= rh <= _GROWTH_KNOTS[-1]:
        raise ValueError(f"rh {rh!r} outside the pinned WASO knot range "
                         f"[{_GROWTH_KNOTS[0]}, {_GROWTH_KNOTS[-1]}]")
    ext550 = _linear_interp(rh, _GROWTH_KNOTS,
                            tuple(WASO_EXT550_PER_PARTICLE_KM[knot]
                                  for knot in _GROWTH_KNOTS))
    return (WASO_NUMBER_MIXING_RATIO*FT_SEA_LEVEL_NUMBER_DENSITY_PER_CM3
            *ext550*FT_GEOMETRY_FACTOR_KM)


# The non-growing soot+insoluble remainder of the published FT reference AOD:
# q_rest = 0.013 - E_waso(50%). Their RH-independence rests on the dataset
# structure (00-only files) and the paper's construction, so the subtraction
# pins the remainder without the unrecovered component files. The gate below
# fails loud at import on any transcription drift (plan risk 1).
Q_REST = FT_AOD550_PUBLISHED - _waso_column_aod550(RH_MIN)
if not 0.0 < Q_REST < FT_AOD550_PUBLISHED:
    raise RuntimeError(
        f"q_rest gate failed: {Q_REST!r} outside (0, {FT_AOD550_PUBLISHED}) - "
        f"the FT composition transcription or the WASO knots drifted")


def ft_growth_factor(rh: float) -> float:
    """Free-troposphere column-AOD growth factor f(RH) relative to 50% RH.

    ``f(RH) = (E_waso(RH) + Q_REST)/0.013`` with the pinned WASO knots and the
    linear-in-RH interpolation between them. Monotone nondecreasing on the
    pinned domain ``[0.5, 0.99]``; outside it the function fails loud - the
    sub-reference direction is deliberately not implemented (no discriminator
    consumes it). ``f(0.5)`` is exactly 1.0, which makes the default factory
    output bit-exact unchanged.
    """
    if not RH_MIN <= rh <= RH_MAX:
        raise ValueError(f"rh {rh!r} outside the pinned humidity domain "
                         f"[{RH_MIN}, {RH_MAX}]")
    e_waso = _waso_column_aod550(rh)
    return (e_waso+Q_REST)/FT_AOD550_PUBLISHED


def ft_humidity_provenance(rh: float, ft_aod550_base: float) -> dict:
    """Humidity-axis provenance fields for a free-troposphere metadata record.

    Single shape consumed by the CLI metadata; the growth factor is recomputed
    through ``ft_growth_factor`` (the authoritative implementation), so the
    metadata cannot disagree with the constructed layer's effective AOD.
    """
    return {"rh": rh,
            "ft_growth_factor": ft_growth_factor(rh),
            "ft_aod550_reference_base": ft_aod550_base,
            "humidity_growth_json": HUMIDITY_GROWTH_JSON,
            "humidity_axis_priors_json": HUMIDITY_AXIS_PRIORS_JSON,
            "waso_ext550_sha256": {f"{knot:g}": sha256
                                   for knot, sha256 in WASO_EXT550_SHA256.items()}}


# The bit-exactness contract of the default: at the 50% reference state the
# growth factor is exactly 1.0. Fail loud at import if the arithmetic identity
# ever breaks (a reference-value change must re-establish it deliberately).
assert ft_growth_factor(RH_MIN) == 1.0, \
    "f(0.5) must be exactly 1.0 (the default output is bit-exact unchanged)"


@dataclass(frozen=True)
class AerosolProfile:
    """Finite exponential layer with a fixed 550 nm vertical optical depth.

    Distances are geometric metres and extinction is m^-1. The defaults come
    from OPAC Tables 3 and 5 (Hess, Koepke & Schult, 1998): free troposphere
    2–12 km, scale height 8 km, AOD550 0.013, and spectral exponents 1.58
    over 500–800 nm and 1.21 over 350–500 nm.
    """

    aod550: float = 0.013
    bottom_m: float = 2_000.0
    top_m: float = 12_000.0
    scale_height_m: float = 8_000.0
    angstrom_visible: float = 1.58
    angstrom_blue: float = 1.21
    break_nm: float = 500.0

    def __post_init__(self):
        values = (self.aod550, self.bottom_m, self.top_m, self.scale_height_m,
                  self.angstrom_visible, self.angstrom_blue, self.break_nm)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("Aerosol parameters must be finite.")
        if self.aod550 < 0:
            raise ValueError("AOD550 must be nonnegative.")
        if self.bottom_m < 0 or self.top_m <= self.bottom_m:
            raise ValueError("Aerosol layer requires 0 <= bottom_m < top_m.")
        if self.scale_height_m <= 0 or self.break_nm <= 0:
            raise ValueError("Scale height and spectral break must be positive.")
        if self.angstrom_visible < 0 or self.angstrom_blue < 0:
            raise ValueError("Angstrom exponents must be nonnegative.")

    @classmethod
    def opac_free_troposphere(cls, aod550: float = 0.013,
                              rh: float = 0.5) -> AerosolProfile:
        """Published OPAC shape and spectral law, with explicit column control.

        Parameters are pinned to ``prior_parameters.json`` rows
        ``free_troposphere_reference`` (aod550_2_to_12km_reference 0.013,
        Angstrom 1.21 over 350-500 nm and 1.58 over 500-800 nm, scale height
        8 km; the 2 km bottom and 12 km top are stated in that row's source
        location). The OPAC PDF sha256 is the module constant
        ``OPAC_PDF_SHA256``.

        ``rh`` applies the pinned FT humidity growth: the effective AOD is
        ``aod550*ft_growth_factor(rh)`` with ``aod550`` always the 50%-RH
        reference base. The default ``rh=0.5`` is the published reference
        state and reproduces the historical factory output bit-exactly
        (``ft_growth_factor(0.5)`` is exactly 1.0); the vertical shape and
        the spectral law never change with RH.
        """
        return cls(aod550=aod550*ft_growth_factor(rh))

    @classmethod
    def opac_stratospheric_background(cls) -> AerosolProfile:
        """OPAC stratospheric-background layer: 12-35 km, near homogeneous.

        Parameters are pinned to ``prior_parameters.json`` rows
        ``stratosphere_background`` (aod550 0.005, Angstrom 0.74 over
        350-500 nm and 1.14 over 500-800 nm, layer 12-35 km, scale height
        99 km). The 99 km scale height over the 23 km thickness is the OPAC
        representation of a homogeneous layer (exp(-23/99) ~ 0.79). OPAC
        fixes this layer independent of volcanic activity.
        """
        return cls(aod550=0.005, bottom_m=12_000.0, top_m=35_000.0,
                   scale_height_m=99_000.0, angstrom_visible=1.14,
                   angstrom_blue=0.74)

    @classmethod
    def none(cls) -> AerosolProfile:
        """Zero extinction, retaining the same profile API for counterfactuals."""
        return cls(aod550=0.0)

    def extinction550(self, altitude_m) -> np.ndarray:
        """Extinction coefficient at 550 nm [m^-1], zero outside the layer.

        The shifted exponential and expm1 denominator preserve the exact
        requested vertical AOD for layers far from sea level or very thin
        compared with their scale height. Layer endpoints have measure zero;
        quadrature must split at both boundaries.
        """
        altitude = np.asarray(altitude_m, dtype=float)
        if np.any(~np.isfinite(altitude)) or np.any(altitude < 0):
            raise ValueError("Geometric altitude must be finite and nonnegative [m].")
        thickness = self.top_m - self.bottom_m
        denominator = self.scale_height_m * -np.expm1(-thickness/self.scale_height_m)
        # Clipping the shifted coordinate avoids evaluating an exponentially
        # large unused value below a high-altitude layer.
        shifted = np.maximum(altitude-self.bottom_m, 0)
        value = self.aod550/denominator * np.exp(-shifted/self.scale_height_m)
        return np.where((altitude >= self.bottom_m) & (altitude <= self.top_m), value, 0.0)

    def spectral_factor(self, wavelength_nm) -> np.ndarray:
        """Continuous broken power law beta(lambda)/beta(550 nm).

        Normalization is tied to the physical 550 nm reference, independent
        of whether 550 nm is present in the caller's wavelength array.
        """
        wavelength = np.asarray(wavelength_nm, dtype=float)
        if np.any(~np.isfinite(wavelength)) or np.any(wavelength <= 0):
            raise ValueError("Wavelengths must be positive and finite [nm].")
        exponent = np.where(wavelength >= self.break_nm,
                            self.angstrom_visible, self.angstrom_blue)
        raw = (wavelength/self.break_nm)**(-exponent)
        reference_exponent = (self.angstrom_visible if 550.0 >= self.break_nm
                              else self.angstrom_blue)
        reference = (550.0/self.break_nm)**(-reference_exponent)
        return raw/reference


@dataclass(frozen=True)
class AerosolColumn:
    """Ordered composition of finite exponential aerosol layers.

    A column's extinction is the sum of its layers' coefficients, for
    arbitrary support layouts: overlapping supports add their optical depths,
    and the vertical AOD of a column is exactly the sum of its layers' AOD550
    (documented invariant, asserted in tests including an overlapping-support
    case). The production default's layers are disjoint; the calibration-lane
    presets stack a third increment layer that may overlap the default bands.
    Each layer carries its own geometry, AOD and broken-power-law spectral
    exponents; no merged spectral factor exists, the only legal expansion is
    the per-layer (tau, profile) pair.
    """

    layers: tuple[AerosolProfile, ...]

    def __post_init__(self):
        if not self.layers:
            raise ValueError("AerosolColumn requires at least one layer.")
        if not all(isinstance(layer, AerosolProfile) for layer in self.layers):
            raise ValueError("AerosolColumn layers must all be AerosolProfile instances.")

    @classmethod
    def opac_default(cls, ft_aod550: float = 0.013, rh: float = 0.5) -> AerosolColumn:
        """Production default: OPAC free troposphere + stratospheric background.

        ``ft_aod550`` overrides only the free-troposphere layer (the CLI
        ``--aod550`` semantics); the stratospheric background stays at its
        published 0.005. ``rh`` delegates the humidity growth to the FT
        factory (the single growth-consumption point): the effective FT AOD
        is ``ft_aod550*ft_growth_factor(rh)`` and ``ft_aod550`` stays the
        50% reference base - callers never pre-multiply the growth.
        """
        return cls((AerosolProfile.opac_free_troposphere(aod550=ft_aod550, rh=rh),
                    AerosolProfile.opac_stratospheric_background()))


# The complete type face of an aerosol prior argument: single layer, layer
# column, or no aerosol at all.
AerosolPrior = AerosolProfile | AerosolColumn | None


def active_aerosol_layers(aerosol: AerosolPrior) -> tuple[AerosolProfile, ...]:
    """Single normalization point for every aerosol prior argument.

    Returns the layers that carry extinction, in order. Zero-AOD layers are
    dropped: they change no physical quantity, and ray tables store tau tables
    only for layers that have one. A column mixing zero and nonzero layers
    therefore normalizes to its active subset, which is exactly what the
    layer-isolation counterfactuals compare.
    """
    if aerosol is None:
        return ()
    if isinstance(aerosol, AerosolColumn):
        layers = aerosol.layers
    elif isinstance(aerosol, AerosolProfile):
        layers = (aerosol,)
    else:
        raise TypeError(f"Unsupported aerosol prior: {type(aerosol).__name__}")
    return tuple(layer for layer in layers if layer.aod550 > 0)
