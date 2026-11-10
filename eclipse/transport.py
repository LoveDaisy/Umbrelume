"""Finite-Sun irradiance from conserved forward annuli, without ray histograms.

All distances are metres. Transport from the Sun through the terrestrial limb
to the receiver is paraxial; the atmosphere's ray bending itself is integrated
through its curved spherical profile. The finite Sun is a uniform disk.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .atmosphere import Atmosphere, EARTH_RADIUS_M
from .aerosol import AerosolColumn, AerosolPrior
from .rays import RayTable, build_ray_table
from .spectra import ozone_cross_section, rayleigh_cross_section, solar_irradiance

AU_M = 149_597_870_700.0
SUN_RADIUS_M = 695_700_000.0
MOON_DISTANCE_M = 384_400_000.0
MOON_RADIUS_M = 1_737_400.0


@dataclass(frozen=True)
class Geometry:
    moon_distance_m: float = MOON_DISTANCE_M
    sun_distance_m: float = AU_M
    sun_radius_m: float = SUN_RADIUS_M

    def __post_init__(self):
        values = (self.moon_distance_m, self.sun_distance_m, self.sun_radius_m)
        if not all(np.isfinite(x) and x > 0 for x in values):
            raise ValueError("Geometry distances must be positive and finite")
        if self.moon_distance_m <= EARTH_RADIUS_M + 120_000:
            raise ValueError("Receiver must be outside Earth's atmosphere")
        if self.sun_radius_m >= self.sun_distance_m:
            raise ValueError("Solar radius must be smaller than solar distance")

    @property
    def magnification(self) -> float:
        return 1 + self.moon_distance_m / self.sun_distance_m

    @property
    def solar_blur_m(self) -> float:
        return self.moon_distance_m * self.sun_radius_m / self.sun_distance_m


def circle_overlap(distance, radius_a: float, radius_b: float) -> np.ndarray:
    """Area common to two circles, including containment and disjoint limits."""
    if not all(np.isfinite(x) and x > 0 for x in (radius_a, radius_b)):
        raise ValueError("Circle radii must be positive and finite")
    d = np.asarray(distance, dtype=float)
    if np.any(~np.isfinite(d)) or np.any(d < 0):
        raise ValueError("Separation must be nonnegative and finite")
    result = np.zeros_like(d)
    inside = d <= abs(radius_a - radius_b)
    result[inside] = np.pi * min(radius_a, radius_b)**2
    partial = (d > abs(radius_a - radius_b)) & (d < radius_a + radius_b)
    dd = d[partial]
    ca = np.clip((dd**2 + radius_a**2 - radius_b**2) / (2*dd*radius_a), -1, 1)
    cb = np.clip((dd**2 + radius_b**2 - radius_a**2) / (2*dd*radius_b), -1, 1)
    radicand = ((-dd + radius_a + radius_b) * (dd + radius_a - radius_b)
                * (dd - radius_a + radius_b) * (dd + radius_a + radius_b))
    result[partial] = (radius_a**2*np.arccos(ca) + radius_b**2*np.arccos(cb)
                       - 0.5*np.sqrt(np.maximum(radicand, 0)))
    return result


def ring_fraction(receiver_radius_m: float, ring_radius_m, blur_radius_m: float):
    """Azimuth fraction of a ring inside the receiver's solar-source disk."""
    s = np.asarray(ring_radius_m, dtype=float)
    r, a = receiver_radius_m, blur_radius_m
    if not np.isfinite(r) or r < 0 or not np.isfinite(a) or a <= 0:
        raise ValueError("Invalid receiver or blur radius")
    if np.any(~np.isfinite(s)) or np.any(s < 0):
        raise ValueError("Ring radii must be nonnegative and finite")
    f = np.zeros_like(s)
    f[r + s <= a] = 1
    partial = (np.abs(r-s) < a) & (r+s > a)
    if r > 0:
        cosine = (r*r + s[partial]**2 - a*a) / (2*r*s[partial])
        f[partial] = np.arccos(np.clip(cosine, -1, 1)) / np.pi
    return f


def geometric_transmission(radius_m, geometry: Geometry = Geometry()):
    """No-atmosphere eclipse in the same finite-distance paraxial geometry."""
    a = geometry.solar_blur_m
    return 1 - circle_overlap(radius_m, a, geometry.magnification*EARTH_RADIUS_M) / (
        np.pi*a*a)


def signed_landing_radius(table: RayTable, geometry: Geometry) -> np.ndarray:
    """First-order finite-source map; negative rays cross the shadow axis."""
    return geometry.magnification*table.impact_m - geometry.moon_distance_m*table.bending_rad


def receiver_quadrature(table: RayTable, radius_m: float, geometry: Geometry,
                        integration_order: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Return incident b samples and conserved finite-source area weights.

    Shared by transport and altitude-contribution diagnostics; no extinction
    or wavelength weighting is applied here. Every landing branch is retained.
    """
    b = table.impact_m
    landing = signed_landing_radius(table, geometry)
    a, m, r = geometry.solar_blur_m, geometry.magnification, radius_m
    nodes, weights = leggauss(integration_order)
    target = np.unique(np.array([r+a, r-a, -r+a, -r-a, 0.0]))
    roots = []
    for boundary in target:
        crossing = ((landing[:-1] - boundary) * (landing[1:] - boundary)) < 0
        roots.extend(b[:-1][crossing] + (boundary - landing[:-1][crossing])
                     / np.diff(landing)[crossing] * np.diff(b)[crossing])
    edges = np.unique(np.concatenate((b, roots)))
    mid = (edges[1:]+edges[:-1])/2
    ymid = np.interp(mid, b, landing)
    support = ring_fraction(r, np.abs(ymid), a) > 0
    lo, hi = edges[:-1][support], edges[1:][support]
    if not len(lo):
        return np.empty(0), np.empty(0)
    bb = ((hi+lo)[:, None]/2 + (hi-lo)[:, None]/2*nodes).ravel()
    db = ((hi-lo)[:, None]/2 * np.broadcast_to(weights, (len(lo), len(weights)))).ravel()
    yy = np.interp(bb, b, landing)
    area = 2*m*m*bb*db/a**2 * ring_fraction(r, np.abs(yy), a)
    return bb, area


def transmission_profile(
    table: RayTable,
    radius_m,
    wavelength_nm,
    *,
    geometry: Geometry = Geometry(),
    ozone: bool = True,
    aerosol: bool = True,
    extinction: bool = True,
    integration_order: int = 8,
) -> np.ndarray:
    """E_lambda/F_lambda at the receiver, preserving absolute annular weights.

    Each ring has incident area 2*pi*b*db. A uniform solar disk convolves its
    landing ring with a disk of radius a. Relative to unocculted irradiance at
    the receiver, the weight is 2*magnification**2*b*db/a**2. No derivative of
    the landing map and no extra focusing factor are required.
    """
    radius = np.asarray(radius_m, dtype=float)
    wl = np.asarray(wavelength_nm, dtype=float)
    if radius.ndim != 1 or np.any(~np.isfinite(radius)) or np.any(radius < 0):
        raise ValueError("Receiver radii must be a finite nonnegative vector")
    if wl.ndim != 1 or wl.size < 1 or np.any(~np.isfinite(wl)) or np.any(wl <= 0):
        raise ValueError("Wavelengths must be a positive finite vector")
    if integration_order < 2:
        raise ValueError("integration_order must be at least 2")
    b = table.impact_m
    if geometry.moon_distance_m <= EARTH_RADIUS_M + table.altitude_m[-1]:
        raise ValueError("Receiver must be outside the modeled atmosphere")
    a, m = geometry.solar_blur_m, geometry.magnification
    outer = 1 - circle_overlap(radius, a, m*b[-1])/(np.pi*a*a)
    result = np.repeat(outer[:, None], len(wl), axis=1)
    sigma_air = rayleigh_cross_section(wl) if extinction else np.zeros_like(wl)
    sigma_o3 = ozone_cross_section(wl) if extinction and ozone else np.zeros_like(wl)
    for i, r in enumerate(radius):
        bb, area = receiver_quadrature(table, float(r), geometry, integration_order)
        if not len(bb):
            continue
        air = np.interp(bb, b, table.air_column_m2)
        o3 = np.interp(bb, b, table.ozone_column_m2)
        tau = air[:, None]*sigma_air + o3[:, None]*sigma_o3
        if extinction and aerosol and table.aerosol_layer_taus is not None:
            # The only extinction hook: each layer contributes its own tau
            # table times its own spectral factor. The (tau, profile) pairs
            # come from the table itself, so they are same-source by
            # construction; a single-layer table reproduces the previous
            # single-term operation order bit-exactly.
            for layer_tau, layer_profile in zip(table.aerosol_layer_taus,
                                                table.aerosol_profiles):
                tau += (np.interp(bb, b, layer_tau)[:, None]
                        * layer_profile.spectral_factor(wl))
        result[i] += area @ np.exp(-tau)
    return result


@dataclass(frozen=True)
class EclipseResult:
    radius_m: np.ndarray
    wavelength_nm: np.ndarray
    irradiance: np.ndarray
    reference_irradiance: np.ndarray
    irradiance_no_ozone: np.ndarray
    ray_table: RayTable
    geometry: Geometry


def simulate(
    radius_m=None,
    wavelength_nm=None,
    *,
    atmosphere: Atmosphere | None = None,
    geometry: Geometry = Geometry(),
    height_step_m: float = 100.0,
    quadrature_order: int = 48,
    integration_order: int = 8,
    aerosol_profile: AerosolPrior = AerosolColumn.opac_default(),
) -> EclipseResult:
    """Precompute forward rays once, then integrate ozone on/off at all receivers.

    The default aerosol prior is the production two-layer OPAC column (free
    troposphere + stratospheric background); pass ``None`` for the
    molecular-only medium, a single ``AerosolProfile`` for one layer, or an
    ``AerosolColumn`` for an explicit layer composition. The parameter keeps
    the historical name ``aerosol_profile`` for every existing call site.
    """
    radius = np.asarray(np.arange(0, 9_000_001, 20_000) if radius_m is None else radius_m,
                        dtype=float)
    wl = np.asarray(np.arange(380, 781, 5) if wavelength_nm is None else wavelength_nm,
                    dtype=float)
    table = build_ray_table(atmosphere or Atmosphere(), height_step_m=height_step_m,
                            quadrature_order=quadrature_order, aerosol_profile=aerosol_profile)
    reference = solar_irradiance(wl)*(AU_M/(geometry.sun_distance_m + geometry.moon_distance_m))**2
    on = transmission_profile(table, radius, wl, geometry=geometry,
                               integration_order=integration_order)
    off = transmission_profile(table, radius, wl, geometry=geometry, ozone=False,
                                integration_order=integration_order)
    return EclipseResult(radius, wl, on*reference, reference, off*reference, table, geometry)
