"""Forward tables of whole-limb refracted rays, before receiver sampling.

The ray label is its vacuum impact parameter b = n(r_t) r_t. Spherical
symmetry reduces bending and optical columns to one-dimensional quadratures.
Geometry is achromatic at the declared reference wavelength; extinction is not.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .atmosphere import (
    Atmosphere, EARTH_RADIUS_M, LAYER_BOUNDARIES_M, REFERENCE_NUMBER_DENSITY,
)
from .spectra import refractive_index_air


@dataclass(frozen=True)
class RayTable:
    altitude_m: np.ndarray
    impact_m: np.ndarray
    bending_rad: np.ndarray
    air_column_m2: np.ndarray
    ozone_column_m2: np.ndarray
    reference_wavelength_nm: float
    refraction: bool


def build_ray_table(
    atmosphere: Atmosphere,
    *,
    height_step_m: float = 100.0,
    quadrature_order: int = 48,
    reference_wavelength_nm: float = 550.0,
    refraction: bool = True,
) -> RayTable:
    """Trace all tangent heights once. Columns include both sides of the limb.

    Integrating in u=sqrt(r-r_t) removes the tangent-point singularity. Each
    US1976 layer is integrated separately, respecting derivative discontinuities.
    The tiny refractivity above the 120 km default top is neglected.
    """
    if not np.isfinite(height_step_m) or height_step_m <= 0:
        raise ValueError("height_step_m must be positive and finite")
    if quadrature_order < 8:
        raise ValueError("quadrature_order must be at least 8")
    if not np.isfinite(reference_wavelength_nm) or reference_wavelength_nm <= 0:
        raise ValueError("reference wavelength must be positive and finite")
    top = atmosphere.top_m
    dense_top = min(60_000.0, top)
    dense = np.linspace(0, dense_top, int(np.ceil(dense_top / height_step_m)) + 1)
    upper = np.linspace(dense_top, top, int(np.ceil((top - dense_top) / max(
        height_step_m, 500.0))) + 1)
    boundaries = np.asarray(LAYER_BOUNDARIES_M)
    boundaries = boundaries[(boundaries > 0) & (boundaries < top)]
    h = np.unique(np.concatenate((dense, upper, boundaries)))
    radius = EARTH_RADIUS_M + h
    refractivity = (float(refractive_index_air(reference_wavelength_nm)) - 1.0
                    if refraction else 0.0)
    coefficient = refractivity / REFERENCE_NUMBER_DENSITY
    if coefficient * float(atmosphere.number_density(top)) > 1e-9:
        raise ValueError("Atmosphere top refractivity must be <= 1e-9; raise top_m "
                         "because a refractive top interface is not modeled")
    n_tangent = 1 + coefficient * atmosphere.number_density(h)
    nprime_tangent = coefficient * atmosphere.density_gradient(h)
    qprime = n_tangent + radius * nprime_tangent
    if np.any(qprime <= 0):
        raise ValueError("Ducting/nonmonotone n*r is outside this ray-table model")
    b = n_tangent * radius
    if np.any(np.diff(b) <= 0):
        raise ValueError("Impact parameters must increase with tangent height")
    bending = np.zeros_like(h)
    air = np.zeros_like(h)
    ozone = np.zeros_like(h)
    nodes, weights = leggauss(quadrature_order)
    layers = np.concatenate(([0.0], boundaries, [top]))
    for lower, upper_bound in zip(layers[:-1], layers[1:]):
        mask = h < upper_bound
        ids = np.flatnonzero(mask)
        ht = h[ids, None]
        u0 = np.sqrt(np.maximum(lower - ht, 0))
        u1 = np.sqrt(upper_bound - ht)
        u = (u1 + u0) / 2 + (u1 - u0) / 2 * nodes
        delta = u * u
        altitude = ht + delta
        rr = EARTH_RADIUS_M + altitude
        density = atmosphere.number_density(altitude)
        n = 1 + coefficient * density
        q = n * rr
        # Near the turning point q-b loses digits in subtraction. Its local
        # analytic limit is q'(r_t)*(r-r_t), including the vacuum limit.
        q_minus_b = np.where(delta < 0.05, qprime[ids, None] * delta,
                             q - b[ids, None])
        if np.any(q_minus_b <= 0):
            raise ValueError("Nonpositive radial ray discriminant")
        denominator = np.sqrt(q_minus_b * (q + b[ids, None]))
        measure = 2 * u * (u1 - u0) / 2 * weights
        path = 2 * q / denominator * measure
        air[ids] += np.sum(density * path, axis=1)
        ozone[ids] += np.sum(atmosphere.ozone_density(altitude) * path, axis=1)
        nprime = coefficient * atmosphere.density_gradient(altitude)
        bending[ids] += np.sum(-2 * b[ids, None] * nprime / n
                               / denominator * measure, axis=1)
    return RayTable(h, b, bending, air, ozone, reference_wavelength_nm, refraction)
