# Free-tropospheric aerosol extinction

The first aerosol increment adds extinction along the already refracted gas
paths. Extinction includes absorption and scattering out of the direct beam.
It does not include aerosol refractivity or light scattered into the receiver.

## Fixed independent prior

The default is the free-tropospheric component of Hess, Koepke and Schult
(1998), *Optical Properties of Aerosols and Clouds: The Software Package OPAC*,
[primary paper](https://lweb.cfa.harvard.edu/HITRAN/HITRAN2012/Aerosols/papers/hess_opac_1998.pdf),
[DOI](https://doi.org/10.1175/1520-0477(1998)079%3C0831:OPOAAC%3E2.0.CO;2).

| Quantity | Fixed value | Source |
|---|---:|---|
| Layer bounds | 2–12 km | Table 3 |
| Vertical optical depth at 550 nm | 0.013 | Table 3 |
| Exponential scale height | 8 km | Table 5, free-troposphere row |
| Spectral exponent, 500–800 nm | 1.58 | Table 3 |
| Spectral exponent, 350–500 nm | 1.21 | Table 3 |

These are published background-model/climatology parameters. They are
**not measurements of the June 2011 eclipse atmosphere** and are not inferred
from that eclipse's brightness curve. The OPAC mixing-layer, transported-dust
and stratospheric components are separate; none is silently added here.

## Definition and units

For lower and upper layer boundaries `h0`, `h1`, scale height H and requested
vertical optical depth tau550,

```text
beta550(h) = tau550 exp[-(h-h0)/H] / {H [1-exp(-(h1-h0)/H)]}, h0 <= h <= h1
           = 0, otherwise.
```

Public altitudes are geometric metres and beta550 is m^-1. Its vertical
integral is exactly tau550. At 11 km the default is approximately
7.39e-4 km^-1 (7.39e-7 m^-1). A shifted exponential and an `expm1` denominator
preserve normalization for thin or high-altitude custom layers.

The spectral factor is a continuous broken power law, normalized at 550 nm,
with a break at 500 nm. The two stated slopes apply on their respective
published intervals; the 380–780 nm baseline lies entirely within those
intervals. The reference does not depend on the wavelength query grid.
Custom wavelengths outside 350–800 nm use the same power laws as explicit
extrapolations and should not be presented as OPAC tabulated optical data.

The path optical depth is `tau_a(lambda) = integral beta550(h(s)) ds * f(lambda)`.
Both limb halves must be included. Path integration must split at **2 and
12 km** as well as the gas layer boundaries, because the prescribed profile
has jumps there. No single-scattering albedo or delta-M factor multiplies
this direct-beam extinction.

## Evidence and interpretation

The first observational comparison found that the gas-only model was
progressively too bright inward. A free-tropospheric component acts strongly
on lower grazing paths while leaving higher paths less affected. A typical
surface-only profile with 1.5 km scale height and AOD550 0.05 has only about
2.2e-5 km^-1 extinction at 11 km; simply importing a boundary-layer preset
would therefore miss an important part of the relevant altitude range.

The default OPAC component is an independently chosen forward-model prior.
Any subsequent AOD calibration against the 2011 curves must be separately
labeled, retain the fixed default comparison and distinguish calibration
bands from held-out bands. Neither success against a calibrated curve nor
the same event's aerosol inversion constitutes independent prediction.
Residual errors can still reflect the solar intensity distribution,
atmospheric temperature/ozone profiles, cloud extinction and missing
instrumental bandpass information. This component does not establish a
unique cause or validate full visible color and the turquoise fringe.
