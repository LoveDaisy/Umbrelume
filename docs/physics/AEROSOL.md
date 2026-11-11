# Aerosol extinction: free troposphere and stratospheric background

The aerosol increments add extinction along the already refracted gas paths.
Extinction includes absorption and scattering out of the direct beam. It does
not include aerosol refractivity or light scattered into the receiver.

## Fixed independent priors

The default composes two published OPAC layers (Hess, Koepke and Schult, 1998,
*Optical Properties of Aerosols and Clouds: The Software Package OPAC*,
[primary paper](https://lweb.cfa.harvard.edu/HITRAN/HITRAN2012/Aerosols/papers/hess_opac_1998.pdf),
[DOI](https://doi.org/10.1175/1520-0477(1998)079%3C0831:OPOAAC%3E2.0.CO;2)):

| Layer | Quantity | Fixed value | Source |
|---|---|---:|---|
| Free troposphere | Layer bounds | 2–12 km | Table 3 |
| Free troposphere | Vertical optical depth at 550 nm | 0.013 | Table 3 |
| Free troposphere | Exponential scale height | 8 km | Table 5, free-troposphere row |
| Free troposphere | Spectral exponent, 500–800 nm | 1.58 | Table 3 |
| Free troposphere | Spectral exponent, 350–500 nm | 1.21 | Table 3 |
| Stratospheric background | Layer bounds | 12–35 km | Table 3/5; "always assumed between 12 and 35 km" (Sec. 5b) |
| Stratospheric background | Vertical optical depth at 550 nm | 0.005 | Table 3 |
| Stratospheric background | Scale height | 99 km | Table 5; homogeneous-layer representation |

Every numeric value is anchored to the pinned rows of
`validation/results/extinction_gap_diagnosis/prior_parameters.json`
(`free_troposphere_reference`, `stratosphere_background`) under the recorded
OPAC PDF sha256; pytest asserts factory parameter == row value, so the source
link is mechanical, not documentary.

These are published background-model/climatology parameters. They are
**not measurements of any eclipse atmosphere** and are not inferred from
eclipse brightness curves. The free-troposphere optical values refer to the
paper's 50% relative-humidity case, which is also the model default (the
`--rh` reference state below); the humidified-scaling prior was unavailable
when this file was first written and is now supplied by the recovered OPAC
component data (see the humidity axis section).
The OPAC mixing-layer and transported-dust components are separate; none is
added here. The wide low-layer increment (antarctic-type proxy, 0–10 km) was
discriminated **against** as a default component: it is a menu-selection
effect that strongly worsens the shape axis (see the standing gap below).

## Humidity axis (free troposphere)

The free-troposphere layer carries an optional, pinned humidity growth
(`eclipse/aerosol.py` `ft_growth_factor`, task-humidity-axis). The OPAC 1998
WASO component's 550 nm extinction per particle is pinned at three
relative-humidity knots from the recovered component files
(`validation/results/prior_source_survey/humidity_growth.json`, sha256 per
knot): ext550 = 6.366e-6 (50%), 7.913e-6 (70%), 4.740e-5 (99%) km⁻¹ per
(particle cm⁻³). The growth composes into a column-AOD factor

```text
f(RH) = (E_waso(RH) + q_rest) / 0.013,
E_waso(RH) = 0.6 * 730 * e_waso(RH) * 8*(exp(-2/8) - exp(-12/8)) km,
q_rest    = 0.013 - E_waso(50%) = 6.0495e-4 (4.65% of the reference),
```

with the paper's FT composition (Sec 5a, printed p. 842: WASO number mixing
ratio 0.6 of a 730 particles cm⁻³ sea-level density) and vertical geometry
(Table 5 FT row H = 8 km; Table 3 FT row 2–12 km). The non-growing
soot+insoluble remainder is pinned by subtraction (`q_rest`): those
components are hydrophobic (pinned only as 00 files), so their extinction
does not grow with RH. f(70%) = 1.2317 and f(99%) = 7.1459. Every value,
quote and ruling is transcribed in
`validation/results/prior_source_survey/humidity_axis_priors.json` (sha256-
anchored; pytest recomputes the chain).

Semantics and scope:

- **Reference state is bit-exact**: f(50%) is exactly 1.0, so the default
  output (and every committed snapshot) is unchanged; the growth is a
  column-AOD rescaling only — the layer's vertical shape and its broken
  power law never change with RH.
- **Domain and interpolation**: RH ∈ [0.5, 0.99], fail loud outside; between
  knots the pinned ext550 values are interpolated linearly in RH (the
  preregistered semantics; a log-RH alternative is recorded as a sensitivity
  bracket in the A/B re-inversion product). The sub-reference direction is
  deliberately not implemented — no discriminator consumes it.
- **Interface**: `--rh` (default 0.5) acts on the free-troposphere layer of
  `--aerosol opac` / `opac-ft` only, in composition with `--aod550` (the
  override stays the 50% reference base; the effective FT AOD is
  `base × f(RH)`, composed once at the factory). It is rejected with
  `calibrated:<name>` and `none`. Metadata records `rh`, the growth factor,
  the reference base and the knot sha256s in the FT layer's source record.
- **Layer scope rulings**: the stratospheric background carries no humidity
  parameter (Table 3 declares RH states only for boundary aerosols and the
  FT; the pinned SUSO data has no knots below 70%, so low-RH growth would
  need unpinned fabrication). Calibration-lane increment layers stay pinned
  at the 50% reference state so a share measurement attributes humidity to
  the prior, never to the calibrated increment.
- **Limitation**: this is a 550 nm scalar column-AOD axis. Real humidification
  also changes particle size and thus flattens the spectral slope; that
  wavelength-resolved response is not simulated (it would entangle the pinned
  power-law shape) and is recorded as the shape-axis cost in the A/B
  re-inversion's sensitivity appendix.

## Definition and units

For lower and upper layer boundaries `h0`, `h1`, scale height H and requested
vertical optical depth tau550,

```text
beta550(h) = tau550 exp[-(h-h0)/H] / {H [1-exp(-(h1-h0)/H)]}, h0 <= h <= h1
           = 0, otherwise.
```

Public altitudes are geometric metres and beta550 is m^-1. Its vertical
integral is exactly tau550. At 11 km the free-troposphere layer is
approximately 7.39e-4 km^-1 (7.39e-7 m^-1). A shifted exponential and an
`expm1` denominator preserve normalization for thin or high-altitude custom
layers. The 99 km scale height over the 23 km stratospheric thickness is the
OPAC representation of a nearly homogeneous layer (exp(-23/99) ~ 0.79).

The spectral factor is a continuous broken power law, normalized at 550 nm,
with a break at 500 nm, applied **per layer**: the free troposphere uses
1.58/1.21 and the stratospheric background uses 1.14/0.74. The two stated
slopes apply on their respective published intervals; the 380–780 nm baseline
lies entirely within those intervals. The reference does not depend on the
wavelength query grid. Custom wavelengths outside 350–800 nm use the same
power laws as explicit extrapolations and should not be presented as OPAC
tabulated optical data.

Each layer's path optical depth is
`tau_l(lambda) = integral beta_l(h(s)) ds * f_l(lambda)` on exactly the same
curved gas ray. The transmission hook sums the layers:
`tau_a = sum_l tau_l(lambda)`. Both limb halves must be included. Path
integration must split at **2, 12 and 35 km** as well as the gas layer
boundaries, because both prescribed profiles have jumps there (the 12 km edge
is shared: free-troposphere top and stratosphere bottom). No
single-scattering albedo or delta-M factor multiplies this direct-beam
extinction.

## Layer composition and interface semantics

`eclipse.aerosol.AerosolColumn` is the ordered composition type: a frozen
tuple of `AerosolProfile` layers whose total vertical AOD is exactly the sum
of the layers' AOD550. Extinction sums per layer for arbitrary support
layouts — overlapping supports add their optical depths and the sum invariant
holds for any layout (asserted in tests including an overlapping-support
case); the production default's two layers are disjoint, while the
calibration-lane presets (`--aerosol calibrated:<name>`) stack a third
increment layer that may overlap the default altitude bands. Each layer
integrates on its own support with the per-tangent quadrature splitting at
that layer's own boundaries (`eclipse/rays.py`), so overlapping supports
compose without interaction. Layers never merge into one effective profile —
the layers' spectral exponents differ, so the only legal expansion is the
per-layer (tau table, spectral factor) pair; the ray table stores one 550 nm
tau table per layer together with that layer's own profile.

Interface mapping (names and semantics stable across the change):

- `--aerosol opac` (new default): the two-layer `AerosolColumn.opac_default()`.
- `--aerosol opac-ft`: free-troposphere layer only — the historical preset,
  unchanged in meaning and reproducible as before.
- `--aerosol none`: molecular-only counterfactual, unchanged.
- `--aod550 X` overrides **only the free-troposphere layer's** AOD; the
  stratospheric layer stays at its published 0.005. `--aod550 0` disables the
  **whole** aerosol component (molecular-only, bit-exact with `--aerosol
  none`); a zero-AOD layer inside an explicit Python-API `AerosolColumn` is
  legal and normalizes to its active subset. In that `--aod550 0` case the
  metadata self-describes asymmetrically: `aerosol` is null and
  `aerosol_preset` keeps the requested choice string, so the actual medium is
  identified by `aerosol_aod_source == "Disabled"` — read all three fields,
  never `aerosol_preset` alone.
- `--rh R` applies the pinned free-troposphere humidity growth (see the
  humidity axis section above); the default 0.5 is the published reference
  state and leaves every output bit-exact unchanged.

Saved ray tables (`rays.npz`) carry the new key `aerosol_layer_taus` (a 2-D
`(n_layers, n_tangent)` stack) whenever the run has an aerosol component, and
keep the historical zero-array key `aerosol_tau550` only for the molecular-only
runs. Metadata records both layers' parameters with their
`prior_parameters.json` row references and the OPAC sha256; the
free-troposphere layer's source record additionally carries the humidity-axis
provenance (`rh`, the growth factor, the 50% reference base and the knot
sha256s).

## Evidence and interpretation

The first observational comparison found that the gas-only model was
progressively too bright inward. A free-tropospheric component acts strongly
on lower grazing paths while leaving higher paths less affected. A typical
surface-only profile with 1.5 km scale height and AOD550 0.05 has only about
2.2e-5 km^-1 extinction at 11 km; simply importing a boundary-layer preset
would therefore miss an important part of the relevant altitude range.

The model's zero extinction above 12 km was a known structural gap: OPAC
itself fixes a 12–35 km background stratospheric layer independent of volcanic
activity. The four-aperture discrimination (task-extinction-gap-diagnosis)
identified this stratospheric-background prior as the only candidate among
twelve that improves the shape axis (June 503/677 ratio 1.1102 → about 1.0231)
with tension-level brightness improvement — and the component entered the
production default as a from-coarse-to-fine completion of that known zero
region, with the discrimination table's numbers serving as the preregistered
prediction (checked item by item in
`validation/results/stratosphere_background/preregistered_check.json`).

**Standing gap (owner ruling 2026-11-11):** the brightness residual is NOT
closed by this component. All four events' brightness closures remain in the
tension band (+0.19 to +0.44 of the committed log-residual). The wide low-layer
increment (antarctic-type proxy, 0–10 km, AOD550 about 0.066) closes the
brightness axis on three of four events but strongly worsens the shape axis
and carries a menu-selection effect; it stays in the calibration-diagnosis
lane and out of the default. Closing the brightness residual therefore needs
a separately discriminated direction (event-dependent profiles, cloud, source
structure), not a widening of this component.

The default OPAC components are independently chosen forward-model priors.
Any subsequent AOD calibration against eclipse curves must be separately
labeled, retain the fixed default comparison and distinguish calibration
bands from held-out bands. Neither success against a calibrated curve nor
the same event's aerosol inversion constitutes independent prediction.
Residual errors can still reflect the solar intensity distribution,
atmospheric temperature/ozone profiles, cloud extinction and missing
instrumental bandpass information. These components do not establish a
unique cause or validate full visible color and the turquoise fringe.

An independent prior-source survey of the low-layer increment (2026-11,
task-prior-source-survey) collected literature evidence for the three
evaluation lanes without changing any default; see
[docs/PRIOR_SOURCE_SURVEY.md](../PRIOR_SOURCE_SURVEY.md) and the products
under `validation/results/prior_source_survey/`.
