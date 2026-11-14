# Numerical validation snapshot

Generated from the first-version implementation, using the locked environment:

```bash
uv run pytest -q
uv run python -m validation.check_physics --out validation/results/physics.json
uv run python -m validation.check_convergence --out validation/results/convergence.json
```

The test suite at this snapshot passes 55 tests. Independent Cartesian ODE rays
at 5, 11, 11.019, 15, 25 and 40 km agree with radial quadrature to a maximum
relative difference of about 2.5e-8 across bending and molecular columns.
The receiver-plane/pupil missing-energy checks close to about 1.2e-9 for
extinction and 3.5e-9 for a transparent refracting atmosphere.

Ray-table height spacing 100→50 m changes the sampled relative photopic
irradiance by at most 0.0791%; 50→25 m changes it by 0.0233%. Spectral spacing
5→2.5 nm changes it by 0.401%. The maximum equal-Y linear RGB channel change
for the spectral refinement is 0.0176 (before gamut clipping). Values and
acceptance criteria are recorded in the JSON files. Above 60 km, the height
table uses 500 m spacing in all these runs.

These checks validate the numerical implementation and selected resolutions
of the declared idealized model. They do not validate real eclipse brightness,
the true ozone profile, atmospheric aerosols, solar limb darkening, gray
refraction, small-angle propagation, lunar reflectance or visual adaptation.

## Observational reference

The [June 2011 comparison](observations2011/comparison.json) is separate from
the numerical pass/fail checks above. It compares the unchanged model with
explicitly digitized published 503/677 nm observed curves, without fitting.
It shows substantial inner-umbra overprediction. Unknown bandpasses and
pointwise observational uncertainties preclude a formal statistical verdict.
See [the analysis and source inventory](../../docs/records/OBSERVATIONAL_VALIDATION.md).

## Round 2: aerosol extinction

`round2/` contains the new analysis, molecular-ray checks, independent aerosol
checks and convergence results with explicit aerosol configuration. The default
uses the fixed OPAC free-tropospheric prior; same-event 677 nm AOD calibration
and withheld 503 nm checks are reported separately. Numerical closure is not
observational/profile validation. Read [ROUND2.md](../../docs/records/ROUND2.md) for the
selection evidence and remaining mismatch.

The original numerical snapshot above is molecular-only. Reproduce that medium
with `validation.check_convergence --aerosol none`; current convergence defaults
to `opac`. Historical JSON files are retained with their original conditions.

## Stratospheric-background rerun

`stratosphere_background/` holds the four-aperture three-variant rerun with the
production two-layer default (`molecular` / `fixed_opac` / `fixed_opac_strat`)
and the preregistered comparison of the committed `stratosphere_background`
probe deltas against that rerun. Generate with:

```bash
uv run python -m validation.extinction_gap_diagnosis --strat-rerun \
    --out validation/results/stratosphere_background
```

`predictions.json` records per-readout predictions for all three variants with
the committed variants' moment-reproduction gaps; `preregistered_check.json`
compares each preregistered item (probe value read from
`extinction_gap_diagnosis/discrimination.json`, never retyped) against the
production rerun under `max(1e-6 x |probe|, 1e-9)`. The brightness axis stays
in the tension band by preregistered expectation — the residual is not closed
by this component. The committed `extinction_gap_diagnosis/` products remain
byte-identical.
