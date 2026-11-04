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
