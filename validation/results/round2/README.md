# Round 2 evidence snapshot

Model v0.2 adds only free-tropospheric aerosol extinction to production. The
default is the independent OPAC prior (AOD550=0.013), not the calibrated AOD.

| File | Scope |
|---|---|
| `analysis.json` / `comparison.csv` | Fixed-prior predictions, molecular counterfactual, 677 nm-only AOD calibration, withheld503nm check, altitude attribution, source/density diagnostics and code hashes |
| `physics.json` | Existing independent molecular-ray and no-aerosol energy checks, retained unchanged in scope |
| `aerosol.json` | Exact straight-ray Bessel check, independent curved Cartesian path check and energy closure with aerosol |
| `convergence.json` | Current OPAC-prior height and wavelength refinements, with explicit model/profile metadata |

The complete test suite passes 98 tests at this snapshot. Numerical checks and
source-mixture checks concern implementation, not whether this atmosphere is
the actual2011atmosphere. The observational readouts have unknown filter
responses and uncertainties and are samples of smoothed published curves.
Their RMS log ratios are descriptive, not independent-sample statistics.

```bash
uv sync --locked
uv run pytest -q
uv run python -m validation.round2_analysis --out outputs/round2
uv run python -m validation.check_physics --out outputs/validation/physics-round2.json
uv run python -m validation.check_aerosol --out outputs/validation/aerosol-round2.json
uv run python -m validation.check_convergence --out outputs/validation/convergence-round2.json
```

The report files were copied unchanged from these outputs. See
[ROUND2.md](../../../docs/ROUND2.md) for interpretation and
[ROADMAP.md](../../../docs/ROADMAP.md) for deferred directions. Original v1
snapshots remain outside this folder with their original conditions.
