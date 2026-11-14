# Validation products index

One entry per subdirectory: what it is, how to regenerate it, and which fields
are authoritative. Numbers live in the committed JSON products (zero hand
copies); this index only routes. Historical snapshots are preserved with their
original conditions — the molecular-only numerical snapshot below is
reproducible with `validation.check_convergence --aerosol none` (the current
convergence default is `opac`).

## Top-level numerical closure (L0)

- `physics.json` — Cartesian ODE rays, geometric and no-atmosphere limits,
  whole-plane energy budget. Regenerate:
  `uv run python -m validation.check_physics --out outputs/validation/physics.json`
- `aerosol.json` — Bessel straight-line reference, independent curved-path ODE,
  layered flux checks. Regenerate:
  `uv run python -m validation.check_aerosol --out outputs/validation/aerosol.json`
- `convergence.json` — height and spectral refinement errors, photopic and
  chromatic errors kept separate. Regenerate:
  `uv run python -m validation.check_convergence --out outputs/validation/convergence.json`

These validate the numerical implementation of the declared idealized model.
They do not validate real eclipse brightness, the true ozone profile, aerosol
identity, solar limb darkening, or visual adaptation.

## Observational, discrimination and calibration products

| Directory | Contents | Entry | Authoritative fields |
|---|---|---|---|
| `observations2011/` | June 2011 digitized-curve fixed-model comparison | `validation.compare_observations` | `comparison.json` per-band ratios |
| `round2/` | aerosol-round analysis, molecular/aerosol checks, convergence with explicit aerosol configuration | `validation.round2_analysis` | `analysis.json` selection and identification_scope |
| `dec2011/` | December 2011 exposure-averaged spectral comparison | `validation.compare_dec2011` | `comparison.json` per-readout records |
| `multiple_events/` | 1990 whole-disk photopic and 2007 extrapolated background comparisons | `validation.multiple_events` | `comparison.csv`/`comparison.json` per-event rows |
| `filter_response_2011/` | assumed-bandpass probe products | `validation.filter_response_2011` | `analysis.json` band analyses |
| `hst2019/` | HST/STIS starborne fixed-model comparison (raw-FITS extraction is the isolated exception) | `validation.compare_hst2019` | `comparison.json` |
| `extinction_gap_diagnosis/` | four-protocol kernel-height statistics, cross-event correlation, literature candidate discrimination | `validation.extinction_gap_diagnosis` | `discrimination.json`, `height_stats.json`, `correlation.json` |
| `stratosphere_background/` | four-aperture three-variant rerun with the production two-layer default, plus the preregistered check of the committed probe | `validation.extinction_gap_diagnosis --strat-rerun` | `predictions.json`, `preregistered_check.json` verdict |
| `atmosphere_inversion/` | per-event low-layer increment inversion (calibration lane) | `validation.atmosphere_inversion` | `inversion.json:cross_event`, `inversion.json:preregistered` |
| `humidity_axis/` | preregistered A/B re-inversion, R_RH verdict, share tables | `validation.humidity_axis_inversion` | `inversion_ab.json` |
| `prior_source_survey/` | source inventory SHA ledger, coverage matrix, humidity-prior transcription | produced by the prior-source-survey task tooling (no auto-regeneration entry; see [survey record](../../docs/records/PRIOR_SOURCE_SURVEY.md)) | `sources.json`, `coverage_matrix.json`, `humidity_axis_priors.json` |

Authoritative numbers are read from these products and never retyped into
prose; [MAINLINE](../../docs/MAINLINE.zh.md) cites fields mechanically.
