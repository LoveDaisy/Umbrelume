# Lunar Eclipse Color

Spectral lunar-eclipse illumination simulator. Python 3.12;
dependencies and environment are managed with `uv sync` and `uv.lock`.

- New implementation lives on `main`; `legacy/*` branches preserve earlier work.
  Do not modify or remove legacy branches unless the user requests it.
- Use topic branches and isolated git worktrees for implementation. Keep tool cwd
  and imports explicitly anchored to the intended worktree. Task management runs
  through the ignored `scratchpad/` directory (a local, git-ignored working area:
  it does not exist in a fresh clone and is not part of the repository contract).
- Run `uv run pytest -q`. For changes to physics or integration, also run
  `uv run python -m validation.check_physics --out outputs/validation/physics.json`
  and `uv run python -m validation.check_aerosol --out outputs/validation/aerosol.json`;
  also run affected resolution checks in `validation/check_convergence.py`.
- `docs/MODEL.md` owns the physical contract and limitations. `eclipse/data/SOURCES.md`
  owns provenance. Keep these and generated baseline metadata consistent with changes.
- Preserve SI units and spectral arrays before display. Include both limb halves,
  all landing-map branches and conserved annular energy weights. Ozone counterfactuals
  change absorption only. Do not tune physical parameters solely to obtain a desired color.
- `figures/` and `validation/results/` hold reviewable results and historical snapshots;
  disposable simulations belong in ignored `outputs/`.
- `docs/ROADMAP.zh.md` owns milestones and deferred directions; `docs/records/ROUND2.md`
  records the aerosol selection and its evidence. Preserve historical snapshots.
  Live comparison snapshots (e.g. `figures/mainline/`, recent `validation/results/`
  entries) are refreshed when their inputs change; frozen archives keep their
  original form. Fixed literature priors, fitted same-event parameters and independent
  prediction are distinct. The default OPAC AOD must not silently become a calibrated value.
- Docs are layered: reader-facing material sits in `docs/` (index at `docs/README.md`);
  internal development records live in `docs/records/`. Narrative docs are bilingual
  (`.zh.md`/`.en.md` pairs, per-pair authority declared in `docs/README.md`, enforced
  by the pairing gate in `tests/test_mainline_docs.py`); contracts and records stay
  single-language.
- Coarse-to-fine organization: `eclipse/` carries only the current default physics,
  layered as spectra → atmosphere → aerosol → rays → transport → color → report → cli.
  Observational comparisons stay in `validation/` (run on demand, not in the pytest
  gate); data acquisition stays in `tools/`.
- Each physics increment enters as a removable component: its own module, exact
  zero-value restoration of the previous state (the aerosol pattern), and
  counterfactual tests. Before adding a third counterfactual switch to
  `transmission_profile`, review composing attenuation terms into a stack.
- Keep the whole pytest run within ~60 s; the independent numerical cross-checks
  are the last tests to cut. Observational comparisons are validation scripts,
  not pytest tests.
