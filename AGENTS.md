# Lunar Eclipse Color

First-version spectral lunar-eclipse illumination simulator. Python 3.12;
dependencies and environment are managed with `uv sync` and `uv.lock`.

- New implementation lives on `main`; `legacy/*` branches preserve earlier work.
  Do not modify or remove legacy branches unless the user requests it.
- Use topic branches and isolated git worktrees for implementation. Keep tool cwd
  and imports explicitly anchored to the intended worktree. No scratchpad workflow
  is required for this project.
- Run `uv run pytest -q`. For changes to physics or integration, also run
  `uv run python -m validation.check_physics --out outputs/validation/physics.json`
  and the affected resolution checks in `validation/check_convergence.py`.
- `docs/MODEL.md` owns the physical contract and limitations. `eclipse/data/SOURCES.md`
  owns provenance. Keep these and generated baseline metadata consistent with changes.
- Preserve SI units and spectral arrays before display. Include both limb halves,
  all landing-map branches and conserved annular energy weights. Ozone counterfactuals
  change absorption only. Do not tune physical parameters solely to obtain a desired color.
- `figures/baseline/` and `validation/results/` hold reviewable first-version results;
  disposable simulations belong in ignored `outputs/`.
