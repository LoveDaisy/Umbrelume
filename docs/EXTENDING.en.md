# How to add a physics increment (removable component)

A human-oriented checklist, distilled from the pattern every task in this
repository has repeated. The agent-facing hard rules live in
[../AGENTS.md](../AGENTS.md), the physical contract in [MODEL.md](MODEL.md),
and the methodology narrative in [MAINLINE.zh.md](MAINLINE.zh.md).

Chinese original: [EXTENDING.zh.md](EXTENDING.zh.md) — **authoritative version**.

## Before entering: discrimination first

- **Probes before physics**: use a same-protocol counterfactual probe to
  confirm the candidate can actually touch the current residual before
  implementing a component (see MAINLINE, "probes before physics").
- **Preregistration**: pin the prediction target, criteria and tolerances
  first, then let another execution path check them; the production
  implementation must never chase the expected result.
- **Three evidence identities**: independent prior / same-event calibration /
  independent prediction. New parameters must state their identity; calibrated
  values never enter defaults (standing ruling in AGENTS.md).

## Implementation checklist

1. **Independent module or function**: never merge spectral laws or action
   paths with existing physics terms; removability is a hard requirement.
2. **Exact zero-value restoration**: switching the new component off must
   restore the previous state bit-exactly (the aerosol pattern:
   `aerosol=False` keeps the counterfactual on the identical ray table;
   humidity f(50%) ≡ 1.0 exactly).
3. **Counterfactual tests**: on/off comparisons plus fail-loud validation
   (parameter domains, source anchoring, structural pairing invariants).
4. **L0 independent gates**: `check_physics` / `check_aerosol` / affected
   `check_convergence` compared against committed snapshots — a pure
   equivalence change must be BYTE-IDENTICAL; a physics change updates the
   snapshots explicitly with a stated reason.
5. **CLI wiring**: fail-loud argument validation; metadata records provenance
   (source, identity, hashes).
6. **Single-sourced numerical defaults**: promote defaults to module constants
   and reference the same constant from signatures and CLI metadata (the
   QUADRATURE_ORDER pattern, established in the 2026-11 consolidation).
7. **Documentation coupling**: MODEL.md (contract), physics/AEROSOL.md
   (formulas and sources), a ROADMAP row (milestone or entry condition),
   SOURCES.md (if new data).
8. **Snapshot discipline**: live comparison snapshots refresh with their input
   fingerprints (`figures/mainline` is mechanically enforced by
   `tests/test_mainline_figures`); frozen archives keep their original form.
9. **Bilingual sync**: README and ROADMAP bilingual pairs update together;
   running pytest triggers the link, field and pairing gates.

## Display and plotting rules

- The matplotlib base style is single-sourced in `eclipse/style.py`:
  `BASE_RC` (7 repo-shared keys) and `HOUSE_STYLE` (report-specific
  extension). New plotting scripts import
  `from eclipse.style import BASE_RC` (or `HOUSE_STYLE`) and layer any extra
  semantic keys locally on top of `BASE_RC` — **never redefine shared keys**.
- `PLOT["font"]` in `tools/render_mainline.py` duplicates
  `BASE_RC["font.family"]` by value (provenance metadata copy); a font change
  must update both.
- In-figure warnings repeat verbatim where they matter (e.g. the equal-Y strip
  "chromaticity diagnostic, not a brightness prediction") — one authoritative
  semantic, several reminders; that repetition is allowed.

## Gate quick reference

```bash
uv run pytest -q                      # ~60 s budget; the independent numerical
                                      # cross-checks are the last tests to cut
uv run python -m validation.check_physics --out outputs/validation/physics.json
uv run python -m validation.check_aerosol --out outputs/validation/aerosol.json
uv run python -m validation.check_convergence --out outputs/validation/convergence.json
```

Observational comparison scripts are validation products, not pytest tests;
after reruns, compare under the snapshot discipline.
