# Independent prior-source survey for the low-layer aerosol increment

Survey of 2026-11-11 (task-prior-source-survey). Purpose: collect
**independent literature/data sources** for the low-troposphere aerosol
increment and assess, per lane of the repository's evaluation discipline,
how strong the evidence is. **This survey collected evidence only; it made
no decision and no default changed.** The lane choice belongs to the owner
and the user.

## What was surveyed

The cross-event inversion (see
`docs/OBSERVATIONAL_VALIDATION.md`, section 跨事件所需大气反演) left a
standing question: the four events need low-layer AOD550 increments spread
over a factor R = 7.110 (pool 0.0259–0.1839, judged `partially_consistent`
under the preregistered criterion), and the humidity-variant axis of the
pinned OPAC prior had been `unavailable` because the OPAC 1998 data files
were unreachable (dead LMU URL, empty LRZ FTP mirror). The survey attacked
both, following the committed retrieval spec
(`validation/results/atmosphere_inversion/inversion.json` →
`retrieval_spec`).

### Headline results

1. **The OPAC 1998 humidity axis is recovered.** 26 of 38 component data
   files of the original 1998 distribution were pinned via GEISA/ESPRI
   copies recovered from the Internet Archive (per-file sha256 + capture
   timestamp in `sources.json`), after version identification against the
   pinned paper (61-wavelength grid 0.25–40 µm; dry-state parameters match
   Table 1c verbatim). 12 files (mid-RH points of WASO, three of SUSO/SSCM,
   and three dry-only components) are not archived anywhere reachable —
   trail in `search_log.md` W1.3–W1.7. The corrected premise: GEISA hosts
   the *classic* BAMS98 dataset alongside a separate OPAC 2015 entry; only
   the classic files were used.
2. **The humidity growth is quantified** (`humidity_growth.json`): at fixed
   number density, the full axis (dry → 99% RH) spans ×12.1–×18.6 across the
   four hygroscopic components; 50% → 99% spans ×7.4–×7.6, the same order as
   R = 7.110; moderate spans are much smaller (50 → 80%: ×1.5; 50 → 90%:
   ×2.1; 70 → 98%: ×3.9–4.0). Magnitude comparison only — no conclusion.
3. **A mid-latitude site distribution brackets the pool band.** GSFC (39.0°N)
   Level-2.0 AOD500, calendar 2011, n = 6677: min 0.0222, median 0.1163,
   P75 0.2056; 71.3% of retrievals fall inside the pool band and the site's
   P90/P10 dynamic range is 6.97 — the same order as R (comparison, not
   explanation).
4. **Clean-marine low end is documented independently** (Dubovik et al. 2002,
   Table 1, Lanai: τ(1020) 0.01–0.2, mean 0.04): the June inverted value
   0.0259 sits at the documented clean tail, below the OPAC boundary-layer
   family floor 0.038 but not below all independent observations.
5. **The 0–10 km exponential layer with H = 8 km has no independent support
   in what was reachable.** The CALIPSO H63 climatology (Winker et al. 2013,
   figure-read) shows extinction scale heights of ~1–3 km typical, at most
   ~5 km; no collected source documents 8 km. The classic standard-atmosphere
   reports (McClatchey 1972, Elterman 1968) are access-blocked; trail in
   `search_log.md` W3.1–W3.3.

## Where the evidence lives

| Product | Content |
|---|---|
| `validation/results/prior_source_survey/sources.json` | Pinned sources: per-file sha256, bytes, retrieval URL, printed location; `unavailable` rows carry search-log pointers |
| `validation/results/prior_source_survey/search_log.md` | Every search action, dead link and budget-gated stop (W1.x–W3.x) |
| `validation/results/prior_source_survey/humidity_growth.json` | f(RH) growth factors of the pinned component classes with per-value source pointers |
| `validation/results/prior_source_survey/coverage_matrix.json` | Retrieval-spec dimensions × sources, verdicts `covered`/`partial`/`unavailable` only (coverage ≠ justification) |
| `validation/results/prior_source_survey/three_paths.json` | The three lanes (fixed higher prior / documented range / keep current), each with evidence, gaps, and what it would still need — no recommendation column |

Downloaded bytes are cached under `outputs/prior_source_survey/raw/`
(git-ignored); every pinned file can be re-downloaded and verified through
its recorded sha256 and retrieval URL.

## Evaluation discipline

The three lanes and their vocabulary are unchanged by this survey:
fixed prior / documented range / calibrated lane stay distinct
(see `eclipse/data/SOURCES.md`, calibration-lane section, and the
discipline paragraph at the end of `docs/physics/AEROSOL.md`). Nothing in
the collected evidence justifies a specific default value; a site or type
climatology covering the pool band is not a justification of any member of
it. No default parameter was touched; `eclipse/` production code is
unchanged by the survey.
