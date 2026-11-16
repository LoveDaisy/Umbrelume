# Browser eclipse demo

[中文权威版本](WEB_DEMO.zh.md) · English translation

The browser is an **interpretation outlet** of the production Python model, not
another atmosphere or refraction solver. It displays illumination on a neutral
lunar receiving screen, not an event prediction, accuracy certification, or
dark-adapted appearance model.

## Controls

- Drag the large lunar screen or the small shadow navigator to move the lunar
  center; arrow keys work after focusing, with Shift for larger steps. Dashed
  navigator circles mark geometric umbra/penumbra boundaries, the central cross
  marks the shadow axis, and distances are in km.
- Ozone changes absorption only. AOD is the FT vertical reference value at 50%
  RH. Humidity scales only the FT amplitude; the two-layer prior's independent
  stratospheric background remains 0.005.
- Exact zero AOD disables **all** aerosol and restores the molecular-only grid
  family. RH is disabled at zero, and its previous value returns when AOD is
  positive. `none` disables both AOD and RH; ozone stays independently adjustable.
- `calibrated:` entries come from the production registry and are explicitly
  labeled fitted, not physical priors. They reject AOD/RH overrides and never
  replace the default.
- Exposure affects display only, not the SI spectral readout. One scene-wide
  solar white is Bradford-adapted to D65. Radial interpolation is in linear RGB,
  before global exposure, negative-channel clipping, common-channel highlight
  compression, and sRGB encoding.

The UI uses fixed geometry/standard molecular atmosphere, 380–780 nm in 5 nm
steps, ozone 0–600 DU, FT reference AOD 0–0.2 (minimum positive step 0.001), RH
50–99%, and exposure −4 to +20 EV. The complete lunar disk must stay inside the
exported 0–9000 km radius domain; out-of-domain values are never clamped to the
last color row. The center moves in two dimensions, up to 9000 km minus the lunar
radius from the shadow axis.

## Data and implementation boundary

The [exporter](../tools/export_web_demo.py) calls production `build_ray_table`,
`receiver_quadrature`, spectral/color APIs, and the calibrated loader. Each grid
family contains:

1. Exactly deduplicated quadrature `b` nodes, CSR radial indices, and conserved
   nonnegative annular weights, retaining every landing branch.
2. Molecular columns interpolated in `b`, ozone columns per DU, and per-layer tau
   per vertical AOD, each paired with its own spectral factor. The analytic
   vacuum/direct term prevents the uneclipsed disk from incorrectly going dark.
3. The uneclipsed SI solar spectrum, CIE integration weights, and linear color
   matrix with one shared solar white.

The browser performs **sum tau → exp → area contraction → SI spectra/CIE → shared
display**. It interpolates production-exported RH growth knots rather than
copying humidity physics, cross-section formulae, or the refraction integrator.
AOD/DU changes do not change W within a family; different active layer supports
can have different nodes and are not claimed to share one grid. Exact zero uses
`none`; ozone-off does not change refraction.

`manifest.json` records schema/software versions, input SHA256 hashes, the
canonical generation command, units/axis order, dtype/little-endian layout,
shape/offset/byteLength, gzip checksums, grid/preset/layer pairings, numerical
parameters, geometry, defaults, and supported domain. Families are shared only
when the entire array payload is byte-identical. Explicit gzip decoding does not
depend on server compression settings. Missing assets, unknown schemas, bad
hashes, and malformed shapes produce visible errors and stop computation; there
is no stale-field fallback.

The loading state remains visible until every preset is preloaded, decoded, and
validated. The measured size allows all exactly shared families to be decoded
once, then compressed copies released, avoiding repeated decompression on
preset changes. The Worker queue holds at most one executing request and one
latest pending request; old results cannot overwrite new parameters. Canvas
backing size stays 560×560 instead of growing without bound with mobile DPR.

## Local reproduction and mechanical gate

Run from the intended worktree root, with absent or empty output directories:

```bash
uv sync --locked --group web
uv run --group web python -m playwright install chromium
uv run python tools/export_web_demo.py --out outputs/web-demo-a
uv run python tools/export_web_demo.py --out outputs/web-demo-b
diff -qr outputs/web-demo-a outputs/web-demo-b
uv run --group web python -m validation.check_web_demo --site outputs/web-demo-a --out outputs/validation/web-demo.json
uv run pytest -q
```

For local browsing, run `uv run python -m http.server 8000 --directory outputs`
and open `http://localhost:8000/web-demo-a/`. Modules/Workers cannot be loaded via
`file://`. There is no npm, backend, Pyodide, or mandatory GPU path. Linux CI adds
`--with-deps` when installing Chromium.

The [independent gate](../validation/check_web_demo.py) builds fresh references
from current production APIs, serves a project subpath, and executes the actual
app→Worker→Canvas path, not just a helper unused by the UI. Coverage includes
parameter endpoints/interior values/seeded random samples, every radial node and
wavelength, real CLI CSV/NPZ files, every radial cell midpoint and geometric
boundary neighborhoods, pixel masks, rapidly interleaved input, invalid
parameters, and fault injection. Where a geometric penumbra center would be
outside the legal center domain, pixel coverage places the outer lunar edge on
that boundary instead of clamping samples.

Fixed outlet tolerances apply pointwise, not to averages:

| Quantity | Criterion |
|---|---|
| Same-node T, relativeY, XYZ / solar reference Y | atol=2e-12, rtol=2e-8 |
| Same-node unclipped linear RGB | atol=2e-11, rtol=2e-8 |
| SI uneclipsed reference spectrum | atol=1e-14 W m⁻² nm⁻¹, rtol=2e-12 |
| Off-node relativeY | atol=5e-9, rtol=0.005 |
| Off-node display RGB at EV=−4/0/8/12/20 | absolute channel error≤3/255 |
| Canvas against the same display convention | ≤1 eight-bit code per channel; separate mask check |

Reports retain maximum errors/worst coordinates, parameters, fingerprints, all
timing samples, buffer inventory, and desktop/narrow-screen screenshots.
Deliberately changing W, aerosol spectral factors, or the color matrix must fail
the numerical gate even after updating asset hashes. Missing browsers/assets or
zero samples cause nonzero exit, never a skipped green result.

After preloading, 30 warm-up and 200 prescribed parameter/preset changes are
measured from application input to complete field/CIE/560×560 drawing. Both p95
and every sample must be ≤100 ms; drag/exposure p95 must be ≤33 ms. Published data
for every supported preset must total ≤32 MiB. The conservative peak inventory
of explicit numerical/image/transfer buffers must be ≤192 MiB, not confused with
whole-browser RSS. A mobile-width screenshot proves layout, not real-device
performance. These are reproduction/responsiveness gates, not physical accuracy
certification.

## CI and publishing

The [web workflow](../.github/workflows/web-demo.yml) builds on every PR and main
push, without path filters that could miss core changes. It regenerates assets
and fresh references, then runs Python, physics/aerosol/convergence checks, and
the complete browser gate. Failure prevents publication. Normal pytest does not
require an installed browser.

Only `v*` tags whose commits belong to main history, or explicit main-branch
workflow_dispatch runs, deploy. A run retains the exact validated site, report,
source SHA, and manifest hash. Deployment downloads and verifies that same
artifact without rebuilding. Build permissions are `contents:read`; deployment
adds `pages:write`/`id-token:write`, the `github-pages` environment, and deployment
concurrency. An authorized owner must trigger the first deployment after merging
the task into main. This implementation does not merge, push, or create tags.
Additional Pages/environment permissions are escalated, never bypassed.

**Actual Actions/deployment records determine acceptance and online status; this
document is not evidence that the site is live.** Real-device interaction and
performance, plus human aesthetic acceptance, remain separate from local green
checks.

## Sources and limitations

The [per-source provenance and terms](../eclipse/data/SOURCES.md), including
partially verified statuses, are copied into the site. MIT source-code licensing
does not make every third-party dataset MIT. Only production numerical derivative
tables are published, not source papers, observational downloads, or a wheel.
Nothing is downloaded implicitly.

The page retains three boundaries: the browser is an interpretation outlet and
the repository CI governs the core solver; this demonstrates phenomena, not
certification/prediction; the model excludes effects such as in-scattered light.
The neutral lunar screen excludes lunar texture/BRDF, observer atmosphere, and
scotopic vision; highlight-compressed pixels are not irradiance measurements.
The authoritative physical contract and outlet registry are in [MODEL](MODEL.md).
