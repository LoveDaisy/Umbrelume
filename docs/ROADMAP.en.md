# The coarse-to-fine project roadmap

> Chinese original: [ROADMAP.zh.md](ROADMAP.zh.md) — **authoritative version**; this
> translation may lag behind it.

## Goal and current boundary

The goal is to compute the spectral irradiance at the Moon's position from the
solar spectrum and its propagation through Earth's atmosphere, then simulate
the colors of the eclipsed lunar surface, including the cold band at the outer
umbral edge. Each round first uses observations or independent computation to
locate one major gap, then adds one piece of physics and compares spectra,
brightness and spatial distribution under identical conditions; tuning to
produce a desired color is never used as an improvement criterion.

The current main output is still the spectral irradiance on the plane through
the Moon, perpendicular to the Sun–Earth axis. Lunar images are
neutral-screen renderings and do not yet include real lunar soil or an
observer. The existing two-band 2011 observed curves do not cover the
turquoise band itself, so no verified claim of full visible colors, band width
or saturation can be made from them.

[MODEL.md](MODEL.md) owns the physical and numerical contract;
[OBSERVATIONAL_VALIDATION.md](records/OBSERVATIONAL_VALIDATION.md) owns the
first-round comparison; [MULTI_EVENT_OBSERVATIONS.md](records/MULTI_EVENT_OBSERVATIONS.md)
records multi-event data, usage conditions and new comparisons. This file
records where the project stands and the entry conditions for future directions.

## Implemented milestones

| Milestone | Output and evidence | Status and boundary |
|---|---|---|
| M0: a new starting point | Old branches preserved as `legacy/*`; the new `main` starts from an independent empty initial commit | Done. Old implementation kept for explicit reference needs |
| M1: end-to-end refracted-transmission model | Sun-side precomputed bent light paths, air/ozone columns; finite solar disk and distance; conserved annular-incidence integration; spectra to color | Done. Includes Rayleigh extinction and ozone absorption; refracted paths use a 550 nm reference; solar radiance uniform |
| M1: independent numerical verification | Cartesian ray ODEs, whole-plane energy budget, no-atmosphere geometric limit, resolution checks | Done. Proves the declared model's implementation and numerical closure; does not replace observational validation |
| M2: first public observational reference | Digitized 503/677 nm observed curves of the western umbra, 2011-06-15; event distance adopted; sources, extraction and comparison scripts saved | Done. First version is over-bright deeper into the umbra; readouts come from published curves, filter response and per-point errors still missing |
| M3: free-tropospheric aerosol extinction | Fixed OPAC prior added: 2–12 km, scale height 8 km, AOD550 = 0.013, segmented spectral exponent 1.58/1.21 | Done; results and improvement magnitude in [ROUND2.md](records/ROUND2.md). The prior is not a measured profile of 2011 |
| M4: broader independent observations | Published photometry tables for 1990/2007, December 2011 observed spectra and the 2019 GF-4 source conditions organized; field-of-view and exposure-time integration completed | Fixed-model comparisons of 1990 brightness and December 2011 continuum done (the latter with independent ephemeris, figure cross-checks and position sensitivity); the 2007 minimum is extrapolated background only; still not a full-color validation |
| M5: stratospheric-background aerosol component | Discrimination table (12 candidates × four protocols) — only `stratosphere_background` (12–35 km, AOD550 = 0.005, α = 0.74/1.14, same OPAC source) improved the shape axis; it entered the production default via `AerosolColumn`; table values checked as a preregistered prediction (`validation/results/stratosphere_background/preregistered_check.json`, 8 rows within_tolerance) | Done. Brightness-axis closure +0.19~+0.44 remains in the tension band — this component does not close the brightness residual (preregistered expectation; standing gap in [physics/AEROSOL.md](physics/AEROSOL.md)); broad low-layer increments (antarctic proxy) ruled out of the default |
| M6: cross-event required-atmosphere inversion (calibration lane) | Per-event low-layer increment inversion across four protocols + preregistered split/discrimination/retrieval spec + per-event `calibrated:` preset registry (`eclipse/data/calibrated_presets.json`); evidence in `validation/results/atmosphere_inversion/` and [OBSERVATIONAL_VALIDATION.md](records/OBSERVATIONAL_VALIDATION.md) | Done. R (`inversion.json:cross_event.R`) falls in the tension band (partially_consistent, preregistered mechanical judgment); all values are relative to the production-default baseline and never enter defaults, default output bit-unchanged; June's inverted value below the pinned BL lineage lower bound is recorded as a finding |
| M7: humidity axis + preregistered A/B re-inversion | Recovered OPAC 1998 WASO component humidity-band values (task-prior-source-survey product, SHA-pinned per band) entered `eclipse/aerosol.py` as optional FT-layer humidity growth (`ft_growth_factor`, f(50%) ≡ 1.0 bit-exact, q_rest pinned by subtraction, domain [0.5, 0.99] fail-loud); CLI `--rh` (opac/opac-ft only); transcription product `validation/results/prior_source_survey/humidity_axis_priors.json` (citations + geometry cross-checks + rulings); preregistered A/B re-inversion (A = RH-only root-finding, B = pinned-RH ladder {50, 70, 99}% × verbatim reuse of the task-5 machine, B(50%) bit-exact pool reproduction as a mechanical gate); (B−A) share table + R_RH verdict in `validation/results/humidity_axis/` | Done. Default RH = 0.5 output bit-unchanged (snapshot cmp); stratospheric and calibrated-increment layers ruled non-humidified (J4/J3, `humidity_axis_priors.json`); humidity spread read as a share measurement, not a causal claim |

First-version figures are kept in [figures/baseline](../figures/baseline) and
the first-round observational comparison in
[figures/observations2011](../figures/observations2011). These are historical
snapshots with fixed conditions and must not be silently redrawn as defaults
change. Later results are stored separately with explicit model and parameter
identification.

## The previous round's physical choice

M3 added only **aerosol-caused extra extinction** to the propagation model. At
the time, the residual's height distribution pointed to missing extinction on
lower-grazing paths; the solar-limb-darkening diagnostic did not remove the
deep-umbra red over-brightness, so the free troposphere was treated first.
This is an evidence-supported choice, not a unique causal identification of
the aerosol.

The default profile and spectral laws come from a background model in
independent literature; see [physics/AEROSOL.md](physics/AEROSOL.md). Default
comparisons fit neither brightness ratios nor spatial offsets, and do not add
boundary layer, stratospheric aerosol or cloud at the same time. Calibrating a
single AOD at 677 nm must always be listed separately as a diagnostic: keep
the independent-prior comparison, use 503 nm as the withheld band, and state
that a withheld band of the same eclipse is still not independent validation
for a new event.

Residuals from this round should continue to be preserved. Fitting a set of
curves by adding parameters cannot substitute for published observational
errors, instrument response and independent events.

## Current priority: broaden events before judging physical refinements

No more fitting degrees of freedom are being added for June 2011 for now. Real
atmospheres differ across events; a generic background model's capability must
be tested over multiple eclipses and different observable types. Aperture
differences from filters, fields of view, time averaging and lunar reflection
should be excluded first. Under the added 1990 whole-disk photometry, the fixed
OPAC model is about 57% too bright — less than the molecular model's ~90% —
but this residual cannot yet be uniquely attributed to an atmospheric process.

For December 2011, the nominal target trajectory and half-hour exposure time
averaging are complete; see [DECEMBER2011.md](records/DECEMBER2011.md): the
default model remains 41%–94% too bright, and the examined position
sensitivity cannot remove the brightness residual alone; aerosol improves
absolute brightness but the nominal blue/red spectral ratio does not follow.
Future physical refinements must explain brightness and shape across events
simultaneously before the necessity of event-specific profiles, ozone terms
and other physics is judged. The 2007 extrapolated minimum does not count as
an independent validation event. Data status, open conditions and reproduction
entries: the [multi-event reference](records/MULTI_EVENT_OBSERVATIONS.md).

## Future directions and entry conditions

The candidate directions below are organized by evidence, not a promise to
implement everything. Each round picks one discriminable direction from the
current residuals.

| Direction | Gap to solve and entry condition | Useful observables or checks |
|---|---|---|
| Complete filter and observational protocol | When residuals approach the influence of unknown bandwidths, sites or lunar geometry: obtain response curves, timings, pointing and errors first; monochromatic predictions must not be treated as strict bandpass predictions | In/eclipse bandpass ratios for the same lunar region; spectra integrated over responses; geometric-condition sensitivity |
| Add independent eclipse events | M4 added 1990 brightness and December 2011 exposure-averaged spectra; before extending further, weigh the processing cost of HST-class data; any calibrated result needs new observational events for testing | Events, positions and bands that did not participate in calibration; measured, extrapolated and figure-readout values recorded separately; keep absolute normalization |
| Extract public HST spectra | Public archive confirmed usable; first recover the actual lunar sampling positions and handle pointing drift, lunar albedo and the publication's correction protocol | Blue-end spectral slopes and in/eclipse spectral ratios. Its bands cannot alone verify the full ~600 nm absorption valley |
| Real temperature and ozone height profiles | Gas-density gradients affect both refraction and extinction; the same ozone column can have different height distributions. Check when position shifts or ozone-absorption shape residuals appear | Lowest path heights, bending, absorption valley, radial spectral ratios; one profile condition changed at a time |
| Solar limb darkening | Currently uniform Sun. It did not solve the main deep-red over-brightness in this diagnosis, but it changes light-path weights across solar regions and is not permanently negligible | Spatial and spectral changes under a fixed total solar flux; no-occultation limit and near-umbral-edge comparisons |
| Observation-constrained aerosol profiles | The default is now the OPAC two-layer column (free troposphere 2–12 km + stratospheric background 12–35 km); the tension-band closure of the brightness residual and cross-event differences suggest the remaining gap is not "one more layer" but event-specific profiles, clouds or source structure — discrimination has falsified the broad-low-layer-increment path into the default | Profile extinction, vertical AOD and slant-path transmission; fixed priors, calibration and independent prediction reported separately |
| Wavelength-dependent refraction | All wavelengths currently follow the reference path. Quantify dispersion if band position or narrow spatial structure becomes the dominant error | Light paths and landing points per wavelength, band position and width; per-wavelength conservation and ray checks |
| Single scattering into the receiver direction | Only unscattered light and losses out of the beam are computed; scattered-in light may matter at short wavelengths, weak-transmission positions or special aerosol conditions | Direct and scattered terms stored separately; spectra with instrument response, angular distributions and independent reference computations |
| Multiple scattering | Enters only when single-scattering results, observations or optical-depth evidence show it is insufficient; avoid taking on the computational cost prematurely | Order-by-order contributions and convergence, independent radiative-transfer references, cross-event spectra; brighter images alone do not prove correctness |
| 3-D atmosphere, clouds and limb weather | A spherically symmetric model cannot explain clear east/west asymmetry or local dark regions; needs atmosphere data with position and time along Earth's limb | 2-D lunar light fields, directional differences, time variation. Clouds and aerosols need separate constraints to avoid compensating for each other |
| Real lunar reflection and observer | The current screen omits lunar soil spectra, BRDF, spherical incidence angles and observer atmosphere; complete before matching lunar photographs or ground-based absolute photometry | Reflectance spectra and radiances of real lunar regions; exoatmospheric vs ground-based distinctions, keeping the illumination field as a separate output |
| Lunar texture maps and realistic rendering | Use ready-made lat/lon maps to restore maria, craters and true orientation; grayscale albedo first, then evaluate color maps and elevation | Neutral-screen vs textured rendering under the same illumination field, texture projection and libration direction; texture work must not mask atmospheric colors, see below |
| Vision and display | Photopic CIE colorimetry with common exposure cannot predict dark-adapted experience; introduce the corresponding models only when discussing naked-eye appearance or camera imaging | Raw spectra and photometric quantities, common-exposure images and chromaticity diagnostics kept separately; annotate camera response, chromatic adaptation and photopic/mesopic/scotopic conditions |

## Evidence the turquoise band requires

Verifying the phenomenon requires observations covering the relevant visible
bands, near the umbral edge, with recorded position and instrument response.
The cool-toned umbral edge seen so far cannot yet be equated with a reproduced
turquoise band. The original two bands and the inner-weighted curves constrain
only part of the transmission behavior; the seven December 2011 points do not
sample the outer edge directly either. Future comparisons should separate:

- spectral shape changes from ozone absorption and other extinction;
- the position and width of cool-colored regions and their relation to the
  brightness gradient;
- visibility under explicit exposure and display rules.

Similar image colors do not close these questions. Changing exposure, white
point, lunar reflection or observing conditions must remain traceable to the
same physical spectra.

[GF-4 2019 observations](https://doi.org/10.3390/rs16224181) provide a concrete
candidate: published blue-band widths of about 120–190 km and spatial maps of
red/blue radiance ratios at four epochs. Real geometry, band responses, lunar
composition and reflection must be matched before comparing position, width
and spectra; radiance ratios must not be read as pure atmospheric ratios or a
specified hue. This round records the evidence only; the gap is not declared
solved by raising saturation.

## Implementation direction for lunar texture

Ready-made assets first: the [NASA CGI Moon Kit](https://svs.gsfc.nasa.gov/4720/)
provides global lat/lon color and elevation maps suitable for adding surface
detail to the current neutral screen. Texture projection and lunar orientation
come first; spherical geometry, libration and elevation follow as needed —
each step verified separately, without requiring a full BRDF at once.

NASA explicitly notes that the color map has been processed (exposure, white
balance, hole filling) for artistic rendering and is not lunar soil spectral
albedo. Grayscale brightness textures can be used first to preserve the
model-computed atmospheric color; color textures must be marked as display
approximations. The newer EXR set is already linear color; sRGB TIFFs must be
decoded to linear space before multiplying with illumination. Scientific lunar
photometry should instead use the linked raw LROC WAC reflectance data, with
soil spectra and directional response handled explicitly. However refined the
rendering becomes, the raw spectral irradiance at the Moon is always kept as a
separate output.

## Per-round advancement criteria

Preserve the current model, parameters and reference conditions first; then
propose candidate causes that could explain the residual distribution; use
low-cost real computations to discriminate candidates and pick one to
implement. After any change, check numerical correctness, the effect on
existing observations, and explicitly unresolved residuals. If an addition
brings no terminal-perceivable benefit or verifiable physical improvement,
re-examine its direction and reason to persist. Unresolved physical gaps,
unclosed observational protocols and pure display choices are recorded
separately, so one class of problem never compensates for another.
