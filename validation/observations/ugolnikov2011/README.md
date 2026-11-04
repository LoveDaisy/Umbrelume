# June 15, 2011: published western-umbra photometry

Source: O. S. Ugolnikov, I. A. Maslov and S. A. Korotkiy (2011),
[arXiv:1106.6178v2](https://arxiv.org/abs/1106.6178v2), Figure 2, PDF page 4.

`points.csv` is **our digital readout of the published observed solid curves**,
not a numerical data release, original CCD exposures or independent measurement
samples. The paper's dashed gas-only model curves are excluded. The fixed PDF
hash, curve selection and coordinate transform are in `source.json`.

## Observable and comparison contract

- The ordinate is the brightness of a lunar surface element during the eclipse
  divided by that same element outside both umbra and penumbra. This substantially
  cancels lunar albedo. Compare the model's E_lambda/F_lambda directly, with no
  arbitrary brightness scaling or fitted offset.
- The abscissa is angular distance from the shadow center in degrees. We adopt
  `r=D*tan(theta)` on the model's plane at representative lunar distance D.
- 503 and 677 nm are the western profile, 30 degrees north of the equator. The
  paper reports standard-star and eclipse-free lunar calibration, sky-background
  correction and correction of the 503 nm filter's weak red leak using the other
  bands. These are published reduction claims; raw reductions are unavailable.
- 867 nm is outside the first model's 380–780 nm scope and is not extrapolated.

The first fifteen minutes of totality supplied the visible observations. We
adopt 2011-06-15 19:30 UTC. `geometry.json` records JPL Horizons DE441 geocentric
Sun and Moon apparent ranges at 19:20, 19:30 and 19:40 UTC, full query URLs and
response hashes. At 19:30 the adopted D is about 374394.877 km and S is
1.015759458 AU. These are range-based plane distances, not a reconstruction of
each lunar surface point, observing site or exposure.

## Graph extraction

```bash
# Optional one-time extraction dependency; not required for ordinary comparisons
uv run --with pymupdf==2.1 python tools/extract_ugolnikov2011.py
uv run python tools/fetch_eclipse2011_geometry.py

# Entirely offline with the committed CSV and geometry
uv run python -m validation.compare_observations --out outputs/observations2011
```

The extraction tool downloads the versioned PDF if absent and rejects a different
SHA-256. The source is a vector PDF. Its filled solid-line outlines are intersected
at chosen angles; the center of the upper/lower outline supplies log10(brightness).
The axes are 0–0.8 degrees horizontally and 1e-5–1e-2 logarithmically vertically.
Coordinates in PDF points are retained in the CSV for audit. Six anchors were
also checked against the rendered figure independently of the vector extraction.

The recorded graph-readout allowance is half the line thickness plus one PDF
point vertically, and one point horizontally. **It is not a measurement error,
standard deviation, confidence interval or a complete digitization-error model.**
Axis placement and source curve smoothing are additional limitations. Readout
precision should not be confused with the much longer decimal values needed for
deterministic reproduction. Do not calculate chi-square significance from these
allowances or treat the eighteen readouts as independent exposures.

## Missing contracts and limits

Filter responses/FWHM, pointwise observational errors and per-readout times are
not given in this paper. We compare at the stated effective wavelengths, without
claiming filter-integrated agreement. The author's gas-only calculation used an
equatorial temperature profile and 250 DU; our unchanged baseline is US1976 and
300 DU. A separate 250 DU run is a declared sensitivity, not a fitted correction
or a reproduction of that author's model.

The selected curve covers 0.400–0.575 degrees at 503 nm and 0.150–0.575 degrees
at 677 nm. Its outermost samples are still inside the umbra; **these data do not
sample the turquoise fringe itself**, and two bands do not determine full visible
color. The paper's separate eastern/western 867 nm asymmetry is not tested by a
spherical atmosphere. No empirical pass/fail threshold is assigned.
