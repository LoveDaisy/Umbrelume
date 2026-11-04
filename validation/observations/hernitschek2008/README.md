# 1990 and 2007: published photopic photometry

[Hernitschek, Schmidt & Vollmer (2008)](https://doi.org/10.1364/AO.47.000H62),
Tables 1–2. `points.csv` transcribes numerical tables, not curves. The publisher
HTML and its hash are identified in `source.json`; the downloaded HTML is an
ignored cache. Instrument details and Eqs. 4–10 were checked in the linked public
author-uploaded full-text HTML. No PDF access is claimed.

The 1990 CS-100 measurement covers the entire Moon (1° aperture, corrected for
lunar solid angle). The 2007 LS-110 covers a centered 1/3° aperture; its minimum
was extrapolated across clouds using the 1990 curve shape. **2007 is context,
not another independent observed minimum.** Its pre-eclipse baseline also carries
the paper's extra-aerosol warning. The printed ± values have unspecified coverage
and covariance; retain them without a Gaussian likelihood.

## Comparison contract

Use `minimum_luminance_cd_m2 / full_luminance_cd_m2`, not the reciprocal of the
rounded published dynamic range. For a neutral, uniform lunar screen this maps
to a ratio of spatially integrated photopic irradiances over the corresponding
aperture. Numerator and denominator must use the same aperture, wavelength grid
and response. A monochromatic value or a profile evaluated only at the Moon's
center is a different observable.

The source reduces observations to **one terrestrial air mass**, not vacuum.
Its correction uses an effective Rayleigh optical depth of 0.098. Eq. 9 scales
luminance to mean Sun–Moon distance; that scale cancels in the same-event ratio.
Eq. 10 also removes lunar-distance scaling from magnitudes. A model integral
using V(lambda) alone approximates this observation: residual observer extinction,
actual detector response and lunar reflectance remain unresolved. Neither these
systematics nor lunar surface structure cancel exactly for a reddened eclipse.

Table 1 geometry is geocentric. Its greatest-eclipse time does not identify the
exact exposure or establish when the measured luminance minimum occurred. Using
that geometry for the minimum is an explicit approximation. The source's
`smallest phase angle` equals shadow-center offset in this geometry, consistent
with `U = (rho_umbra + rho_moon - offset) / (2*rho_moon)` within tabular rounding.

Verify the table transcription against the pinned HTML cache with:

```bash
uv run python tools/extract_hernitschek2008.py
```

This checks raw numerical transcription; it does not validate observing conditions
or turn the reported error ranges into complete prediction tolerances.
