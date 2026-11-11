# Spectral data and atmosphere provenance

The v0.2 aerosol component has its own [formula and provenance](../../docs/physics/AEROSOL.md).
It uses fixed OPAC background-model parameters, not a measured 2011 profile.
The two source exponents approximate OPAC interval slopes; no Mie library or
humidity evolution is imported from the Sky Color Simulation project.

These three tables were copied without modification from the **Sky Color
Simulation** project's `skycolor/data/`, source commit
`96b021eaaf96e0406b5303b90093fc888b76cf95`. This project packages its own copies
and has no runtime dependency on that repository. Values below describe the
data, not a claim that the other project's transport solver was validated for
lunar eclipses.

| File | Quantity | Original source |
|---|---|---|
| `solar_am0.csv` | Extraterrestrial solar irradiance at 1 au, W m⁻² nm⁻¹, 300–830 nm | ASTM G173-03 **extraterrestrial** column, distributed in [pvlib's reference table](https://github.com/pvlib/pvlib-python/blob/main/pvlib/data/ASTMG173.csv). This is the spectrum before atmospheric attenuation, not either AM1.5 surface column. |
| `cie1931_2deg.csv` | CIE 1931 2° XYZ colour matching functions, dimensionless, 300–830 nm | [CVRL, UCL](https://cvrl.ioo.ucl.ac.uk/database/data/cmfs/ciexyz31_1.csv), with zero padding below the original 360 nm lower boundary. |
| `o3_cross_section.csv` | O₃ absorption cross section, **cm² molecule⁻¹**, 300–830 nm | JPL 2006 evaluation as distributed with [NCAR TUV-x](https://github.com/NCAR/tuv-x), `data/cross_sections/O3_4.nc`. |

All files have a header and a 1 nm output grid. That grid is a resampling grid;
it does not imply that the original ozone measurements resolved 1 nm detail.
`ozone_cross_section()` converts cm² to m² by multiplying by `1e-4`. The API
rejects requests outside the packaged range instead of extrapolating. The
extraterrestrial spectrum integrates to 529.965 W m⁻² between 400 and 700 nm.

The visible ozone columns in the upstream file are identical at 218 and 295 K.
These copies therefore contain no visible temperature dependence. Real
Chappuis absorption has some temperature dependence. High-resolution gas
lines, other absorbing species, and wavelength-dependent solar limb darkening
are outside this first version.

SHA-256 checksums of the copied files:

```text
3f8de45150219f165c9a57c20346e3a7f04e183c3a1cdd2029ca3e8d3e9bff10  cie1931_2deg.csv
4606941959561aeeb9bb6f967661888cb7f82ba92e4009b739106ac091ac2462  o3_cross_section.csv
a532555ff4d5bd8071b7aabed80aa0c167ddf6fd6e2d51c3cfe1916d265e793c  solar_am0.csv
```

## Analytic atmosphere and molecular optics

- **Air density:** the [U.S. Standard Atmosphere, 1976](https://ntrs.nasa.gov/citations/19770009539)
  hydrostatic equations and geopotential lapse-rate layers up to 84.852 km
  geopotential altitude (approximately 86 km geometric). The standard's
  effective geopotential radius is 6,356,766 m. The separate spherical Earth
  radius used for ray geometry is 6,371,000 m. Above 86 km an explicitly
  approximate isothermal, constant-composition hydrostatic tail continues to
  120 km; the thermosphere part of US1976 is not implemented.
- **Reference number density:** the ideal-gas value at 101325 Pa and 288.15 K,
  using the exact SI Boltzmann constant, 2.5469164933 × 10²⁵ molecules m⁻³.
- **Refractivity:** [Edlén (1966), *The Refractive Index of Air*](https://doi.org/10.1088/0026-1394/2/2/002),
  for historical standard dry air, scaled linearly with molecular number
  density by the ray tracer. The visible-band baseline does not adjust humidity
  or CO₂ composition.
- **Rayleigh cross section:** the first-principles expression described by
  [Bodhaine et al. (1999)](https://doi.org/10.1175/1520-0426(1999)016%3C1854:ORODC%3E2.0.CO;2),
  with Edlén refractivity and a constant depolarisation ratio of 0.0279. At
  550 nm this yields 4.50698 × 10⁻³¹ m² per molecule. This is an approximation
  to the paper's recommended composition-dependent refractivity/King-factor
  prescription, not its full precision implementation.
- **Ozone vertical distribution:** an idealised Gaussian, centred at 25 km
  with standard deviation 8 km by default, truncated and exactly normalised
  over 0–120 km to the requested Dobson-unit column. One DU is
  2.687 × 10²⁰ molecules m⁻². This is an adjustable illustrative stratospheric
  profile; it is not a measured profile for a particular eclipse. Molecular
  density scaling and ozone column scaling are independent controls.

The model is dry, horizontally homogeneous, and cloud free. A finite top
leaves a tiny residual air density immediately below the vacuum interface;
the continuous interior gradient does not include that cutoff's jump.

## Calibration-lane preset registry

`calibrated_presets.json` is a data-layer registry of rendering presets built
from the cross-event atmospheric inversion (calibration lane; see
`docs/OBSERVATIONAL_VALIDATION.md`, section 跨事件所需大气反演). Each entry is
a complete atmosphere: the production two-layer OPAC default plus one
inverted low-layer increment whose AOD is fitted to that event's brightness
residual. The values are mechanically anchored to the committed inversion
product (`validation/results/atmosphere_inversion/inversion.json`, whose
sha256 each entry carries), and the pytest gate re-derives the anchor.
These are **calibrated values, not physical priors**: they never enter the
default model, and every consumer must carry the forced
`calibrated_not_physical_prior` provenance marker. The name-to-event mapping
lives solely in each entry's `event` field.

## Independent prior-source survey

The literature and data sources evaluated for the low-layer aerosol
increment's justification lanes (fixed prior / documented range / keep
current) are pinned in
`validation/results/prior_source_survey/sources.json` — per-file SHA-256,
byte count, retrieval URL and printed location, with `search_log.md`
recording every dead end. The pinned artefacts include the OPAC 1998
component data files recovered via GEISA/ESPRI Internet-Archive captures
(26 of 38 files, version-identified as the original 1998 distribution
against the pinned Hess, Koepke and Schult paper), an AERONET Level-2.0
site-year export, and open-access climatology papers. The survey changed
no default; its conclusion section is
[docs/PRIOR_SOURCE_SURVEY.md](../../docs/PRIOR_SOURCE_SURVEY.md).
