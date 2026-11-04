# First-version physical and numerical contract

## Quantity being computed

The primary result is spectral irradiance `E_lambda(r)` [W m^-2 nm^-1] on a
plane at lunar distance perpendicular to the Sun–Earth axis. A spherical,
horizontally uniform atmosphere makes this field radial. It includes
unscattered rays, both unobscured and refracted, after extinction. It excludes
light scattered into the receiving direction.

The visualized lunar disk is a uniform neutral screen sampling that plane.
It is an illumination diagnostic, not a lunar BRDF or human night-vision model.

## Forward curved-ray table

For tangent radius `r_t`, label the ray by the vacuum impact parameter

`b = n(r_t) r_t`.

With `q(r) = n(r) r`, the whole-limb bending and molecular columns are

```text
alpha(b) = -2 b ∫[r_t, r_top] (n'(r)/n(r)) / sqrt(q(r)^2-b^2) dr
C_i(b)   =  2   ∫[r_t, r_top] N_i(r) q(r) / sqrt(q(r)^2-b^2) dr
```

Both limbs are included. Substitution `r = r_t + u^2` removes the turning-point
singularity. Gauss–Legendre integration is split at the US1976 layer boundaries.
The density gradient is analytic within each layer. A local analytic limit of
`q(r)-b` avoids subtracting nearly equal numbers at the tangent point.

Refractivity is `(n_550-1) N_air/N_reference`. The reference density and Edlén
refractivity refer to dry air at 288.15 K and 101325 Pa. **This ray geometry is
achromatic**; extinction retains its wavelength dependence:

`T_lambda(b) = exp[-sigma_R(lambda) C_air(b) - sigma_O3(lambda) C_O3(b)]`.

An opaque spherical surface blocks any ray whose turning point is below sea
level. Thus the admitted interval starts at `b_ground = n(R) R`, rather than R.
The residual refractivity at the default 120 km atmosphere top is neglected;
tops with refractivity above 1e-9 are rejected because their vacuum interface
would require an additional refractive term.
No ducting or nonmonotone `n(r) r` is supported.

## Finite-distance Sun and conserved annular energy

Let S be Earth–Sun distance, D Earth–receiver distance and R_s the solar radius:

```text
m = 1 + D/S
a = D R_s/S
y(b) = m b - D alpha(b)
F_lambda = F_AM0(lambda) [AU/(S+D)]^2
```

These are consistent first-order small-angle (thin-lens) relations. Spherical
parallel-beam ray tables supply the local deflection law. Finite source
distance expands Earth's projected radius by m, retaining the approximately
16 km correction at mean lunar distance. The remaining finite-angle error is
not included in the numerical convergence claim.

A source point shifts each landing ring. A uniform solar disk therefore
convolves point-source illumination with a normalized disk of radius a.
For a receiver at radius r and landing ring radius `s=abs(y)`, the azimuth
fraction of the ring inside that source disk is

```text
f(r,s,a) = 1                                  if r+s <= a
         = 0                                  if |r-s| >= a
         = acos[(r^2+s^2-a^2)/(2rs)] / pi      otherwise.
```

Tangencies are measure-zero boundaries. The center r=0 uses the containment
rule directly. A negative y denotes rays crossing the shadow axis and is
retained through `abs(y)`, not discarded.

An incident annulus carries area `2 pi b db`. Relative to the unobscured solar
irradiance at the receiver, its contribution is

`d(E_lambda/F_lambda) = 2 m^2 b/a^2 f(r,abs(y(b)),a) T_lambda(b) db`.

Outside the atmospheric pupil, the vacuum contribution is analytic:

`1 - overlap_area(r, a, m*b_top)/(pi*a^2)`.

This avoids histograms, derivatives of the focusing map and separate focusing
corrections. The finite Sun regularizes the axial point-source caustic.
Irradiance is not forcibly capped at the unobscured value.

The landing map need not be monotone. The sharp lapse-rate change at the
idealized US1976 tropopause produces a small fold, independently confirmed
with a Cartesian ODE. Quadrature retains every table segment and splits at
every source-disk support/containment crossing of each segment. Sorting the
landing map or globally inverting it would lose physical ray branches.

## Numerical checks and their scope

1. **Geometric eclipse:** zero air/ozone recovers circle-overlap occultation
   with projected Earth radius mR, umbra radius `mR-a` and penumbra radius `mR+a`.
2. **Independent ray equation:** a Cartesian path-length ODE integrates
   `du/ds = grad(log n) - u (u dot grad(log n))`, plus both molecular columns.
   It checks bending, full path columns, the Bouguer invariant and direction norm.
3. **Global energy closure:** receiver-plane missing flux must equal

   `2 pi ∫ r[1-E_lambda/F_lambda]dr = pi(m*b_top)^2 - 2 pi m^2 ∫ b T_lambda db`.

   This is independent of where refracted rays land. Both extinguishing and
   lossless refracting atmospheres are checked.
4. **Resolution:** ray-table height spacing and spectral sampling are separately
   refined; photopic irradiance and normalized linear RGB differences are saved.

The independent ray solver shares the declared molecular medium with production
code, but uses a different equation, parameterization and integrator. It validates
transport, not whether the chosen atmosphere describes a particular eclipse.
See `validation/results/` for the committed checks and their parameters.

## Color and display

Spectra integrate against CIE 1931 2-degree CMFs. Relative Y retains physical
illumination ratios. One unobscured solar white is Bradford-adapted to D65 for
all positions and both ozone conditions. A shared exposure scales the entire
scene. Negative sRGB components are clipped; highlights outside display range
undergo common-channel compression, so they no longer retain brightness ratios.
Raw spectra and relative Y remain available independently of the image.

The equal-Y image is a clearly labeled diagnostic: local brightness is removed
to reveal spectral chromaticity. It is not a photometric rendering. Photopic
colorimetry does not account for mesopic/scotopic vision during dark totality.

## References and next discriminating experiments

- Vollmer & Gedzelman (2008), [spherically symmetric lunar eclipse irradiance](https://doi.org/10.1364/AO.47.000H52).
- Gedzelman & Vollmer (2008), [irradiance and color with satellite data](https://doi.org/10.1364/AO.47.00H149).
- García Muñoz & Pallé (2011), [direct and diffuse eclipse illumination](https://doi.org/10.1016/j.jqsrt.2011.03.017).
- Spectral inputs, atmospheric standard, refractivity and optical approximations:
  [data provenance](../eclipse/data/SOURCES.md).

The baseline produces a red/brown umbra and a pale cool edge that disappears
when ozone is removed. It does not yet establish the observed turquoise band's
width, saturation or absolute brightness. The first [observational comparison](OBSERVATIONAL_VALIDATION.md)
now shows substantial inner-umbra overprediction against June 2011 published
503/677 nm photometry. It motivates investigating additional low-altitude
extinction, with background aerosols as one candidate, alongside matched
atmospheric profiles and solar limb darkening. The digitized curves and unknown
filter responses are explicitly provisional; the production physics is unchanged.
Wavelength-dependent refraction and diffuse illumination remain further increments.
Each should be compared using raw spectra and common display settings.
