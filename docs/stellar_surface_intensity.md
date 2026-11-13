# Stellar continuum intensity in direct transmission

## Implementation

run_xds_model now accepts optional stellar_surface=StellarSurfaceConfig(...).
Default None preserves the old uniform fully-backed reference. StellarSurfaceConfig
is exported from xds. This is continuum background weighting, not an irradiance
penumbra or a stellar spectral/Rossiter-McLaughlin model.

```python
from xds import StellarSurfaceConfig, run_xds_model

surface = StellarSurfaceConfig(
    radius_rp=8.8465079373, orbital_separation_rp=359.44239577,
    orbital_period_hours=720., impact_parameter=0., u1=.4, u2=.25,
)
result = run_xds_model(planet, moon, sim, stellar_surface=surface)
```

Quadratic specific intensity:

I(mu)/Icenter = 1-u1(1-mu)-u2(1-mu)^2,
Imean/Icenter = 1-u1/3-u2/6.

Local weights are w0 * I(mu)/Imean for rays backed by the stellar disk, and zero
outside it. No renormalization of the visible annulus or surviving cells is
applied. The existing atmospheric amplitude, proportional to area/Rstar^2,
therefore remains normalized to OUT-OF-TRANSIT total stellar flux, not center
intensity. Both moonless reference and moon-affected map use the same local
stellar weights. The reference can now vary with time; the residual must be
A(t,v)-A0(t,v), not subtraction of a fixed A0(v).

In the observer sky basis (u,v), the planet center is placed at
(a sin(2 pi t/P), bstar Rstar cos(2 pi t/P)). Stellar radius and a are in Rp;
period/time are hours. bstar is the PLANETARY impact parameter in stellar radii,
distinct from lunar minimum projected separation in Rp. An explicit circular
orbit and fixed sky trajectory are assumed; the module does not infer them
from the moon configuration.

Stellar brightness is evaluated at each representative annular cell position.
Fractional lunar blocking remains area-integrated, with brightness and velocity
held constant per cell. Sharp stellar-limb crossings are NOT fractionally
integrated by this implementation, so ingress/egress needs further validation.
The reported illumination/observer/combined fractions retain their geometrical
area meaning and are not brightness-weighted intensity fractions.

## Pilot experiment

Reproduce: `PYTHONPATH=src python scripts/audit_stellar_surface.py`.
Requires the tracked pilot_phase_convergence.csv from the previous pilot.
Output table: stellar_surface_comparison.csv; figure:
outputs/stellar_surface_comparison.png/pdf.

Same physical host, atmosphere-linked annulus, moon radii .03,.06,.09 Rp,
central/grazing trajectories, and observing assumptions as pilot_phase_convergence.md.
Stellar radius 8.84651 Rp, separation 359.44240 Rp, period 720 h, planetary
impact parameter zero. Adopt illustrative u1=.4,u2=.25, uniform across wavelength.
These are chosen coefficients, not a prediction for a specific star/bandpass.
Lunar phases .25 (observer-side) and .75 (upstream) and four exposure origins
are compared with the original uniform-star pilot at equal phase. Original
zero-contact phases remain zero in this fully-backed direct-ray model and are
not rerun. No phase marginalization or population probability is estimated.

Both baseline and jointly refined production configurations are rerun with
the stellar surface: 480x16/order16,1201x401 time/velocity versus
960x32/order32,2401x801. Fluxes use identical atmospheric normalization. Noise
is held at the same out-of-transit continuum level to isolate weighting effects.
Photon-noise changes from opaque planetary/lunar transit photometry and noise
from constructing a finite reference are not modeled.

## Numerical comparison results

Across all three radii and four exposure starts, refined S/N ratios are:

| Trajectory | Lunar phase | Limb-darkened / uniform S/N |
| --- | ---: | ---: |
| Central | .25 | 1.1334--1.1346 |
| Central | .75 | 1.1895--1.1901 |
| Grazing | .25 | 1.2088--1.2089 |
| Grazing | .75 | 1.2089--1.20894 |

The two central phase groups sample different stellar background locations as
the planet moves along the chosen chord; the difference is geometric weighting,
not an assumed difference in atmospheric response upstream/downstream. At fixed
total stellar flux, central regions are brighter than the disk mean, so ratios
above one are expected. A different planetary chord or coefficients can change
both sign and magnitude of this comparison.

Maximum joint base/refined S/N difference is 0.07771%, below the 1% working
tolerance. Highest refined S/N is 0.11881: no case approaches the adopted
detection threshold. The stellar continuum correction does not reverse the
pilot's low-S/N conclusion.

## Validation status

Tests verify projected disk-mean intensity equals one under independent radial
area integration, zero background outside the stellar disk, exact recovery of
the previous residual with u1=u2=0 in a fully-backed window, time dependence of
the moonless reference, and A=A0+delta consistency. Coefficients producing
negative intensity on any part of mu in [0,1] are rejected analytically.

## Referee implications

Finite stellar size enters direct transmission through the specific intensity
behind each observer-directed atmospheric ray. In a fully-backed uniform star,
the previous direct-ray signal is recovered. Limb darkening redistributes the
stellar flux spatially and can increase central-transit signals relative to a
uniform star at FIXED TOTAL stellar flux. This does not create extra photons.
The effect differs from suppression of total incident irradiance in a penumbra.

This is an extension of continuum geometry, not yet realistic stellar line
formation. Stellar rotation, wavelength-dependent limb darkening, center-to-limb
spectral variations, spots, tellurics and atmospheric molecular templates remain
outside the experiment. New numerical convergence results should not be promoted
to a complete observational-feasibility claim.
