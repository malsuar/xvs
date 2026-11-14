# Ingress and egress: joint stellar/lunar cell coverage

## Geometry and numerical method

Enable `SimulationConfig(stellar_cell_integration=True)` together with
projected_annulus and fractional blocking, and supply StellarSurfaceConfig to
run_xds_model. The former representative-point stellar weighting remains the
default for reproducibility. `time_center` (hours from planetary mid-transit)
permits short explicitly placed observing windows; its default zero preserves
previous time grids. Lunar phase is still defined at absolute t=0.

For each cell we now integrate two quantities separately:

Vcell = cell-average[I/Imean * indicator(stellar disk)],
Jcell = cell-average[I/Imean * indicator(stellar disk AND lunar disk)].

Reference weights are w0 Vcell. Missing absorption uses w0 Jcell when the
corresponding upstream/downstream blocking is enabled. Multiplying stellar
coverage by lunar coverage would be wrong in general because the intersections
are correlated within the cell. There is no renormalization of covered cells.
The reference spectrum remains A0(t,v), including stellar ingress/egress.

At fixed radius, circle boundaries become angular intervals. Their intersection
with the angular cell and with each other is clipped exactly, including angular
wrapping. Radial integration uses Gauss-Legendre quadrature in r squared, split
at stellar/lunar radial tangencies. Angular integration of limb darkening uses
eight Gaussian nodes within each retained interval (uniform intensity is exact
in angle). Local velocity remains evaluated at the representative cell point.
This is area/intensity integration, not exact intra-cell velocity integration.

The area arrays illumination_fraction/observer_fraction/combined_fraction retain
their historical geometric lunar-mask meaning and do not include stellar
backing or brightness. Use reference_absorption and delta_absorption for the
observable; do not interpret those fraction arrays as observable coverage during
ingress/egress. Fractional lunar mode still requires equal parallel ray directions.

## Analytic controls

Uniform stellar backing is compared against the difference of stellar-disk
overlap with Rout and Rin, divided by annulus area. Fully backed uniform stars
recover the existing lunar fractional integrator. Empty backing yields zero
weight. A coincident stellar/lunar disk test verifies Jcell=Vcell for partially
covered cells, whereas multiplying fractions would incorrectly yield Vcell^2.
Intensity and spectral normalization controls from previous tests still apply.

## Designed limb encounter experiment

Run `PYTHONPATH=src python scripts/audit_ingress_egress.py`.
Output table: ingress_egress_comparison.csv. Figure:
outputs/ingress_egress_comparison.png/pdf.

Same host, atmosphere-linked annulus h=0.01025833, veq=12 km/s, sigma=1.5 km/s
and .03 Rp moon as the pilot. Stellar surface uses illustrative u1=.4,u2=.25,
planetary impact zero, Rstar=8.84651 Rp, aplanet=359.44240 Rp, period720 h.
Short .8-hour windows are centered where a ray at sky coordinates (u=0,
v=1+h/2) lies on the stellar boundary. Lunar inclination=90 deg, node=0 gives
a polar orbit for the revised +z spin. Each lunar phase is deliberately selected
to put the moon center at that ray at window center, so the lunar and stellar
boundaries intersect the same annular region. This is a stress control, not uniform phase selection
or a population prediction. It is not the previous mid-transit pilot window.

Three calculations per window: old point-sampled stellar weighting on the base
grid, joint fractional stellar/lunar weighting on the same grid, and jointly
refined production settings. Base uses reference_simulation with time_center;
refined doubles spatial, temporal and spectral intervals and radial order.
Actual exposures remain five minutes, four start phases; noise and physical
amplitude remain identical. Window-edge partial exposures are discarded; these
are intentionally limited observing windows, not whole-transit S/N values.

The uniform-star analytic coverage test uses 31 separate epochs per window.
Its maximum absolute fraction error is not a relative error near zero backing.

## Results

| Window | Maximum joint S/N refinement change | Refined S/N range over four exposure starts |
| --- | ---: | ---: |
| Ingress | 0.1726% | 0.002348--0.003310 |
| Egress | 0.1801% | 0.002719--0.003341 |

Both cases pass the 1% working numerical tolerance. Maximum absolute error in
uniform stellar coverage versus independent circle overlap is 1.37e-12 across
the sampled epochs. These small area errors are not claims about limb-darkened
spectral accuracy or continuous-time contact timing.

On the SAME base grid, point-sampled stellar weighting differs from joint
integration in S/N by at most 0.2641% in these cases. Thus the new method removes
a conceptual/edge-integration limitation but does not imply a large correction
to every already-resolved calculation. The eight phase-specific comparisons are
preserved in the table. `scripts/plot_limb_snr_comparison.py` regenerates a second
figure, outputs/limb_snr_comparison.png/pdf, directly from that table.

The reference absorption rises smoothly from zero as the annulus gains stellar
backing and vanishes during egress. Lunar residuals are evaluated only over
the common stellar/lunar region. The very low S/N is conditional on these
designed polar .03 Rp lunar encounters and short observing windows; it should
not be compared as a detection forecast to the longer mid-transit pilot with
different phases and alignments.

## Limits for the referee

This explicitly models the atmospheric annulus entering/leaving the stellar
background, not total opaque-planet/lunar photometric light curves. Photon-noise
changes during photometric transit, finite reference noise, stellar spectra,
stellar rotation/CLV, refraction and scattering remain absent. The paraxial
observer-directed rays and sky-plane effective atmosphere assumptions persist.
Independent convergence of angular limb-intensity quadrature is still required
for an arbitrary limb law; joint mesh refinement alone does not certify it.

An additional fixed-geometry check compares 8 versus 16 angular intensity nodes
at nine epochs around a stellar-limb crossing with the same quadratic law.
Maximum absolute cell-mean reference-intensity change is 8.26e-6; jointly blocked
intensity change is below 5e-16 for that particular check. This is a local
quadrature diagnostic, not an independent full S/N convergence proof.

Results should support the revised numerical method only within the tested
limb encounters. Do not extrapolate to every grazing planetary transit or
describe these designed phases as favorable population fractions.
