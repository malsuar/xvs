# Spectral convergence, faster rotation and narrow local profiles

Reproduce: `PYTHONPATH=src python scripts/audit_spectral_convergence.py`.
Data: spectral_convergence.csv. Figure: outputs/spectral_convergence.png/pdf.

## Controlled experiment

Fractional blocking, annulus [1,1.03] Rp, moon radius .03 Rp, central b=0 and
nearly tangent b=1.057. Other physical/noise assumptions match the preceding
fractional-cadence audit. Four combinations (veq,sigma), in km/s:
(12,1.5), (25,1.5), (25,.5), (25,.2). These are numerical stress controls,
not a new population or measured atmospheric widths. Gaussian amplitude stays
fixed, so narrower lines have smaller equivalent width; absolute S/N comparisons
between widths cannot isolate Doppler resolution alone.

Input velocities span [-40,40] km/s with 101,201,401,801,1601 points, giving
steps .8,.4,.2,.1,.05 km/s. The observing operator retains R=100000 (FWHM
about 3 km/s), five-minute exposures and four exposure phases. Physical
observing bins and noise are unchanged during spectral refinement.

Spatial grid is 480x16, quadrature order 16, midpoint velocities and zero
angular-origin shift. The spectral reference is 1601 points ON THAT GRID.
The spatial control compares 960x32/order32 at the same 1601 velocities, so
it includes both spatial refinement and radial quadrature-order change.
All reported errors are maxima across the four observing phases at fixed phase.
Time sampling remains 1201 points; previous temporal checks are not repeated
for every faster-spin/narrow-line combination. Velocity-window padding and
arbitrary velocity-grid origins are not varied.

## Results

| veq | sigma | Max spectral S/N error at step .8 | At step .4 | S/N change doubling spatial grid |
| ---: | ---: | ---: | ---: | ---: |
| 12 | 1.5 | <1e-10% | <1e-10% | .00884% |
| 25 | 1.5 | <1e-10% | <1e-10% | .03862% |
| 25 | .5 | .00406% | <1e-10% | .06837% |
| 25 | .2 | 1.109% | .1469% | .07523% |

At step .2, even sigma=.2 gives a maximum relative S/N difference below
2.5e-10 (dimensionless), against the .05 reference. This very small discrepancy
concerns smooth, Gaussian, instrument-broadened templates at fixed spatial/time
sampling; it is not total physical/numerical accuracy. Plot differences below
1e-12 are displayed at a floor of 1e-12 for readability.

For sigma=.2 and step .8, central and grazing errors are 1.109% and .903%,
respectively. Thus a coarse velocity grid can fail the 1% criterion even when
the instrumental width is much larger than the intrinsic profile. Local
features must be integrated before/while applying the instrumental response;
instrumental broadening does not recover undersampled input information.

## Provisional numerical policy

For these tested widths/resolution, fractional 480x16 and input step <=.2 km/s
pass the spatial and spectral S/N portions of the 1% criterion. A conservative
starting rule for other Gaussian profiles is delta_v <= min(sigma/2,(c/R)/8),
followed by explicit factor-two refinement. This rule is an engineering margin,
not a theorem or empirically validated bound outside the tested cases.
Keep spectral padding large enough for rotation at Rout plus intrinsic/LSF
wings, and verify it separately for the final configurations.

The representative-cell approximation remains: this audit supports it for
these integrated-S/N controls at R=100000, not for resolving all intrinsic
profile moments. Higher resolving power, different sky angles, winds, narrower
spectra or spatially heterogeneous opacity may need intra-cell spectral
quadrature. Time/cadence and grid-origin tests for the stress combinations
should precede final population guarantees. None of these checks demonstrates
recovery of a lunar signal against nuisance atmospheric variability.

## Referee material

We explicitly tested independent velocity sampling at fixed instrumental
resolution and exposure duration, including faster rotation and narrower
Gaussian features. The observed S/N converges under spectral and spatial
refinement in these controls; coarse sampling can fail for narrow intrinsic
features. These checks justify numerical choices for the controlled Gaussian
application, but do not substitute for realistic molecular templates and
injection-recovery experiments requested by the referee.
