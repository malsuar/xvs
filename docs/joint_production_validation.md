# Joint refinement and production/audit consistency

Reproduce: `PYTHONPATH=src python scripts/audit_joint_production.py`.
Table: joint_production_convergence.csv. Figure:
outputs/joint_production_convergence.png/pdf (regenerable, not tracked).

## Scope and parameters

Four controls: (veq,sigma)=(12,1.5),(25,.2) km/s; moon radius .03 Rp,
central b=0 and grazing b=1.057; annulus [1,1.03] Rp. Circular a=3 Rp,
P=21.238292071543107 h; orbital inclination asin(b/3), node=90 deg,
phase=.25. Observer spin convention i=90, lambda=0. Parallel stellar/observer
directions -x. Same uniform atmosphere and observing assumptions as preceding
audits (0.7 Mj,1.1 Rj,500 K, stellar radius one solar, five scale heights;
five-minute exposures, R=100000, reference continuum S/N 300 at one minute,
200 equivalent independent lines, one transit and independent bin noise floor
10 ppm). Four observing phases 0,.25,.5,.75 exposure.

The actual production factory reference_simulation is used, not a manually
copied configuration. Base: spatial 480x16, order16, 1201 time points over
four hours (12 s), 401 velocity points for sigma=1.5 (.2 km/s) and 801 for
sigma=.2 (.1 km/s), both over [-40,40]. Joint refined: spatial 960x32,
order32, 2401 times (6 s), 801/1601 velocities (.1/.05 km/s). Physical
observing bins and phases are identical. No population search is run.

## Production consistency

For each grid and case run_xds_model is compared against the audit's direct
blocked-contribution sum, with identical orbital/configuration parameters.
Time/velocity grids are asserted identical. Checks compare the full intrinsic
residual map, blocked-fraction time series, absorption=reference+residual,
and S/N from detectability_report versus observe_residual plus the documented
matched-filter norm. Maximum observed map/fraction/S/N discrepancies are zero
at floating-point precision for these cases.

Both paths share the low-level geometry, orbit, fractional overlap and observing
operator. This is an orchestration/parameter-consistency check, NOT independent
physical validation. Independent circle-overlap geometry tests and Gaussian
observing-integral tests remain the separate analytic controls.

Fractional production synthesis now computes the negative blocked contribution
directly and forms absorption=reference+residual, avoiding cancellation from
subtracting nearly equal profiles. This algebraically equivalent change also
avoids repeatedly summing all unblocked profiles. The point/legacy path retains
its historical calculation. A regression checks that integrated deficit equals
d*sigma*sqrt(2*pi)*blocked_fraction and that illumination/observer contributions
sum to the combined map in the tested parallel geometry.

## Joint refinement results

Maximum absolute relative S/N change |SNR_base/SNR_refined-1| over four observing
phases, in percent:

| veq km/s | sigma km/s | Impact Rp | Joint change |
| ---: | ---: | ---: | ---: |
| 12 | 1.5 | 0 | .05230% |
| 12 | 1.5 | 1.057 | .02568% |
| 25 | .2 | 0 | .05157% |
| 25 | .2 | 1.057 | .09318% |

Largest analytic integrated-area error across all base/refined runs is 0.001027%.
The combined S/N change passes the 1% working tolerance, rather than inferring
a total tolerance from separate axes. The refined grid is still numerical and
not an exact continuous-time/spectral solution.

## Decision and referee wording

The numerical reference factory passes this joint check and agrees with audit
orchestration for the tested small central/grazing controls. It is suitable as
a starting point for a limited pilot sample with per-candidate refinement and
phase-aware reporting. It is not blanket certification of the full search
parameter space. No thresholds or percentages from the historical population
are validated by this experiment.

Before a full search: freeze the physical/observing assumptions, explicit
orbit-spin geometry, population/phase selection and a candidate acceptance
policy. Preserve numerical convergence flags, refine near-threshold cases,
and profile computational cost on a small pilot. Molecular templates,
atmospheric alternatives, correlated noise and injection/recovery remain
necessary for the referee's scientific feasibility concerns.
