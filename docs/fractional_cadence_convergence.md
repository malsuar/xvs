# Fractional blocking: angular origins and observing cadence

Reproduce: `PYTHONPATH=src python scripts/audit_fractional_cadence.py`.
Data: fractional_cadence_convergence.csv. Figure:
outputs/fractional_cadence_convergence.png/pdf. Source scripts and table are
tracked; generated figures are ignored by Git.

## Experiment

Four controls: lunar radii 0.06 and 0.03 Rp, central b=0 and shallow grazing
b=1.03+0.9 Rm. Same orbit, spin, atmosphere and noise assumptions as
fractional_cell_integration.md. Uniform annulus [1,1.03] Rp; parallel ray
directions. Five-minute exposures, R=100000; reference continuum S/N 300 for
one minute at the same R, 200 equivalent independent lines, one transit,
independent bin floor 10 ppm. These remain controlled observing assumptions,
not instrument-specific detection predictions.

Grids: 240x8, 480x16, 960x32, fractional radial quadrature order 16.
Azimuth samples AND angular cell boundaries rotate consistently by 0,0.25,0.5
cell spacings; the physical orbit/spin do not rotate. Radial representative
positions remain midpoints. Exposure offsets: 0,0.25,0.5,0.75 times exposure
duration. Window [-2,2] h; 1201 temporal samples and 401 velocities over
[-30,30] km/s. Signal is zero at the window edges, so discarded incomplete
boundary exposures do not clip a signal.

Reference: fractional 1920x64, quadrature order 32, angular shift zero, evaluated
separately for each exposure phase. Spatial errors compare the SAME observing
phase to that reference. They include cell-velocity and radial-quadrature
approximations; this is not an analytic S/N reference. Independent exact
circle-overlap areas validate the geometric fractions.

## Results

Maximum relative S/N differences, in percent, over all three angular origins
and four observing phases:

| Rm/Rp | Impact Rp | 240x8 | 480x16 | 960x32 |
| --- | ---: | ---: | ---: | ---: |
| .06 | 0 | .001682 | .000382 | .000115 |
| .06 | 1.084 | .04267 | .009859 | .002215 |
| .03 | 0 | .001498 | .000357 | .0000875 |
| .03 | 1.057 | .07701 | .02686 | .01035 |

At 480x16, the maximum integrated-area error across all cases and angular
origins is 0.001027% against analytic overlap at identical epochs.

The observed S/N range over the four exposure phases is:

| Rm/Rp | Impact Rp | Min S/N | Max S/N | (max-min)/mean |
| --- | ---: | ---: | ---: | ---: |
| .06 | 0 | .0386167 | .0394219 | 2.059% |
| .06 | 1.084 | .00193141 | .00194253 | .5740% |
| .03 | 0 | .00967445 | .01324943 | 31.083% |
| .03 | 1.057 | .000392772 | .000397690 | 1.245% |

These are extrema on FOUR sampled phases, not bounds over every possible start.
The 31.083% range normalized to the mean corresponds to a 36.95% increase
from the lowest-S/N phase to the highest. No case approaches a detection under
these observing/atmospheric assumptions.

An independent temporal check at 960x32, zero angular shift, doubles time
samples to 2401. Maximum S/N changes across the four exposure phases are
0.0150%, 0.00889%, 0.0525%, 0.01698% in table order. Thus the large cadence
dependence in the small central case is not explained by the measured temporal
discretization error. It comes from dividing short contacts among finite
exposures; spatial refinement does not eliminate it.

## Resolution decision and remaining limits

Fractional 480x16 with order 16 passes the proposed 1% spatial-area and S/N
tolerance for these four controls and angular shifts. It is now a defensible
INITIAL working grid for these assumptions, unlike point masks on the same
grid. Still compare selected candidates with a doubled grid before accepting
population statistics. No default settings were changed by this experiment.

Keep observing start phase as an explicit parameter. Report a phase range or
a declared phase marginalization, not an arbitrarily selected favorable phase.
If a detection classification changes with cadence, flag it as cadence-sensitive.

This audit does not repeat radial-origin shifts at all angular/exposure
combinations, nor establish convergence for faster spins, narrower local lines,
nonuniform opacity, other annular thicknesses or nonparallel ray directions.
Velocity sampling remains fixed: repeat spectral refinement in the final
showcase and before population recalculation. The local spectrum still uses
one representative velocity per cell. Finite-reference uncertainty, correlated
noise and injection/recovery remain separate tasks.

## Referee-response material

The revised numerical validation separates discretization from physical
exposure averaging. Fractional area integration yields sub-percent spatial
S/N convergence across small central/grazing controls, while short events retain
measurable sensitivity to observing cadence. This supports explicit exposure
integration and reporting cadence-dependent ranges; it does not by itself
establish observational feasibility.
