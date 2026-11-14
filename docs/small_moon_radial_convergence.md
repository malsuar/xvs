# Radial cell positions and smaller moons

Reproduce from repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/tmp/xvs-mpl OPENBLAS_NUM_THREADS=1 python scripts/audit_small_moons.py
```

Numeric results: small_moon_radial_convergence.csv. Figure:
outputs/small_moon_radial_convergence.png/pdf. Figures are regenerated, not
tracked. The underlying baseline is retained; defaults are unchanged.

## Design

Moon radii 0.06 and 0.03 Rp, circular a=3 Rp, P=21.238292071543107 h;
central b=0 and nearly tangent b=1.03+0.9 Rm. Thus grazing penetration below
outer tangency is only 0.1 Rm (0.006 or 0.003 Rp). Annulus [1,1.03] Rp.
Spin, spectrum, physical amplitude, noise and observing assumptions match
grazing_convergence.md. Inclination is chosen to achieve each b; these cases
do not hold the mutual spin-orbit angle fixed.

Spatial grids: 480x16, 960x32, 1920x64. Within each equal-area radial cell,
sample r^2 = r_in^2 + f (r_out^2-r_in^2), for f=0.25,0.5,0.75. Cell areas
and normalized weights are unchanged. Only f=0.5 is the midpoint quadrature
used by the forward-model default; off-center rules deliberately test spatial
sensitivity and generally introduce first-order radial bias. These are not
equally accurate alternative realizations of the midpoint rule.

The annulus builder exposes radial_sample_fraction for direct controls. It is
not added as a new physical SimulationConfig parameter. Angular origin is
held fixed in this experiment; joint angular/radial-origin combinations remain
untested. Each spatial grid is evaluated at four exposure start phases.

For each exposure phase the S/N reference is the mean of the three radial
positions on the finest grid. This is a diagnostic ensemble, not a proven
converged reference or an uncertainty distribution. Exact area overlap provides
an independent check, but is not an exact velocity-resolved S/N calculation.

## Results

Maximum S/N difference from that fine-grid ensemble, over radial positions and
four exposure phases:

| Moon radius Rp | Encounter | 480x16 | 960x32 | 1920x64 |
| --- | --- | ---: | ---: | ---: |
| 0.06 | Central | 2.812% | 0.313% | 0.00457% |
| 0.06 | Nearly tangent | 13.573% | 7.369% | 3.476% |
| 0.03 | Central | 3.948% | 2.655% | 0.0748% |
| 0.03 | Nearly tangent | 38.780% | 16.947% | 7.538% |

Maximum integrated-area error against analytic overlap:

| Moon radius Rp | Encounter | 480x16 | 960x32 | 1920x64 |
| --- | --- | ---: | ---: | ---: |
| 0.06 | Central | 1.933% | 0.315% | 0.0646% |
| 0.06 | Nearly tangent | 14.796% | 8.239% | 3.650% |
| 0.03 | Central | 3.033% | 2.092% | 0.381% |
| 0.03 | Nearly tangent | 40.585% | 17.389% | 7.389% |

For the MIDPOINT fine-grid rule alone, integrated-area errors are 0.0620%,
0.1207%, 0.3196%, 1.766% in table order. Thus radial-origin sensitivity must
not be presented as the exact midpoint S/N error. Nevertheless, even the
midpoint area test fails the 1% tolerance in the smallest grazing case.

At the finest midpoint grid, doubling time samples from 1201 to 2401 changes
S/N by at most 0.0141%, 0.1187%, 0.0255%, 0.6348%, respectively, over observing
phases. Temporal error is becoming relevant in the smallest grazing case.
Spectral sampling is fixed at 401 points over [-30,30] km/s; a separate spectral
refinement for these cases is still needed.

## Updated resolution decision

The previous 960x32 provisional grid is NOT sufficient across these smaller
moons under the requested 1% robustness criterion. Even 1920x64 cannot be
certified for the grazing cases using these controls. Do not launch the
population search using a blanket resolution guarantee.

Retain midpoint evaluation for ordinary runs. Refine based on moon radius,
annulus thickness and grazing penetration, and require BOTH analytic area and
observed-S/N checks. For central 0.03 Rp cases 1920x64 is a candidate starting
grid, not a universal guarantee. For near tangency use further adaptive
refinement or fractional-cell integration; the latter is preferable to assigning
an entire cell from one point. Contact timing needs its own analytic control.

The angular scale must resolve the lunar footprint; the radial scale must also
resolve shallow overlap, which can be much smaller than either Rm or annular
thickness. Numerical origin sensitivity is not physical atmospheric variation.
Cadence sensitivity must continue to be reported separately.

## Referee-response material

We extended the spatial audit to smaller satellites and deliberately displaced
radial quadrature points while preserving projected areas. The checks identify
limitations of point-mask integration in shallow grazing contacts, rather than
supporting a universal grid. We therefore withhold new population detection
claims until those configurations pass adaptive or fractional-area controls.
These tests concern numerical accuracy; they do not establish astrophysical
detectability or atmospheric-model realism.
