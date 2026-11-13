# Fractional-cell lunar blocking

## Implementation and scope

Select explicitly:

```python
sim = SimulationConfig(geometry="projected_annulus", spin_convention="observer",
                       blocking_integration="fractional",
                       blocking_quadrature_order=16, n_azimuths=480, n_radii=16)
```

Default point blocking remains unchanged for reproducibility. Fractional mode
currently requires equal parallel star-to-planet and observer directions, with
the effective annulus in the observer sky plane. Other geometries raise an
error instead of silently applying this approximation. Lunar upstream/downstream
position selects illumination or observer blocking. At exactly zero axial
separation the legacy strict ordering convention gives neither mask.

For a radial location r and a moon at projected distance d, occulted azimuths
lie within alpha=acos((r^2+d^2-Rm^2)/(2rd)) of the moon direction. We intersect
this interval analytically with each angular cell, including periodic wrapping,
then integrate its length in r squared using Gauss-Legendre quadrature. Radial
integration is split at disk tangencies r=|d-Rm| and r=d+Rm. Concentric cases
use an exact area formula. The blocked fraction f lies in [0,1] and the effective
weight is w0*(1-f); no renormalization of surviving weights is applied.

This integrates blocked AREA within cells; each cell still has one local
Gaussian profile evaluated at its representative velocity. It is not full
intra-cell spectral quadrature or vertical radiative transfer. The geometry
and physical amplitude thickness remain separately configured as before.

## Experiments

Run `PYTHONPATH=src python scripts/audit_fractional_cells.py`.
Numeric table: fractional_cell_comparison.csv. Figure:
outputs/fractional_cell_comparison.png/pdf. Same four smaller-moon cases and
observing assumptions as small_moon_radial_convergence.md. Exposure start is
zero here; earlier cadence sweeps are not implicitly repeated.

Practical grids 240x8, 480x16, 960x32 use order 16; reference 1920x64 uses
order 32. An additional order-32 run at 480x16 shifts the representative radial
position through fractions 0.25,0.5,0.75. Fractions of blocked area no longer
depend on that representative position; only cell velocities do. Order changes
and spatial refinements provide separate controls.

At 480x16, errors (percent):

| Rm/Rp | Impact Rp | Point area | Fractional area | Point S/N | Fractional S/N |
| --- | ---: | ---: | ---: | ---: | ---: |
| .06 | 0 | 1.892 | 0.00000804 | 1.873 | 0.000122 |
| .06 | 1.084 | 0.926 | 0.000343 | 0.0551 | 0.00986 |
| .03 | 0 | 2.542 | 0.0000191 | 2.543 | 0.000122 |
| .03 | 1.057 | 10.025 | 0.001027 | 8.989 | 0.00983 |

Area errors compare time-integrated fractions to analytic circle overlap at
the same epochs; S/N errors compare to the fine fractional grid, not an exact
S/N solution. Very small area errors do not measure continuous-time accuracy.
The numerical reference retains time/spectral and representative-velocity
errors. Current trials support fractional area integration, not universal
precision claims for arbitrary opacities, rotations or orbital geometries.

## Referee-response implications

Point masks can assign a whole cell to the wrong side of a lunar edge. We now
integrate fractional overlap and validate it against independent disk-overlap
areas, especially for small satellites and shallow grazing encounters. This
removes the dominant boundary-sampling error in the tested cases without
requiring enormous meshes. The original implementation remains available.
Before population recalculation, repeat cadence and angular-origin sweeps with
fractional integration, confirm temporal/spectral tolerances, and test faster
rotation/narrower profiles for intra-cell velocity errors. A 480x16 fractional
grid is a promising starting resolution for these controls, not an unconditional
default for the new population.
