# Baseline geometry and S/N audit

Baseline: fe70080. Development branch: codex/v0.2-scientific-revision.
This is an exploratory audit, not a corrected-model release. No forward-model
equations or historical results have been changed in this audit.

## Reproduction

From the repository root, with NumPy installed:

```sh
PYTHONPATH=src python scripts/audit_baseline.py
```

## Confirmed findings

### G1: The implemented grid is not a thin transmission annulus

The grid uses x=sin(longitude)cos(latitude),
y=cos(longitude)cos(latitude), z=sin(latitude), with longitude in [-pi/2,pi/2].
It covers a spherical hemisphere with y>=0. With the default observer along
-x, its sky projection occupies half a disk, including projected radii below
Rp. In the 81 by 55 grid, 0.414459 of the normalized weight lies at projected
radius below 0.9 Rp. Weights sum to one, but normalization does not establish
the physical validity of the emitting/absorbing geometry.

The mirror option fills the other projected half; it does not construct a
transmission annulus. Earlier explanations treating this option as a physical
annulus correction were insufficient. A separate annulus geometry is needed,
with projected radius at or above Rp and explicit area/opacity weighting.
PlanetConfig.atmosphere_half_width is not used by the current grid builder.

### G2: Spin parameters do not implement the stated observer inclination

The axis is Rz(lambda) Rx(i) zhat. For lambda=0 and the default observer -x,
both i=0 and i=90 degrees produce an axis 90 degrees from the observer.
Thus i cannot currently be described as inclination to the line of sight.
The vector cross-product velocity formula itself remains consistent with its
supplied axis. Define an observer-based basis and a documented sky reference
before changing this parameterization. Also document whether velocities are
positive toward or away from the observer.

### G3: Basic ray ordering passes

A controlled point and upstream/downstream moon correctly distinguish
illumination blocking from observer blocking. Both default direction vectors
being -x is consistent: one denotes photon propagation from the star, the
other propagation to the observer. This control does not validate finite-star
illumination or rays through an extended atmosphere.

### S1: Reported matched-filter S/N depends on numerical sampling

At fixed duration, exposure time, continuum S/N, instrument and geometry:

| Time samples | Velocity samples | S/N | Ratio |
| --- | --- | --- | --- |
| 61 | 121 | 0.02387383 | 1.00000 |
| 121 | 121 | 0.03348150 | 1.40244 |
| 61 | 241 | 0.03376270 | 1.41421 |
| 121 | 241 | 0.04735000 | 1.98334 |

The routine sums signal squared over grid points with constant noise and no
time or velocity quadrature weights. Numerical samples are therefore counted
as independent measurements. These are sampling experiments, not changes to
observational information. This prevents interpreting the baseline S/N as a
sampling-independent observing forecast.

### S2: Instrument and exposure are not fully simulated

Instrument resolution enters drift limits, but spectra are not convolved with
an instrumental line-spread function. Exposure time enters a heuristic
smearing penalty, but does not set the actual time sampling or perform exposure
integration. A systematic floor is added in quadrature as independent noise;
it is not a correlated-noise floor. The square-root gain from line count and
transit count assumes independent, equivalent contributions.

## Consequences for existing results

The archived maps remain outputs of the baseline mathematical model, but their
interpretation as a physical transmission annulus needs revision. Absolute
S/N values, detectable counts, parameter percentiles selected by S/N, and
instrument forecasts require recalculation after the corrections. Duration
statistics conditioned on those detectable subsets inherit that selection
issue. These findings alone do not establish whether the physical effect is
detectable or undetectable.

## Next implementation steps

1. Add an explicit projected annulus and analytic area/velocity controls;
   preserve the legacy geometry for reproducibility with explicit naming.
2. Define spin orientation relative to the observer and a stated sky axis.
3. Construct independent observing bins and exposure-integrated signals;
   apply the instrumental line-spread function before sampling.
4. Define noise on those bins, including reference-spectrum uncertainty and
   explicit assumptions for correlations and line stacking.
5. Require convergence under internal spatial/time/spectral refinement before
   injection-recovery or a new population search.

## Material for the referee response

Completed: baseline audit and reproducible sampling experiment.
Pending: physical geometry replacement, observing model and new results.
Do not describe pending work as completed in the response letter. Once
validated, explain that the revised application uses a projected atmospheric
annulus and independent observing samples, and replace the original numerical
detectability estimates with the recalculated values.
