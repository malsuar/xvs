# Temporal convergence and consolidated numerical reference

Reproduce: `PYTHONPATH=src python scripts/audit_temporal_convergence.py`.
Data: temporal_convergence.csv. Figure: outputs/temporal_convergence.png/pdf.

## Temporal controls

Moon radius .03 Rp, a=3 Rp, P=21.238292071543107 h, central b=0 and nearly
tangent b=1.057. Fractional annulus integration on 480x16 cells, radial order
16, observer spin convention, uniform annulus [1,1.03] Rp. Velocity step .1
km/s over [-40,40]. Four five-minute exposure origins: 0,.25,.5,.75 exposure.
Other physical/noise assumptions match fractional_cadence_convergence.md.

Temporal samples 601,1201,2401,4801 over four hours correspond to 24,12,6,3
seconds. Observing bins, photon noise and exposures remain fixed. Each S/N is
compared to the SAME exposure phase of the 3-second calculation, at fixed
spatial and spectral grids. The reference is numerical, not analytic.

Maximum S/N errors across exposure phases, in percent:

| veq km/s | sigma km/s | Encounter | 24 s | 12 s | 6 s |
| ---: | ---: | --- | ---: | ---: | ---: |
| 12 | 1.5 | Central | .3226 | .06640 | .01394 |
| 12 | 1.5 | Grazing | .08667 | .02087 | .00419 |
| 25 | .5 | Central | .3224 | .06638 | .01394 |
| 25 | .5 | Grazing | .09986 | .02305 | .00467 |
| 25 | .2 | Central | .3224 | .06638 | .01393 |
| 25 | .2 | Grazing | .1022 | .02357 | .00478 |

Twelve seconds passes the temporal 1% criterion with margin in these controls.
Cadence-dependent S/N variation remains at the fine sampling; it must not be
treated as numerical noise. No intrinsic-profile or arbitrary-phase convergence
claim follows from convergence of this instrument-broadened, integrated S/N.

## Reference configuration available in code

`xds.reference.reference_simulation` constructs explicit NEW settings without
changing historical defaults. Example:

```python
from xds import PlanetConfig, MoonConfig, run_xds_model
from xds.config import InstrumentConfig, ObservationConfig
from xds.reference import reference_simulation, NumericalTolerances

planet = PlanetConfig(v_eq=25., line_sigma=.2)
instrument = InstrumentConfig(spectral_resolution=100000.)
observation = ObservationConfig(exposure_time=5/60)
sim = reference_simulation(planet, 4., observation, instrument)
# Explicit aligned orbit for the revised +z spin and observer -x.
moon = MoonConfig(radius=.03, semi_major_axis=3.,
                  period=21.238292071543107, inclination=0., phase=.08)
result = run_xds_model(planet, moon, sim)
tolerances = NumericalTolerances()
```

This illustrates numerical selection, not a new physical showcase or the exact
orbit of the audit. The caller must set physical parameters and orbital/spin
alignment consistently. No atmosphere/noise normalization is inferred.

Starting settings:

- projected annulus, observer-relative spin, fractional blocking;
- 480 angular x 16 radial cells, Gaussian radial order 16;
- target time step exposure/25 (12 seconds for five minutes); the actual
  uniform grid uses enough samples that its step is no larger;
- input velocity step <= min(.2 km/s, sigma/2, (c/R)/8);
- half-window >= max(40 km/s, veq*(1+h)+8 combined Gaussian widths), rounded
  upward in km/s; the combined width includes intrinsic and instrumental sigma.

The window rule supplies generous padding, but has not been separately swept
in all controls. More restrictive settings than strictly necessary are chosen
to provide a margin; they are not derived universal error bounds. The sample
count can be large for long windows/short exposures, so do not launch a large
search without profiling and convergence checks.

## Acceptance policy

NumericalTolerances records relative 1% tolerances for spatial S/N, temporal
S/N, spectral S/N and integrated area. It also records the four exposure-phase
fractions. This object is a policy record, NOT an automatic acceptance engine.

Require per-candidate comparisons to refined grids, including angular and radial
origin checks where needed. Compare equal observing phase, and report physical
cadence variation separately. Validate geometric area against an independent
reference where available. Contact-time accuracy requires a separate criterion.
Each tolerance is per-axis; three errors each below 1% do NOT guarantee total
error below 1%. Use a jointly refined configuration if a total-error claim or
a near-threshold detection classification is needed. The audit family supports
these starting settings only for the uniform Gaussian parallel-ray model.

## Consolidated referee status

Completed: baseline geometry/S/N audit; projected annulus and spin convention;
physical observing bins, LSF and finite exposures; fractional cell blocking;
analytic overlap checks; spatial, sampling-origin, exposure-phase, spectral and
temporal controls, including small grazing moons and rapid rotation/narrow lines.

Still required before a new population result: final showcase geometry and
physical amplitude choice, production-versus-audit consistency, joint refined
controls for selected candidates, and phase-aware reporting. Realistic molecular
spectra, nuisance atmosphere, noise covariance and injection-recovery remain
scientific tasks distinct from numerical convergence. The new numerical checks
cannot validate the original detectable fractions or establish feasibility.
