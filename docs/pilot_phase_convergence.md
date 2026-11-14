# Deterministic pilot: physical assumptions, phases and joint convergence

Reproduce: `PYTHONPATH=src python scripts/pilot_phase_convergence.py`.
Machine-readable assumptions: pilot_assumptions.json. Results:
pilot_phase_convergence.csv. Figure: outputs/pilot_phase_convergence.png/pdf.
Generated figures are ignored by Git; scripts, assumptions and tables are tracked.

## Fixed physical scenario

One host, not a population: Mp=0.7 Mj, Rp=1.1 Rj; solar mass/radius star with
Teff=5772 K; circular planetary period 30 days. Zero-albedo full-redistribution
temperature estimate is 640.30 K; adopt isothermal 640 K, mean molecular weight
2.3, five effective scale heights and unit line contrast. Rigid equatorial
speed 12 km/s, intrinsic Gaussian sigma=1.5 km/s; observer spin axis +z.

Unlike earlier fixed-thickness controls, geometric annular thickness is tied
to the same atmosphere used for physical amplitude: h=5H/Rp=0.0102583.
The amplitude still uses the thin-annulus approximation 2 Rp(5H)/Rstar^2;
the omitted quadratic area term is about h/2=0.513% relative to that scale.
This is a specified approximation, not evidence of equal physical accuracy.

Moon radii .03,.06,.09 Rp with assumed bulk density 5.5 g/cm^3. Circular a=3 Rp;
Keplerian period 21.23829 h (test-particle approximation). Roche radius 1.23385 Rp,
adopted prograde outer limit .35 Hill=7.62520 Rp. Spin period 11.4379 h passes
the breakup check. These are simple instantaneous dynamical filters, not tidal
survival, formation or atmospheric-evolution models.

For each moon radius, central and shallow grazing projected trajectories are
chosen: b=0 and b=1+h+0.9 Rm. Orbital inclination asin(b/3), node=90 deg.
Thus mutual spin-orbit angle varies between configurations. These are designed
test geometries, not draws from an inclination distribution. The planetary
transit impact parameter is separate from this lunar projected minimum.

## Window, phases and observing assumptions

Observe a fixed four-hour interval centered on planetary mid-transit. Approximate
central planet-center transit duration is 5.641 h (ignoring planet radius), so
this is a mid-transit subwindow, not the entire transit or an extension outside
transit. Stellar ingress/egress/limb darkening are not synthesized. Ray directions
remain parallel/point-source; finite stellar angular size and penumbra are not
included and require assessment before observational claims, especially for
small moons at compact separations.

Uniform deterministic lunar phases at t=0: 0,.25,.5,.75 cycles, all retained
with equal status and no peak centering. Four exposure start offsets: 0,.25,.5,.75
of five minutes. These are sparse controlled phase grids, not converged orbital
phase marginalization or occurrence probabilities. No heuristic score is used.

R=100000, reference continuum S/N 300 at one minute and R=100000, actual exposure
five minutes, 200 equivalent independent lines, one transit, independent-bin
noise floor 10 ppm. The forecast assumes a noiseless known template and reference;
no molecular spectra, correlated noise, contaminating signals or recovery test.

## Per-candidate convergence

Six designed candidates x four lunar phases =24 model configurations, each at
base and jointly refined grids, and four exposure phases =96 observing rows.
Base is reference_simulation: fractional 480x16/order16, 1201 times,401 velocities
over [-40,40]. Refined:960x32/order32,2401 times,801 velocities. Equal observing
phase is compared using |SNRbase/SNRrefined-1|.

All 48 NONZERO observing rows pass the 1% joint criterion. Maximum change is
0.07764%, in the smallest central moon. The other 48 rows have exactly zero
signal on both grids; relative errors are recorded as zero by convention and
flagged no_signal=True. They are consistent nondetections, not relative-accuracy
certificates. Half of the tested lunar phases have no contacts in this window;
this is not a population detection fraction.

S/N ranges over the four exposure offsets, conditional on contacting lunar
phase .25 (phase .75 gives the symmetric result):

| Rm/Rp | Geometry | Minimum refined S/N | Maximum refined S/N | Cadence range / mean |
| --- | --- | ---: | ---: | ---: |
| .03 | Central | .013627 | .017018 | 22.06% |
| .03 | Grazing | .001479 | .001506 | 1.80% |
| .06 | Central | .051387 | .051748 | .701% |
| .06 | Grazing | .007280 | .007333 | .716% |
| .09 | Central | .097954 | .099874 | 1.94% |
| .09 | Grazing | .018261 | .018353 | .500% |

No candidate reaches S/N=5. This pilot demonstrates convergence and phase
sensitivity, not feasibility for Earth-sized moons or validation of the old
sweet-spot results. No favorable host was selected to manufacture detection.

## Decision before a full search

The reference setup is numerically usable for a broader PILOT under these
uniform Gaussian assumptions, with candidate refinement flags. A full search
needs a frozen physical parameter domain and more densely sampled lunar phases;
otherwise phase aliases can miss contacts. Report a declared phase statistic
(including zeros) and cadence envelope rather than best phase alone. Near a
detection threshold, refine observing and orbital phases as well as numerical
grids. Add spatial-origin/analytic checks for new footprints where appropriate.

The current low S/N should prompt realistic spectral and instrumental feasibility
experiments before spending resources on a large optimization. Real atmospheric
templates, finite-star effects and injection/recovery remain outstanding
scientific checks; numerical convergence does not settle them.

## Material for the referee

We fixed a physically specified pilot with Keplerian satellite periods, basic
Roche/Hill/breakup checks, atmosphere-linked annular thickness and uniform lunar
phases. We retained zero-signal cases, evaluated observing-phase ranges and
checked joint numerical convergence by configuration. The resulting low S/N
and strong cadence sensitivity of the smallest central satellite are reported
without interpreting the deliberately sparse pilot as a population prediction.
