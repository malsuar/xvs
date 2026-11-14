# Direction review and first phase/period-blind protocol

## Are we moving in the right direction?

Yes for model validation, but not yet for an observational detection claim.
The projected annulus, fractional cell blocking, instrumental bins, exposures,
direct-transmission stellar weighting and molecular opacity have passed
controlled checks. The selected physical system remains far below useful S/N
at the original noise. Recovery near the five-unit information scale has
required deliberately reduced noise, not an established instrument budget.

The atmospheric tests have also changed the scientific interpretation:
EXOVEL is not just a line-depth change, but flexible velocity/opacity structure
can remove lunar information or bias an incomplete model. Temporal coherence
helps only conditionally. A moon preference can coexist with a poor atmosphere
fit. Numerical agreement, model distinguishability, discovery and orbital
measurement are different milestones.

The next useful step is therefore not another population sweep. It is a
conditional blind search, empirical false-alarm calibration and an assessment
of orbital identifiability. After that, realistic stellar/telluric contaminants,
throughput/noise and atmospheric-model uncertainty are needed before revising
the population detectability claims. Stop interpreting high scores as evidence
if atmospheric controls produce comparable scores; stop reporting a measured
period if aliases or search boundaries dominate the recovered periods.

## Scope fixed before calibration

- Same host, moon radius, aligned circular geometry, four-hour window and CO
  model as the preceding bounded nonlinear pilot. Radius, orientation and host
  properties are held fixed; only lunar phase and period are unknown here.
- Period and separation are coupled by Kepler's law at 0.7 Mj and 1.1 Rj.
  With a prescribed lunar density 3 g/cm3, the fluid-Roche lower bound is
  1.5101 Rp. The adopted 0.35-Hill outer bound is 7.6252 Rp. These yield
  7.5849--86.0626 h. This is a dynamical screen, not a lifetime calculation.
- Phase covers [0,1) with no injected coordinates added to the bank. Initial
  grid: 73 logarithmic periods and 1024 phases. Top three separated grid seeds
  receive continuous local score refinement within their grid neighborhoods.
- Compare the three- and five-control bounded atmospheric families. The
  nonlinear atmosphere-only fit is repeated on EVERY simulated data set with
  the same fixed starting points. Failures are retried, never silently removed.
- Physical lunar blocking is between zero and one. No signal amplification.
  Main calibration retains the explicitly idealized 11.39465 ppm noise from
  the preceding control; it is not a new claim about the native instrument.
- Exposure integration uses sixteen Gaussian nodes after checking the
  eight-node approximation against independent production spectra.

## What statistic is actually calibrated?

The detection stage is a bounded efficient-score scan after a nonlinear
atmosphere-only fit. At each trial orbit, the molecular lunar template is
computed at that fitted atmospheric history. Its projection onto the local
atmospheric Jacobian is removed, retaining the nuisance correction in the full
spectral space. If q is that projected template and r the whitened residual,
c=<r,q>, d=||q||^2, f_hat=clip(c/d,0,1), and

T(P,phase)=max[2 f_hat c - f_hat^2 d,0].

The statistic is the maximum over the specified grid PLUS its fixed local
refinement procedure. It is NOT a chi-square variable, a Gaussian sigma, or
the globally maximized nonlinear likelihood ratio. The tangent nuisance space
is unbounded even when atmospheric fitting has active bounds; this convention
is reproduced in every calibration and validation draw.

The best candidates also receive joint nonlinear refinement of atmosphere,
blocking, period and phase, with the molecular lunar spectrum changing with
the atmosphere. That is a physical consistency and localization check. Its
likelihood improvement is NOT assigned the score search's false-alarm p-value.
Calibrating a full nonlinear global search would require bootstrapping that
entire different algorithm.

## Independent null and validation ensembles

Three explicitly specified null generators contain no moon:

1. Reference atmosphere with Gaussian noise.
2. Smooth two-sector atmosphere: independent normalized control coefficients
   uniform in [-0.3,0.3] for each physical parameter and sector.
3. The smooth atmosphere plus a signed uniform LOS-wind pulse: amplitude
   uniform 0.3--0.9 km/s, center uniform -1.6--1.6 h, Gaussian width uniform
   0.2--0.6 h. The pulse changes physical wind before nonlinear synthesis,
   remains inside the atmospheric bounds, and is outside the fitted low-order
   temporal family in general.

Calibration and held-out validation use distinct random streams. The first
pilot requests 128 calibration and 64 validation data sets PER generator,
plus 32 injections per comparison group. The nominal diagnostic false-alarm
level is 5%, not a journal-discovery threshold. The minimum finite Monte Carlo
tail estimate is 1/129; five-sigma claims are impossible with these ensembles.

Each family supplies a tail estimate (1 + count[T_null >= T_obs])/(N+1).
Use the MAXIMUM tail estimate across the three families for the candidate
flag. This guards against the represented null distributions rather than
depending on an unspecified mixture weight. It does NOT take a supremum over
all possible atmospheric histories or prove a fixed-nuisance-parameter error
bound. Report held-out false-alarm fractions with binomial intervals and also
show how a Gaussian-only threshold behaves on atmospheric controls.
Report both predeclared temporal models separately and a conservative consensus
that takes the maximum of their calibrated p-values. Do not choose whichever
temporal model gives the more favorable p-value after inspecting results.

## Recovery and limits

The main comparison injects the previous physical moon at an off-grid phase
0.7473 and its Keplerian period, without supplying either to the search.
Separate phase-uniform injections assess phase selection; distinguish periods
of maximum score from orbital confidence intervals. Report continuous period
errors, circular phase errors, boundaries and aliases rather than choosing a
post-hoc tolerance that labels favorable cases as recovered.

Compare an oracle fixed-orbit score and its own calibration for the showcase,
but never insert that orbit into the blind candidate selection. A larger best
score is expected after searching; only the search-maximum calibration accounts
for that selection. Injections and atmospheric ensembles using the accelerated
spectral grid are controlled-family tests; independent direct line-by-line
spectra validate the moving-shadow approximation before calibration. They do
not establish uniqueness against unrepresented atmospheric physics.
