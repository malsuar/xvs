# EXOVEL: first phase- and period-blind recovery pilot

## Direction and scope

We are making progress on a reproducible forward model and a tested detection
procedure, not yet demonstrating observational feasibility. The native-noise
benchmark remains too weak. The earlier bounded nonlinear recovery assumed
a known orbit and did not include the search penalty. The present experiment
addresses that gap and tests atmospheric false alarms with independent data.

The numerical precision is now adequate for this pilot. The dominant issues
are null-model adequacy, observing noise and phase coverage, not adding ever
denser geometry grids. A significant candidate, correct orbital recovery and
unique identification as a moon remain different claims. The original
population search and its selection were NOT changed or repeated.

The predeclared design and strategic review are in
`research_direction_blind_protocol.md`. This report adds the completed results.

## What is searched and calibrated

The host, radius of the moon, circular prograde aligned geometry, spin,
instrumental resolution, four-hour fully backed window and CO line model are
fixed. Lunar period and phase are unknown. Period determines semimajor axis
through Kepler's law for the fixed host, neglecting satellite mass; it is NOT
varied at fixed separation. The adopted fluid-Roche and 0.35-Hill screen gives
7.5849--86.0626 h (1.5101--7.6252 Rp), using lunar density 3 g/cm3. These are
instantaneous dynamical screens, not a tidal lifetime calculation.

The bank has 73 logarithmic periods and 1024 uniform phases (74752 trial
orbits), including orbits with no visible signal. The three best separated
grid seeds receive fixed local score refinement. No injected coordinates are
added. The showcase phase is 0.7473, off-grid, and period is 21.238292 h.

Each data set receives a new bounded nonlinear atmosphere-only fit, using
either three or five temporal controls (18 or 30 atmospheric parameters).
Moon spectra are evaluated at that fitted atmosphere for every orbit. The
efficient-score statistic projects out the local atmospheric Jacobian and
profiles a physical blocking fraction 0<=f<=1. The maximum of that specified
finite search/refinement procedure is calibrated empirically. The nuisance
tangent projection itself is unbounded, including when the null fit has an
active bound; the same convention is used in every simulation.

This is NOT a globally maximized nonlinear likelihood-ratio search. Selected
candidates also receive joint nonlinear atmosphere/blocking/period/phase fits,
but their likelihood improvements do not inherit the score calibration.
We report them as model-consistency and localization checks only.

## Physical and numerical implementation

An aligned projected lunar disk is represented by area-weighted polar-cell
overlaps and four spatial moments. The stellar intensity is expanded to
second order over the small occulted patch. All radial CO opacity, saturation,
rotation, line shifts and instrumental bins are retained through a separate
moving-shadow spectral library. Its 160 geometric modes and 97 spectral modes
accelerate arbitrary orbital phases without shifting/scaling a single
fixed-orbit time map. Independent ray/cell integration and direct line-by-line
spectra validate the approximation before calibration.

Five production comparisons spanning the period interval give relative
spectral errors 0.012--0.056%, or at most 0.00271 control-noise units over a
complete map. Sixteen Gaussian nodes integrate each five-minute exposure.
Moving-atmosphere spectral derivatives agree with numerical derivatives to
about 2e-9. The denser 97x1536 search changes the tested maxima by at most
0.260%; this is a finite audit, not a proof of a continuum global maximum.

One noise example changes its best period from 7.626 to 10.042 h while the
score changes only 0.260%. Stable evidence can coexist with unstable orbital
coordinates. This is why score convergence is not a period-error estimate.

There were 228 seed refinements without reported convergence among 4224
Monte Carlo seed refinements (5.40%). No trial was removed: the best computed
score and a grid fallback were retained, with the identical rule in calibration
and validation. This limits claims about an exact continuous maximum, not
the definition of the calibrated finite algorithm. The two nonlinear null
starts agree to at most 1.14e-6 in chi2, and all 60 selected-candidate nonlinear
starts report convergence. Unit tests check phase periodicity, sparse/full
spectral projection, bound profiling and separation of calibration/validation.

## Ensembles and decision rule

The noise is the preceding idealized 11.39465 ppm per resolution bin. Lunar
signals are never amplified. This is NOT an instrument forecast and does not
preserve an independently fixed systematic floor. No native-noise blind
false-alarm calibration is claimed here.

Three no-moon generators are used: reference atmosphere plus Gaussian noise;
bounded smooth two-sector atmospheric histories; and those histories plus
signed wind pulses with random amplitude 0.3--0.9 km/s, center -1.6--1.6 h and
Gaussian width 0.2--0.6 h. Temperature, abundance and wind affect actual
nonlinear spectra. Pulses generally lie outside the fitted temporal family.
The protocol gives the precise distributions; they are not population priors.

There are 384 calibration null data sets (128 per generator), 192 independent
validation nulls (64 per generator), and 128 injected systems (32 per group).
The same data are evaluated by both temporal models, giving 1408 blind
searches. Models are paired, not independent duplicate observations. Each
calibration, validation and injection stream has distinct random seeds.

For each generator, p_MC=(1+count[T_null>=T_obs])/(N+1). The atmospheric
decision uses the maximum p across generators. A separate conservative
consensus uses the maximum across both temporal models as well; neither the
better model nor a smaller p is selected after inspecting a candidate.

The diagnostic level is 5%. The minimum reported p is 1/129=0.007752, not
zero. These ensembles cannot establish a five-sigma discovery threshold.
Calibration applies to the specified generative distributions, NOT every
fixed physical atmospheric history or all possible contaminants.

## False alarms on independent no-moon data

| Validation generator | 3 controls, Gaussian-only calibration | 3 controls, atmospheric calibration | 5 controls, Gaussian-only calibration | 5 controls, atmospheric calibration | Consensus |
| --- | ---: | ---: | ---: | ---: | ---: |
| Reference | 1/64 | 0/64 | 5/64 | 0/64 | 0/64 |
| Smooth atmosphere | 2/64 | 0/64 | 3/64 | 1/64 | 0/64 |
| Wind pulse | 35/64 | 4/64 | 13/64 | 2/64 | 2/64 |

For wind controls, Gaussian-only calibration gives 54.69% false candidates
with three controls and 20.31% with five controls. Atmospheric calibration
reduces these to 6.25% and 3.13%. Their exact 95% binomial intervals are
1.73--15.24% and 0.38--10.84%, respectively; the small validation set does not
verify a precisely known 5% rate. Consensus also gives 2/64 wind false alarms.
Zero counts in other groups have a 95% upper bound of 5.60%, not zero risk.

The strict critical sqrt(T) levels implied by the finite Monte Carlo rule are:

| Fitted temporal controls | Gaussian-only global threshold | Worst represented atmospheric global threshold |
| --- | ---: | ---: |
| 3 | 4.0804 | 5.4168 |
| 5 | 4.0465 | 4.8111 |

These are score thresholds, NOT Gaussian sigmas. A more flexible temporal
model loses some lunar information but can substantially reduce atmospheric
null tails. Maximizing S/N alone is therefore not an adequate model-selection
criterion.

## Recovery and the unknown-orbit penalty

For 32 reference-atmosphere injections at the favorable off-grid phase, a
fixed-orbit oracle flags 32/32 with atmospheric calibration, whereas the blind
search flags 24/32 for each temporal model. Consensus flags 22/32. The oracle
is only a benchmark: its known coordinates are never used to seed the search.

| Injected family | Atmospheric flags, 3 controls | Atmospheric flags, 5 controls | Consensus flags |
| --- | ---: | ---: | ---: |
| Reference, fixed phase | 24/32 | 24/32 | 22/32 |
| Smooth atmosphere, fixed phase | 24/32 | 25/32 | 24/32 |
| Wind pulse + moon, fixed phase | 23/32 | 23/32 | 20/32 |
| Smooth atmosphere, uniform phase | 4/32 | 7/32 | 4/32 |

These are candidate-flag fractions, not probabilities of recovering the
correct orbit. Consensus fixed-phase smooth recovery is 75.0% (95% interval
56.6--88.5%), versus 12.5% (3.5--29.0%) for uniform phase. Eleven of the
32 random-phase moons have zero occulted atmospheric signal in the available
window. They are retained in all primary rates; none of those eleven was
consensus-flagged. This is a phase-coverage penalty, not an arbitrary cut.

For reference fixed-phase injections, period maxima have q16--q50--q84 of
20.04--20.94--21.93 h (three controls) and 20.16--21.12--22.38 h (five
controls), around the injected 21.2383 h. These are Monte Carlo estimator
distributions in a favorable fixed geometry, not confidence intervals for a
particular observation or a general exomoon period precision.

With wind pulses, candidate fractions are similar between temporal models,
but orbital recovery is not: the all-injection period q16--q84 expands to
9.44--52.21 h with three controls, versus 18.05--22.43 h with five. Among
consensus-flagged wind injections, the 84th percentile of the absolute log
period disagreement between models is 1.02, a factor about 2.8. Consensus
evidence is not consensus about an orbit. Random-phase cases show strong
aliases; reporting just the maximum would conceal this lack of identification.

## Noiseless and nonlinear checks

The independent direct line-by-line off-grid injection gives a score maximum
near P=21.2395 h, phase=0.7472996; subsequent joint nonlinear fitting returns
P=21.23816 h, phase=0.7473001. This is a numerical recovery control, NOT an
observational precision claim. The atmosphere-plus-moon residual norm is
only 0.00284 control-noise units. All selected nonlinear starts converge.

For this example, the oracle p reaches the Monte Carlo floor, but the blind
atmospheric p values are 0.04651 (three controls) and 0.03101 (five controls),
with consensus 0.04651. This is at most a candidate at the diagnostic 5%
level in an idealized experiment, emphatically not a five-sigma detection.

The illustrative random wind-pulse-plus-moon draw has p=0.08527 with three
controls and the Monte Carlo floor with five, hence consensus p=0.08527.
Choosing the more favorable temporal model after seeing the candidate would
change the decision. Its nonlinear residuals remain appreciable; the higher
likelihood score is not evidence of a uniquely adequate moon model.

The illustrative draws use different nuisance histories and should not be
subtracted as paired atmospheric experiments. The noise ensembles and the
explicit oracle comparisons, rather than those individual drawings, define
the rate comparisons. No p-value in this report is transferred to the
uncalibrated nonlinear-fit likelihood ratio.

## What this means for the manuscript and next work

We can support a conditional spectroscopic diagnostic and a reproducible
finite-search test. We cannot identify S/N>=5 with five-sigma discovery, use
the highest-scoring system to define a favorable population, or treat one
best period as an orbital measurement. Atmospheric calibration is essential:
Gaussian-only controls substantially underestimate false alarms here.

Next priorities are predictive atmospheric adequacy, realistic stellar/telluric
contamination and a photon/instrument budget for concrete targets. A fully
nonlinear blind likelihood-ratio search would need its OWN end-to-end
calibration. Inclination, eccentricity, radius, spin and host uncertainty are
still fixed. The current model does not establish uniqueness against arbitrary
localized atmospheric patches. Do not repeat the population sweep before
these statistical and observational assumptions are settled.

For population work, the appropriate objective is recovery probability at a
controlled false-alarm level, marginalized over unknown phase, not maximum
S/N at a favorable phase. This experiment establishes the infrastructure and
illustrates why that distinction matters; it is not a final selection function.

## Reproduction and assets

See the README command sequence. `build_blind_co_library.py` requires the
previous nonlinear emulator. `validate_blind_co.py` gates the main run.
`blind_co_recovery.py --tag pilot` resumes only identical procedure settings;
raw trials are checkpointed. `audit_blind_grid.py` audits resolution and
derivatives independently. `summarize_blind_comparisons.py` and
`plot_blind_co.py` export comparisons and figures without rerunning trials.

Tracked outputs include `co_blind_pilot_calibration.csv`, `decisions.csv`,
`rates.csv`, `orbit_errors.csv`, `examples.csv`, `summary.json`,
`candidate_fits.json`, `oracle_comparison.csv`, `model_agreement.csv`,
`phase_coverage.csv` and `extra_checks.json` (all with the co_blind_pilot_
prefix). Production and grid checks are in `co_blind_forward_validation.json`
and `co_blind_grid_audit.json`.

Regenerable figures in outputs/:

- `co_blind_pilot_calibration.png/pdf`: fixed-orbit and search-maximum null
  tails, including atmospheric controls and the independent injection.
- `co_blind_pilot_validation.png/pdf`: held-out false alarms and injection
  candidate fractions with binomial intervals.
- `co_blind_pilot_orbits.png/pdf`: score maps and noisy orbital maxima;
  these maps are not posterior distributions.
- `co_blind_pilot_orbital_errors.png/pdf`: period/phase errors, retaining
  unflagged and invisible injections; filled markers denote consensus flags.

The molecular data provenance and separate ExoMol data license remain in the
CO data manifest and preceding reports. Existing manuscript assets, original
population samples and earlier experiment results were preserved.
