# Atmospheric alternatives to a known-trajectory lunar CO signal

## Scope and reproduction

This is a follow-up to `co_template_injection_recovery.md`, not a population
search or a demonstration of observational feasibility. The native lunar signal
in that experiment had S/N well below unity. A shape may be mathematically
distinct even when it cannot be measured at the adopted noise level.

```sh
python -m pip install -e '.[molecular]'
python scripts/fetch_co_template_data.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/atmospheric_nuisance_recovery.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/audit_atmospheric_nuisance.py
```

Source spectroscopy remains ExoMol 12C16O Li2015:
https://exomol.com/data/molecules/CO/12C-16O/Li2015/.
Data provenance, units, line-strength equations and separate CC BY-SA4.0 data
license are described in the preceding report and data manifest. The default
previous synthesis was checked to reproduce its saved reference and lunar
maps exactly after extracting the reusable geometric kernels.
With identical sector parameters, splitting the annulus into two velocity
halves preserves the full reference to floating-point precision (maximum
absolute difference 4.1e-20 in the small regression control) and the lunar map
exactly. Sector splitting does not introduce an additional flux normalization.

## What changes physically

The planet, moon, stellar continuum, instrumental resolution and observing
sequence are unchanged: 0.7 Mj, 1.1 Rj, equatorial rotation 12 km/s; a moon of
0.09 Rp at 3 Rp with period 21.23829 h and phase .75; R=100000; 48 five-minute
exposures in a four-hour fully backed central window. Independent noise is
1490.75 ppm per observed resolution bin. There is no extra line-count gain.

The reference is a 640 K, isothermal CO atmosphere with X_CO=1e-4 and pressure
1e-3 bar at the fixed opaque radius. The new domain extends to 12 H(640 K)
with 40 radial cells, versus 10H/32 in the preceding pilot. No opaque-radius
adjustment or arbitrary peak renormalization is made when parameters change.

Temperature changes LTE strengths and partition function, thermal broadening,
hydrostatic scale height, radial pressure and tangent column. CO abundance
changes the column BEFORE applying 1-exp(-tau), preserving saturation. Wind
adds a signed line-of-sight velocity to each local rotational velocity before
instrumental convolution and binning. These are prescribed effective winds,
not a circulation solution or a mass-conserving three-dimensional wind field.

Two types of perturbation are tested for each parameter:

- A uniform change across the projected annulus.
- An antisymmetric change between the negative- and positive-rotation-velocity
  halves: theta_minus=theta0-delta, theta_plus=theta0+delta. For CO the
  antisymmetry is in ln(X_CO), not linear abundance. These are fixed spatial
  sectors, NOT sectors moving along the lunar trajectory.

Positive wind means addition to the model's signed velocity coordinate. The
sector labels refer to that coordinate; they do not introduce a new sign
convention for observed radial velocity.

## Fitted nuisance models

Let B(t,nu) be the moonless spectrum and S(t,nu) its native lunar residual.
At every exposure we fit an independent linear combination of spectral modes.
The model family comprises depth; depth plus each uniform parameter separately;
depth plus all uniform derivatives; depth plus all uniform and antisymmetric
derivatives; and the last model plus separate single-mode curvatures.

The tangent modes are central finite differences at the reference atmosphere.
Nominal steps are 10 K, 0.1 in ln(X_CO) and 0.2 km/s; reported fits use half
these steps, with the nominal steps retained for a derivative-convergence test.
Curvatures are (B_plus+B_minus-2B)/step^2. Mixed second derivatives are NOT
included. The curvature model is an expanded linear spectral basis, not a
fully nonlinear fit of temperature, abundance and wind.

Small wind curvatures exposed a problem with piecewise-linear Doppler shifting:
the numerical second derivative changed about 20% when its step was halved.
This experiment instead integrates a cubic-spline antiderivative of the
instrument-broadened local profile directly between Doppler-shifted detector
bin edges. This is both smooth and avoids resampling an entire fine spectrum
for every surface velocity. The previous pilot retains its linear default.
Exact polynomial-bin and wind-curvature tests check the new integration. The
experiment records the change in the native lunar map between the methods.

Each mode is unit-normalized before an exposure-by-exposure singular-value
decomposition. The relative rank tolerance is 1e-8. Zero singular modes are
discarded. Coefficients are unrestricted and have no temporal smoothness or
physical priors. Such freedom is a conservative local spectral-degeneracy
test, NOT a claim that atmospheric parameters can change arbitrarily between
five-minute exposures. Large fitted coefficients need not describe any physical
atmosphere, especially in artificially amplified injection controls.

With P_N the orthogonal projector onto the nuisance modes, the retained lunar
component is Q=(I-P_N)S. Under the adopted constant independent noise sigma,
expected known-trajectory S/N=||Q||/sigma, and the signed recovery score is
Z=<Y,Q>/(sigma ||Q||). All exposures and spectral bins enter the inner product.
There is no phase/orbit search or look-elsewhere correction in this test.

## Noise and finite-amplitude controls

2000 actual Gaussian noise maps, seed 63118, are projected onto all retained
lunar templates. The same noise is used for paired comparisons. A single
fixed laboratory lunar multiplier is chosen to give expected Z=5 with only
depth fitted, then held fixed as the atmospheric model gains flexibility.
Separate multipliers needed for Z=5 are diagnostic numbers, not physically
validated moons or instrument improvements. Null draws cannot establish a
five-sigma false-alarm rate with this sample size.

Held-out finite alternatives are T=695 K, X_CO=2.3e-4, wind=1.3 km/s,
antisymmetric T offsets +/-35 K, ln(CO) offsets +/-0.45, and wind offsets
+/-0.7 km/s. Combined cases use uniform (675 K, 1.6e-4, 0.6 km/s), and
sector values (615/665 K, X0 exp[-.3/+.3], -0.4/+0.4 km/s). These points
are not the finite-difference training endpoints.

Each alternative is injected WITHOUT a moon, first as a static atmospheric
change and then as a two-column spectral mixture with covering fraction
f(t)=exp[-0.5((t+1.1 h)/0.35 h)^2]. This timing deliberately overlaps one
lunar contact and is a stress test, not a probability model for variability.
The mixture residual f(t)[B_alt-B] is NOT the same as evaluating a uniform
atmosphere at an interpolated T(t), X(t) or wind(t). No dynamical explanation
of that covering-fraction evolution is assumed.

For each finite control we report residual S/N, postfit residual S/N, signed
overlap with Q, and expected false lunar Z. We also report the expected score
when that control contaminates the same laboratory lunar injection. A small
false score at this noise level does NOT prove that its shape is orthogonal;
the overlap and residual norms must also be considered. Negative scores bias
or suppress recovery rather than trigger a positive lunar detection.

## Results of the specified pilot

The full native lunar S/N is 0.04431. The table below reports the information
remaining after each nuisance fit. The last column uses the SAME amplified
lunar residual, with multiplier 130.8285, throughout; it is not retuned for
each row.

| Nuisance model | Native S/N | S/N retained relative to depth-only | Expected laboratory Z |
| --- | ---: | ---: | ---: |
| Depth only | 0.038218 | 100% | 5.000 |
| Depth + uniform temperature | 0.038211 | 99.98% | 4.999 |
| Depth + uniform CO abundance | 0.038207 | 99.97% | 4.999 |
| Depth + uniform wind | 0.030989 | 81.08% | 4.054 |
| Depth + all uniform changes | 0.030935 | 80.94% | 4.047 |
| Depth + uniform and two-limb changes | 0.013458 | 35.21% | 1.761 |
| Previous row + single-mode curvatures | 0.012124 | 31.72% | 1.586 |

Uniform temperature and abundance changes add little degeneracy to an already
fitted depth. Uniform wind matters more, and antisymmetric structure matters
substantially. These percentages refer to retained S/N (a spectral norm), not
retained atmospheric area or the fraction of blocked photons. The two-limb
linear fit removes about 87.6% of the squared depth-profiled spectral norm.

The ranks are 1, 2, 2, 2, 4, 7 and 10 across the rows, constant in time. In
this additive two-sector construction, the antisymmetric single-parameter
curvatures duplicate their uniform counterparts; there are not thirteen
independent modes. In 2000 paired realizations, the same laboratory moon
exceeds signed Z=3 in 97.5%, 84.15%, 11.45% and 8.5% of draws for depth-only,
all-uniform, two-limb-linear and curvature models, respectively. These are
controlled recovery fractions at a chosen threshold, not calibrated discovery
rates for an unknown-orbit search. Individual null-score standard deviations
range from 1.016 to 1.030, compatible with the fixed-linear-filter control.
Explicit nested least-squares examples satisfy delta chi2=Z^2 in every model.

No held-out finite atmospheric control reaches a significant lunar score at
the native noise level. Across both static and transient tests, maximum
absolute false scores are 0.0700 (depth), 0.0828 (all uniform), 0.0219
(two-limb linear), and 0.00157 (curvatures). This does NOT establish uniqueness:
the atmospheric residuals themselves are weak under these noise assumptions.
For example, the static +1.3 km/s wind leaves a false score 0.02192 under the
two-limb linear fit, equivalent to an expected fitted multiplier 1.63 of the
NATIVE lunar template, despite being observationally insignificant. An
instrumental noise reduction would strengthen atmospheric contaminants too;
the moon-only laboratory amplification does not simulate that situation.

After fitting, some finite controls have signed spectral overlap around
0.20--0.25 with the lunar component. Curvature fitting can remove most of
their amplitudes without forcing the normalized overlap to vanish. Both
overlap and absolute residual strength are therefore reported, not just the
absence of high-Z false positives.

Halving derivative steps changes each reported S/N by at most 0.0844% and the
retained-component maps by at most 0.230%. The wind curvature itself now
changes by about 0.118%, versus about 20% with the earlier nonsmooth shift.
Switching spectral integration changes the native lunar map by 0.0684%.
Extending the domain from 12H to 14H changes the hotter and more abundant
finite atmospheric residual maps by 0.0527% and 0.0136%, respectively. The
separate joint-refinement audit is recorded in co_atmospheric_convergence.json.
Doubling angular, radial and temporal sampling, halving the fine spectral and
velocity-channel steps, and doubling the line-wing cutoff changes all retained
S/N values by 0.048--0.082%. Projected spectral maps change by at most 0.126%.
Ranks are unchanged. Thus the large loss of information after fitting sector
inhomogeneities is not explained by these numerical sampling choices.

Suggested interpretation for the manuscript/referee: the lunar residual is
not equivalent to a simple line-depth modulation, but its recoverable
information is sensitive to atmospheric velocity structure and hemispheric
inhomogeneity. The tested atmospheric modes do not fully reproduce its spectral
shape, yet substantially weaken its recovery when freely fitted. This selected
native configuration remains undetectable; these experiments do not establish
that all planetary atmospheric alternatives can be excluded.

## Output inventory

- `co_atmospheric_recovery.csv`: native and laboratory lunar recovery for each
  nuisance model, null moments, ranks and finite-step sensitivity.
- `co_atmospheric_false_positives.csv`: static/transient finite controls,
  atmospheric residual norms, false scores and contaminated lunar scores.
- `co_atmospheric_summary.json`: assumptions, finite-difference changes and
  outer-boundary checks for the hotter and more abundant alternatives.
- `co_atmospheric_convergence.json`: independent joint spatial/time/spectral/
  velocity-channel/wing refinement of every projected lunar component.
- `outputs/co_atmospheric_degeneracy.png/pdf`: retained lunar information,
  time dependence, paired recovery distributions and finite-control overlaps.
- `outputs/co_atmospheric_line_shapes.png/pdf`: a native spectral line and
  its atmospheric approximations, alongside the remaining lunar component.
- `outputs/co_atmospheric_convergence.png/pdf`: joint-refinement changes in
  retained S/N and spectral maps.
- `outputs/co_atmospheric_nuisance.npz`: native maps, retained components,
  null scores and finite transient maps, with the original physical amplitude.

## Limits for the referee response

This can test whether these particular uniform/two-sector CO modes remove
lunar information and whether specified atmospheric changes bias a fixed-orbit
score. It cannot establish a unique lunar interpretation against arbitrary
localized clouds, abundance/temperature structure, or a velocity patch that
follows the same projected trajectory. Such flexible local atmospheric
suppression can be intrinsically degenerate with geometrical removal of the
same contributions. External constraints and multi-line/multi-species or
repeated-event information would be needed to break those degeneracies.

Other exclusions remain: non-isothermal radiative transfer, equilibrium
chemistry, multi-species spectra, ingress/egress in this pilot, telluric and
stellar spectral contamination, correlated noise, trajectory uncertainty and
blind-search false-alarm calibration. Do not repeat the full sweet-spot search
or claim that static planetary effects are universally distinguishable on the
basis of this one selected system.
