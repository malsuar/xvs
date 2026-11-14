# Bounded nonlinear CO recovery with temporal coherence

## Scientific question

The previous linear nuisance experiment allowed an unrestricted coefficient
for every spectral mode at every exposure. It showed substantial degeneracy
with two-limb atmospheric structure, but did not establish that a bounded,
temporally coherent atmosphere could realize those spectral combinations.
Here we fit the actual nonlinear CO spectral family and compare an
atmosphere-only hypothesis with the same atmosphere plus a known-orbit moon.

This remains a controlled, single-system experiment. Its parameter bounds and
temporal basis are stated hypotheses, not independently measured constraints
or universal restrictions on giant-planet atmospheres.

## Reproduction

```sh
python -m pip install -e '.[molecular]'
python scripts/fetch_co_template_data.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/nonlinear_co_recovery.py --build-only
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/nonlinear_co_recovery.py --n-mc 16
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/audit_nonlinear_co.py
```

The molecular source remains ExoMol 12C16O Li2015, with the provenance and
separate data license in `data/exomol/CO_Li2015/manifest.json` and its README.
The earlier Gaussian, molecular and linear-nuisance pilots are preserved.
Relevant algorithm documentation: bounded least squares,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html;
interpolating B-splines,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_interp_spline.html.

## Physical forward model and bounds

The host, moon orbit, stellar weighting, wavelength window, R=100000 and
48 five-minute exposures are unchanged from the preceding molecular pilot.
The two sectors are the negative- and positive-rotation-velocity halves of
the projected annulus. Temperature, CO mixing ratio and effective LOS wind
are independent in each sector. Their bounded ranges are:

| Parameter | Allowed range | Meaning |
| --- | --- | --- |
| T | 560--720 K | LTE excitation, thermal width, pressure and scale height change |
| X_CO | 1e-4/3--3e-4 | Abundance enters the slant optical depth before saturation |
| Wind | -1.5--+1.5 km/s | Additive local signed LOS velocity |
| Lunar blocking f | 0--1 | Transparent-to-opaque known-radius lunar mask |

CO abundance is prescribed, not derived from chemical equilibrium. The wind
is an effective velocity, not a circulation solution. Fixed Rp, P0=1e-3 bar,
gravity, mean molecular weight and isothermal radial structure follow the
preceding experiment. A 16H(640 K) computational annulus with 52 radial cells
extends the domain to accommodate the hotter/more abundant permitted states.
There is no free gray line-depth multiplier and no unconstrained signed
spectral amplitude: atmospheric absorption is determined by the opacity and
geometry. Consequently the comparison with the earlier tangent-space fit
changes both the physical family and its temporal freedom, not only the
optimizer.

If B(theta) is the moonless atmosphere and S(theta) the lunar residual, the
fitted differential prediction is

Y_model = B(theta(t)) - B(theta0) + f S(theta(t)).

The lunar spectrum is recomputed/interpolated with the SAME atmosphere as B;
it is not a fixed spectral residual added to an unrelated atmosphere. The
union of the original blocking masks supplies S. With fixed theta, f between
zero and one interpolates between the unblocked and fully blocked spectra.
Opaque geometric blocking has f=1; f is a known-trajectory recovery diagnostic,
not a retrieved moon radius or a proof of a partially transparent real moon.

## Coherence in time

For each sector and parameter, normalized values are represented by clamped,
nonnegative B-spline basis functions over [-2,2] h. Control coefficients lie
in [-1,1], and partition of unity guarantees the entire trajectory remains
within the parameter bounds, not just the exposure centers. We compare:

- One control per parameter: constant atmosphere, six parameters.
- Three controls: a quadratic Bernstein trajectory, eighteen parameters.
- Five controls: a C2 cubic spline with an interior knot at zero,
  thirty parameters.

These families disallow independent five-minute atmospheric changes, but do
not derive their coherence from dynamics, heat transport or chemistry. A
shorter-lived, localized atmospheric perturbation could fall outside them.
There is no extra derivative penalty or hidden probability prior. States
are frozen within each five-minute exposure after integrating the geometric
kernels; a direct 2.5-minute calculation tests that approximation separately.

## Spectral acceleration and independent checks

Direct line-by-line spectra train a tensor-product nonlinear interpolation
grid in normalized T, ln(X_CO) and wind, with 5x7x9 nodes. Both the reference
and missing-light spectra are retained separately for each sector and time.
The local line shapes are Voigt profiles with saturation, hydrostatic slant
columns, instrumental convolution and exact smooth shifted-bin integration.

A common 128-column orthonormal wavelength basis compresses these spectra.
Its range is built from normalized physical reference and lunar snapshots,
not from the observed noise. The same basis is used for data and every
hypothesis. Independent Gaussian coefficients then have the exact projected
noise covariance; omitted orthogonal noise contributes an identical constant
to chi2 differences only to the extent that model compression is accurate.

Independent non-grid physical spectra test interpolation AND compression in
units of the control noise. Inference is blocked if their full spectral error
norm exceeds 0.1. The check is not a substitute for an atmospheric-model
uncertainty budget. Key fitted trajectories are also re-evaluated directly,
without the grid or spectral basis. The cache key records source-code/data
fingerprints and grid/basis choices. Cached nodes are regenerable outputs.

Analytic spline/interpolation derivatives drive SciPy bounded trust-region
least squares. Tests compare the Jacobian with numerical derivatives and
verify bound preservation, moon/atmosphere coupling and nested recovery on
an exactly specified nonlinear toy spectrum. Multiple starting atmospheres,
optimizer status, active bounds and objective spreads are retained for every
noiseless fit; a smaller multi-start set is used in the noise pilot. Multi-start
agreement is a local robustness check, not a proof of the global optimum.

## Injections and fair comparison

The native opaque moon is NEVER amplified in this experiment. We use the
original independent noise (about 1490.75 ppm per observed bin) and, separately,
an idealized noise divided by 130.8285 (about 11.39 ppm). The latter matches
the earlier depth-profiled laboratory information scale while keeping the
absorption and blocking physically bounded. ALL atmospheric contaminants and
the lunar signal see the same reduced noise. This is not an instrument
forecast: the original systematic floor is not independently held fixed and
no collecting-area, throughput or telluric model justifies this noise level.

Line-by-line injections, not grid interpolation, supply these cases:

- Native moon, reference atmosphere, original noise.
- Same physical moon and atmosphere, idealized control noise.
- A held-out smooth, bounded two-sector atmosphere without a moon.
- The same atmosphere with an opaque moon.
- The smooth atmosphere plus a nonlinear LOS wind pulse without a moon.
- The wind-pulse atmosphere with a moon.

The smooth atmosphere is a three-control trajectory with different T, CO
and wind histories in each sector. Its values are not training-grid nodes.
It is a well-specified recovery control, not a test against arbitrary
atmospheric histories. The additional wind pulse has amplitude 0.9 km/s,
center -1.1 h and Gaussian sigma 0.35 h, overlapping one lunar contact. It
changes the actual wind before spectral synthesis, unlike the preceding
experiment's mixture of endpoint spectra. The pulse is smooth and bounded
but outside the fitted low-order trajectory family; it explicitly tests
temporal-model misspecification. No circulation model is claimed for it.

Noiseless (Asimov) comparisons measure the distance between model families.
The paired Gaussian-noise pilot refits both hypotheses for each draw, rather
than reusing a linear-filter variance. Known orbit/phase, noiseless reference,
no telluric/stellar lines and independent Gaussian noise remain assumptions.

## Interpretation of recovery statistics

We report fitted f, chi2 for both hypotheses and sqrt(delta chi2). This last
quantity is a model-separation diagnostic, NOT a calibrated Gaussian sigma:
f=0 is a boundary, f=1 is another boundary for the injected opaque moon,
nonlinear nuisance parameters can be weakly identified, and the trajectory
has not been searched. The small ensemble cannot determine a discovery
false-alarm rate. Boundary pile-up in recovered f is expected and must not
be described as an unbiased amplitude measurement.

Profile fits at fixed f test the shape of the likelihood for selected cases.
Local Fisher projections using the SAME temporal controls provide an
unbounded tangent-space comparison at the reference atmosphere and f=0.
Disagreement between these controls and
the nonlinear fit may reflect physical bounds, finite curvature, or failed
optimization; it is not automatically evidence for a unique lunar signature.

## Results of the completed pilot

The following are NOISELESS model separations at the idealized control noise,
not detection significances. The atmosphere-only and atmosphere-plus-moon
hypotheses have the same temporal freedom in each column.

| Injected system | Constant atmosphere (6 parameters) | Three controls (18 parameters) | Five controls (30 parameters) |
| --- | ---: | ---: | ---: |
| Physical moon, reference atmosphere | 5.649 | 5.502 | 5.131 |
| Smooth atmosphere, no moon | 3.509 | 0.000 | 0.000 |
| Smooth atmosphere + moon | 8.290 | 5.253 | 4.905 |
| Smooth atmosphere + wind pulse, no moon | 7.581 | 2.429 | 0.877 |
| Smooth atmosphere + wind pulse + moon | 10.861 | 7.267 | 5.727 |

Entries are sqrt(delta chi2). With the ORIGINAL noise, the moon-only entries
are only 0.04318, 0.04205 and 0.03922. The original system remains undetectable;
reducing the noise for a laboratory control does not change that conclusion.

The much smaller information loss than in the previous per-exposure fit is
conditional on the new temporal/physical family. In the moon-only case, local
Fisher controls predict 5.6496, 5.5022 and 5.1329, agreeing with the nonlinear
fits within 0.035%. For the three-control atmosphere-only fit, maximum control
offsets are about 1.05 K, 0.0238 in ln(X_CO), and 0.036 km/s. They lie far from
the adopted atmospheric bounds. Thus, for THIS well-specified control, the
restored separation is not a consequence of clipping the atmospheric fit at
its boundaries. Coherence and removal of unrestricted spectral freedom matter
more than finite nonlinear curvature. A wider-domain global search was not
performed, so this is not a universal exclusion of other physical atmospheres.

Model misspecification is visible rather than hidden. A constant atmosphere
fitted to the smooth no-moon injection favors f=0.617 with separation 3.509,
but leaves large noiseless residuals (joint chi2=666.4). The three- and
five-control models correctly prefer f=0 for this matched-family control.

The out-of-family wind pulse remains partially confounding. With no moon,
it produces f=0.463 and separation 2.429 for three controls, falling to f=0.179
and separation 0.877 with five controls. Even after adding a moon, the
three-control pulse fit leaves joint chi2=196.1 without a real moon and 198.8
with one. These are residual model distances, not calibrated goodness-of-fit
p-values. A preference for the moon is not enough: both hypotheses can omit
important atmospheric structure. The apparent separation 7.267 for a real moon
plus pulse must not be interpreted as clean lunar information or an improved
instrumental sensitivity.

For the primary three-control, smooth-atmosphere experiment, sixteen paired
noise draws give:

| Injection | sqrt(delta chi2), q16--q50--q84 | Recovered f, q16--q50--q84 |
| --- | --- | --- |
| No moon | 0.00005--0.302--0.797 | 0.000--0.057--0.152 |
| Opaque moon | 5.047--5.546--5.992 | 0.963--1.000--1.000 |

None of the sixteen no-moon controls exceeds sqrt(delta chi2)=3 (maximum
1.598), whereas all sixteen moon controls do. This is a descriptive pilot,
NOT a 100% recovery probability or a calibrated false-positive rate. Eleven
of sixteen moon fits pile up at f=1, illustrating the physical upper boundary;
the fitted f distribution is not an unbiased Gaussian amplitude estimator.

All recorded optimizer starts terminate successfully. Objective spreads are
below about 6.4e-8 for the noiseless joint fits and 1.6e-8 for the noise fits;
null-fit spreads are smaller. Fixed-f profile starts also succeed. Multiple
starts support local robustness but do not prove the global optimum.

## Independent numerical validation

The maximum non-grid spectral approximation error is 0.07484 control-noise
units over the complete map, with pointwise errors below 0.019 ppm in these
checks. Re-evaluating key fitted trajectories directly gives approximation
errors 0.00045--0.01436 in the same whole-map noise norm. In particular:

- Moon-only atmosphere-only fit: direct chi2=30.268531 versus30.268499 on
  the interpolated/compressed model.
- Wind-pulse atmosphere-only fit: direct chi2=202.048148 versus202.046614.
- Matched atmosphere-plus-moon fit: direct chi2 is below 1e-7.

The finer exposure treatment changes the matched atmospheric-plus-moon map
by 0.00231 control-noise units. Extending the annulus 16H to 20H changes the
reference by at most 0.00543 and the lunar residual by less than 0.00010 control-noise
units in the reference and upper-T/CO checks. The original molecular pilot's
saved reference and residual are reproduced EXACTLY with the default API.

The cache fingerprint includes the geometric kernels, grid/basis choices,
data manifest and generation code. During finalization the already verified
nodes were relabeled under this expanded fingerprint; no spectra or fit
results changed. The independent validation is regenerated under the final
key recorded in the summary. Raw/derived ExoMol spectral outputs retain their
separate data-license attribution.

## Implication for the referee response

The previously severe degeneracy is substantially reduced when fitting this
bounded nonlinear, temporally coherent two-sector family. It does not disappear
in general: underfitting atmospheric time dependence can create spurious lunar
preference or inflate the apparent separation. The appropriate statement is
conditional distinguishability with explicit temporal/model assumptions, not
an inability of planetary atmospheres to mimic EXOVEL. The noise-reduced,
known-orbit recovery is a method validation; the native system is not a viable
detection under the original observing assumptions.

## Outputs

- `co_nonlinear_emulator_checks.json`: independent non-grid approximation checks.
- `co_nonlinear_recovery.csv`: noiseless separation, blocking, bounds and starts.
- `co_nonlinear_noise_recovery.csv`: paired nonlinear noise refits.
- `co_nonlinear_profiles.csv`: fixed-blocking profile objectives.
- `co_nonlinear_summary.json`: bounds, coherence, noise, direct-fit and exposure checks.
- `co_nonlinear_fits.json`: full fitted coefficients and every starting-point result.
- `co_nonlinear_geometry_checks.json`: outer-boundary and backward-compatibility checks.
- `outputs/co_nonlinear_recovery.png/pdf`: recovery and atmospheric trajectories.
- `outputs/co_nonlinear_profiles.png/pdf`: fixed-blocking profiles and residual
  model distances, including the misspecified wind-pulse case.
- `outputs/co_nonlinear_trajectories.png/pdf`: temperature, abundance and wind
  histories against the declared bounds and independent injection.
- `outputs/co_nonlinear_emulator.npz`: regenerable spectral grid and orthonormal basis.
- `outputs/co_nonlinear_truth.npz`: independent physical atmospheric injections.

The scientific conclusions must remain conditional on these bounds, temporal
bases and two-sector geometry. Localized clouds, opacity or velocity patches,
non-isothermal/multi-species atmospheres, blind orbital recovery and correlated
observational contaminants are not excluded by success in this experiment.
