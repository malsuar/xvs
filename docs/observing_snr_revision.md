# Observing-bin S/N revision

detectability_report now defaults to method="observed". Explicit method="legacy"
retains the historical sum over numerical grid points and heuristic smearing.
Historical search scripts must request legacy to reproduce baseline S/N; new
default runs intentionally produce revised estimates, not the archived values.

The observing operator integrates a piecewise-linear time series over complete
exposures anchored at the start of the model window. It convolves the residual
with a unit-area Gaussian LSF of FWHM c/R and averages over nonoverlapping bins
of width c/R. Spectral integration uses trapezoidal quadrature of the input
grid and analytic Gaussian bin probabilities. Residual outside the supplied
velocity window is assumed zero; provide generous spectral margins. Incomplete
final time and velocity bins are discarded. The input grids must resolve the
intrinsic spectrum and temporal motion; independence of sampling means
convergence under refinement, not recovery of unresolved structure.

For fixed photon flux per unit velocity, continuum noise is

sigma_ph = SNR_ref^-1 sqrt(t_ref/t_exp * R/R_ref).

ObservationConfig.snr_reference_exposure_time is in hours (default one minute),
and snr_reference_resolution defaults to 100000. Thus the configured continuum
S/N refers to these reference conditions, not an arbitrary simulation sample.
Noise per observed bin is sqrt(sigma_ph^2 + systematic_floor^2).
The matched-filter forecast is sqrt(sum(signal/noise)^2) times
sqrt(n_lines*n_transits), with signal=delta_atm times the observed residual.
No additional heuristic smearing penalty is applied: exposure averaging already
incorporates drift. The report's smearing_penalty is 1 by convention and its
cadence_limited_exposure is NaN (not evaluated by this estimator).
n_effective_pixels now counts physical observing bins, including zero-signal
bins. active_threshold is retained for legacy compatibility only.

## Scope and limitations for the referee

This corrects the artificial gain from counting internal samples as independent
data. It includes instrumental broadening and finite exposures. It remains a
known-template, independent Gaussian-noise forecast, not an injection-recovery
significance. The reference is assumed noiseless; finite-reference uncertainty,
correlated systematics, telluric/stellar contamination, dead time and throughput
are not modeled. The configured floor is independent bin noise, not a correlated
systematic floor. Lines and transits are assumed equivalent and independent;
square-root gains are conditional on that assumption. A fixed photon spectral
density across resolving powers is an approximation, not an instrument budget.

Do not reuse the old detectable counts or percentiles without recalculation.
Spatial geometry, physical atmosphere and phase selection also require separate
validation. Analytic tests cover Gaussian broadening plus bin integration,
linear time averaging, numerical-grid convergence and exposure invariance for
a stationary signal in the photon-noise-only limit.
