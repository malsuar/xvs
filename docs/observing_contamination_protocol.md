# Predeclared contamination/count-budget pilot

This extends the phase/period-blind pilot, without changing its physical moon,
validated geometry, 73x1024 orbit bank or bounded efficient-score statistic.
We fix three temporal atmospheric controls, rather than selecting between
models after seeing a candidate. Every realization refits the nonlinear null.

## Observing assumptions

- R=100000, same CO interval and four-hour window, 48 five-minute cadences.
- Eight30-second integrations give240 seconds live time in the8m case. The
  39m/K6 case uses84x1.427s=119.868s live time to avoid nonlinearity. These
  duty cycles are assumptions, not observatory overhead forecasts. The39m/K3
  budget fails the adopted linearity limit and is only a photon upper bound.
- 49 or 978 m2 effective geometrical collecting area; efficiency 8%, excluding
  atmospheric transmission. The 39 m scenario is hypothetical, NOT a certified
  instrument performance claim. K=6, plus a K=3 count-budget comparison.
- 666.7 Jy Ks zero point, approximate monochromatic normalization at 2.159 um,
  5700 K continuum slope and PHOENIX line spectrum (5700K/logg4.5/solar).
- ESO SkyCalc Paranal, airmass1.3, PWV2.5mm, physical emission and transmission;
  noise uses this fixed mean condition, not a full changing-airmass ETC.
- Extraction0.08arcsec2, 8 pixels/resolution bin, 12e read noise, 0.003e/s dark,
  telescope thermal surrogate280K/emissivity0.1, sky subtraction factor2.
- A conservative peak/mean pixel factor2 and12900e linearity limit are adopted;
  real spatial profiles, nonlinearity and ramp/overhead performance need an ETC.
- A one-hour, 12-exposure shared reference has finite noise. Its time covariance
  I+11^T/12 is used EXACTLY in fits, scores and simulated draws; it is not
  treated as independent noise in every exposure.
- Fixed transmission cut >=0.3 and stellar normalized flux>=0.2, set before any
  outcomes. No data-driven masks or selection of detected injections.

## Residuals and controls

All generators include the preceding bounded smooth atmosphere plus wind-pulse
distribution. Four families add no contamination, stellar residuals, telluric
residuals, or both. This is a sensitivity model, not a calibrated distribution
of residuals from a real reduction pipeline.

Stellar residuals retain a random3--15% of a disk-minus-occulted-profile RM
surrogate, using PHOENIX surface flux as a local spectrum, vsini2km/s, the
actual central transit chord and a phenomenological20%(1-mu) depth variation.
This is NOT a specific-intensity CLV calculation, spot model or full stellar
surface integration. The moon's own photometric stellar occultation is absent.

Telluric residuals have optical-depth errors uniform0.001--0.005 and wavelength
errors uniform0.003--0.015km/s, each with independent zero-mean unit-variance
Gaussian-process histories, squared-exponential timescale0.7h. Stellar retention
also has a10% correlated modulation. The injected patterns are shifted into the
planetary frame using the planet's 30-day orbital acceleration; an adopted
0.1km/s/h terrestrial frame drift is included. Systemic offset is fixed at zero
for a potentially overlapping case, not marginalized over target velocities.

Two predeclared reductions are compared: unprotected and protected. The latter
projects six fixed nuisance spectra (tau, tau derivative, tau squared, stellar
contrast and its first/second derivatives) with free coefficients PER exposure.
Their shapes are known exactly here; neither amplitudes nor history are taken
from the injected truth. Actual shifted nonlinear contaminants are not forced
to lie in their span. Removed moon information is explicitly measured.

Two noise levels use the SAME photon-noise spectral shape: the 39m/K6 count
budget, and an unphysical median11.39465ppm laboratory control. No lunar signal
amplification, equal-noise assumption or extra independent spectral lines.
Only the count budget is realistic in origin; efficiencies/backgrounds and
correction accuracy remain scenario assumptions, not target-specific forecasts.

## Statistics

Seed63150; 64 calibration nulls and 32 independent held-out nulls per family,
24 injections per group (uncontaminated fixed phase, combined fixed phase,
combined uniform phase). All groups include atmospheric wind variability.
The noise levels and two reductions use paired data; their rates are NOT
independent observations. Seed streams separate calibration/validation/injection.

Each reduction/noise level has its own global search-max calibration. Compare
the atmosphere-only null tail to the maximum tail over all four generators.
p=(1+number of null maxima>=observed)/(N+1), diagnostic alpha0.05; minimum
p=1/65. Held-out rates receive exact95% binomial intervals. These finite
ensembles cannot demonstrate discovery significance or protect against all
possible contaminants. The nonlinear candidate LRT remains uncalibrated.

Normalized planetary spectra are already instrument-binned. Combining their
signal with stellar/telluric residuals at that level is a small-signal observation
surrogate, NOT a simulation of multiplication before the LSF, detector ramps,
order gaps, wavelength calibration or a real detrending pipeline. The current
planetary spectral subspace discards external spectral information that could
otherwise constrain the contaminants. We must state these limits in results.

Sources: ESO CRIRES performance page (detector quantities), ESO SkyCalc
documentation, Göttingen PHOENIX archive, and IPAC 2MASS calibration. Request
parameters, source URLs and downloaded-file hashes accompany the template cuts.
