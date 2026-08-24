# Scientific Landing Checklist for the XDS Manuscript

The current manuscript establishes the kinematic mechanism of the Exomoon
Doppler Shadow. To make the result quantitatively testable, it needs a second
layer: an observability and model-comparison budget with explicit assumptions.

## Required Additions

1. Define the two satellite masks separately:
   - illumination shadowing, when the satellite lies between star and
     atmospheric terminator;
   - observer-line occultation, when the satellite lies between terminator and
     observer.

2. Replace the conceptual efficiency curve with a physical detectability metric:

   ```text
   SNR_XDS = sqrt(sum[(Delta F_XDS / sigma)^2]) * sqrt(N_lines * N_transits)
   ```

   with penalties for spectral smearing, finite exposure time, and systematic
   floors.

3. State the atmospheric amplitude model:

   ```text
   delta_atm ~= 2 Rp N_H H / Rstar^2
   H = k_B T / (mu m_H g)
   ```

   and identify which molecular bands provide the line list.

4. Quantify cadence and integration limits:
   - contact duration through the terminator;
   - velocity drift of the shadow;
   - maximum exposure before the signal is smeared by more than a chosen
     fraction of one resolution element.

5. Include injection-recovery tests:
   - generate synthetic time-velocity maps;
   - inject photon noise and systematic floors;
   - recover moon parameters with and without the dual-mask model.

6. Compare against false positives:
   - winds;
   - evolving clouds/hazes;
   - oblateness;
   - rings;
   - stellar residuals;
   - static rotation mis-modelling.

7. Provide target-specific forecasts:
   - Kepler-167e-like cold giant;
   - warm long-duration giant;
   - hot giant control case where moon survival and spin are unfavourable.

## Minimum Quantities to Report

- peak XDS depth in ppm;
- matched-filter SNR;
- required number of spectral lines for SNR 5 in one transit;
- required number of transits for a fixed line list;
- exposure-time limit from velocity smearing;
- contact windows for illumination and observer occultation separately.

## Interpretation

A scientifically grounded claim should not be "XDS exists" only. It should be:

```text
Given planet, atmosphere, instrument, and observing sequence assumptions,
the XDS signal is detectable/non-detectable at SNR X, with required cadence Y
and integration time Z.
```
