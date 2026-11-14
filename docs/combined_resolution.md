# Combined-grid and exposure-phase audit

Reproduce with `PYTHONPATH=src python scripts/audit_combined_resolution.py`.
Table: combined_resolution.csv. Figure: outputs/combined_resolution_audit.png/pdf.

## Experiment

Uses the three circular encounters, physical normalization and observing
assumptions in grazing_convergence.md. Combined grids: 480x16, 960x32,
1920x64. Azimuthal samples are rotated by 0, 0.25 and 0.5 cell spacings,
with orbit and spin held fixed. Radial samples retain equal-area midpoints;
radial counts change but radial quadrature origins are not shifted here.
Exposure starts shift by 0, 0.25, 0.5 and 0.75 times five minutes.

The new ObservationConfig.exposure_start_offset is in hours, defaults to zero,
and must lie in [0,exposure_time). Only complete exposures within the supplied
window are retained. Shifted runs lose a partial boundary exposure; in these
experiments the signal is zero at both boundaries, so the comparison is not
driven by clipping a contact. This must be rechecked in other windows.

For each exposure phase the reference S/N is the average of the three shifted
1920x64 grids. Reference rows therefore have nonzero errors: they measure
finite spatial-phase scatter, not disagreement with an exact S/N solution.
The reference itself retains discretization error. All errors below are maxima
over tested cell shifts and observing phases, not formal uncertainty bounds.

## Results

| Impact Rp | Worst S/N error 480x16 | Worst S/N error 960x32 | Fine-grid phase scatter about reference |
| --- | ---: | ---: | ---: |
| 0 | 0.370% | 0.0465% | 0.0089% |
| 1.08 | 0.0422% | 0.00947% | 0.00277% |
| 1.14 | 0.432% | 0.302% | 0.101% |

Worst relative integrated-area errors versus analytic overlap at 960x32 are
0.0437%, 0.00467%, 0.1178% respectively. At 480x16 they are 0.306%,
0.0619%, 0.663%. Near tangency, even the fine-grid integrated-area error varies
with sample origin (up to 0.148%); a single favorable phase can hide error.

Changing observing phase on the fine grid yields peak-to-peak S/N ranges
divided by mean S/N of approximately 2.11--2.12% (central), 0.0198--0.0219%
(grazing), and 0.0699--0.1159% (nearly tangent). This is a cadence-dependent
physical effect plus residual numerical error, not an error that spatial
refinement must drive to zero. Near tangency the cell-phase scatter is comparable
to this small observing-phase effect, so its precise size is not yet robust.

## Operational resolution criterion before a population search

Use 960x32 as a provisional initial grid for moons around 0.12 Rp and annuli
of thickness 0.03 Rp, NOT as a universal resolution guarantee. For a candidate:

1. Compare to a grid doubled in BOTH directions; require corrected S/N changes
   below 1% at EACH tested observing phase, not only in a phase average.
2. Shift angular samples by quarter and half cells and require a numerical S/N
   envelope below 1%. Do not demand observing-phase variation below 1%; report
   its range separately. If classification at S/N=5 changes across this range,
   label it cadence-sensitive rather than assigning a unique detection status.
3. Where analytic uniform-annulus overlap applies, require integrated-area error
   below 1% across phases, and examine pointwise curves to avoid cancellation.
4. Independently require temporal/spectral refinement below 1%. Previous temporal
   controls support this only for the tested geometry and unshifted observing
   sequence. Repeat checks for shifted sequences near a decision threshold.
5. Contact timing needs a separate tolerance tied to exposure duration: use
   analytic boundaries when available; otherwise verify a timing tolerance
   before claiming contact durations. S/N convergence alone is insufficient.
6. Refine nearly tangent or smaller-moon systems adaptively. Validate grid criteria
   across new radii, annular thicknesses, orbital orientations and spatial opacity
   structure before using them as a blanket population guarantee.

The 1% threshold is an explicit working numerical tolerance, not a physical
precision claim. These cases pass the measured S/N and integrated-area portions
at 960x32. Full population readiness additionally requires the other checks.
None of these controls changes the atmospheric physics or demonstrates recovery
against stellar/atmospheric alternatives.

## Referee-response material

We have now separated spatial quadrature, cell-origin sensitivity and observing
cadence. A reproducible grid ensemble demonstrates sub-percent corrected S/N
differences for the tested central and grazing encounters. We retain explicit
qualifications for contact timing, finite reference resolution and broader
parameter coverage. Proposed population criteria are not yet a validation of
the old selected population and its detectable fractions.
