# Controlled comparison of legacy and observing-bin S/N

Estimator correction committed as 8210a54. Reproduce from repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/tmp/xvs-mpl python scripts/compare_snr_revision.py
```

Outputs: docs/snr_comparison.csv and outputs/snr_revision_diagnostics.png/pdf.
Generated figures are ignored by Git; the script and numeric table are tracked.

## Explicit assumptions

All comparisons use the SAME residual map for both estimators, with projected
annulus geometry and observer-relative spin. This isolates estimator changes;
it is NOT a reproduction of the old geometry's published S/N.
Host: 0.7 Mj, 1.1 Rj, 500 K, stellar radius 1 solar; molecular weight 2.3,
five scale heights and unit contrast. Geometric annulus thickness remains 0.03
Rp, independently of the physical amplitude approximation. Reference continuum
S/N 300 per resolution element at R=100000 for a one-minute exposure; actual
exposures five minutes, 200 equivalent independent lines, one transit, independent
bin floor 10 ppm. Window [-8,8] hours. No full stellar transit envelope is applied.

The script paper-showcase has Rm/Rp=0.0892, whereas the historical notebook used
0.12. Both are compared explicitly. Base veq=12 km/s, sigma=1.5 km/s,
am/Rp=3, Pm=21.238292 h, phase=0.08. Alternative cases change only spin speed
or lunar radius, not a self-consistent satellite mass/evolution model.
The new observer spin default points along +z; old default axis was -y.
The moon orbit is left unchanged. Its alignment with the planetary spin must
be selected explicitly before treating these as final manuscript configurations.

## Results on the finest internal grid

| Case | Legacy S/N | Corrected S/N | Corrected peak ppm |
| --- | ---: | ---: | ---: |
| Script paper radius 0.0892 | 0.82874 | 0.08662 | 3.372 |
| Notebook paper radius 0.12 | 1.29516 | 0.13991 | 4.688 |
| Slow spin 3 km/s | 0.85135 | 0.08784 | 3.467 |
| Larger moon radius 0.2 | 2.63201 | 0.29821 | 7.329 |

Legacy/corrected ratios are grid-dependent, not universal correction factors.
For the first case, internal (time,velocity) grids of (481,201), (961,401),
(1921,801) give legacy S/N 0.20722, 0.41347, 0.82874 versus corrected
0.08587, 0.08622, 0.08662. Corrected variation across this range is below 0.9%.
This is temporal/spectral convergence at fixed spatial sampling (240 azimuths,
four radial cells), not a completed spatial convergence study.

## Interpretation for the referee

These controls demonstrate removal of spurious numerical sampling gain and
show finite-exposure/instrumental effects. They do not establish a detection:
none reaches S/N 5 under these assumptions. The slightly higher S/N of the slow
rotator reflects total missing absorption, not evidence of a resolved drifting
signature. The current matched filter can detect amplitude deficits as well as
kinematic structure. A future test must distinguish a moving shadow from nuisance
line-depth variations and atmospheric alternatives.

Next: validate spatial contact convergence and orbital/spin alignment; compare
alternative physical normalizations and observing conditions; then perform
injection-recovery before rerunning the selected population.
