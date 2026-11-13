# First molecular template and known-trajectory injection/recovery

## Reproduction and data provenance

```sh
python -m pip install -e '.[molecular]'
python scripts/fetch_co_template_data.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/molecular_injection_recovery.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/audit_co_template.py
```

Data source: ExoMol 12C16O Li2015,
https://exomol.com/data/molecules/CO/12C-16O/Li2015/;
Li et al. 2015, ApJS216,15, https://doi.org/10.1088/0067-0049/216/1/15.
Line-intensity expression follows ExoCross:
https://doi.org/10.1051/0004-6361/201732531.
H2 broadening: https://exomol.com/data/data-types/broadening_coefficients/CO/.
Pinned source URLs, retrieval date and SHA256 checksums are in
data/exomol/CO_Li2015/manifest.json. Fetching validates them. Source/derived
molecular spectra use CC BY-SA4.0 attribution, separately from the MIT code;
see the data folder README. Raw files are excluded from Git.

Tracked results: co_injection_summary.json, co_injection_recovery.csv,
co_template_convergence.json. Regenerable outputs: co_injection_recovery.png/pdf
and co_molecular_injection.npz. The latter stores native observed maps and
local molecular profiles, not amplified physical predictions.

## Atmospheric model: explicit advance and remaining assumptions

This replaces identical synthetic lines with a line-by-line LTE, single-species
isothermal transmission calculation. It is NOT a radiative-convective or
chemical-equilibrium atmosphere, or a complete replacement for the Gaussian
run_xds_model API. Molecular synthesis is presently a separate pilot using the
same orbit, annulus, fractional-blocking, stellar-intensity and rotation helpers.

Host: 0.7 Mj,1.1 Rj,640 K; solar star on the pilot's 30-day planetary orbit.
Continuum/opaque radius Rp is assigned pressure P0=1e-3 bar (a model choice).
Constant CO mixing fraction=1e-4; mean molecular weight 2.3. CO abundance is
assumed, not established as equilibrium chemistry at 640 K. Constant gravity
and scale height give P(z)=P0 exp(-z/H). At tangent radius r,
N_CO=X_CO [P(z)/(kT)] sqrt(2 pi r H), with units converted to cm^-2.
Absorption is 1-exp[-N_CO sigma_nu(T,P)] in each radial layer.

The annulus extends to 10H and uses 32 equal-area radial cells. Upper/lower
sample columns are about 5.41e16 and 8.63e20 cm^-2. Its exact projected area,
Rp^2[(1+10H/Rp)^2-1]/Rstar^2, normalizes the geometric sum. The associated
area scale is 529.69 ppm; it is NOT the atmospheric line peak or lunar residual.
No arbitrary peak normalization is applied to CO opacity. The physical
normalization is now exact annular area rather than the earlier thin-annulus
approximation, and the atmospheric extent differs from the five-H pilot.

Line positions use upper-minus-lower energies; Einstein A and state degeneracy
give integrated strengths in cm/molecule:
S=gu A/(8 pi c nu^2 Q) exp(-c2 El/T) [1-exp(-c2 nu/T)].
Partition function at 640 K is Q=233.9567 from the source table. Voigt profiles
use SciPy's implementation; Gaussian width combines CO thermal motion (28 amu)
and an assumed 1 km/s microturbulent sigma. Lorentz HWHM uses the ExoMol a0
lower-J H2 coefficients at 296 K and one bar, scaled to H2 partial pressure
0.85 P; nearest endpoint is used above the tabulated J=80. He broadening,
pressure shifts, isotopologues, clouds, CIA, scattering and winds are absent.
The opaque lower boundary and CO abundance need physical justification in a
final atmospheric application. This experiment does not infer them from data.

Window 4250--4400 cm^-1 (about 2.27--2.35 micron) includes padding in synthesis.
145 transitions survive a strength threshold of 1e-8 times the strongest in
the padded 4245--4405 range; retained summed strength is 0.9999999974 of THAT
range at 640 K, not of all CO transitions. Wings are truncated at the larger
of 100 Gaussian sigmas or 100 Lorentz HWHM. The refinement doubles the cutoff.

## Geometry and observing sequence

Aligned circular moon .09 Rp at 3 Rp, Keplerian period21.23829 h; lunar phase
.75 at mid-transit. Four-hour central window, fully backed by the stellar disk;
quadratic continuum u1=.4,u2=.25. Representative-cell stellar weighting is used
in this smooth fully-backed window; ingress/egress is not part of this pilot.
The cell velocity uses solid-body rotation with surface veq12 km/s.

Fractional blocking on 480 angular x32 radial cells, order16; 1201 times. Cell
velocities are linearly deposited into .25 km/s velocity channels, preserving
weights and first moments. Both reference and missing-contribution histograms
are averaged over exposures, then combined with shifted radial CO templates.
Input log-wavenumber step is .2 km/s, with 50 km/s padding. Gaussian instrumental
FWHM c/R at R=100000 is applied, then piecewise-linear integration over fixed
c/R bins. Linear exposure/spectral operations commute here; pressure opacity
and exp(-tau) are evaluated BEFORE rotational/instrumental mixing.

48 five-minute exposures and 3468 spectral bins are used. Noise=1490.75 ppm
per bin, from continuum S/N300 at one minute, with independent 10 ppm floor.
NO extra sqrt(N_lines) is applied: the molecular spectrum already contains the
line information. Reference atmosphere and stellar spectrum are assumed known
without noise. Photon noise is independent Gaussian, with no telluric or stellar
line contamination and no correlated systematic error model.

## Depth-only versus kinematic models

Let S be the observed native lunar residual and B the moonless atmospheric
template. At EACH exposure, allow an independent unconstrained depth change
beta_t B_t. Its best approximation to S is D_t=beta_t B_t. The remaining
kinematic component is Q=S-D, orthogonal to B in the adopted equal-noise metric.
This is more flexible than one global depth change, but it is not a model of
all atmospheric variability. Temperature/abundance changes, winds and stellar
spectral variations can change shapes and are not removed by this projection.

For a FIXED known orbit/trajectory, signed Z=<Y,Q>/(sigma ||Q||); expected
kinematic S/N=||Q||/sigma. The fitted lunar amplitude is Z/SNR_kin. Adding a
signed lunar amplitude to the depth-only model gives delta chi2=Z^2; this is
explicitly checked on one injected noisy spectral sequence. Positive and
negative amplitudes are allowed in this diagnostic. Unknown-orbit searches,
positivity constraints and template mismatch need different calibration.

## Results and injection controls

Native full-template S/N=0.04428; after depth nuisance fitting,0.03819.
The kinematic norm is 86.25% of the full norm (about74.4% of its squared norm).
Thus the selected spectrum is not equivalent to a simple depth change, but
the native perturbation is UNDETECTABLE under these observing assumptions.
Do not confuse mathematical distinguishability with measured significance.

Seed63117, 2000 actual independent Gaussian noise maps projected onto the
templates give null mean Z=.0297 and standard deviation .9999. Two null draws
exceed Z=3. This is a descriptive low-tail control; 2000 draws do not calibrate
a 5-sigma false-alarm probability or an orbital-template search.

A controlled residual multiplier130.919 gives EXPECTED kinematic Z=5. This
deliberately amplified injection is a numerical laboratory test, not a physical
moon prediction, larger radius, or validated instrument improvement. It may
violate the atmospheric amplitude budget and cannot establish feasibility.
Its recovered Z ensemble has the expected unit scatter; 98.15% exceed Z=3.
One explicit noisy realization has Z=4.177, fitted multiplier109.38 and
delta chi2=17.451 versus the depth-only fit, consistent with the known linear
Gaussian estimator. The native signal remains buried in noise.

A separate amplified DEPTH-ONLY template is scaled to give expected Z=5 in
the unprofiled lunar filter. Its expected depth-profiled score is zero (within
floating-point precision). This illustrates why the full residual S/N alone
does not establish kinematic information. The projected null distributions are
identical by linear algebra; they are not separate independent noise ensembles.

## New-model convergence, including atmospheric extent

The joint refined molecular run doubles angular/radial/time/spectral sampling,
uses .125 km/s velocity channels, and doubles the wing cutoff. Its S/N changes
by 0.0756% (full) and0.0797% (depth-profiled); residual-map norm difference is
0.1071%. The fine values are0.04431 and0.03822. These checks concern this
molecular pilot, rather than inheriting the earlier Gaussian-map convergence.

Five-H versus ten-H models differ in kinematic S/N by3.97%; truncation at5H
was not sufficient for the adopted CO abundance/column. Extending10H to12H,
at comparable radial sampling, changes kinematic S/N by0.0255%. The10H domain
is therefore a checked starting choice for this atmosphere, not a universal
height cutoff. The reference-pressure/radius convention still carries physical
model uncertainty beyond these numerical comparisons.

## Referee status and next experiment

Completed: traceable molecular lines with radial slant opacity, physical
amplitude rather than peak rescaling, instrument/exposure integration, a
depth-only nuisance model, native nondetection, and a known-trajectory recovery
laboratory control. Tests include optically thin integrated strength, exact
linear bin averaging and null projection of pure depth variations.

Not yet demonstrated: blind recovery, atmospheric model uniqueness, realistic
instrument/telluric systematics or astrophysical detectability. Next add
temperature/composition/wind alternatives and an orbit/phase search with
false-alarm calibration. A realistic chemistry/multi-species atmosphere and
instrument throughput/noise budget are still needed before population claims.
Do not replace the old sweet-spot figures with an inferred detection region
from this single, phase-selected experiment.
