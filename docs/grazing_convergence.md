# Grazing contacts: separate spatial convergence and observed S/N

Reproduce from repository root:

```sh
PYTHONPATH=src MPLCONFIGDIR=/tmp/xvs-mpl OPENBLAS_NUM_THREADS=1 python scripts/audit_grazing_convergence.py
```

Numeric output: docs/grazing_convergence.csv. Figure:
outputs/grazing_spatial_snr_convergence.png/pdf (regenerable, ignored by Git).

## Geometry and controls

Circular moon: a=3 Rp, radius=0.12 Rp, period=21.238292071543107 h.
Uniform annulus: [1,1.03] Rp. Observer and stellar photon directions: -x.
Spin: observer convention, i=90, lambda=0, axis +z, veq=12 km/s.
Intrinsic Gaussian sigma=1.5 km/s. Orbit node=90 deg,
inclination=asin(b/a), phase=0.25. This gives a true projected minimum
separation b at t=0; lunar inclination is global, not inclination to observer.
These are different orbital alignments as well as different impact parameters;
they are not an impact-only comparison at fixed mutual angle.

Cases: b=0 (central, two separate annular crossings), b=1.08 (grazing),
b=1.14 (nearly tangent; outer tangency is b=1.15). Window [-2,2] h, 1201
time samples and 401 velocity samples over [-30,30] km/s. Combined blocking
equals the projected lunar disk for these parallel directions and nonzero
moon x. A test checks this equivalence against compute_shadow_masks.
The residual is evaluated directly as the negative blocked weighted Gaussian
sum, algebraically equivalent to subtracting the reference without cancellation.

Blocked area has an independent analytic reference: difference of two
circle-overlap areas divided by annular area. The outer contact span is
P/pi * asin(sqrt(((Rout+Rm)^2-b^2)/(a^2-b^2))). This span includes the
unblocked gap in the central case and must not be called total active duration.

## Refinement design

Angular sweep: Nphi=120,240,480,960 at Nr=64.
Radial sweep: Nr=2,4,8,16,32 at Nphi=1920.
S/N reference: Nphi=1920, Nr=64. This is a numerical reference, not an exact
S/N solution. Area errors use the analytic overlap. Separate sweeps do not
establish the combined error of, for example, a 480 x 32 grid.

Observed S/N: R=100000, five-minute integrations, continuum S/N=300 at
reference one minute/R=100000, 200 independent equivalent lines, one transit,
independent bin floor=10 ppm. Amplitude atmosphere: 0.7 Mj,1.1 Rj,500 K,
stellar radius one solar, molecular weight 2.3, five scale heights.
Geometric thickness and amplitude remain independent parameters. No stellar
transit envelope, atmospheric inhomogeneity or recovery experiment is included.

## Results

| b/Rp | Reference corrected S/N | Area integral error of reference | S/N change doubling time samples |
| --- | ---: | ---: | ---: |
| 0 | 0.117445 | 0.0118% | 0.00647% |
| 1.08 | 0.119831 | 0.00701% | 0.00466% |
| 1.14 | 0.0067473 | 0.000682% | 0.0636% |

The very small integrated-area error in the last case partly reflects error
cancellation: maximum instantaneous fraction error is 2.72e-5, not zero.

| Sweep | b=0 S/N error | b=1.08 S/N error | b=1.14 S/N error |
| --- | ---: | ---: | ---: |
| Angular 120, radial 64 | 4.180% | 0.0118% | 0.298% |
| Angular 480, radial 64 | 0.305% | 0.00673% | 0.183% |
| Angular 960, radial 64 | 0.0162% | 0.000882% | 0.209% |
| Radial 4, angular 1920 | 0.0455% | 0.0894% | 7.092% |
| Radial 16, angular 1920 | 0.00429% | 0.0171% | 0.388% |
| Radial 32, angular 1920 | 0.0275% | 0.00816% | 0.0505% |

Near tangency, radial 4 also biases integrated area by 6.97%; radial 2 by
22.58%. At radial 32 this area error is 0.102%. Small S/N error does not prove
pointwise convergence or correct contact times. For example, the near-tangent
outer span is underestimated by 4.93 minutes with radial 4 (at angular 1920),
and by 0.531 minutes in the 1920 x 64 reference. Doubling time samples reduces
the latter to 0.331 minutes. Midpoint point masks do not resolve cell boundaries.
The temporal control does not replace an independent velocity-grid convergence
test for each geometry.

## Conclusions for the referee

The new observing estimator remains sensitive to spatial integration errors,
as it should: correcting sample-count normalization does not cure unresolved
geometry. Central encounters need adequate angular sampling; nearly tangent
encounters require more radial sampling. Errors can oscillate under refinement
because the moon crosses discrete cell centers. Report both area and S/N
convergence and compare contacts to analytic boundaries, rather than choosing
resolution from S/N alone.

These results support case-by-case tolerance checks, not a universal default.
Before a population search, test combined grids and smaller moons, vary azimuthal
grid phase and exposure alignment, and refine any nearly tangent configurations.
For precise contact times prefer analytic geometry or fractional-cell integration.
No change to default resolution or archived population results is made here.
