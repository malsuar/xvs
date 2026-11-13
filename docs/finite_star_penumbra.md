# Finite star: direct transmission versus irradiance penumbra

Reproduce: `PYTHONPATH=src python scripts/audit_penumbra_pilot.py`.
Data: penumbra_pilot.csv; source geometry: penumbra_assumptions.json.
Figure: outputs/penumbra_pilot.png/pdf.

## Physical distinction that must be resolved in the manuscript

An extended star generates a penumbra in TOTAL incident irradiance at an
atmospheric point. Pure direct transmission toward a distant observer instead
follows the specific stellar ray through that point. With no scattering or
refraction coupling directions, its intensity is I_star(ray) exp(-tau(ray));
the moon removes that ray if intercepted upstream or downstream. Integrating
irradiance over all stellar directions and multiplying the transmission opacity
by that fraction is a DIFFERENT physical model, requiring an atmospheric response
or redistribution mechanism. It is not automatically a finite-star correction
to the direct transmitted spectrum.

The earlier wording suggesting finite-star penumbra would necessarily suppress
the direct transmission signal was too broad. We preserve the direct-ray
production model and add a SEPARATE uniform-disk irradiance diagnostic. No
finite-star mode is silently substituted into run_xds_model or the pilot.

## Direct-ray control

Pilot star radius=8.84651 Rp, distance=359.44240 Rp, angular radius about
0.024612 rad. Within the central four-hour window, the sky displacement of
the planet plus annulus radius stays inside the uniform stellar disk; the
script asserts this sufficient bound using circular planetary motion.
Therefore uniform I_star weights remain constant and the ideal direct-ray
mask/S/N are unchanged. This is a geometric result under the stated uniform,
fully backed, pure-absorption assumptions, not a full limb-darkened stellar
spectral simulation. Stellar ingress/egress, rotation, limb darkening, scattering
and refraction can change the conclusion and are not included here.

## Irradiance diagnostic implementation

Source plane at distance D, uniform circular stellar radius Rs. An upstream
moon at axial distance L has g=1-L/D. In the effective atmosphere sky plane,
shadow-center position is moon_transverse/g, lunar shadow radius Rm/g and
source blur radius B=L Rs/(D g). The fractional irradiance blocked at a point
is the overlap area of disks B and Rm/g divided by pi B^2. This includes finite
distance magnification in the paraxial thin-occulter approximation; it does not
trace an exact spherical stellar photosphere or tilted lunar silhouette.
The B=0 branch reduces to fractional point-source disk overlap, and tends to
the parallel-ray mask as D goes to infinity. Downstream observer blocking is
unchanged.

Cell averages integrate this probability by tensor Gauss-Legendre quadrature
in r squared and azimuth, orders 8 and 16. Local spectra still use representative
cell velocities. The diagnostic assumes immediate LOCAL LINE-AMPLITUDE
suppression proportional to missing irradiance, with no chemistry, thermal
response or scattering calculation. This assumption is deliberately labelled
irradiance_toy in outputs: its spectral residual/S/N is hypothetical.

## Results relative to the physical pilot

Same host, atmosphere-linked annulus, radii .03,.06,.09 Rp, central and shallow
grazing trajectories. Lunar phase .75 gives upstream illumination; .25 gives
observer occultation. Other two pilot phases have no contacts and are not
rerun here. Four exposure origins and the same noise assumptions are used.

At L around 3 Rp, the source blur scale is about .074 Rp, comparable to or
larger than the smaller moons. Central irradiance peaks are therefore diluted
and spread over time. Shallow grazing geometries can gain irradiance-shadow
overlap that was absent in a collimated beam. The hypothetical S/N ratio is:

| Rm/Rp | Central toy/direct range | Grazing toy/direct range |
| --- | ---: | ---: |
| .03 | .702--.876 | 8.630--8.787 |
| .06 | .857--.864 | 3.356--3.380 |
| .09 | .911--.928 | 1.884--1.894 |

Ranges cover four observing phases, not all possible phases. Observer-side
ratios are exactly one, by construction, because only the upstream irradiance
model changes. Large grazing ratios compare against very small direct-ray S/N;
they are not evidence for detectability or amplified direct transmission.
Maximum relative S/N change from quadrature order 8 to 16 is 1.19e-6
(0.000119%). This validates only this quadrature-order control at fixed
spatial/time/spectral sampling, not a full new convergence envelope. Existing
direct-ray grid tests do not automatically validate every irradiance model.

## Tests and referee implications

Tests cover circle-overlap limits, central penumbra probability, point-source
limit, absence of downstream irradiance shadow and cell-quadrature convergence.
The production direct-ray treatment remains the default. The irradiance module
is exploratory and must not replace the pilot's S/N in a new population search.

For the referee, explicitly state which observable is modeled: missing
observer-directed transmission or atmosphere changes driven by irradiation.
The latter requires an atmospheric response model and may have different sign,
amplitude and timescale. For a direct-transmission paper, next implement
stellar surface intensity/background weighting and realistic local spectra,
rather than claiming that an irradiance penumbra is a validated transmission
correction. The two historical labels illumination/observer should be explained
as ray ordering unless additional atmospheric-response physics is introduced.
