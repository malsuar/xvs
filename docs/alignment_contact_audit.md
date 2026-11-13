# Orbital alignment and spatial contact audit

Reproduce: `PYTHONPATH=src python scripts/audit_alignment_contacts.py`.
Numeric results: contact_convergence.csv. Generated figures in outputs:
alignment_contact_audit.png/pdf and alignment_residual_maps.png/pdf.

## Orbital convention

The lunar inclination is measured from global +z, not from the observer or
planetary spin. Its angular-momentum normal is Rz(node) Rx(inclination) zhat.
The argument of pericenter does not change this normal. The new orbital_normal
and spin_orbit_angle helpers expose this convention and compute the true mutual
angle acos(normal dot spin). Tests check the plane and prograde direction.

For observer -x and revised spin i=90, lambda=0, the spin is +z. An equatorial
prograde orbit therefore uses lunar inclination=0, node=0. The previous lunar
inclination=90, node=0 gives normal -y and a polar mutual angle of 90 degrees.
Under the LEGACY spin default (-y), that same orbit was aligned. Thus changing
the spin convention alone changes the physical spin-orbit configuration.
No orbital defaults or historical presets are silently changed in this audit.

The maps compare aligned and polar orbits at fixed radius 0.12 Rp, separation
3 Rp, period 21.238292 h and phase 0.08. Both are edge-on to observer -x and
have the same projected separation history. Aligned contacts remove high-speed
limbs of the +z rotating annulus; polar contacts remove material near zero
projected velocity. The earlier polar S/N comparison is a valid control of the
estimator but not an aligned showcase. It should not be used as the final
paper configuration without explicitly choosing its mutual inclination.

## Independent analytic contact reference

Uniform annulus Rin=1, Rout=1.03, lunar radius 0.12, circular separation 3,
period 21.238292 h. Around conjunction d(t)=a|sin(2 pi t/P)|.
On one increasing-separation passage the nonzero-overlap boundaries are
t_in=P/(2 pi) asin((Rin-Rm)/a)=1.00631887 h and
t_out=P/(2 pi) asin((Rout+Rm)/a)=1.32977138 h.
One passage lasts 19.4072 minutes; this is not a universal event duration.

Blocked fraction is exactly [C(Rout,Rm,d)-C(Rin,Rm,d)] /
[pi(Rout^2-Rin^2)], where C is the two-circle overlap area. The spatial test
compares point-mask quadrature directly with that formula, independently of
the spectral model. It uses a fixed temporal step of about 0.133 minutes.

| Azimuth x radial cells | Maximum absolute fraction error | Integrated fraction relative error | Entry / exit error [min] |
| --- | ---: | ---: | ---: |
| 120 x 4 | 0.005105 | 3.044% | +0.548 / -0.541 |
| 240 x 8 | 0.004039 | 1.798% | +0.283 / -0.275 |
| 480 x 16 | 0.001085 | 0.301% | +0.150 / -0.143 |
| 960 x 32 | 0.000465 | 0.040% | +0.150 / -0.143 |

Maximum errors are fractions of the entire annulus, not relative percentages
of the blocked signal. Integrated error measures time-integrated blocked area,
not matched-filter S/N. Contact estimates plateau partly because of fixed time
sampling and midpoint cells; exact boundaries should not be inferred from a
single active grid element. This experiment refines both spatial directions
together; it does not establish separate radial and angular error budgets.

## Implications and referee material

The circular aligned case approaches an independent analytic overlap reference.
480 x 16 cells is a provisional working resolution for THIS geometry at roughly
0.3% integrated-area error, not a universal justified default. Smaller moons,
grazing encounters, oblique orbits and spectral moments/S/N require their own
convergence controls. Default 360 x 8 remains unchanged pending broader tests.
Figures are intrinsic geometric diagnostics; they do not include the observing
operator or establish detection significance. The full stellar transit envelope
is still absent. The nearly empty time-velocity maps reflect localized annular
contacts, not the extended spherical-surface trails of the baseline model.

Before population recalculation: choose/document the showcase mutual angle,
verify separate radial/angular refinement including grazing cases, and propagate
the converged geometry through the observing operator and recovery experiment.
