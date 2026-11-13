# Projected annulus: first geometry revision

The legacy builder and default configuration are preserved for historical
reproducibility. Select the revised geometry explicitly:

```python
planet = PlanetConfig(atmosphere_half_width=0.03)
sim = SimulationConfig(geometry="projected_annulus", n_azimuths=360, n_radii=8)
```

For this geometry, the historical atmosphere_half_width field specifies the
FULL outward thickness h in Rp units: Rin=1, Rout=1+h. It is not a half-width
about Rp. All positions and orbital distances remain in Rp units. The unused
PlanetConfig.radius is not a new physical scaling factor.

The basis (u,v) spans the plane perpendicular to the observer. Elements are
r(cos(phi)u+sin(phi)v). Radial cells are uniform in r squared and azimuths
use periodic midpoints. Cell areas are exactly
pi*(r_outer^2-r_inner^2)/Nphi. Weights are cell areas divided by total area;
there is no spherical cosine illumination factor and no post-mask normalization.
The total projected area is pi*((1+h)^2-1). No element lies inside the opaque
planetary disk. Both projected limbs are sampled without reflection.

This is an effective sky-plane transmission model with uniform absorption per
projected area. It is not vertical radiative transfer. Representative points
lie at zero line-of-sight coordinate; their use in ray ordering is an effective
thin-plane approximation. Finite-star illumination, stellar limb darkening,
planetary ingress/egress and non-collinear star-observer geometries have not
been validated by the present controls. Rotation remains rigid with altitude,
so local speed can exceed the equatorial surface value at radius above Rp.
The legacy spin parameter convention is unchanged and remains a pending issue.

The independent detectability scale 2 Rp N_H H/Rstar^2 is still a thin-annulus
approximation. Thickness h is currently an explicit geometric parameter, not
automatically N_H H/Rp. The tests do not validate absolute physical amplitudes
or S/N; those require the next observing-model revision.

## Analytic controls

Tests check exact total area, normalized weights, radial bounds, sky-plane
orthogonality for a general observer direction, zero first spatial moment,
and the uniform-area second moment <y^2>=(Rout^2+Rin^2)/4. For rotation with
vlos proportional to y this also fixes its unbroadened second moment.
A concentric lunar disk cutting the annulus has blocked fraction
(Rm^2-Rin^2)/(Rout^2-Rin^2), tested in both ray orderings to radial-cell accuracy.
Forward synthesis checks a symmetric moonless profile, exact reference recovery
with masks disabled, and nonpositive residuals when blocking is active.

## Referee-response material

Implemented: a distinct projected-annulus geometry, replacing the spherical
surface interpretation for new explicitly selected runs; analytic controls of
area and spatial moments and a known occultation fraction. Preserved: legacy
geometry and published-baseline configuration. Pending: convergence for
arbitrary moon contacts, spin convention, observing/noise model and recalculated
population results. Do not claim revised detectability from these geometry tests.
