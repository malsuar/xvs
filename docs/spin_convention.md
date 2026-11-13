# Spin convention revision

Select SimulationConfig(spin_convention="observer") for the revised convention.
The default "legacy" preserves historical results. Geometry and spin options
are independent; new physical runs should explicitly select both
geometry="projected_annulus" and spin_convention="observer".

Let n point toward the observer and (u,v,n) be the orthonormal basis returned
by sky_basis. We define

Omega_hat = cos(i) n + sin(i) [cos(lambda) v + sin(lambda) u].

Thus i is the angle to the observer, and lambda is the projected spin-axis
angle measured from v toward u. For the default observer n=-x, v=+z and u=-y.
For other observers the reference is set by sky_basis, including its fallback
near the global z axis; it is not automatically a projected orbital axis.
The historical spin_obliquity field stores lambda but should be described as
a sky-plane spin angle, not an orbital obliquity. At i=0 or 180 degrees the
projected angle is undefined and the projected rotation vanishes.

The velocity formula remains v_eq (Omega_hat cross r) dot n. Its positive
direction is TOWARD the observer (approaching/blueshift); conventional positive
recession velocity would be its negative. Existing figures/captions associating
negative velocity with approach must be revised or explicitly use a sign-flipped
observational velocity coordinate. No automatic sign flip is made here.

The previous Rz(lambda) Rx(i) z axis is retained as legacy. It does not in
general make i an observer inclination. Under the revised convention the default
i=90, lambda=0 axis changes from -y to +z; old numeric angles cannot be reused
assuming they describe the same physical orientation. For example, the old
default physical axis corresponds to revised i=90, lambda=90 for observer -x.

## Proposed manuscript formulation

The spin inclination i_rot is measured from the direction toward the observer.
The projected spin angle lambda_rot is measured from a specified sky-plane
reference v toward u, with Omega_hat given above. It is a projected obliquity
only when that reference is explicitly identified with the projected orbital
angular-momentum direction. The cross-product formula uses velocities positive
toward the observer.

## Validation and referee status

Analytic tests verify Omega_hat dot n=cos(i) for two observer directions and
several inclinations/angles, sky-axis orientation, opposite limb signs, and
v_eq sin(i) on opposite unit-radius limbs. They also preserve the legacy axis.
This validates angle interpretation; it does not resolve atmospheric winds,
instrumental response, or the baseline detectability issues.
