"""Rigid-rotation velocity fields."""

from __future__ import annotations

import numpy as np

from .geometry import rotation_x, rotation_z, unit_vector, sky_basis


def spin_axis(inclination_deg: float, obliquity_deg: float,
              observer_direction=(-1.0, 0.0, 0.0),
              convention: str = "legacy") -> np.ndarray:
    """Construct the legacy axis or an observer-relative inclination/sky angle.

    In observer convention lambda is measured from sky_basis v toward u.
    For observer -x these reference directions are +z and -y, respectively.
    """

    inc = np.deg2rad(inclination_deg)
    obl = np.deg2rad(obliquity_deg)
    if convention == "observer":
        u, v, n = sky_basis(observer_direction)
        return unit_vector(np.cos(inc)*n + np.sin(inc)*(np.cos(obl)*v + np.sin(obl)*u))
    if convention != "legacy":
        raise ValueError(f"Unknown spin convention: {convention}")
    axis = rotation_z(obl) @ rotation_x(inc) @ np.array([0.0, 0.0, 1.0])
    return unit_vector(axis)


def line_of_sight_velocity(
    positions: np.ndarray,
    v_eq: float,
    inclination_deg: float,
    obliquity_deg: float,
    observer_direction: np.ndarray | tuple[float, float, float],
    convention: str = "legacy",
) -> np.ndarray:
    """Compute v_los = v_eq (Omega x r) . n_obs."""

    omega = spin_axis(inclination_deg, obliquity_deg, observer_direction, convention)
    n_obs = unit_vector(observer_direction)
    return v_eq * (np.cross(omega, positions) @ n_obs)
