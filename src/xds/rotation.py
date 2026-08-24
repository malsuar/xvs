"""Rigid-rotation velocity fields."""

from __future__ import annotations

import numpy as np

from .geometry import rotation_x, rotation_z, unit_vector


def spin_axis(inclination_deg: float, obliquity_deg: float) -> np.ndarray:
    """Return a unit spin axis using RM-like inclination and sky angle."""

    inc = np.deg2rad(inclination_deg)
    obl = np.deg2rad(obliquity_deg)
    axis = rotation_z(obl) @ rotation_x(inc) @ np.array([0.0, 0.0, 1.0])
    return unit_vector(axis)


def line_of_sight_velocity(
    positions: np.ndarray,
    v_eq: float,
    inclination_deg: float,
    obliquity_deg: float,
    observer_direction: np.ndarray | tuple[float, float, float],
) -> np.ndarray:
    """Compute v_los = v_eq (Omega x r) . n_obs."""

    omega = spin_axis(inclination_deg, obliquity_deg)
    n_obs = unit_vector(observer_direction)
    return v_eq * (np.cross(omega, positions) @ n_obs)
