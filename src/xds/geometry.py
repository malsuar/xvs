"""Vector geometry helpers."""

from __future__ import annotations

import numpy as np


def unit_vector(vec: np.ndarray | tuple[float, float, float]) -> np.ndarray:
    arr = np.asarray(vec, dtype=float)
    norm = np.linalg.norm(arr)
    if norm == 0.0:
        raise ValueError("zero vector cannot be normalised")
    return arr / norm


def rotation_x(angle_rad: float) -> np.ndarray:
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def rotation_z(angle_rad: float) -> np.ndarray:
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def sky_basis(axis: np.ndarray | tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return two perpendicular sky axes plus the supplied propagation axis."""

    w = unit_vector(axis)
    helper = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(w, helper)) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, w)
    u = unit_vector(u)
    v = np.cross(w, u)
    return u, v, w


def projected_offsets(points: np.ndarray, center: np.ndarray, axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Project point-center offsets onto the plane perpendicular to axis."""

    u, v, _ = sky_basis(axis)
    delta = np.asarray(points, dtype=float) - np.asarray(center, dtype=float)
    return delta @ u, delta @ v


def distance_to_ray(points: np.ndarray, ray_origin: np.ndarray, ray_axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return perpendicular distance to a ray and signed coordinate along it."""

    axis = unit_vector(ray_axis)
    delta = np.asarray(points, dtype=float) - np.asarray(ray_origin, dtype=float)
    along = delta @ axis
    perp = delta - np.outer(along, axis)
    return np.linalg.norm(perp, axis=1), along
