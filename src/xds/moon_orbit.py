"""Keplerian moon positions in planet-radius units."""

from __future__ import annotations

import numpy as np

from .config import MoonConfig
from .geometry import rotation_x, rotation_z


def solve_kepler(mean_anomaly: np.ndarray, eccentricity: float, iterations: int = 8) -> np.ndarray:
    mean_anomaly = np.asarray(mean_anomaly, dtype=float)
    if eccentricity == 0.0:
        return mean_anomaly
    eccentric_anomaly = mean_anomaly.copy()
    for _ in range(iterations):
        eccentric_anomaly -= (
            eccentric_anomaly
            - eccentricity * np.sin(eccentric_anomaly)
            - mean_anomaly
        ) / (1.0 - eccentricity * np.cos(eccentric_anomaly))
    return eccentric_anomaly


def moon_positions(times: np.ndarray, config: MoonConfig) -> np.ndarray:
    """Return moon positions in the planet-centred frame."""

    mean_anomaly = 2.0 * np.pi * (np.asarray(times) / config.period + config.phase)
    ecc = config.eccentricity
    e_anom = solve_kepler(mean_anomaly, ecc)
    true_anom = 2.0 * np.arctan2(
        np.sqrt(1.0 + ecc) * np.sin(0.5 * e_anom),
        np.sqrt(1.0 - ecc) * np.cos(0.5 * e_anom),
    )
    radius = config.semi_major_axis * (1.0 - ecc * np.cos(e_anom))
    orbital = np.column_stack(
        [
            radius * np.cos(true_anom),
            radius * np.sin(true_anom),
            np.zeros_like(radius),
        ]
    )
    rot = (
        rotation_z(np.deg2rad(config.longitude_node))
        @ rotation_x(np.deg2rad(config.inclination))
        @ rotation_z(np.deg2rad(config.argument_pericenter))
    )
    return orbital @ rot.T
