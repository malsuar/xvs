"""Discretised atmospheric terminator."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TerminatorGrid:
    longitude: np.ndarray
    latitude: np.ndarray
    positions: np.ndarray
    weights: np.ndarray
    cell_area: np.ndarray


def build_terminator(
    n_longitudes: int,
    n_latitudes: int,
    include_mirror: bool = False,
) -> TerminatorGrid:
    """Build a thin spherical terminator grid in units of planetary radius."""

    lon = np.linspace(-0.5 * np.pi, 0.5 * np.pi, n_longitudes)
    lat = np.linspace(-0.5 * np.pi, 0.5 * np.pi, n_latitudes)
    lon2, lat2 = np.meshgrid(lon, lat, indexing="xy")

    x = np.sin(lon2) * np.cos(lat2)
    y = np.cos(lon2) * np.cos(lat2)
    z = np.sin(lat2)
    positions = np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    dlon = lon[1] - lon[0] if n_longitudes > 1 else np.pi
    dlat = lat[1] - lat[0] if n_latitudes > 1 else np.pi
    cell_area = np.cos(lat2).clip(min=0.0).ravel() * dlon * dlat
    mu = np.cos(lon2).clip(min=0.0).ravel()
    weights = mu * cell_area

    if include_mirror:
        positions = np.vstack([positions, positions * np.array([1.0, -1.0, 1.0])])
        lon2 = np.concatenate([lon2.ravel(), lon2.ravel()])
        lat2 = np.concatenate([lat2.ravel(), lat2.ravel()])
        cell_area = np.concatenate([cell_area, cell_area])
        weights = np.concatenate([weights, weights])
        weights *= 0.5
        weights /= np.sum(weights)
        return TerminatorGrid(lon2, lat2, positions, weights, cell_area)

    weights /= np.sum(weights)

    return TerminatorGrid(lon2.ravel(), lat2.ravel(), positions, weights, cell_area)
