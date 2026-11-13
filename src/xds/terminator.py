"""Discretised atmospheric terminator."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import sky_basis


@dataclass(frozen=True)
class TerminatorGrid:
    longitude: np.ndarray
    latitude: np.ndarray
    positions: np.ndarray
    weights: np.ndarray
    cell_area: np.ndarray


def build_projected_annulus(
    n_azimuths: int,
    n_radii: int,
    thickness: float,
    observer_direction=(-1.0, 0.0, 0.0),
) -> TerminatorGrid:
    """Uniform effective annulus at projected radii [1, 1+thickness].

    Positions are sky-plane representatives, not 3D ray formation locations.
    Equal-area radial cells and midpoint azimuths avoid duplicated endpoints.
    Angular coordinates are stored in the legacy longitude/latitude fields.
    """
    if n_azimuths < 4 or n_radii < 1:
        raise ValueError("annulus requires at least four azimuths and one radial cell")
    if not np.isfinite(thickness) or thickness <= 0:
        raise ValueError("annulus thickness must be finite and positive")
    u, v, _ = sky_basis(observer_direction)
    edges_squared = np.linspace(1.0, (1.0 + thickness)**2, n_radii + 1)
    radii = np.sqrt(0.5 * (edges_squared[:-1] + edges_squared[1:]))
    angles = (np.arange(n_azimuths) + 0.5) * 2.0 * np.pi / n_azimuths
    radius, angle = np.meshgrid(radii, angles, indexing="ij")
    positions = (radius[..., None] *
                 (np.cos(angle)[..., None] * u + np.sin(angle)[..., None] * v))
    area = np.repeat(np.pi * np.diff(edges_squared) / n_azimuths, n_azimuths)
    weights = area / area.sum()
    return TerminatorGrid(angle.ravel(), np.zeros(angle.size),
                          positions.reshape(-1, 3), weights, area)


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
