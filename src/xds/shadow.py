"""Illumination and observer-line occultation masks from an orbiting moon."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import distance_to_ray, unit_vector


@dataclass(frozen=True)
class ShadowMasks:
    illumination: np.ndarray
    observer: np.ndarray

    @property
    def combined(self) -> np.ndarray:
        return self.illumination | self.observer


def illumination_shadow_mask(
    terminator_positions: np.ndarray,
    moon_position: np.ndarray,
    moon_radius: float,
    star_to_planet_direction: np.ndarray | tuple[float, float, float],
) -> np.ndarray:
    """Mask elements whose incoming stellar ray intersects the moon."""

    ray_axis = unit_vector(star_to_planet_direction)
    distances, along = distance_to_ray(terminator_positions, moon_position, ray_axis)
    moon_is_upstream = along > 0.0
    return (distances <= moon_radius) & moon_is_upstream


def observer_occultation_mask(
    terminator_positions: np.ndarray,
    moon_position: np.ndarray,
    moon_radius: float,
    terminator_to_observer_direction: np.ndarray | tuple[float, float, float],
) -> np.ndarray:
    """Mask elements whose line of sight to the observer intersects the moon."""

    ray_axis = unit_vector(terminator_to_observer_direction)
    distances, along = distance_to_ray(terminator_positions, moon_position, ray_axis)
    moon_is_downstream = along < 0.0
    return (distances <= moon_radius) & moon_is_downstream


def compute_shadow_masks(
    terminator_positions: np.ndarray,
    moon_position: np.ndarray,
    moon_radius: float,
    star_direction: np.ndarray | tuple[float, float, float],
    observer_direction: np.ndarray | tuple[float, float, float],
) -> ShadowMasks:
    return ShadowMasks(
        illumination=illumination_shadow_mask(
            terminator_positions, moon_position, moon_radius, star_direction
        ),
        observer=observer_occultation_mask(
            terminator_positions, moon_position, moon_radius, observer_direction
        ),
    )
