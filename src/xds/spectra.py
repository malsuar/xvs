"""Forward synthesis of Exomoon Doppler Shadow maps."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MoonConfig, PlanetConfig, SimulationConfig
from .moon_orbit import moon_positions
from .rotation import line_of_sight_velocity
from .shadow import compute_shadow_masks
from .terminator import TerminatorGrid, build_terminator, build_projected_annulus


@dataclass(frozen=True)
class XDSResult:
    time: np.ndarray
    velocity: np.ndarray
    reference_absorption: np.ndarray
    absorption: np.ndarray
    delta_absorption: np.ndarray
    moon_positions: np.ndarray
    illumination_fraction: np.ndarray
    observer_fraction: np.ndarray
    combined_fraction: np.ndarray
    terminator: TerminatorGrid
    v_los: np.ndarray


def gaussian_profiles(velocity: np.ndarray, centers: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((velocity[None, :] - centers[:, None]) / sigma) ** 2)


def absorption_profile(
    velocity: np.ndarray,
    v_los: np.ndarray,
    weights: np.ndarray,
    line_depth: float,
    line_sigma: float,
) -> np.ndarray:
    kernels = gaussian_profiles(velocity, v_los, line_sigma)
    return weighted_absorption(kernels, weights, line_depth)


def weighted_absorption(kernels: np.ndarray, weights: np.ndarray, line_depth: float) -> np.ndarray:
    return line_depth * (np.asarray(weights, dtype=float) @ kernels)


def run_xds_model(
    planet: PlanetConfig,
    moon: MoonConfig,
    sim: SimulationConfig,
) -> XDSResult:
    """Run the velocity-resolved XDS forward model."""

    if sim.geometry == "legacy":
        terminator = build_terminator(
            sim.n_longitudes, sim.n_latitudes,
            include_mirror=sim.include_mirrored_terminator,
        )
    elif sim.geometry == "projected_annulus":
        if sim.include_mirrored_terminator:
            raise ValueError("mirroring applies only to legacy geometry")
        terminator = build_projected_annulus(
            sim.n_azimuths, sim.n_radii, planet.atmosphere_half_width,
            sim.observer_direction,
        )
    else:
        raise ValueError(f"Unknown geometry: {sim.geometry}")
    time = np.linspace(-0.5 * sim.transit_duration, 0.5 * sim.transit_duration, sim.n_times)
    velocity = np.linspace(sim.velocity_min, sim.velocity_max, sim.n_velocities)
    v_los = line_of_sight_velocity(
        terminator.positions,
        planet.v_eq,
        planet.spin_inclination,
        planet.spin_obliquity,
        sim.observer_direction,
    )

    kernels = gaussian_profiles(velocity, v_los, planet.line_sigma)
    reference = weighted_absorption(kernels, terminator.weights, planet.line_depth)
    positions = moon_positions(time, moon)

    absorption = np.empty((time.size, velocity.size))
    illumination_fraction = np.empty(time.size)
    observer_fraction = np.empty(time.size)
    combined_fraction = np.empty(time.size)

    for idx, moon_pos in enumerate(positions):
        masks = compute_shadow_masks(
            terminator.positions,
            moon_pos,
            moon.radius,
            sim.star_direction,
            sim.observer_direction,
        )
        active = np.zeros(terminator.weights.size, dtype=bool)
        if sim.use_illumination_shadow:
            active |= masks.illumination
        if sim.use_observer_occultation:
            active |= masks.observer

        weights = np.where(active, 0.0, terminator.weights)
        absorption[idx] = weighted_absorption(kernels, weights, planet.line_depth)
        illumination_fraction[idx] = np.sum(terminator.weights[masks.illumination])
        observer_fraction[idx] = np.sum(terminator.weights[masks.observer])
        combined_fraction[idx] = np.sum(terminator.weights[active])

    return XDSResult(
        time=time,
        velocity=velocity,
        reference_absorption=np.repeat(reference[None, :], time.size, axis=0),
        absorption=absorption,
        delta_absorption=absorption - reference[None, :],
        moon_positions=positions,
        illumination_fraction=illumination_fraction,
        observer_fraction=observer_fraction,
        combined_fraction=combined_fraction,
        terminator=terminator,
        v_los=v_los,
    )
