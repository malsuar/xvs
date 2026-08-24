"""Simple dynamical scalings for XDS parameter searches."""

from __future__ import annotations

import numpy as np

from .detectability import (
    GRAVITATIONAL_CONSTANT,
    JUPITER_MASS_KG,
    JUPITER_RADIUS_KM,
    SOLAR_RADIUS_KM,
)
from .efficiency import hill_radius_over_planet_radius

SECONDS_PER_HOUR = 3_600.0


def equatorial_velocity_km_s(radius_jupiter: np.ndarray | float, rotation_period_hours: np.ndarray | float) -> np.ndarray:
    """Return equatorial rotation speed in km/s."""

    radius_km = np.asarray(radius_jupiter, dtype=float) * JUPITER_RADIUS_KM
    period_s = np.asarray(rotation_period_hours, dtype=float) * SECONDS_PER_HOUR
    return 2.0 * np.pi * radius_km / period_s


def breakup_rotation_period_hours(mass_jupiter: np.ndarray | float, radius_jupiter: np.ndarray | float) -> np.ndarray:
    """Critical spin period where centrifugal and gravitational accelerations match."""

    mass = np.asarray(mass_jupiter, dtype=float) * JUPITER_MASS_KG
    radius_m = np.asarray(radius_jupiter, dtype=float) * JUPITER_RADIUS_KM * 1_000.0
    period_s = 2.0 * np.pi * np.sqrt(radius_m**3 / (GRAVITATIONAL_CONSTANT * mass))
    return period_s / SECONDS_PER_HOUR


def moon_period_hours(
    semi_major_axis_rp: np.ndarray | float,
    planet_mass_jupiter: np.ndarray | float,
    planet_radius_jupiter: np.ndarray | float,
) -> np.ndarray:
    """Circular satellite period around the planet."""

    axis_m = (
        np.asarray(semi_major_axis_rp, dtype=float)
        * np.asarray(planet_radius_jupiter, dtype=float)
        * JUPITER_RADIUS_KM
        * 1_000.0
    )
    mass = np.asarray(planet_mass_jupiter, dtype=float) * JUPITER_MASS_KG
    period_s = 2.0 * np.pi * np.sqrt(axis_m**3 / (GRAVITATIONAL_CONSTANT * mass))
    return period_s / SECONDS_PER_HOUR


def fluid_roche_limit_rp(planet_density_g_cm3: np.ndarray | float, moon_density_g_cm3: float = 3.0) -> np.ndarray:
    """Fluid Roche limit in planetary radii."""

    return 2.456 * (np.asarray(planet_density_g_cm3, dtype=float) / moon_density_g_cm3) ** (1.0 / 3.0)


def planet_density_g_cm3(mass_jupiter: np.ndarray | float, radius_jupiter: np.ndarray | float) -> np.ndarray:
    """Mean density in g/cm^3."""

    jupiter_density = 1.326
    return jupiter_density * np.asarray(mass_jupiter, dtype=float) / np.asarray(radius_jupiter, dtype=float) ** 3


def stable_outer_moon_limit_rp(
    planet_period_days: np.ndarray | float,
    planet_mass_jupiter: np.ndarray | float,
    star_mass_solar: float = 1.0,
    planet_radius_jupiter: np.ndarray | float = 1.0,
    stability_fraction: float = 0.35,
) -> np.ndarray:
    """Approximate prograde stability limit as a fraction of Hill radius."""

    hill = hill_radius_over_planet_radius(
        np.asarray(planet_period_days, dtype=float),
        planet_mass_jupiter=float(np.asarray(planet_mass_jupiter).flat[0])
        if np.asarray(planet_mass_jupiter).size == 1
        else 1.0,
        star_mass_solar=star_mass_solar,
        planet_radius_jupiter=1.0,
    )
    # Recompute explicitly to support array masses/radii.
    period_years = np.asarray(planet_period_days, dtype=float) / 365.25
    a_au = (period_years**2 * star_mass_solar) ** (1.0 / 3.0)
    mass_ratio = np.asarray(planet_mass_jupiter, dtype=float) * 0.0009543 / (3.0 * star_mass_solar)
    hill_au = a_au * mass_ratio ** (1.0 / 3.0)
    rj_au = 0.0004778945 * np.asarray(planet_radius_jupiter, dtype=float)
    hill = hill_au / rj_au
    return stability_fraction * hill


def transit_duration_hours(
    planet_period_days: np.ndarray | float,
    stellar_radius_solar: np.ndarray | float,
    star_mass_solar: float = 1.0,
    impact_parameter: float = 0.0,
) -> np.ndarray:
    """Central-transit duration for a circular orbit, ignoring planet radius."""

    period_days = np.asarray(planet_period_days, dtype=float)
    period_years = period_days / 365.25
    a_au = (period_years**2 * star_mass_solar) ** (1.0 / 3.0)
    rstar_au = np.asarray(stellar_radius_solar, dtype=float) * SOLAR_RADIUS_KM / 149_597_870.7
    chord = 2.0 * rstar_au * np.sqrt(max(0.0, 1.0 - impact_parameter**2))
    duration_days = period_days / np.pi * np.arcsin(np.clip(chord / (2.0 * a_au), 0.0, 1.0))
    return 24.0 * duration_days
