"""Heuristic XDS efficiency metrics."""

from __future__ import annotations

import numpy as np


def hill_radius_over_planet_radius(
    planet_period_days: np.ndarray,
    planet_mass_jupiter: float = 1.0,
    star_mass_solar: float = 1.0,
    planet_radius_jupiter: float = 1.0,
) -> np.ndarray:
    """Estimate Hill radius in planet radii for circular star-planet orbits."""

    period_years = np.asarray(planet_period_days, dtype=float) / 365.25
    a_au = (period_years**2 * star_mass_solar) ** (1.0 / 3.0)
    mass_ratio = planet_mass_jupiter * 0.0009543 / (3.0 * star_mass_solar)
    hill_au = a_au * mass_ratio ** (1.0 / 3.0)
    rj_au = 0.0004778945 * planet_radius_jupiter
    return hill_au / rj_au


def xds_efficiency_metric(
    planet_period_days: np.ndarray,
    spin_period_days: np.ndarray | float = 0.42,
    refraction_period_scale_days: float = 1000.0,
) -> np.ndarray:
    """Conceptual normalised efficiency curve for exploratory plots."""

    period = np.asarray(planet_period_days, dtype=float)
    hill = hill_radius_over_planet_radius(period)
    f_hill = 1.0 - np.exp(-hill / 30.0)
    f_spin = 1.0 / (1.0 + (np.asarray(spin_period_days) / 2.0) ** 2)
    f_atm = np.exp(-period / refraction_period_scale_days)
    metric = f_hill * f_spin * f_atm
    return metric / np.nanmax(metric)
