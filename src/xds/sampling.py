"""Sampling utilities for XDS parameter searches."""

from __future__ import annotations

import numpy as np


PARAMETER_RANGES = {
    "planet_mass_jupiter": ("log", 0.08, 2.0),
    "planet_radius_jupiter": ("linear", 0.75, 1.75),
    "stellar_radius_solar": ("linear", 0.45, 1.05),
    "planet_period_days": ("log", 3.0, 300.0),
    "rotation_period_hours": ("log", 3.0, 24.0),
    "moon_separation_rp": ("linear", 1.8, 12.0),
    "moon_radius_rp": ("linear", 0.12, 0.35),
}


def latin_hypercube_unit(n: int, dimensions: int, seed: int = 7) -> np.ndarray:
    """Generate a simple Latin hypercube in the unit cube."""

    rng = np.random.default_rng(seed)
    sample = np.empty((n, dimensions), dtype=float)
    base = (np.arange(n, dtype=float) + rng.random((dimensions, n))) / n
    for dim in range(dimensions):
        sample[:, dim] = rng.permutation(base[dim])
    return sample


def scale_unit_interval(unit_values: np.ndarray, mode: str, lower: float, upper: float) -> np.ndarray:
    if mode == "log":
        return 10.0 ** (np.log10(lower) + unit_values * (np.log10(upper) - np.log10(lower)))
    if mode == "linear":
        return lower + unit_values * (upper - lower)
    raise ValueError(f"unknown scaling mode: {mode}")


def latin_hypercube_candidates(n: int, seed: int = 7) -> dict[str, np.ndarray]:
    names = list(PARAMETER_RANGES)
    unit = latin_hypercube_unit(n, len(names), seed=seed)
    return {
        name: scale_unit_interval(unit[:, idx], *PARAMETER_RANGES[name])
        for idx, name in enumerate(names)
    }
