"""Time-series diagnostics extracted from XDS maps."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class XDSDiagnostics:
    centroid: np.ndarray
    skewness: np.ndarray
    blue_signal: np.ndarray
    red_signal: np.ndarray
    shadow_velocity: np.ndarray
    shadow_strength: np.ndarray


def map_diagnostics(time: np.ndarray, velocity: np.ndarray, delta_absorption: np.ndarray) -> XDSDiagnostics:
    """Compute compact diagnostics from the differential map."""

    deficit = np.clip(-delta_absorption, 0.0, None)
    total = np.sum(deficit, axis=1)
    safe_total = np.where(total > 0.0, total, np.nan)

    centroid = np.sum(deficit * velocity[None, :], axis=1) / safe_total
    variance = np.sum(deficit * (velocity[None, :] - centroid[:, None]) ** 2, axis=1) / safe_total
    sigma = np.sqrt(variance)
    skewness = (
        np.sum(deficit * (velocity[None, :] - centroid[:, None]) ** 3, axis=1)
        / safe_total
        / np.where(sigma > 0.0, sigma**3, np.nan)
    )

    blue = velocity < 0.0
    red = velocity > 0.0
    blue_signal = np.trapz(delta_absorption[:, blue], velocity[blue], axis=1)
    red_signal = np.trapz(delta_absorption[:, red], velocity[red], axis=1)

    idx = np.argmax(deficit, axis=1)
    shadow_velocity = velocity[idx]
    shadow_strength = deficit[np.arange(deficit.shape[0]), idx]
    shadow_velocity = np.where(total > 0.0, shadow_velocity, np.nan)
    shadow_strength = np.where(total > 0.0, shadow_strength, 0.0)

    return XDSDiagnostics(
        centroid=centroid,
        skewness=skewness,
        blue_signal=blue_signal,
        red_signal=red_signal,
        shadow_velocity=shadow_velocity,
        shadow_strength=shadow_strength,
    )


def finite_difference_rate(time: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Return d(values)/dt, preserving NaNs outside active intervals."""

    values = np.asarray(values, dtype=float)
    rates = np.full_like(values, np.nan)
    valid = np.isfinite(values)
    if np.count_nonzero(valid) < 2:
        return rates

    valid_idx = np.flatnonzero(valid)
    dt = np.median(np.diff(time)) if time.size > 1 else np.inf
    breaks = np.where(np.diff(valid_idx) > 1)[0] + 1
    blocks = np.split(valid_idx, breaks)
    for block in blocks:
        if block.size < 2:
            continue
        block_time = time[block]
        if np.any(np.diff(block_time) > 1.5 * dt):
            sub_breaks = np.where(np.diff(block_time) > 1.5 * dt)[0] + 1
            for sub_block in np.split(block, sub_breaks):
                if sub_block.size >= 2:
                    rates[sub_block] = np.gradient(values[sub_block], time[sub_block])
            continue
        rates[block] = np.gradient(values[block], block_time)
    return rates
