"""Cadence and exposure-time estimates for XDS observability."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import InstrumentConfig
from .diagnostics import finite_difference_rate, map_diagnostics
from .spectra import XDSResult

SPEED_OF_LIGHT_KM_S = 299_792.458


@dataclass(frozen=True)
class ContactInterval:
    kind: str
    start: float
    end: float
    duration: float
    peak_fraction: float


@dataclass(frozen=True)
class CadenceReport:
    intervals: list[ContactInterval]
    velocity_resolution: float
    min_contact_duration: float
    max_shadow_drift_rate: float
    cadence_from_contacts: float
    max_integration_from_smearing: float
    recommended_cadence: float
    recommended_integration_time: float


def contact_intervals(
    time: np.ndarray,
    fraction: np.ndarray,
    kind: str,
    threshold: float = 1.0e-5,
) -> list[ContactInterval]:
    """Convert a shadowed-fraction time series into contiguous intervals."""

    active = np.asarray(fraction) > threshold
    intervals: list[ContactInterval] = []
    if not np.any(active):
        return intervals

    starts: list[int] = []
    ends: list[int] = []
    in_event = False
    for idx, is_active in enumerate(active):
        if is_active and not in_event:
            starts.append(idx)
            in_event = True
        if in_event and (not is_active or idx == active.size - 1):
            end_idx = idx - 1 if not is_active else idx
            ends.append(end_idx)
            in_event = False

    for start_idx, end_idx in zip(starts, ends):
        start = float(time[start_idx])
        end = float(time[end_idx])
        intervals.append(
            ContactInterval(
                kind=kind,
                start=start,
                end=end,
                duration=end - start,
                peak_fraction=float(np.max(fraction[start_idx : end_idx + 1])),
            )
        )
    return intervals


def cadence_report(
    result: XDSResult,
    instrument: InstrumentConfig,
    threshold: float = 1.0e-5,
) -> CadenceReport:
    """Estimate cadence and exposure times needed to resolve the XDS signal."""

    diagnostics = map_diagnostics(result.time, result.velocity, result.delta_absorption)
    # The deepest-pixel velocity is useful for plotting but can jump between
    # neighbouring velocity bins. Exposure smearing is better tied to the
    # deficit centroid, which tracks the bulk shadow motion continuously.
    measurable = result.combined_fraction > max(threshold, 0.1 * np.nanmax(result.combined_fraction))
    centroid_for_drift = np.where(measurable, diagnostics.centroid, np.nan)
    drift = finite_difference_rate(result.time, centroid_for_drift)
    max_drift = float(np.nanmax(np.abs(drift))) if np.any(np.isfinite(drift)) else 0.0
    velocity_resolution = SPEED_OF_LIGHT_KM_S / instrument.spectral_resolution

    intervals = []
    intervals.extend(contact_intervals(result.time, result.illumination_fraction, "illumination", threshold))
    intervals.extend(contact_intervals(result.time, result.observer_fraction, "observer", threshold))
    intervals.extend(contact_intervals(result.time, result.combined_fraction, "combined", threshold))

    durations = [interval.duration for interval in intervals if interval.duration > 0.0]
    min_contact_duration = float(min(durations)) if durations else float("nan")
    cadence_from_contacts = (
        min_contact_duration / instrument.samples_per_contact if durations else float("nan")
    )

    if max_drift > 0.0:
        max_integration_from_smearing = (
            instrument.smear_fraction * velocity_resolution / max_drift
        )
    else:
        max_integration_from_smearing = float("inf")

    candidates = [
        value
        for value in [cadence_from_contacts, max_integration_from_smearing]
        if np.isfinite(value) and value > 0.0
    ]
    recommended = float(min(candidates)) if candidates else float("nan")

    return CadenceReport(
        intervals=intervals,
        velocity_resolution=velocity_resolution,
        min_contact_duration=min_contact_duration,
        max_shadow_drift_rate=max_drift,
        cadence_from_contacts=cadence_from_contacts,
        max_integration_from_smearing=max_integration_from_smearing,
        recommended_cadence=recommended,
        recommended_integration_time=recommended,
    )
