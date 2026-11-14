"""Physical amplitude and SNR estimates for Exomoon Doppler Shadows."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import AtmosphereConfig, InstrumentConfig, ObservationConfig
from .diagnostics import finite_difference_rate, map_diagnostics
from .observability import SPEED_OF_LIGHT_KM_S, cadence_report
from .spectra import XDSResult
from .observing import observe_residual

BOLTZMANN = 1.380649e-23
GRAVITATIONAL_CONSTANT = 6.67430e-11
ATOMIC_MASS_UNIT = 1.66053906660e-27
JUPITER_MASS_KG = 1.89813e27
JUPITER_RADIUS_KM = 71_492.0
SOLAR_RADIUS_KM = 695_700.0


@dataclass(frozen=True)
class DetectabilityReport:
    scale_height_km: float
    atmospheric_line_depth: float
    peak_xds_depth: float
    matched_filter_snr: float
    peak_pixel_snr: float
    n_effective_pixels: int
    exposure_time: float
    velocity_resolution: float
    smearing_penalty: float
    cadence_limited_exposure: float
    notes: tuple[str, ...]


@dataclass(frozen=True)
class DetectionRequirements:
    target_snr: float
    required_lines_one_transit: float
    required_transits: float
    current_snr: float


def surface_gravity_m_s2(atmosphere: AtmosphereConfig) -> float:
    mass = atmosphere.planet_mass_jupiter * JUPITER_MASS_KG
    radius_m = atmosphere.planet_radius_jupiter * JUPITER_RADIUS_KM * 1_000.0
    return GRAVITATIONAL_CONSTANT * mass / radius_m**2


def scale_height_km(atmosphere: AtmosphereConfig) -> float:
    gravity = surface_gravity_m_s2(atmosphere)
    particle_mass = atmosphere.mean_molecular_weight * ATOMIC_MASS_UNIT
    height_m = BOLTZMANN * atmosphere.temperature / (particle_mass * gravity)
    return height_m / 1_000.0


def atmospheric_line_depth(atmosphere: AtmosphereConfig) -> float:
    """Approximate line depth from an effective annulus in stellar-flux units."""

    rp = atmosphere.planet_radius_jupiter * JUPITER_RADIUS_KM
    rstar = atmosphere.stellar_radius_solar * SOLAR_RADIUS_KM
    height = atmosphere.effective_scale_heights * scale_height_km(atmosphere)
    return atmosphere.line_contrast * 2.0 * rp * height / rstar**2


def exposure_smearing_penalty(
    result: XDSResult,
    exposure_time: float,
    instrument: InstrumentConfig,
) -> tuple[float, float]:
    """Estimate amplitude loss from velocity drift during one exposure."""

    diagnostics = map_diagnostics(result.time, result.velocity, result.delta_absorption)
    active = result.combined_fraction > 0.1 * np.nanmax(result.combined_fraction)
    drift = finite_difference_rate(result.time, np.where(active, diagnostics.centroid, np.nan))
    max_drift = float(np.nanmax(np.abs(drift))) if np.any(np.isfinite(drift)) else 0.0
    dv = SPEED_OF_LIGHT_KM_S / instrument.spectral_resolution
    if max_drift == 0.0 or not np.isfinite(max_drift):
        return 1.0, float("inf")
    limiting_exposure = instrument.smear_fraction * dv / max_drift
    penalty = min(1.0, limiting_exposure / exposure_time)
    return float(penalty), float(limiting_exposure)


def detectability_report(
    result: XDSResult,
    atmosphere: AtmosphereConfig,
    observation: ObservationConfig,
    instrument: InstrumentConfig,
    active_threshold: float = 1.0e-5,
    method: str = "observed",
) -> DetectabilityReport:
    """Use independent observing bins, or explicitly request method='legacy'."""

    if method == "observed":
        depth = atmospheric_line_depth(atmosphere)
        observed = observe_residual(result, depth, observation, instrument)
        if observation.n_lines < 1 or observation.n_transits < 1:
            raise ValueError("line and transit counts must be positive")
        snr = np.linalg.norm(observed.signal/observed.noise) * np.sqrt(observation.n_lines*observation.n_transits)
        sampling_notes = []
        if np.max(np.diff(result.time)) > observation.exposure_time/2:
            sampling_notes.append("input_time_grid_coarse_relative_to_exposure_check_convergence")
        if np.max(np.diff(result.velocity)) > observed.velocity_step/4:
            sampling_notes.append("input_velocity_grid_coarse_relative_to_lsf_check_convergence")
        return DetectabilityReport(
            scale_height_km(atmosphere), depth, float(np.max(np.abs(observed.signal))),
            float(snr), float(np.max(np.abs(observed.signal))/observed.noise),
            observed.signal.size, observation.exposure_time, observed.velocity_step,
            1.0, float("nan"),
            ("exposure_and_lsf_integrated_no_extra_smearing_penalty",
             "independent_noise_noiseless_reference_equal_independent_lines",
             "partial_final_bins_discarded", "legacy_cadence_limit_not_evaluated") + tuple(sampling_notes),
        )
    if method != "legacy":
        raise ValueError(f"Unknown detectability method: {method}")

    line_depth = atmospheric_line_depth(atmosphere)
    signal = line_depth * result.delta_absorption
    active = result.combined_fraction > active_threshold
    active_map = active[:, None] & np.isfinite(signal)
    n_pixels = int(np.count_nonzero(active_map))

    base_noise = 1.0 / observation.continuum_snr_per_resolution_element
    if observation.systematic_floor > 0.0:
        noise = np.sqrt(base_noise**2 + observation.systematic_floor**2)
    else:
        noise = base_noise

    if n_pixels > 0:
        snr_single_line = np.sqrt(np.sum((signal[active_map] / noise) ** 2))
        peak_pixel_snr = float(np.nanmax(np.abs(signal[active_map])) / noise)
    else:
        snr_single_line = 0.0
        peak_pixel_snr = 0.0

    smear_penalty, limiting_exposure = exposure_smearing_penalty(
        result, observation.exposure_time, instrument
    )
    snr = (
        snr_single_line
        * np.sqrt(observation.n_lines)
        * np.sqrt(observation.n_transits)
        * smear_penalty
    )

    cadence = cadence_report(result, instrument, threshold=active_threshold)
    notes: list[str] = []
    if observation.exposure_time > limiting_exposure:
        notes.append("exposure_time_exceeds_smearing_limit")
    if np.isfinite(cadence.cadence_from_contacts) and observation.exposure_time > cadence.cadence_from_contacts:
        notes.append("exposure_time_under_samples_contact")
    if observation.systematic_floor > 0.0 and observation.systematic_floor > base_noise:
        notes.append("systematics_dominate_over_photon_noise")

    return DetectabilityReport(
        scale_height_km=scale_height_km(atmosphere),
        atmospheric_line_depth=line_depth,
        peak_xds_depth=float(np.nanmax(np.abs(signal))) if signal.size else 0.0,
        matched_filter_snr=float(snr),
        peak_pixel_snr=peak_pixel_snr,
        n_effective_pixels=n_pixels,
        exposure_time=observation.exposure_time,
        velocity_resolution=SPEED_OF_LIGHT_KM_S / instrument.spectral_resolution,
        smearing_penalty=smear_penalty,
        cadence_limited_exposure=limiting_exposure,
        notes=tuple(notes),
    )


def detection_requirements(
    report: DetectabilityReport,
    observation: ObservationConfig,
    target_snr: float = 5.0,
) -> DetectionRequirements:
    """Invert the SNR scaling to estimate required lines or transits."""

    if report.matched_filter_snr <= 0.0:
        return DetectionRequirements(
            target_snr=target_snr,
            required_lines_one_transit=float("inf"),
            required_transits=float("inf"),
            current_snr=report.matched_filter_snr,
        )

    scale = (target_snr / report.matched_filter_snr) ** 2
    required_lines = observation.n_lines * observation.n_transits * scale
    required_transits = observation.n_transits * scale
    return DetectionRequirements(
        target_snr=target_snr,
        required_lines_one_transit=float(required_lines),
        required_transits=float(required_transits),
        current_snr=report.matched_filter_snr,
    )
