"""Configuration dataclasses for XDS simulations."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PlanetConfig:
    """Planet and atmospheric-line parameters in planet-radius units."""

    radius: float = 1.0
    v_eq: float = 12.0
    spin_inclination: float = 90.0
    spin_obliquity: float = 0.0
    line_depth: float = 1.0
    line_sigma: float = 2.0
    atmosphere_half_width: float = 0.03


@dataclass(frozen=True)
class MoonConfig:
    """Satellite orbital and size parameters in planet-radius units."""

    radius: float = 0.25
    semi_major_axis: float = 3.0
    period: float = 5.0
    phase: float = 0.0
    inclination: float = 90.0
    longitude_node: float = 0.0
    argument_pericenter: float = 0.0
    eccentricity: float = 0.0


@dataclass(frozen=True)
class SimulationConfig:
    """Numerical grid and geometry configuration."""

    transit_duration: float = 5.0
    n_times: int = 181
    n_longitudes: int = 161
    n_latitudes: int = 121
    velocity_min: float = -25.0
    velocity_max: float = 25.0
    n_velocities: int = 401
    star_direction: tuple[float, float, float] = (-1.0, 0.0, 0.0)
    observer_direction: tuple[float, float, float] = (-1.0, 0.0, 0.0)
    use_illumination_shadow: bool = True
    use_observer_occultation: bool = True
    include_mirrored_terminator: bool = False
    geometry: str = "legacy"
    n_azimuths: int = 360
    n_radii: int = 8
    spin_convention: str = "legacy"


@dataclass(frozen=True)
class InstrumentConfig:
    """Observation parameters used for cadence and integration estimates."""

    spectral_resolution: float = 100_000.0
    samples_per_contact: int = 5
    smear_fraction: float = 0.5
    min_exposures_per_event: int = 3


@dataclass(frozen=True)
class AtmosphereConfig:
    """Bulk atmospheric parameters used to scale XDS amplitudes."""

    temperature: float = 500.0
    mean_molecular_weight: float = 2.3
    effective_scale_heights: float = 5.0
    line_contrast: float = 1.0
    planet_mass_jupiter: float = 1.0
    planet_radius_jupiter: float = 1.0
    stellar_radius_solar: float = 1.0


@dataclass(frozen=True)
class ObservationConfig:
    """Noise model for a high-resolution transmission sequence."""

    continuum_snr_per_resolution_element: float = 200.0
    exposure_time: float = 1.0 / 60.0
    n_transits: int = 1
    n_lines: int = 30
    systematic_floor: float = 0.0
