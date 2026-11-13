"""Explicit numerical starting settings; acceptance still requires refinement."""
from dataclasses import dataclass
import numpy as np
from .config import PlanetConfig, SimulationConfig, InstrumentConfig, ObservationConfig
from .observability import SPEED_OF_LIGHT_KM_S


@dataclass(frozen=True)
class NumericalTolerances:
    spatial_snr_relative: float = .01
    temporal_snr_relative: float = .01
    spectral_snr_relative: float = .01
    integrated_area_relative: float = .01
    exposure_phase_fractions: tuple = (0., .25, .5, .75)


def reference_simulation(planet: PlanetConfig, duration_hours: float,
                         observation: ObservationConfig,
                         instrument: InstrumentConfig) -> SimulationConfig:
    """Construct initial grids for a uniform annulus; does not certify accuracy.

    Recheck convergence for each new footprint/velocity field. Spin and orbital
    angles must be explicitly chosen by the caller in the observer convention.
    """
    if duration_hours<=0 or planet.line_sigma<=0 or planet.v_eq<0:
        raise ValueError("positive duration/width and nonnegative speed required")
    if observation.exposure_time<=0 or instrument.spectral_resolution<=0:
        raise ValueError("positive exposure and resolution required")
    dv_inst=SPEED_OF_LIGHT_KM_S/instrument.spectral_resolution
    step=min(.2,planet.line_sigma/2,dv_inst/8)
    # Eight combined Gaussian widths beyond the maximum rotating radius.
    width=np.hypot(planet.line_sigma,dv_inst/np.sqrt(8*np.log(2)))
    limit=max(40.,float(np.ceil(planet.v_eq*(1+planet.atmosphere_half_width)+8*width)))
    target_dt=observation.exposure_time/25
    return SimulationConfig(geometry='projected_annulus',spin_convention='observer',
        blocking_integration='fractional',blocking_quadrature_order=16,
        n_azimuths=480,n_radii=16,transit_duration=duration_hours,
        n_times=int(np.ceil(duration_hours/target_dt))+1,
        velocity_min=-limit,velocity_max=limit,
        n_velocities=int(np.ceil(2*limit/step))+1)
