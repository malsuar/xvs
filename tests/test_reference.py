import numpy as np
from xds.config import PlanetConfig, ObservationConfig, InstrumentConfig
from xds.reference import reference_simulation, NumericalTolerances


def test_reference_grid_spacing_and_modes():
    sim=reference_simulation(PlanetConfig(v_eq=25,line_sigma=.2),4.,
                             ObservationConfig(exposure_time=5/60),InstrumentConfig())
    assert sim.blocking_integration=='fractional'
    assert sim.geometry=='projected_annulus' and sim.spin_convention=='observer'
    assert sim.transit_duration/(sim.n_times-1)<=5/60/25+1e-14
    assert (sim.velocity_max-sim.velocity_min)/(sim.n_velocities-1)<=.1+1e-14
    assert sim.velocity_max>=40
    assert NumericalTolerances().temporal_snr_relative==.01
