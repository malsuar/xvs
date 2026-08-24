import numpy as np

from xds.config import MoonConfig, PlanetConfig, SimulationConfig
from xds.spectra import run_xds_model


def test_disabling_both_masks_recovers_reference_profile():
    result = run_xds_model(
        PlanetConfig(),
        MoonConfig(),
        SimulationConfig(
            n_times=5,
            n_longitudes=21,
            n_latitudes=15,
            n_velocities=41,
            use_illumination_shadow=False,
            use_observer_occultation=False,
        ),
    )
    assert np.allclose(result.absorption, result.reference_absorption)


def test_line_depth_scales_absorption_and_residual_linearly():
    common_sim = SimulationConfig(
        n_times=7,
        n_longitudes=21,
        n_latitudes=15,
        n_velocities=41,
    )
    moon = MoonConfig(radius=0.3, semi_major_axis=2.5, period=5.0, phase=0.1)

    unit = run_xds_model(PlanetConfig(line_depth=1.0), moon, common_sim)
    scaled = run_xds_model(PlanetConfig(line_depth=2.5), moon, common_sim)

    assert np.allclose(scaled.reference_absorption, 2.5 * unit.reference_absorption)
    assert np.allclose(scaled.absorption, 2.5 * unit.absorption)
    assert np.allclose(scaled.delta_absorption, 2.5 * unit.delta_absorption)
