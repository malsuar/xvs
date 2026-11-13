import numpy as np

from xds import MoonConfig, PlanetConfig, SimulationConfig, run_xds_model
from xds.shadow import compute_shadow_masks
from xds.terminator import build_projected_annulus


def test_radial_sampling_preserves_cell_areas_and_bounds():
    base=build_projected_annulus(120,8,.03)
    for phase in [.25,.5,.75]:
        grid=build_projected_annulus(120,8,.03,radial_sample_fraction=phase)
        assert np.array_equal(grid.weights,base.weights)
        assert np.array_equal(grid.cell_area,base.cell_area)
        r2=np.sum(grid.positions**2,axis=1).reshape(8,120)
        edges=np.linspace(1.,1.03**2,9)
        assert np.all(r2>edges[:-1,None]) and np.all(r2<edges[1:,None])


def test_area_radii_and_velocity_second_moment():
    h = 0.1
    grid = build_projected_annulus(360, 16, h)
    r2 = np.sum(grid.positions**2, axis=1)
    assert np.all(r2 > 1) and np.all(r2 < (1+h)**2)
    assert np.allclose(grid.positions[:, 0], 0)
    assert np.isclose(grid.cell_area.sum(), np.pi*((1+h)**2-1))
    assert np.isclose(grid.weights.sum(), 1)
    assert np.allclose(grid.weights @ grid.positions, 0, atol=1e-14)
    # For uniform projected area, <y^2> = (Rout^2+Rin^2)/4.
    assert np.isclose(grid.weights @ grid.positions[:, 1]**2, ((1+h)**2+1)/4)


def test_observer_basis_and_known_occulted_area():
    grid = build_projected_annulus(720, 256, 0.1, (1, 2, 3))
    assert np.allclose(grid.positions @ np.array([1, 2, 3]), 0, atol=1e-14)
    grid = build_projected_annulus(720, 256, 0.1)
    # A concentric disk with radius inside the annulus removes an exact area.
    expected = (1.05**2-1)/(1.1**2-1)
    for x, kind in [(3., "illumination"), (-3., "observer")]:
        masks = compute_shadow_masks(grid.positions, np.array([x, 0, 0]),
                                     1.05, (-1, 0, 0), (-1, 0, 0))
        assert abs(grid.weights @ getattr(masks, kind) - expected) < 1/256


def test_annulus_forward_model_reference_and_deficit():
    sim = SimulationConfig(geometry="projected_annulus", n_azimuths=120,
                           n_radii=4, n_times=7, n_velocities=101,
                           use_illumination_shadow=False, use_observer_occultation=False)
    result = run_xds_model(PlanetConfig(), MoonConfig(), sim)
    assert np.allclose(result.absorption, result.reference_absorption)
    assert np.allclose(result.reference_absorption, result.reference_absorption[:, ::-1])
    sim = SimulationConfig(geometry="projected_annulus", n_azimuths=120,
                           n_radii=4, n_times=21, n_velocities=101)
    result = run_xds_model(PlanetConfig(), MoonConfig(radius=0.3), sim)
    assert np.max(result.combined_fraction) > 0
    assert np.all(result.delta_absorption <= 1e-14)
