import numpy as np

from xds.rotation import line_of_sight_velocity, spin_axis


def test_observer_inclination_and_projected_angle():
    for observer in [(-1., 0., 0.), (1., 2., 3.)]:
        n = np.array(observer)/np.linalg.norm(observer)
        for inc in [0., 30., 90., 180.]:
            for angle in [0., 45., 90.]:
                axis = spin_axis(inc, angle, observer, "observer")
                assert np.isclose(axis @ n, np.cos(np.deg2rad(inc)))
    assert np.allclose(spin_axis(90, 0, convention="observer"), [0, 0, 1])
    assert np.allclose(spin_axis(90, 90, convention="observer"), [0, -1, 0])


def test_velocity_inclination_scaling_and_sign():
    points = np.array([[0., 1., 0.], [0., -1., 0.]])
    for inc in [0., 30., 90., 180.]:
        vel = line_of_sight_velocity(points, 12., inc, 0., (-1, 0, 0), "observer")
        assert np.allclose(vel, 12*np.sin(np.deg2rad(inc))*np.array([1, -1]))
    # Legacy axis is preserved exactly for historical configurations.
    assert np.allclose(spin_axis(90, 0), [0, -1, 0])


def test_velocity_field_changes_sign_across_terminator():
    positions = np.array([[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    velocity = line_of_sight_velocity(positions, 10.0, 90.0, 0.0, (-1.0, 0.0, 0.0))
    assert np.isclose(velocity[0], -velocity[1])
