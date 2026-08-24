import numpy as np

from xds.rotation import line_of_sight_velocity


def test_velocity_field_changes_sign_across_terminator():
    positions = np.array([[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    velocity = line_of_sight_velocity(positions, 10.0, 90.0, 0.0, (-1.0, 0.0, 0.0))
    assert np.isclose(velocity[0], -velocity[1])
