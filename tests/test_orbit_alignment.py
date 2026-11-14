import numpy as np
from xds import MoonConfig
from xds.moon_orbit import orbital_normal, spin_orbit_angle, moon_positions


def test_orbital_plane_normal_and_direction():
    for inc, node in [(0., 0.), (90., 0.), (43., 27.), (180., 0.)]:
        moon = MoonConfig(inclination=inc, longitude_node=node)
        normal = orbital_normal(moon)
        p = moon_positions(np.array([0., 1e-4]), moon)
        assert np.allclose(p @ normal, 0, atol=1e-14)
        assert np.cross(p[0], p[1]) @ normal > 0
    assert np.isclose(spin_orbit_angle(MoonConfig(inclination=0), [0,0,1]), 0)
    assert np.isclose(spin_orbit_angle(MoonConfig(inclination=90), [0,0,1]), 90)
    assert np.isclose(spin_orbit_angle(MoonConfig(inclination=180), [0,0,1]), 180)
