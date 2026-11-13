import numpy as np
from xds import MoonConfig
from xds.moon_orbit import moon_positions
from xds.terminator import build_projected_annulus
from xds.shadow import compute_shadow_masks


def test_constructed_impact_and_projected_union():
    grid = build_projected_annulus(240,8,.03)
    for b in [0.,1.08,1.14]:
        moon=MoonConfig(semi_major_axis=3.,radius=.12,phase=.25,
                        inclination=np.rad2deg(np.arcsin(b/3)),longitude_node=90.)
        points=moon_positions(np.array([0.,.3,-.3]),moon)
        assert np.isclose(np.linalg.norm(points[0,1:]),b)
        for p in points:
            masks=compute_shadow_masks(grid.positions,p,.12,(-1,0,0),(-1,0,0))
            projected=np.sum((grid.positions[:,1:]-p[1:])**2,axis=1)<=.12**2
            assert np.array_equal(masks.combined,projected)
