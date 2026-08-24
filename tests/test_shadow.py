import numpy as np

from xds.shadow import illumination_shadow_mask, observer_occultation_mask


def test_illumination_and_observer_masks_are_distinct():
    points = np.array([[0.0, 0.0, 0.0]])
    star_direction = (-1.0, 0.0, 0.0)
    observer_direction = (-1.0, 0.0, 0.0)

    upstream_moon = np.array([0.5, 0.0, 0.0])
    downstream_moon = np.array([-0.5, 0.0, 0.0])

    assert illumination_shadow_mask(points, upstream_moon, 0.1, star_direction)[0]
    assert not observer_occultation_mask(points, upstream_moon, 0.1, observer_direction)[0]
    assert not illumination_shadow_mask(points, downstream_moon, 0.1, star_direction)[0]
    assert observer_occultation_mask(points, downstream_moon, 0.1, observer_direction)[0]
