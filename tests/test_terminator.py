import numpy as np

from xds.terminator import build_terminator


def test_mirrored_terminator_has_two_projected_branches():
    single = build_terminator(11, 9)
    mirrored = build_terminator(11, 9, include_mirror=True)

    assert mirrored.positions.shape[0] == 2 * single.positions.shape[0]
    assert np.any(mirrored.positions[:, 1] < 0.0)
    assert np.any(mirrored.positions[:, 1] > 0.0)
    assert np.isclose(np.sum(mirrored.weights), 1.0)
