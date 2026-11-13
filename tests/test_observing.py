from types import SimpleNamespace

import numpy as np

from xds.config import InstrumentConfig, ObservationConfig
from xds.observing import observe_residual


def synthetic(nt, nv):
    t = np.linspace(0, 1, nt)
    v = np.linspace(-30, 30, nv)
    return SimpleNamespace(time=t, velocity=v,
                           delta_absorption=-(1+t[:, None])*np.exp(-v[None, :]**2/8))


def test_sampling_convergence():
    obs = ObservationConfig(exposure_time=0.1)
    scores = []
    for nt, nv in [(101, 601), (201, 1201), (401, 2401)]:
        measured = observe_residual(synthetic(nt, nv), 1., obs, InstrumentConfig())
        scores.append(np.linalg.norm(measured.signal)/measured.noise)
        assert measured.signal.shape == (10, 20)
    assert np.allclose(scores, scores[-1], rtol=1e-4)


def test_gaussian_lsf_bin_integral_and_linear_time_average():
    from math import erf
    result = synthetic(101, 2401)
    measured = observe_residual(result, 1., ObservationConfig(exposure_time=.1), InstrumentConfig())
    sigma = np.sqrt(4 + (measured.velocity_step/np.sqrt(8*np.log(2)))**2)
    a = measured.velocity-measured.velocity_step/2
    b = measured.velocity+measured.velocity_step/2
    exact = 2*np.sqrt(np.pi/2)/measured.velocity_step * np.array([
        erf(y/(np.sqrt(2)*sigma))-erf(x/(np.sqrt(2)*sigma)) for x,y in zip(a,b)])
    assert np.allclose(measured.signal, -(1+measured.time[:, None])*exact, atol=1e-7)


def test_photon_noise_fixed_duration_exposure_invariance():
    result = synthetic(101, 601)
    result.delta_absorption[:] = result.delta_absorption[0]
    scores = []
    for exp in [.05, .1, .2]:
        measured = observe_residual(result, 1., ObservationConfig(exposure_time=exp), InstrumentConfig())
        scores.append(np.linalg.norm(measured.signal)/measured.noise)
    assert np.allclose(scores, scores[0])
