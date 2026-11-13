"""Exposure integration and Gaussian instrumental response on physical bins."""

from dataclasses import dataclass
from math import erf

import numpy as np

from .config import InstrumentConfig, ObservationConfig
from .observability import SPEED_OF_LIGHT_KM_S


@dataclass(frozen=True)
class ObservedResidual:
    time: np.ndarray
    velocity: np.ndarray
    signal: np.ndarray
    noise: float
    time_step: float
    velocity_step: float


def observe_residual(result, amplitude: float, observation: ObservationConfig,
                     instrument: InstrumentConfig) -> ObservedResidual:
    """Bin a continuous residual; noise assumes fixed flux per velocity and time.

    Continuum S/N is specified at the reference exposure and resolution.
    The template is noiseless and lines/transits are independent equivalents.
    """
    t, v = np.asarray(result.time), np.asarray(result.velocity)
    exp = observation.exposure_time
    resolution = instrument.spectral_resolution
    if (exp <= 0 or resolution <= 0 or observation.snr_reference_exposure_time <= 0
            or observation.snr_reference_resolution <= 0
            or observation.continuum_snr_per_resolution_element <= 0
            or observation.systematic_floor < 0):
        raise ValueError("exposure, resolution and reference S/N must be positive; floor nonnegative")
    if len(t) < 2 or len(v) < 2 or np.any(np.diff(t) <= 0) or np.any(np.diff(v) <= 0):
        raise ValueError("time and velocity require strictly increasing grids")
    data = np.asarray(result.delta_absorption)
    if data.shape != (len(t), len(v)) or not np.all(np.isfinite(data)):
        raise ValueError("residual must be finite and match the grids")
    dv = SPEED_OF_LIGHT_KM_S / resolution
    nt = int(np.floor((t[-1]-t[0])/exp + 1e-10))
    nv = int(np.floor((v[-1]-v[0])/dv + 1e-10))
    if nt < 1 or nv < 1:
        raise ValueError("window must contain at least one full exposure and resolution element")
    # Anchor bins to the window start; discard incomplete final bins.
    te = t[0] + np.arange(nt+1)*exp
    ve = v[0] + np.arange(nv+1)*dv
    time_weights = np.zeros((nt, len(t)))
    for i in range(nt):
        nodes = np.r_[te[i], t[(t > te[i]) & (t < te[i+1])], te[i+1]]
        segments = np.diff(nodes)
        weights = np.r_[segments, 0] / 2 + np.r_[0, segments] / 2
        left = np.clip(np.searchsorted(t, nodes, side="right")-1, 0, len(t)-2)
        frac = (nodes-t[left])/(t[left+1]-t[left])
        np.add.at(time_weights[i], left, weights*(1-frac)/exp)
        np.add.at(time_weights[i], left+1, weights*frac/exp)
    sigma = dv / np.sqrt(8*np.log(2))
    cdf = np.vectorize(erf, otypes=[float])((ve[:, None]-v[None, :])/(np.sqrt(2)*sigma))/2
    quadrature = np.r_[np.diff(v), 0]/2 + np.r_[0, np.diff(v)]/2
    spectral_weights = np.diff(cdf, axis=0)*quadrature[None, :]/dv
    signal = amplitude * (time_weights @ data @ spectral_weights.T)
    reference_dv = SPEED_OF_LIGHT_KM_S / observation.snr_reference_resolution
    photon_noise = (1/observation.continuum_snr_per_resolution_element
                    * np.sqrt(observation.snr_reference_exposure_time/exp * reference_dv/dv))
    noise = float(np.hypot(photon_noise, observation.systematic_floor))
    return ObservedResidual((te[:-1]+te[1:])/2, (ve[:-1]+ve[1:])/2,
                            signal, noise, exp, dv)
