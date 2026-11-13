"""LTE CO slant-column templates from ExoMol with SciPy Voigt profiles.

This is a specified isothermal, single-species atmosphere, not a retrieval or
an equilibrium-chemistry/radiative-convective calculation.
"""
import bz2
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.special import voigt_profile
from scipy.interpolate import CubicSpline
from .observability import SPEED_OF_LIGHT_KM_S as C


@dataclass(frozen=True)
class COLines:
    wavenumber: np.ndarray
    strength: np.ndarray
    gamma_ref: np.ndarray
    temperature_exponent: np.ndarray
    partition: float
    retained_strength_fraction: float


def load_co_lines(folder: Path,temperature: float,nu_min=4245.,nu_max=4405.):
    states=np.loadtxt(bz2.open(folder/'12C-16O__Li2015.states.bz2','rt'),usecols=(0,1,2,3))
    transitions=np.loadtxt(bz2.open(folder/'12C-16O__Li2015.trans.bz2','rt'),usecols=(0,1,2))
    pf=np.loadtxt(folder/'12C-16O__Li2015.pf')
    if not pf[0,0]<=temperature<=pf[-1,0]:raise ValueError('temperature outside partition table')
    q=np.interp(temperature,pf[:,0],pf[:,1])
    lookup=np.zeros((int(states[:,0].max())+1,3));lookup[states[:,0].astype(int)]=states[:,1:]
    upper=transitions[:,0].astype(int);lower=transitions[:,1].astype(int)
    nu=lookup[upper,0]-lookup[lower,0]
    mask=(nu>=nu_min)&(nu<=nu_max)
    nu=nu[mask];lower=lower[mask];upper=upper[mask];a=transitions[mask,2]
    c2=1.438776877
    strength=lookup[upper,1]*a/(8*np.pi*2.99792458e10*nu**2)*np.exp(-c2*lookup[lower,0]/temperature)*(-np.expm1(-c2*nu/temperature))/q
    broad=np.loadtxt(folder/'12C-16O__H2.broad',usecols=(1,2,3))
    # a0 recipe depends on lower J; nearest endpoint for J above supplied table.
    j=lookup[lower,2]
    gamma=np.interp(j,broad[:,2],broad[:,0]);exponent=np.interp(j,broad[:,2],broad[:,1])
    keep=strength>=strength.max()*1e-8
    return COLines(nu[keep],strength[keep],gamma[keep],exponent[keep],float(q),float(strength[keep].sum()/strength.sum()))


def co_absorption(nu,lines,temperature,pressure_bar,column_cm2,microturbulence_km_s=1.,wing_sigma=100.):
    """1-exp(-N*sigma_nu), H2 pressure broadening, 28 amu CO Doppler mass."""
    speed=np.hypot(np.sqrt(1.380649e-23*temperature/(28*1.66053906660e-27))/1000,microturbulence_km_s)
    sigma=lines.wavenumber*speed/C
    gamma=lines.gamma_ref*.85*pressure_bar*(296/temperature)**lines.temperature_exponent
    opacity=np.zeros_like(nu)
    for center,strength,s,g in zip(lines.wavenumber,lines.strength,sigma,gamma):
        cutoff=max(wing_sigma*s,100*g)
        lo,hi=np.searchsorted(nu,[center-cutoff,center+cutoff])
        opacity[lo:hi]+=strength*voigt_profile(nu[lo:hi]-center,s,g)
    return -np.expm1(-column_cm2*opacity)


def bin_average(grid,values,edges):
    """Exact piecewise-linear integration along the last dimension."""
    if (len(grid)<2 or np.any(np.diff(grid)<=0) or np.any(np.diff(edges)<=0)
            or edges[0]<grid[0]-1e-10 or edges[-1]>grid[-1]+1e-10):
        raise ValueError('bin edges must increase and lie inside the integration grid')
    dx=np.diff(grid)
    cumulative=np.concatenate([np.zeros(values.shape[:-1]+(1,)),np.cumsum((values[...,1:]+values[...,:-1])*dx/2,axis=-1)],axis=-1)
    i=np.clip(np.searchsorted(grid,edges,side='right')-1,0,len(grid)-2)
    offset=edges-grid[i]
    integral=cumulative[...,i]+values[...,i]*offset+(values[...,i+1]-values[...,i])/dx[i]*offset**2/2
    return np.diff(integral,axis=-1)/np.diff(edges)


def cubic_shifted_bin_average(grid,profile,edges,shifts):
    """Integrate smooth shifted profiles directly over detector bins.

    The antiderivative avoids resampling a full fine spectrum for every local
    velocity and makes small wind derivatives independent of linear-grid kinks.
    """
    shifted_edges=np.asarray(edges)[None,:]-np.asarray(shifts)[:,None]
    if (np.any(np.diff(edges)<=0) or shifted_edges.min()<grid[0]
            or shifted_edges.max()>grid[-1]):
        raise ValueError('shifted detector bins must lie inside padded spectrum')
    integral=CubicSpline(grid,profile,extrapolate=False).antiderivative()
    return np.diff(integral(shifted_edges),axis=1)/np.diff(edges)


def depth_orthogonal_component(signal,reference):
    """Remove one freely fitted depth coefficient per exposure (equal noise)."""
    norm=np.sum(reference**2,axis=1,keepdims=True)
    coefficient=np.divide(np.sum(signal*reference,axis=1,keepdims=True),norm,out=np.zeros_like(norm),where=norm>0)
    return signal-coefficient*reference


def nuisance_basis(templates,rtol=1e-8):
    """Per-exposure orthonormal spectral nuisance rows (equal noise metric).

    Input shape is (n_modes, n_exposures, n_bins). Unit-normalizing each
    nonzero mode makes rank decisions independent of parameter units.
    Zero singular modes are retained as zero rows for a stable array shape.
    """
    a=np.asarray(templates,dtype=float)
    if a.ndim!=3 or not np.all(np.isfinite(a)) or not 0<rtol<1:
        raise ValueError('finite (mode,time,bin) templates and 0<rtol<1 required')
    a=np.swapaxes(a,0,1)
    norms=np.linalg.norm(a,axis=-1,keepdims=True)
    a=np.divide(a,norms,out=np.zeros_like(a),where=norms>0)
    _,singular,vh=np.linalg.svd(a,full_matrices=False)
    active=singular>rtol*singular[:,:1]
    return vh*active[...,None],np.sum(active,axis=1)


def project_nuisance(signal,basis):
    """Remove the fitted per-exposure nuisance subspace from a spectral map."""
    signal=np.asarray(signal)
    if signal.shape!=(basis.shape[0],basis.shape[2]):
        raise ValueError('signal and per-exposure basis shapes do not agree')
    coefficients=np.einsum('tv,tmv->tm',signal,basis)
    return signal-np.einsum('tm,tmv->tv',coefficients,basis)
