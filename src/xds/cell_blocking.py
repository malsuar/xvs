"""Projected disk overlap with polar annular cells."""
from functools import lru_cache
import numpy as np
from .geometry import sky_basis


@lru_cache(maxsize=16)
def quadrature(order):
    return np.polynomial.legendre.leggauss(order)


def blocked_cell_fractions(n_azimuths, n_radii, thickness, moon_position,
                           moon_radius, observer_direction=(-1.,0.,0.),
                           order=16, angular_phase=0.):
    """Area fractions, radial-major ordering; exact angular interval overlap.

    Integrate in r squared by Gaussian quadrature, splitting at radial disk
    tangencies. Local spectra still use the representative cell velocity.
    """
    if order < 2 or moon_radius <= 0 or thickness <= 0:
        raise ValueError("positive radius/thickness and quadrature order >=2 required")
    u,v,_=sky_basis(observer_direction)
    center=np.array([np.dot(moon_position,u),np.dot(moon_position,v)])
    d=np.linalg.norm(center)
    edges=np.linspace(1.,(1+thickness)**2,n_radii+1)
    result=np.zeros((n_radii,n_azimuths))
    if d-moon_radius >= 1+thickness or d+moon_radius <= 1:
        return result.ravel()
    if d < 1e-14:
        covered=np.clip(moon_radius**2-edges[:-1],0,np.diff(edges)) / np.diff(edges)
        return np.repeat(covered,n_azimuths)
    theta=np.arctan2(center[1],center[0])
    width=2*np.pi/n_azimuths
    centers=(np.arange(n_azimuths)+.5+angular_phase)*width
    delta=(centers-theta+np.pi)%(2*np.pi)-np.pi
    max_angle=np.pi if d<=moon_radius else np.arcsin(moon_radius/d)
    selected=np.flatnonzero(abs(delta)<=max_angle+width/2)
    nodes,weights=quadrature(order)
    for j,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
        cuts=sorted(set([lo,hi]+[x*x for x in [abs(d-moon_radius),d+moon_radius] if lo<x*x<hi]))
        for a,b in zip(cuts[:-1],cuts[1:]):
            s=(a+b)/2+(b-a)/2*nodes
            alpha=np.arccos(np.clip((s+d*d-moon_radius**2)/(2*np.sqrt(s)*d),-1,1))
            length=np.zeros((order,len(selected)))
            for wrap in [-2*np.pi,0.,2*np.pi]:
                c=delta[selected]+wrap
                length+=np.maximum(0,np.minimum(c[None,:]+width/2,alpha[:,None])-
                                   np.maximum(c[None,:]-width/2,-alpha[:,None]))
            result[j,selected]+=(weights@length)*(b-a)/(2*(hi-lo)*width)
    return np.clip(result.ravel(),0,1)
