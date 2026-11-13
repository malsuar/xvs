"""Continuum background specific intensity for observer-directed rays."""
import numpy as np
from .geometry import sky_basis


def background_intensity(positions,time,surface,observer_direction=(-1.,0.,0.)):
    """I(mu)/disk-mean I; zero outside the stellar disk. No renormalization."""
    if surface.radius_rp<=0 or surface.orbital_separation_rp<=surface.radius_rp or surface.orbital_period_hours<=0:
        raise ValueError('invalid stellar size or circular orbit')
    if abs(surface.impact_parameter*surface.radius_rp)>surface.orbital_separation_rp:
        raise ValueError('impact parameter incompatible with orbit')
    # Validate the intensity polynomial on its entire physical domain.
    z=[0.,1.]
    if surface.u2!=0:
        vertex=-surface.u1/(2*surface.u2)
        if 0<vertex<1:z.append(vertex)
    if min(1-surface.u1*x-surface.u2*x*x for x in z)<0:
        raise ValueError('negative stellar intensity')
    mean=1-surface.u1/3-surface.u2/6
    if mean<=0:raise ValueError('nonpositive disk mean intensity')
    u,v,_=sky_basis(observer_direction)
    phase=2*np.pi*time/surface.orbital_period_hours
    x=positions@u+surface.orbital_separation_rp*np.sin(phase)
    y=positions@v+surface.impact_parameter*surface.radius_rp*np.cos(phase)
    r2=(x*x+y*y)/surface.radius_rp**2
    mu=np.sqrt(np.maximum(0,1-r2))
    return np.where(r2<=1,(1-surface.u1*(1-mu)-surface.u2*(1-mu)**2)/mean,0.)
