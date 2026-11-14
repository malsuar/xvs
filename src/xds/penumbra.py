"""Uniform-disk IRRADIANCE shadow diagnostic, not direct transmission transfer."""
import numpy as np
from .geometry import sky_basis
from .cell_blocking import quadrature,blocked_cell_fractions


def circle_overlap(r, q, distance):
    d=np.asarray(distance,dtype=float)
    result=np.zeros_like(d)
    inside=d<=abs(r-q)
    result[inside]=np.pi*min(r,q)**2
    cross=(d>abs(r-q))&(d<r+q)
    x=d[cross]
    if x.size:
        a=np.clip((x*x+r*r-q*q)/(2*x*r),-1,1)
        b=np.clip((x*x+q*q-r*r)/(2*x*q),-1,1)
        result[cross]=r*r*np.arccos(a)+q*q*np.arccos(b)-.5*np.sqrt(np.maximum(0,
            (-x+r+q)*(x+r-q)*(x-r+q)*(x+r+q)))
    return result


def irradiance_blocked_fractions(na,nr,thickness,moon_position,moon_radius,
                                 stellar_radius_rp,stellar_distance_rp,
                                 observer_direction=(-1.,0.,0.),order=8):
    """Average blocked incident irradiance in polar cells (uniform source).

    Paraxial stellar plane and thin circular occulter; integrates directions,
    NOT specific intensity along the observer ray. No atmospheric response.
    """
    if stellar_radius_rp<0 or stellar_distance_rp<=0 or order<2:
        raise ValueError('invalid source size, distance or quadrature')
    u,v,n=sky_basis(observer_direction)
    length=-np.dot(moon_position,n)
    if length<=0:return np.zeros(na*nr)
    if length>=stellar_distance_rp:raise ValueError('moon must lie between source plane and atmosphere')
    g=1-length/stellar_distance_rp
    center=np.array([moon_position@u,moon_position@v])/g
    radius=moon_radius/g
    blur=length*stellar_radius_rp/stellar_distance_rp/g
    if blur==0:
        projected=center[0]*u+center[1]*v
        return blocked_cell_fractions(na,nr,thickness,projected,radius,observer_direction,order=max(16,order))
    d=np.linalg.norm(center); reach=radius+blur
    result=np.zeros((nr,na))
    if d-reach>=1+thickness or d+reach<=1:return result.ravel()
    width=2*np.pi/na
    phi=(np.arange(na)+.5)*width
    delta=(phi-np.arctan2(center[1],center[0])+np.pi)%(2*np.pi)-np.pi
    angle=np.pi if d<=reach else np.arcsin(reach/d)
    selected=np.flatnonzero(abs(delta)<=angle+width/2)
    nodes,w=quadrature(order)
    angles=phi[selected,None]+nodes[None,:]*width/2
    edges=np.linspace(1,(1+thickness)**2,nr+1)
    for j,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
        if np.sqrt(hi)<d-reach or np.sqrt(lo)>d+reach:continue
        r=np.sqrt((hi+lo)/2+(hi-lo)*nodes/2)
        x=r[:,None,None]*np.cos(angles)[None,:,:]-center[0]
        y=r[:,None,None]*np.sin(angles)[None,:,:]-center[1]
        probability=circle_overlap(blur,radius,np.hypot(x,y))/(np.pi*blur**2)
        result[j,selected]=np.einsum('i,ijk,k->j',w,probability,w)/4
    return np.clip(result.ravel(),0,1)
