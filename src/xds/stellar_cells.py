"""Joint stellar-disk/lunar-disk integration in polar annular cells."""
import numpy as np
from .geometry import sky_basis
from .cell_blocking import quadrature
from .stellar_surface import background_intensity


def stellar_cell_weights(na,nr,h,time,surface,moon_position,moon_radius,
                         observer_direction=(-1.,0.,0.),order=16,intensity_order=8):
    """Return cell-mean normalized intensity and jointly blocked intensity.

    Angular circle intersections are clipped exactly at each radial node.
    Gaussian quadrature integrates r squared and the smooth limb law.
    Velocity remains represented by one value per cell.
    """
    background_intensity(np.zeros((1,3)),time,surface,observer_direction)
    u,v,_=sky_basis(observer_direction)
    phase=2*np.pi*time/surface.orbital_period_hours
    star=np.array([-surface.orbital_separation_rp*np.sin(phase),
                   -surface.impact_parameter*surface.radius_rp*np.cos(phase)])
    moon=np.array([moon_position@u,moon_position@v])
    mean=1-surface.u1/3-surface.u2/6
    edges=np.linspace(1,(1+h)**2,nr+1);width=2*np.pi/na
    cell_lo=np.arange(na)*width
    nodes,w=quadrature(order)
    an,aw=quadrature(intensity_order)
    visible=np.zeros((nr,na));blocked=np.zeros_like(visible)
    if np.linalg.norm(star)>surface.radius_rp+1+h:return visible.ravel(),blocked.ravel()

    def intervals(center,radius,r):
        d=np.linalg.norm(center)
        if d<1e-14:
            alpha=np.where(r<=radius,np.pi,0.)
            theta=np.pi
        else:
            alpha=np.arccos(np.clip((r*r+d*d-radius*radius)/(2*r*d),-1,1))
            theta=np.arctan2(center[1],center[0])%(2*np.pi)
        return [(theta-alpha[:,None]+wrap,theta+alpha[:,None]+wrap)
                for wrap in [-2*np.pi,0.,2*np.pi]]

    def intensity_integral(lo,hi,r):
        length=np.maximum(0,hi-lo)
        if surface.u1==0 and surface.u2==0:return length
        result=np.zeros_like(length)
        ii,jj=np.nonzero(length>0)
        if len(ii)==0:return result
        phi=(lo[ii,jj]+hi[ii,jj])[:,None]/2+length[ii,jj,None]*an/2
        x=r[ii,None]*np.cos(phi)-star[0]
        y=r[ii,None]*np.sin(phi)-star[1]
        mu=np.sqrt(np.maximum(0,1-(x*x+y*y)/surface.radius_rp**2))
        intensity=(1-surface.u1*(1-mu)-surface.u2*(1-mu)**2)/mean
        result[ii,jj]=length[ii,jj]*(intensity@aw)/2
        return result

    for j,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        cuts=[a,b]
        for center,radius in [(star,surface.radius_rp),(moon,moon_radius)]:
            for boundary in [abs(np.linalg.norm(center)-radius),np.linalg.norm(center)+radius]:
                if a<boundary**2<b:cuts.append(boundary**2)
        cuts=sorted(set(cuts))
        for lo_s,hi_s in zip(cuts[:-1],cuts[1:]):
            r=np.sqrt((lo_s+hi_s)/2+(hi_s-lo_s)*nodes/2)
            factor=(hi_s-lo_s)/(2*(b-a)*width)
            moon_intervals=intervals(moon,moon_radius,r)
            for slo,shi in intervals(star,surface.radius_rp,r):
                lo=np.maximum(cell_lo[None,:],slo)
                hi=np.minimum(cell_lo[None,:]+width,shi)
                visible[j]+=factor*(w@intensity_integral(lo,hi,r))
                for mlo,mhi in moon_intervals:
                    blocked[j]+=factor*(w@intensity_integral(np.maximum(lo,mlo),np.minimum(hi,mhi),r))
    return visible.ravel(),np.minimum(visible,blocked).ravel()
