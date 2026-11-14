"""Aligned moving shadows factored into projected-cell and stellar moments.

Only the circular, equatorial, central-transit pilot is represented here.
Changing period changes orbital separation through Kepler's law.
"""
import numpy as np
from scipy.interpolate import CubicSpline,make_interp_spline
from .cell_blocking import blocked_cell_fractions
from .terminator import build_projected_annulus
from .geometry import sky_basis
from .rotation import line_of_sight_velocity
from .dynamics import moon_period_hours,planet_density_g_cm3,fluid_roche_limit_rp,stable_outer_moon_limit_rp


REFERENCE_PERIOD=float(moon_period_hours(3.,.7,1.1))
MOON_RADIUS=.09


def period_limits():
    density=planet_density_g_cm3(.7,1.1)
    axes=np.array([fluid_roche_limit_rp(density,3.),stable_outer_moon_limit_rp(30.,.7,1.,1.1)])
    return moon_period_hours(axes,.7,1.1)


def separation_from_period(period):
    if np.any(np.asarray(period)<=0):raise ValueError('positive orbital period required')
    return 3.*(np.asarray(period)/REFERENCE_PERIOD)**(2/3)


def stellar_moments(time,center,rp):
    """Quadratic Taylor coefficients for I at (u=center,z=0), normalized moments."""
    radius=695700/rp;axis=(30/365.25)**(2/3)*149597870.7/rp
    x=np.asarray(center)+axis*np.sin(2*np.pi*np.asarray(time)/720.)
    mu=np.sqrt(1-(x/radius)**2)
    if np.any(~np.isfinite(mu)) or np.any(mu<=0):raise ValueError('moment expansion requires fully backed rays')
    a0=1-.4-.25;a1=.4+2*.25;a2=-.25;mean=1-.4/3-.25/6
    common=(a1/mu+2*a2)/mean
    intensity=(a0+a1*mu+a2*mu*mu)/mean
    first=-common*x/radius**2
    second=-common/radius**2-a1*x*x/(mean*radius**4*mu**3)
    zsecond=-common/radius**2
    return np.stack([intensity,first*MOON_RADIUS,.5*second*MOON_RADIUS**2,.5*zsecond*MOON_RADIUS**2],axis=-1)


def moment_geometry(base,n_positions=193,n_modes=160):
    nr=len(base['radii']);na=480;h=np.sqrt(1+2*nr*(base['radii'][0]**2-1))-1
    grid=build_projected_annulus(na,nr,h);u,v,_=sky_basis((-1.,0.,0.))
    sky_u=grid.positions@u;sky_z=grid.positions@v
    velocity=line_of_sight_velocity(grid.positions,12.,90.,0.,(-1.,0.,0.),'observer')
    shifts=base['shifts'];step=shifts[1]-shifts[0]
    index=np.clip(np.floor((velocity-shifts[0])/step).astype(int),0,len(shifts)-2)
    fraction=(velocity-shifts[index])/step
    flat=np.repeat(np.arange(nr),na)*len(shifts)+index
    centers=np.linspace(1-MOON_RADIUS,1+h+MOON_RADIUS,n_positions)
    kernels=[]
    for sign in [-1.,1.]:
        for center in sign*centers:
            blocked=blocked_cell_fractions(na,nr,h,center*u,MOON_RADIUS,order=16)
            delta=(sky_u-center)/MOON_RADIUS
            for moment in [np.ones_like(delta),delta,delta**2,(sky_z/MOON_RADIUS)**2]:
                w=grid.weights*blocked*moment
                kernels.append(np.bincount(flat,weights=w*(1-fraction),minlength=nr*len(shifts))+
                    np.bincount(flat+1,weights=w*fraction,minlength=nr*len(shifts)))
    kernels=np.array(kernels)
    rng=np.random.default_rng(63130);rank=min(n_modes+32,min(kernels.shape))
    q=np.linalg.qr(kernels@rng.normal(size=(kernels.shape[1],rank)),mode='reduced')[0]
    q=np.linalg.qr(kernels@(kernels.T@q),mode='reduced')[0]
    left,s,right=np.linalg.svd(q.T@kernels,full_matrices=False)
    left=(q@left[:,:n_modes]).reshape(2,n_positions,4,n_modes)
    modes=(s[:n_modes,None]*right[:n_modes]).reshape(n_modes,nr,len(shifts))
    approximation=left.reshape(-1,n_modes)@modes.reshape(n_modes,-1)
    error=float(np.linalg.norm(kernels-approximation)/np.linalg.norm(kernels))
    geometry=dict(base,full=np.zeros((1,*modes.shape)),missing=modes[None],time=np.arange(n_modes,dtype=float))
    return centers,left,geometry,error


class MovingShadow:
    def __init__(self,axes,library,left,centers,spectral_rotation,rp):
        self.library=library;self.centers=centers;self.rotation=spectral_rotation;self.rp=rp
        self.axes=axes;self.splines=[make_interp_spline(a,np.eye(len(a)),k=3) for a in axes]
        self.position=[CubicSpline(centers,a,axis=0,extrapolate=False) for a in left]

    def atmospheric_states(self,p,derivatives=False):
        """p: (limb,time,3); output (limb,time,mode,spectral)."""
        p=np.asarray(p)
        if np.any(abs(p)>1+1e-12):raise ValueError('atmosphere outside bounded grid')
        values=[];jac=[]
        cube=self.library.reshape(np.prod([len(a) for a in self.axes]),-1)
        for g in range(2):
            w=[s(p[g,:,i]) for i,s in enumerate(self.splines)]
            def tensor(v):return (v[0][:,:,None,None]*v[1][:,None,:,None]*v[2][:,None,None,:]).reshape(len(p[g]),-1)
            v=[tensor(w)]
            if derivatives:
                for i,s in enumerate(self.splines):
                    d=w.copy();d[i]=s(p[g,:,i],nu=1);v.append(tensor(d))
            result=(np.stack(v)@cube).reshape(len(v),len(p[g]),*self.library.shape[-2:])
            values.append(result[0])
            if derivatives:jac.append(result[1:])
        return np.array(values),np.array(jac) if derivatives else None

    def orbit_weights(self,periods,phases,time,order=8):
        """Sparse active (orbit,time,limb) weights, exposure-averaged over 5 min."""
        periods=np.atleast_1d(periods);phases=np.atleast_1d(phases)
        if len(periods)!=len(phases):raise ValueError('period and phase arrays must match')
        nodes,weights=np.polynomial.legendre.leggauss(order)
        tq=np.asarray(time)[None,:,None]+nodes[None,None,:]*(5/120)
        center=-separation_from_period(periods)[:,None,None]*np.cos(2*np.pi*(tq/periods[:,None,None]+phases[:,None,None]))
        active=(abs(center)>self.centers[0])&(abs(center)<self.centers[-1])
        oi,ti,qi=np.where(active);uc=center[oi,ti,qi]
        branch=(uc>=0).astype(int);limb=1-branch
        mode_count=self.library.shape[-2]
        totals=np.zeros((len(periods)*len(time)*2,mode_count))
        for b in [0,1]:
            selected=branch==b
            if not np.any(selected):continue
            coeff=self.position[b](abs(uc[selected]))
            star=stellar_moments(tq[0,ti[selected],qi[selected]],uc[selected],self.rp)
            combined=np.einsum('nmk,nm->nk',coeff,star)*(weights[qi[selected]]/2)[:,None]
            row=(oi[selected]*len(time)+ti[selected])*2+limb[selected]
            np.add.at(totals,row,combined)
        rows=np.flatnonzero(np.linalg.norm(totals,axis=1)>1e-18)
        return rows//(2*len(time)),(rows//2)%len(time),rows%2,totals[rows]

    def signal(self,p,period,phase,time,derivatives=False,order=8):
        states,jac=self.atmospheric_states(p,derivatives)
        _,ti,g,w=self.orbit_weights([period],[phase],time,order)
        reduced=np.zeros((len(time),self.library.shape[-1]))
        np.add.at(reduced,ti,np.einsum('nk,nkr->nr',w,states[g,ti]))
        if not derivatives:return reduced@self.rotation.T
        derivative=np.zeros((2,3,len(time),self.library.shape[-1]))
        for limb in [0,1]:
            selected=g==limb
            for i in range(3):
                np.add.at(derivative[limb,i],ti[selected],np.einsum('nk,nkr->nr',w[selected],jac[limb,i,ti[selected]]))
        return reduced@self.rotation.T,derivative@self.rotation.T
