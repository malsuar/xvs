"""Bounded, temporally coherent two-limb CO inference on a validated spectral grid.

The lookup is built by the molecular production pilot, not by fitting spectral
derivatives. Compression uses a common orthonormal wavelength basis.
"""
import numpy as np
from scipy.interpolate import BSpline,make_interp_spline
from scipy.optimize import least_squares


def physical_parameters(normalized):
    """Coordinates [...,3] -> T, CO abundance, signed effective LOS wind."""
    p=np.asarray(normalized)
    return 640+80*p[...,0],1e-4*np.exp(np.log(3)*p[...,1]),1.5*p[...,2]


def temporal_basis(time,n_controls):
    """Nonnegative clamped B-splines over [-2,2] h preserve control-point bounds."""
    if n_controls<1:raise ValueError('at least one temporal control required')
    if n_controls==1:return np.ones((len(time),1))
    degree=min(3,n_controls-1)
    interior=np.linspace(-2,2,n_controls-degree+1)[1:-1]
    knots=np.r_[np.repeat(-2.,degree+1),interior,np.repeat(2.,degree+1)]
    if np.any(np.asarray(time)<-2) or np.any(np.asarray(time)>2):raise ValueError('time outside transit window')
    return BSpline(knots,np.eye(n_controls),degree)(time)


class CoherentCOModel:
    """Nonlinear bounded CO interpolation with analytic parameter Jacobian."""

    def __init__(self,axes,library,time,n_controls):
        self.axes=axes;self.library=library;self.time=np.asarray(time)
        self.temporal=temporal_basis(time,n_controls);self.n_controls=n_controls
        self.splines=[make_interp_spline(a,np.eye(len(a)),k=3) for a in axes]
        self.cubes=[np.ascontiguousarray(library[:,:,:,g].transpose(3,0,1,2,4,5)
            .reshape(len(time),-1,2*library.shape[-1])) for g in range(2)]
        self.n_parameters=6*n_controls
        zero=np.zeros(self.n_parameters)
        self.reference=self.evaluate(zero,0.,subtract_reference=False)[0]

    def trajectory(self,coefficients):
        return np.einsum('gck,tk->gtc',np.asarray(coefficients).reshape(2,3,self.n_controls),self.temporal)

    def evaluate(self,coefficients,blocking,subtract_reference=True):
        if not 0<=blocking<=1:raise ValueError('lunar blocking must lie in [0,1]')
        p=self.trajectory(coefficients)
        if np.any(abs(p)>1+1e-12):raise ValueError('atmosphere outside lookup bounds')
        fields=[];derivatives=[]
        for g in range(2):
            weights=[s(p[g,:,i]) for i,s in enumerate(self.splines)]
            def tensor(w):return (w[0][:,:,None,None]*w[1][:,None,:,None]*w[2][:,None,None,:]).reshape(len(self.time),-1)
            interpolation=[tensor(weights)]
            for i,spline in enumerate(self.splines):
                w=weights.copy();w[i]=spline(p[g,:,i],nu=1)
                interpolation.append(tensor(w))
            evaluated=(np.stack(interpolation,axis=1)@self.cubes[g]).reshape(len(self.time),4,2,-1)
            fields.append(evaluated[:,0]);derivatives.append([evaluated[:,i] for i in range(1,4)])
        fields=np.array(fields)
        result=np.sum(fields[:,:,0]+blocking*fields[:,:,1],axis=0)
        if subtract_reference:result-=self.reference
        jac=[]
        for g in range(2):
            for i in range(3):
                derivative=derivatives[g][i][:,0]+blocking*derivatives[g][i][:,1]
                jac.extend(derivative*self.temporal[:,k,None] for k in range(self.n_controls))
        return result,np.stack(jac,axis=-1),np.sum(fields[:,:,1],axis=0)


def fit_coherent(model,data,noise,starts,with_moon=False,max_nfev=120,fixed_blocking=0.):
    """Bounded least squares, returning every start rather than hiding failures."""
    runs=[]
    for initial in starts:
        x0=np.r_[initial,.5] if with_moon else np.asarray(initial)
        lower=np.r_[np.full(model.n_parameters,-1.),0.] if with_moon else -1.
        upper=np.r_[np.ones(model.n_parameters),1.] if with_moon else 1.
        cached={}
        def values(x):
            if 'x' not in cached or not np.array_equal(cached['x'],x):
                a=x[-1] if with_moon else fixed_blocking
                pred,jac,ds=model.evaluate(x[:model.n_parameters],a)
                if with_moon:jac=np.concatenate([jac,ds[...,None]],axis=-1)
                cached.update(x=x.copy(),residual=((pred-data)/noise).ravel(),
                    jac=jac.reshape(-1,len(x))/noise)
            return cached
        fit=least_squares(lambda x:values(x)['residual'],x0,jac=lambda x:values(x)['jac'],
            bounds=(lower,upper),method='trf',ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=max_nfev)
        runs.append(dict(coefficients=fit.x[:model.n_parameters],blocking=float(fit.x[-1]) if with_moon else float(fixed_blocking),
            chi2=float(2*fit.cost),success=bool(fit.success),status=int(fit.status),nfev=int(fit.nfev),
            optimality=float(fit.optimality),active_bounds=int(np.sum(abs(fit.active_mask)>0))))
    best=min(runs,key=lambda r:r['chi2'])
    return best,runs
