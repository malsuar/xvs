"""Count-based K-band budgets and a shared-reference observation metric.

Template residual amplitudes are sensitivity assumptions, not ETC predictions.
The sufficient spectral subspace is fixed by the planetary library, never data.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import qr
from .orbit_search import OrbitBank,bounded_improvement


@dataclass(frozen=True)
class CountBudget:
    area_m2: float = 49.
    k_magnitude: float = 6.
    throughput: float = .08
    live_seconds: float = 240.
    resolving_power: float = 100000.
    extraction_arcsec2: float = .08
    pixels_per_bin: int = 8
    reads_per_exposure: int = 4
    read_noise_e: float = 6.
    dark_e_s_pixel: float = .003
    thermal_temperature_k: float = 280.
    thermal_emissivity: float = .1


def electron_budget(wavelength_um,stellar_relative,transmission,sky_radiance,budget):
    """Counts per resolution element. Sky is ph/s/m2/um/arcsec2.

    Ks zero point 666.7 Jy is a monochromatic normalization at 2.159 um;
    stellar_relative includes the specified continuum slope and line spectrum.
    A paired sky subtraction contributes twice the background variance.
    """
    b=budget
    if min(b.area_m2,b.throughput,b.live_seconds,b.resolving_power)<=0:
        raise ValueError('positive collecting area, efficiency, exposure and resolution required')
    w=np.asarray(wavelength_um);dl=w/b.resolving_power
    trans=np.asarray(transmission);star=np.asarray(stellar_relative)
    if np.any(trans<=0) or np.any(trans>1+1e-10) or np.any(star<=0):
        raise ValueError('positive stellar flux and physical transmission required')
    h=6.62607015e-34;c=299792458.;kb=1.380649e-23
    fnu=666.7e-26*10**(-.4*b.k_magnitude)*star
    source=fnu/(h*b.resolving_power)*trans*b.area_m2*b.throughput*b.live_seconds
    sky=np.asarray(sky_radiance)*dl*b.extraction_arcsec2*b.area_m2*b.throughput*b.live_seconds
    lam=w*1e-6
    thermal_photons=2*c/lam**4/np.expm1(h*c/(lam*kb*b.thermal_temperature_k))*1e-6
    thermal=thermal_photons*b.thermal_emissivity*(b.extraction_arcsec2/206265**2)*dl*b.area_m2*b.throughput*b.live_seconds
    dark=np.full_like(w,b.pixels_per_bin*b.dark_e_s_pixel*b.live_seconds)
    read=np.full_like(w,b.pixels_per_bin*b.reads_per_exposure*b.read_noise_e**2)
    variance=source+2*(sky+thermal+dark+read)
    return dict(source=source,sky=sky,thermal=thermal,dark=dark,read_variance=read,
        variance=variance,sigma=np.sqrt(variance)/source)


def time_whiten(values,n_reference):
    """(I+11^T/Nref)^(-1/2), along the exposure axis (axis 0)."""
    if n_reference<=0:raise ValueError('positive reference exposure count required')
    a=1-np.sqrt(n_reference/(n_reference+len(values)))
    return values-a*np.mean(values,axis=0,keepdims=True)


def spectral_metric(basis,sigma,keep,templates=None):
    """Whiten, optionally remove fixed nuisance spectra, then preserve all
    retained planetary-library directions by QR. Returns Q and coefficient map.
    """
    b=np.asarray(basis)[keep];s=np.asarray(sigma)[keep]
    if np.any(s<=0):raise ValueError('positive uncertainties required')
    weighted=b/s[:,None]
    if templates is not None:
        design=np.asarray(templates)[:,keep].T/s[:,None]
        u,sv,_=np.linalg.svd(design,full_matrices=False)
        u=u[:,sv>sv[0]*1e-9]
        weighted-=u@(u.T@weighted)
    q,r=qr(weighted,mode='economic')
    if np.linalg.matrix_rank(r)!=r.shape[0]:raise ValueError('planetary subspace rank lost')
    return q,r.T


class ObservedCOModel:
    """Exact linear observation operator applied to the bounded CO model."""
    def __init__(self,model,transform,n_reference=12):
        self.model=model;self.transform=transform;self.n_reference=n_reference
        self.n_parameters=model.n_parameters;self.time=model.time
        self.n_controls=model.n_controls

    def trajectory(self,coefficients):return self.model.trajectory(coefficients)

    def evaluate(self,coefficients,blocking):
        value,jac,moon=self.model.evaluate(coefficients,blocking)
        return (time_whiten(value@self.transform,self.n_reference),
            time_whiten(np.einsum('tsp,sr->trp',jac,self.transform),self.n_reference),
            time_whiten(moon@self.transform,self.n_reference))


class ReferenceOrbitBank(OrbitBank):
    """Same fixed orbital search, with exact shared-reference covariance.

    Reuses the validated geometrical weights. The shadow's spectral rotation
    must be orthonormal in the new whitened/projected spectral coordinates.
    """
    def __init__(self,base_bank,shadow,n_reference=12):
        self.__dict__.update(base_bank.__dict__)
        self.shadow=shadow;self.n_reference=n_reference

    def conditional(self,null,noise):
        if noise!=1.:raise ValueError('observation coordinates already whiten spectral noise')
        states,_=self.shadow.atmospheric_states(null['parameters'])
        u=null['u'];rank=u.shape[1];v=self.shadow.rotation
        native_u=time_whiten(u.reshape(len(self.time),v.shape[0],rank),self.n_reference)
        ur=np.einsum('tvr,vs->tsr',native_u,v)
        residual=time_whiten(null['residual'],self.n_reference)@v
        return dict(states=states,ur=ur,residual=residual,
            utr=u.T@null['residual'].ravel(),noise=1.)

    def score_grid(self,condition):
        states=condition['states'];rank=condition['ur'].shape[-1]
        corr=np.zeros(self.n_templates);norm=np.zeros_like(corr)
        projections=np.zeros((self.n_templates,rank))
        sums=np.zeros((self.n_templates,self.shadow.library.shape[-1]))
        for g,t,idx in self.groups:
            if not len(idx):continue
            values=self.weights[idx]@states[g,t];oi=self.oi[idx]
            corr+=np.bincount(oi,weights=values@condition['residual'][t],minlength=self.n_templates)
            norm+=np.bincount(oi,weights=np.sum(values*values,axis=1),minlength=self.n_templates)
            np.add.at(projections,oi,values@condition['ur'][t])
            np.add.at(sums,oi,values)
        norm-=np.sum(sums*sums,axis=1)/(self.n_reference+len(self.time))
        corr-=projections@condition['utr'];norm-=np.sum(projections*projections,axis=1)
        return bounded_improvement(corr,norm)

    def score_point(self,condition,period,phase):
        _,ti,g,w=self.shadow.orbit_weights([period],[phase%1],self.time,self.order)
        values=np.zeros((len(self.time),self.shadow.library.shape[-1]))
        np.add.at(values,ti,np.einsum('nk,nkr->nr',w,condition['states'][g,ti]))
        projection=np.einsum('ts,tsr->r',values,condition['ur'])
        corr=np.sum(values*condition['residual'])-projection@condition['utr']
        norm=np.sum(time_whiten(values,self.n_reference)**2)-projection@projection
        score,f=bounded_improvement(corr,norm)
        return float(score),float(f)
