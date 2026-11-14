"""Empirically calibrated moving-shadow score search, not a nonlinear global LRT."""
import numpy as np
from scipy.optimize import minimize,least_squares
from scipy.stats import beta
from .nonlinear_atmosphere import fit_coherent
from .moving_shadow import period_limits


def bounded_improvement(correlation,norm_squared):
    """Exact scalar least-squares gain for a fixed projected template, 0<=f<=1."""
    c=np.asarray(correlation);d=np.asarray(norm_squared)
    f=np.clip(np.divide(c,d,out=np.zeros_like(c,dtype=float),where=d>1e-20),0,1)
    return np.maximum(2*f*c-f*f*d,0),f


def empirical_tail(calibration,statistic):
    values=np.asarray(calibration)
    if values.size==0 or np.any(~np.isfinite(values)):raise ValueError('finite calibration required')
    return float((1+np.sum(values>=statistic))/(len(values)+1))


def binomial_interval(k,n,confidence=.95):
    if n<1 or not 0<=k<=n:raise ValueError('invalid binomial counts')
    alpha=(1-confidence)/2
    return (0. if k==0 else float(beta.ppf(alpha,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha,k+1,n-k)))


def fit_null(model,data,noise):
    starts=[np.zeros(model.n_parameters),np.resize(np.array([.2,-.2,.1]),model.n_parameters)]
    best,runs=fit_coherent(model,data,noise,starts)
    if not all(r['success'] for r in runs):
        best,runs=fit_coherent(model,data,noise,starts,max_nfev=300)
    if not all(r['success'] for r in runs):raise RuntimeError('null optimizer failed; do not discard this trial')
    prediction,jac,_=model.evaluate(best['coefficients'],0.)
    design=jac.reshape(-1,model.n_parameters)
    norms=np.linalg.norm(design,axis=0)
    normalized=np.divide(design,norms,out=np.zeros_like(design),where=norms>0)
    u,s,_=np.linalg.svd(normalized,full_matrices=False);u=u[:,s>s[0]*1e-8]
    residual=(data-prediction)/noise
    return dict(best=best,residual=residual,u=u,parameters=model.trajectory(best['coefficients']),
        start_spread=max(r['chi2'] for r in runs)-min(r['chi2'] for r in runs))


class OrbitBank:
    def __init__(self,shadow,time,n_periods=73,n_phases=1024,order=16):
        self.shadow=shadow;self.time=time;self.order=order
        limits=period_limits();self.period_axis=np.geomspace(*limits,n_periods)
        self.phase_axis=np.arange(n_phases)/n_phases
        p,ph=np.meshgrid(self.period_axis,self.phase_axis,indexing='ij')
        self.periods=p.ravel();self.phases=ph.ravel();self.n_templates=len(self.periods)
        self.dlog=np.log(self.period_axis[1]/self.period_axis[0]);self.dphase=1/n_phases
        all_rows=[]
        for start in range(0,self.n_templates,256):
            oi,ti,g,w=shadow.orbit_weights(self.periods[start:start+256],self.phases[start:start+256],time,order)
            all_rows.append((oi+start,ti,g,w.astype(np.float32)))
        self.oi=np.concatenate([r[0] for r in all_rows]);self.ti=np.concatenate([r[1] for r in all_rows])
        self.g=np.concatenate([r[2] for r in all_rows]);self.weights=np.concatenate([r[3] for r in all_rows])
        self.groups=[(g,t,np.flatnonzero((self.g==g)&(self.ti==t))) for g in range(2) for t in range(len(time))]

    def conditional(self,null,noise):
        states,_=self.shadow.atmospheric_states(null['parameters'])
        u=null['u'];rank=u.shape[1];v=self.shadow.rotation
        ur=np.einsum('tvr,vs->tsr',u.reshape(len(self.time),v.shape[0],rank),v)
        r=null['residual'];reduced=r@v;utr=u.T@r.ravel()
        return dict(states=states,ur=ur,residual=reduced,utr=utr,noise=noise)

    def score_grid(self,condition):
        states=condition['states'];rank=condition['ur'].shape[-1];noise=condition['noise']
        corr=np.zeros(self.n_templates);norm=np.zeros_like(corr);projections=np.zeros((self.n_templates,rank))
        for g,t,idx in self.groups:
            if not len(idx):continue
            values=self.weights[idx]@states[g,t]/noise
            oi=self.oi[idx]
            corr+=np.bincount(oi,weights=values@condition['residual'][t],minlength=self.n_templates)
            norm+=np.bincount(oi,weights=np.sum(values*values,axis=1),minlength=self.n_templates)
            np.add.at(projections,oi,values@condition['ur'][t])
        corr-=projections@condition['utr'];norm-=np.sum(projections*projections,axis=1)
        score,f=bounded_improvement(corr,norm)
        return score,f

    def score_point(self,condition,period,phase):
        _,ti,g,w=self.shadow.orbit_weights([period],[phase%1],self.time,self.order)
        values=np.zeros((len(self.time),self.shadow.library.shape[-1]))
        np.add.at(values,ti,np.einsum('nk,nkr->nr',w,condition['states'][g,ti])/condition['noise'])
        projection=np.einsum('ts,tsr->r',values,condition['ur'])
        corr=np.sum(values*condition['residual'])-projection@condition['utr']
        norm=np.sum(values*values)-projection@projection
        score,f=bounded_improvement(corr,norm)
        return float(score),float(f)

    def search(self,null,noise,n_seeds=3,keep_grid=False):
        condition=self.conditional(null,noise);score,f=self.score_grid(condition)
        selected=[]
        for index in np.argsort(score)[::-1]:
            p=self.periods[index];ph=self.phases[index]
            if any(abs(np.log(p/self.periods[j]))<2*self.dlog and
                   abs((ph-self.phases[j]+.5)%1-.5)<2*self.dphase for j in selected):continue
            selected.append(index)
            if len(selected)==n_seeds:break
        candidates=[]
        lo,hi=np.log(period_limits())
        for index in selected:
            lp=np.log(self.periods[index]);ph=self.phases[index]
            bounds=[(max(-1.5,(lo-lp)/self.dlog),min(1.5,(hi-lp)/self.dlog)),(-1.5,1.5)]
            def objective(z):return -self.score_point(condition,np.exp(lp+z[0]*self.dlog),ph+z[1]*self.dphase)[0]
            fit=minimize(objective,[0.,0.],method='L-BFGS-B',bounds=bounds,
                options=dict(maxiter=35,maxfun=150,ftol=1e-10,eps=1e-5))
            period=float(np.exp(lp+fit.x[0]*self.dlog));phase=float((ph+fit.x[1]*self.dphase)%1)
            value,fraction=self.score_point(condition,period,phase)
            if value<score[index]:period=float(self.periods[index]);phase=float(self.phases[index]);value=float(score[index]);fraction=float(f[index])
            candidates.append(dict(statistic=value,period=period,phase=phase,blocking=fraction,refinement_success=bool(fit.success)))
        best=max(candidates,key=lambda c:c['statistic'])
        result=dict(best,grid_max=float(score.max()),candidates=candidates)
        if keep_grid:result['grid']=score.reshape(len(self.period_axis),len(self.phase_axis))
        return result,condition


def nonlinear_candidate_fit(model,shadow,data,noise,null,candidate):
    """Local joint nonlinear refinement; NOT the statistic used for false alarms."""
    p0=candidate['period'];phase0=candidate['phase'];lp_bounds=np.log(period_limits())
    starts=[np.r_[null['best']['coefficients'],candidate['blocking'],np.log(p0),0.],
            np.r_[np.zeros(model.n_parameters),.5,np.log(p0),0.]]
    lower=np.r_[np.full(model.n_parameters,-1.),0.,lp_bounds[0],-.03]
    upper=np.r_[np.ones(model.n_parameters),1.,lp_bounds[1],.03]
    runs=[]
    for start in starts:
        cache={}
        def values(x):
            if 'x' not in cache or not np.array_equal(cache['x'],x):
                coeff=x[:model.n_parameters];f=x[-3];period=np.exp(x[-2]);phase=phase0+x[-1]
                reference,jb,_=model.evaluate(coeff,0.)
                p=model.trajectory(coeff);s,ds=shadow.signal(p,period,phase,model.time,derivatives=True,order=16)
                jac=jb.copy();column=0
                for g in range(2):
                    for i in range(3):
                        for k in range(model.n_controls):
                            jac[...,column]+=f*ds[g,i]*model.temporal[:,k,None];column+=1
                delta=1e-5
                dp=(shadow.signal(p,period*np.exp(delta),phase,model.time,order=16)-shadow.signal(p,period*np.exp(-delta),phase,model.time,order=16))/(2*delta)
                dph=(shadow.signal(p,period,phase+delta,model.time,order=16)-shadow.signal(p,period,phase-delta,model.time,order=16))/(2*delta)
                jac=np.concatenate([jac,s[...,None],f*dp[...,None],f*dph[...,None]],axis=-1)
                cache.update(x=x.copy(),r=((reference+f*s-data)/noise).ravel(),j=jac.reshape(-1,len(x))/noise)
            return cache
        fit=least_squares(lambda x:values(x)['r'],start,jac=lambda x:values(x)['j'],bounds=(lower,upper),
            x_scale=np.r_[np.ones(model.n_parameters+1),.05,.002],max_nfev=150,ftol=1e-9,xtol=1e-9,gtol=1e-7)
        runs.append(dict(chi2=float(2*fit.cost),period=float(np.exp(fit.x[-2])),phase=float((phase0+fit.x[-1])%1),
            blocking=float(fit.x[-3]),coefficients=fit.x[:model.n_parameters].tolist(),success=bool(fit.success),nfev=int(fit.nfev)))
    return min(runs,key=lambda r:r['chi2']),runs
