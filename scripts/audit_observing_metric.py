"""Compare the observation-space sparse score to direct full-map algebra."""
import json
from pathlib import Path
import numpy as np
from blind_co_recovery import load_models
from observing_contamination_pilot import templates,setup_analysis
from xds.observing_contamination import CountBudget,electron_budget,time_whiten
from xds.orbit_search import fit_null,bounded_improvement
from xds.moving_shadow import REFERENCE_PERIOD


def main():
    root=Path(__file__).resolve().parents[1]
    old,shadow,models,bank=load_models(root,3,64)
    tpl=templates(root,old);keep=(tpl['transmission']>=.3)&(tpl['disk']>=.2)
    budget=CountBudget(area_m2=978.,live_seconds=84*1.427,reads_per_exposure=84,read_noise_e=12.)
    sigma=electron_budget(tpl['wavelength'],tpl['stellar_relative'],tpl['transmission'],tpl['radiance'],budget)['sigma']
    rng=np.random.default_rng(63152);results=[]
    for protected in [False,True]:
        a=setup_analysis(old,shadow,models[3],bank,sigma,keep,tpl,protected)
        coeff=rng.uniform(-.2,.2,models[3].n_parameters)
        prediction,jac,_=a['model'].evaluate(coeff,0.)
        maximum=0.
        for column in [0,4,8,12,17]:
            delta=np.zeros_like(coeff);delta[column]=1e-5
            finite=(a['model'].evaluate(coeff+delta,0.)[0]-a['model'].evaluate(coeff-delta,0.)[0])/2e-5
            maximum=max(maximum,float(np.linalg.norm(finite-jac[:,:,column])/max(np.linalg.norm(finite),1e-30)))
        raw=rng.normal(size=prediction.shape)+rng.normal(size=(1,prediction.shape[1]))/np.sqrt(12)
        null=fit_null(a['model'],prediction+time_whiten(raw,12),1.)
        condition=a['bank'].conditional(null,1.);grid,f=a['bank'].score_grid(condition)
        indices=np.unique(np.r_[np.linspace(0,len(grid)-1,8,dtype=int),np.argmax(grid)])
        errors=[];spectral_errors=[]
        for i in indices:
            p=bank.periods[i];ph=bank.phases[i]
            original=shadow.signal(null['parameters'],p,ph,old['time'],order=16)
            s=time_whiten(original@a['transform'],12).ravel()
            q=s-null['u']@(null['u'].T@s)
            direct,fraction=bounded_improvement(q@null['residual'].ravel(),q@q)
            actual,amplitude=a['bank'].score_point(condition,p,ph)
            errors.append(abs(actual-direct));spectral_errors.append(abs(float(grid[i])-direct))
            assert abs(amplitude-fraction)<1e-9
        result=dict(protected=protected,maximum_derivative_relative_error=maximum,
            maximum_point_score_absolute_error=max(errors),maximum_grid_score_absolute_error=max(spectral_errors),
            spectral_q_orthogonality_error=float(np.max(abs(a['q'].T@a['q']-np.eye(a['q'].shape[1])))))
        results.append(result)
        assert maximum<1e-6 and max(errors)<1e-9 and max(spectral_errors)<1e-6
    covariance=np.eye(48)+np.ones((48,48))/12
    w=time_whiten(np.eye(48),12)
    report=dict(passed=True,checks=results,
        shared_reference_whitening_error=float(np.max(abs(w@covariance@w.T-np.eye(48)))))
    signal=shadow.signal(np.zeros((2,len(old['time']),3)),REFERENCE_PERIOD,.7473,old['time'],order=16)
    rows=[]
    for name,area,mag,live,reads in [('8m_K6',49.,6.,240.,8),
        ('39m_K6',978.,6.,84*1.427,84),('39m_K3_invalid',978.,3.,84*1.427,84)]:
        spec=CountBudget(area_m2=area,k_magnitude=mag,live_seconds=live,reads_per_exposure=reads,read_noise_e=12.)
        uncertainty=electron_budget(tpl['wavelength'],tpl['stellar_relative'],tpl['transmission'],tpl['radiance'],spec)['sigma']
        for protected in [False,True]:
            a=setup_analysis(old,shadow,models[3],bank,uncertainty,keep,tpl,protected)
            s=time_whiten(signal@a['transform'],12)
            null=fit_null(a['model'],np.zeros_like(s),1.)
            q=s.ravel()-null['u']@(null['u'].T@s.ravel())
            rows.append(dict(scenario=name,protected=protected,
                whitened_full_signal_norm=float(np.linalg.norm(s)),known_orbit_kinematic_norm=float(np.linalg.norm(q))))
    report['information']=rows
    (root/'docs/co_observing_metric_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
