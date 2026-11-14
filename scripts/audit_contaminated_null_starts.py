"""Post-run numerical audit, never substituted into the frozen calibration."""
import csv
import argparse
import json
from pathlib import Path
import numpy as np
from blind_co_recovery import load_models,generate,TRUE_PHASE
from observing_contamination_pilot import templates,contamination,setup_analysis,FAMILIES
from xds.observing_contamination import CountBudget,electron_budget,time_whiten
from xds.nonlinear_atmosphere import fit_coherent
from xds.orbit_search import fit_null


def main(tag):
    root=Path(__file__).resolve().parents[1]
    with (root/f'docs/co_observing_{tag}_calibration.csv').open() as handle:rows=list(csv.DictReader(handle))
    with (root/f'docs/co_observing_{tag}_decisions.csv').open() as handle:rows+=list(csv.DictReader(handle))
    selected=[max([r for r in rows if r['level']==level and r['method']==method],key=lambda r:float(r['null_start_spread']))
        for level in ['control','39m_K6'] for method in ['unprotected','protected']]
    selected.append(next(r for r in rows if r['level']=='control' and r['method']=='unprotected' and r['stream']=='calibration' and r['family']=='atmosphere' and r['draw']=='0'))
    old,shadow,models,bank=load_models(root)
    model=models[3];tpl=templates(root,old);keep=(tpl['transmission']>=.3)&(tpl['disk']>=.2)
    budget=CountBudget(area_m2=978.,live_seconds=84*1.427,reads_per_exposure=84,read_noise_e=12.)
    real=electron_budget(tpl['wavelength'],tpl['stellar_relative'],tpl['transmission'],tpl['radiance'],budget)['sigma']
    results=[]
    for row in selected:
        sigma=real if row['level']=='39m_K6' else real*11.394654178131821e-6/np.median(real[keep])
        analysis=setup_analysis(old,shadow,model,bank,sigma,keep,tpl,row['method']=='protected')
        stream=row['stream'];label=row['family'];draw=int(row['draw'])
        labels=FAMILIES if stream!='injection' else ['clean_fixed','combined_fixed','combined_random']
        rng=np.random.default_rng(np.random.SeedSequence([63150,{'calibration':0,'validation':1,'injection':2}[stream],labels.index(label),draw]))
        phase=float(rng.uniform()) if label=='combined_random' else TRUE_PHASE
        truth,_,_=generate(shadow,model,rng,'wind_pulse',phase,stream=='injection')
        family=('atmosphere' if label=='clean_fixed' else 'combined') if stream=='injection' else label
        contaminants,_,_=contamination(tpl,old['time'],rng,family)
        white=rng.normal(size=(len(old['time']),int(keep.sum())))
        reference=rng.normal(size=int(keep.sum()))/np.sqrt(12)
        data=truth@analysis['transform']+(contaminants[:,keep]/sigma[keep])@analysis['q']
        data+=(white+reference)@analysis['q']
        data=time_whiten(data,12)
        null=fit_null(analysis['model'],data,1.)
        original,_=analysis['bank'].search(null,1.)
        assert abs(original['statistic']-float(row['statistic']))<1e-6
        starts=[np.zeros(model.n_parameters),np.resize([.2,-.2,.1],model.n_parameters)]
        starts+=list(np.random.default_rng(63153).uniform(-.7,.7,(6,model.n_parameters)))
        best,runs=fit_coherent(analysis['model'],data,1.,starts,max_nfev=300)
        if not all(r['success'] for r in runs):best,runs=fit_coherent(analysis['model'],data,1.,starts,max_nfev=1000)
        pred,jac,_=analysis['model'].evaluate(best['coefficients'],0.)
        design=jac.reshape(-1,model.n_parameters);norm=np.linalg.norm(design,axis=0)
        u,s,_=np.linalg.svd(design/norm,full_matrices=False);u=u[:,s>s[0]*1e-8]
        other=dict(best=best,residual=data-pred,u=u,parameters=model.trajectory(best['coefficients']))
        audited,_=analysis['bank'].search(other,1.)
        result=dict(level=row['level'],method=row['method'],stream=stream,family=label,draw=draw,
            original_start_spread=float(row['null_start_spread']),
            null_chi2_two_starts=null['best']['chi2'],null_chi2_eight_starts=best['chi2'],
            improvement= null['best']['chi2']-best['chi2'],all_audit_starts_success=all(r['success'] for r in runs),
            original_score=original['statistic'],audited_score=audited['statistic'],
            original_period=original['period'],audited_period=audited['period'])
        results.append(result);print(result,flush=True)
    report=dict(selection='largest recorded null-start spread per level/reduction, plus a clean calibration control',
        note='Post-run numerical audit only; no modified trial replaces the frozen calibration or receives its p value',checks=results)
    (root/f'docs/co_observing_{tag}_null_start_audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot')
    main(parser.parse_args().tag)
