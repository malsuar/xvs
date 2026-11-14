"""Independent finite-search and moving-spectrum derivative controls."""
import json
from pathlib import Path
import numpy as np
from blind_co_recovery import load_models,generate,direct_example,TRUE_PHASE
from xds.orbit_search import OrbitBank,fit_null
from xds.moving_shadow import REFERENCE_PERIOD


def main():
    root=Path(__file__).resolve().parents[1]
    old,shadow,models,primary=load_models(root)
    fine=OrbitBank(shadow,old['time'],97,1536)
    noise=float(old['noise'])/130.8285010091114;rows=[]
    for i,label in enumerate(['reference','smooth','wind_pulse','direct_offgrid']):
        rng=np.random.default_rng([63160,i])
        if label=='direct_offgrid':data=direct_example(root,old)
        else:
            truth,_,_=generate(shadow,models[3],rng,label)
            data=truth+noise*rng.normal(size=truth.shape)
        for n,model in models.items():
            null=fit_null(model,data,noise)
            a,_=primary.search(null,noise);b,_=fine.search(null,noise)
            row=dict(case=label,n_controls=n,primary_statistic=a['statistic'],fine_statistic=b['statistic'],
                relative_statistic_difference=float(abs(a['statistic']-b['statistic'])/max(1,b['statistic'])),
                primary_period=a['period'],fine_period=b['period'],primary_phase=a['phase'],fine_phase=b['phase'])
            rows.append(row);print(row,flush=True)
    p=np.zeros((2,48,3));signal,ds=shadow.signal(p,REFERENCE_PERIOD,TRUE_PHASE,old['time'],True,16)
    derivatives=[];step=1e-4
    for g in range(2):
        for parameter in range(3):
            plus=p.copy();minus=p.copy();plus[g,:,parameter]+=step;minus[g,:,parameter]-=step
            numerical=(shadow.signal(plus,REFERENCE_PERIOD,TRUE_PHASE,old['time'],order=16)-
                       shadow.signal(minus,REFERENCE_PERIOD,TRUE_PHASE,old['time'],order=16))/(2*step)
            derivatives.append(dict(limb=g,parameter=parameter,
                relative_error=float(np.linalg.norm(numerical-ds[g,parameter])/np.linalg.norm(numerical))))
    result=dict(primary_grid=[73,1024],fine_grid=[97,1536],comparisons=rows,derivative_checks=derivatives,
        all_statistic_changes_below_1percent=all(r['relative_statistic_difference']<.01 for r in rows),
        note='finite predefined cases, not proof of a continuum global maximum; no retuning after calibration')
    (root/'docs/co_blind_grid_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
