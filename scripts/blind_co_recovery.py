"""First period/phase-blind CO score search with held-out atmospheric nulls."""
import argparse
import csv
import hashlib
import json
import time as clock
from pathlib import Path
import numpy as np

from nonlinear_co_recovery import AXES,direct
from molecular_injection_recovery import build_molecular_geometry
from validate_blind_co import independent_geometry
from xds.nonlinear_atmosphere import CoherentCOModel,temporal_basis
from xds.moving_shadow import MovingShadow,REFERENCE_PERIOD,period_limits,separation_from_period
from xds.orbit_search import OrbitBank,fit_null,empirical_tail,binomial_interval,nonlinear_candidate_fit


FAMILIES=['reference','smooth','wind_pulse']
INJECTIONS=['showcase_reference','showcase_smooth','showcase_wind_pulse','uniform_phase']
TRUE_PHASE=.7473


def load_models(root,n_periods=73,n_phases=1024):
    old=np.load(root/'outputs/co_nonlinear_emulator.npz');moving=np.load(root/'outputs/co_blind_shadow.npz')
    if not np.allclose(old['basis'],moving['basis'],rtol=0,atol=1e-14):raise RuntimeError('spectral bases differ')
    shadow=MovingShadow(AXES,moving['library'],moving['left'],moving['centers'],moving['rotation'],float(moving['rp']))
    models={n:CoherentCOModel(AXES,old['library'],old['time'],n) for n in [3,5]}
    started=clock.monotonic();bank=OrbitBank(shadow,old['time'],n_periods,n_phases)
    print('Orbit bank ready:',bank.n_templates,'orbits;',len(bank.oi),'active rows;',clock.monotonic()-started,'s',flush=True)
    return old,shadow,models,bank


def generate(shadow,reference_model,rng,family,phase=None,moon=False):
    time=reference_model.time
    if family=='reference':p=np.zeros((2,len(time),3))
    else:
        controls=rng.uniform(-.3,.3,(2,3,3))
        p=np.einsum('gck,tk->gtc',controls,temporal_basis(time,3))
    if family=='wind_pulse':
        amplitude=float(rng.uniform(.3,.9));sign=float(rng.choice([-1,1]));center=float(rng.uniform(-1.6,1.6));width=float(rng.uniform(.2,.6))
        p[:,:,2]+=sign*(amplitude/1.5)*np.exp(-.5*((time-center)/width)**2)[None,:]
    fields,_=reference_model.fields_at(p,derivatives=False)
    truth=np.sum(fields[:,:,0],axis=0)-reference_model.reference
    signal=np.zeros_like(truth)
    if moon:signal=shadow.signal(p,REFERENCE_PERIOD,phase,time,order=16)
    return truth+signal,p,float(np.linalg.norm(signal))


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def direct_example(root,old):
    path=root/'outputs/co_blind_direct_truth.npz'
    if not path.exists():
        geometry=build_molecular_geometry(nr=52,n_scale_heights=16.,split_limbs=True)
        geometry=independent_geometry(geometry,REFERENCE_PERIOD,TRUE_PHASE)
        p=np.array([[-.23,-.16,.21],[.19,.22,-.27]])
        exact=direct(root,geometry,p)
        np.savez_compressed(path,data=(exact['reference']-old['reference']+exact['signal'])@old['basis'],
            period=REFERENCE_PERIOD,phase=TRUE_PHASE,parameters=p)
    return np.load(path)['data']


def run(root,args):
    if not (root/'docs/co_blind_forward_validation.json').exists():raise RuntimeError('run validate_blind_co.py first')
    validation=json.loads((root/'docs/co_blind_forward_validation.json').read_text())
    if not validation['passed']:raise RuntimeError('moving-shadow validation did not pass')
    old,shadow,models,bank=load_models(root,args.n_periods,args.n_phases)
    noise=float(old['noise'])/130.8285010091114
    config=dict(seed=63140,n_calibration=args.n_calibration,n_validation=args.n_validation,n_injections=args.n_injections,
        n_periods=args.n_periods,n_phases=args.n_phases,models=[3,5],noise_ppm=noise*1e6,
        period_bounds=period_limits().tolist(),phase_bounds=[0.,1.],true_period=REFERENCE_PERIOD,true_phase=TRUE_PHASE,
        detection_statistic='max bounded efficient-score gain after nonlinear null fit, top3 local score refinements',
        generators=FAMILIES,alpha=.05)
    fingerprint=hashlib.sha256(json.dumps(config,sort_keys=True).encode())
    for name in ['scripts/blind_co_recovery.py','src/xds/orbit_search.py','src/xds/moving_shadow.py','src/xds/nonlinear_atmosphere.py']:
        fingerprint.update((root/name).read_bytes())
    fingerprint.update((root/'docs/co_blind_library.json').read_bytes())
    config['procedure_fingerprint']=fingerprint.hexdigest()
    folder=root/'outputs/co_blind_runs'/args.tag;folder.mkdir(parents=True,exist_ok=True)
    metadata=folder/'config.json'
    if metadata.exists() and json.loads(metadata.read_text())!=config:raise RuntimeError('run configuration changed; choose a new --tag')
    metadata.write_text(json.dumps(config,indent=2)+'\n')
    raw=folder/'trials.jsonl';records=[]
    if raw.exists():records=[json.loads(line) for line in raw.read_text().splitlines() if line]
    completed={(r['stream'],r['family'],r['draw'],r['n_controls']) for r in records}
    streams=[('calibration',FAMILIES,args.n_calibration,0),('validation',FAMILIES,args.n_validation,1),
             ('injection',INJECTIONS,args.n_injections,2)]
    started=clock.monotonic()
    with raw.open('a') as handle:
        for stream,families,count,stream_id in streams:
            for family_id,label in enumerate(families):
                for draw in range(count):
                    if all((stream,label,draw,n) in completed for n in models):continue
                    rng=np.random.default_rng(np.random.SeedSequence([63140,stream_id,family_id,draw]))
                    moon=stream=='injection'
                    family=label.removeprefix('showcase_') if moon and label!='uniform_phase' else ('smooth' if label=='uniform_phase' else label)
                    phase=float(rng.uniform(0,1)) if label=='uniform_phase' else TRUE_PHASE
                    truth,p,signal_norm=generate(shadow,models[3],rng,family,phase,moon)
                    data=truth+noise*rng.normal(size=truth.shape)
                    for n,model in models.items():
                        key=(stream,label,draw,n)
                        if key in completed:continue
                        null=fit_null(model,data,noise)
                        search,condition=bank.search(null,noise)
                        oracle,oracle_f=bank.score_point(condition,REFERENCE_PERIOD,TRUE_PHASE)
                        row=dict(stream=stream,family=label,draw=draw,n_controls=n,statistic=search['statistic'],
                            grid_statistic=search['grid_max'],oracle_statistic=oracle,period=search['period'],phase=search['phase'],
                            blocking=search['blocking'],true_period=REFERENCE_PERIOD if moon else None,
                            true_phase=phase if moon else None,physical_signal_norm=signal_norm/noise,
                            null_chi2=null['best']['chi2'],null_start_spread=null['start_spread'],
                            null_active_bounds=null['best']['active_bounds'],
                            refinement_failures=sum(not c['refinement_success'] for c in search['candidates']))
                        handle.write(json.dumps(row)+'\n');handle.flush();records.append(row);completed.add(key)
                    if draw%8==0 or draw==count-1:print(stream,label,draw+1,'/',count,'elapsed',round(clock.monotonic()-started,1),'s',flush=True)
    # Fixed noiseless examples and independent direct-source off-grid injection.
    observed=[];fits=[]
    for i,label in enumerate(['reference_moon','smooth_moon','wind_pulse','wind_pulse_moon','direct_offgrid']):
        rng=np.random.default_rng([63141,i]);family='reference' if label=='reference_moon' else ('smooth' if label=='smooth_moon' else 'wind_pulse')
        truth,_,_=generate(shadow,models[3],rng,family,TRUE_PHASE,label!='wind_pulse')
        if label=='direct_offgrid':truth=direct_example(root,old)
        for n,model in models.items():
            null=fit_null(model,truth,noise);search,condition=bank.search(null,noise,keep_grid=True)
            oracle,_=bank.score_point(condition,REFERENCE_PERIOD,TRUE_PHASE)
            np.savez_compressed(folder/f'grid_{label}_{n}.npz',statistic=search.pop('grid'),periods=bank.period_axis,phases=bank.phase_axis)
            local=[]
            for candidate in search['candidates']:
                best,runs=nonlinear_candidate_fit(model,shadow,truth,noise,null,candidate)
                local.append(dict(best=best,runs=runs))
            best=min(local,key=lambda c:c['best']['chi2'])['best']
            observed.append(dict(case=label,n_controls=n,statistic=search['statistic'],oracle_statistic=oracle,
                period=search['period'],phase=search['phase'],blocking=search['blocking'],
                nonlinear_period=best['period'],nonlinear_phase=best['phase'],nonlinear_blocking=best['blocking'],
                nonlinear_root_delta_chi2=float(np.sqrt(max(null['best']['chi2']-best['chi2'],0))),
                nonlinear_residual_norm=float(np.sqrt(best['chi2'])),null_residual_norm=float(np.sqrt(null['best']['chi2']))))
            fits.append(dict(case=label,n_controls=n,candidates=search['candidates'],nonlinear=local))
            print('Observed',label,n,search['statistic'],search['period'],search['phase'],flush=True)
    (folder/'candidate_fits.json').write_text(json.dumps(fits,indent=2)+'\n')
    summarize(root,folder,records,observed,config,args.tag)


def summarize(root,folder,records,observed,config,tag):
    alpha=config['alpha'];calibration={}
    for n in [3,5]:
        for family in FAMILIES:
            rows=[r for r in records if r['stream']=='calibration' and r['family']==family and r['n_controls']==n]
            calibration[n,family]=np.array([r['statistic'] for r in rows])
    def tails(n,value,oracle=False):
        p={}
        for family in FAMILIES:
            if oracle:
                values=[r['oracle_statistic'] for r in records if r['stream']=='calibration' and r['family']==family and r['n_controls']==n]
            else:values=calibration[n,family]
            p[family]=empirical_tail(values,value)
        return p
    for row in records:
        if row['stream']=='calibration':continue
        p=tails(row['n_controls'],row['statistic'])
        row.update(p_gaussian=p['reference'],p_smooth=p['smooth'],p_wind=p['wind_pulse'],p_robust=max(p.values()))
        row['flag_robust']=row['p_robust']<=alpha;row['flag_gaussian']=row['p_gaussian']<=alpha
        # The oracle is fixed at the showcase coordinates. Do not assign its
        # calibration to uniform-phase injections with a different true orbit.
        row['oracle_p_robust']=max(tails(row['n_controls'],row['oracle_statistic'],True).values()) if row['family']!='uniform_phase' else None
    for row in observed:
        p=tails(row['n_controls'],row['statistic']);row.update(p_gaussian=p['reference'],p_smooth=p['smooth'],p_wind=p['wind_pulse'],p_robust=max(p.values()))
        row['oracle_p_robust']=max(tails(row['n_controls'],row['oracle_statistic'],True).values())
    # Report both predeclared temporal models, and a conservative consensus.
    # Never select the smaller p-value after seeing which model is favorable.
    for row in records:
        if row['stream']=='calibration':continue
        peers=[r for r in records if r['stream']==row['stream'] and r['family']==row['family'] and r['draw']==row['draw']]
        row['p_consensus']=max(r['p_robust'] for r in peers)
        row['flag_consensus']=row['p_consensus']<=alpha
    for row in observed:
        row['p_consensus']=max(r['p_robust'] for r in observed if r['case']==row['case'])
    rates=[]
    for n in [3,5]:
        for stream,labels in [('validation',FAMILIES),('injection',INJECTIONS)]:
            for label in labels:
                rows=[r for r in records if r['stream']==stream and r['family']==label and r['n_controls']==n]
                for rule in ['robust','gaussian']:
                    k=sum(r['flag_'+rule] for r in rows);lo,hi=binomial_interval(k,len(rows))
                    rates.append(dict(stream=stream,family=label,n_controls=n,rule=rule,n=len(rows),flags=k,
                        fraction=k/len(rows),ci95_low=lo,ci95_high=hi))
    for stream,labels in [('validation',FAMILIES),('injection',INJECTIONS)]:
        for label in labels:
            rows=[r for r in records if r['stream']==stream and r['family']==label and r['n_controls']==3]
            k=sum(r['flag_consensus'] for r in rows);lo,hi=binomial_interval(k,len(rows))
            rates.append(dict(stream=stream,family=label,n_controls=0,rule='consensus',n=len(rows),flags=k,
                fraction=k/len(rows),ci95_low=lo,ci95_high=hi))
    inference=[]
    for row in records:
        if row['stream']!='injection':continue
        row['period_relative_error']=float((row['period']-row['true_period'])/row['true_period'])
        row['phase_circular_error']=float(abs((row['phase']-row['true_phase']+.5)%1-.5))
        inference.append(row)
    metrics=[]
    for n in [3,5]:
        for label in INJECTIONS:
            for selected_only in [False,True]:
                rows=[r for r in inference if r['n_controls']==n and r['family']==label and (not selected_only or r['flag_robust'])]
                for parameter in ['period','period_relative_error','phase_circular_error']:
                    q=np.quantile([r[parameter] for r in rows],[.16,.5,.84]).tolist() if rows else [None]*3
                    metrics.append(dict(n_controls=n,family=label,flagged_only=selected_only,parameter=parameter,n=len(rows),q16=q[0],q50=q[1],q84=q[2]))
    calibration_rows=[{k:r[k] for k in ['family','draw','n_controls','statistic','oracle_statistic']} for r in records if r['stream']=='calibration']
    decision_rows=[r for r in records if r['stream']!='calibration']
    # Uniform field names despite nullable physical truth metadata.
    fields=list(dict.fromkeys(k for r in decision_rows for k in r))
    decision_rows=[{k:r.get(k) for k in fields} for r in decision_rows]
    for filename,rows in [('calibration',calibration_rows),('decisions',decision_rows),('rates',rates),('orbit_errors',metrics),('examples',observed)]:
        write_csv(root/f'docs/co_blind_{tag}_{filename}.csv',rows)
    summary=dict(config,tag=tag,n_grid_templates=config['n_periods']*config['n_phases'],
        empirical_p_resolution=1/(config['n_calibration']+1),
        represented_null='maximum tail probability across three specified generative distributions, not all possible atmospheres',
        nonlinear_refit_fap_calibrated=False,score_fap_calibrated=True,
        physical_scope='unknown period and phase; fixed circular aligned orbit family, moon radius and host',
        max_null_start_spread=max(r['null_start_spread'] for r in records),
        refinement_limit_count=sum(r['refinement_failures'] for r in records),
        native_noise_blind_calibration=False)
    (root/f'docs/co_blind_{tag}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (root/f'docs/co_blind_{tag}_candidate_fits.json').write_text((folder/'candidate_fits.json').read_text())


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tag',default='pilot');parser.add_argument('--n-calibration',type=int,default=128)
    parser.add_argument('--n-validation',type=int,default=64);parser.add_argument('--n-injections',type=int,default=32)
    parser.add_argument('--n-periods',type=int,default=73);parser.add_argument('--n-phases',type=int,default=1024)
    args=parser.parse_args()
    if min(args.n_calibration,args.n_validation,args.n_injections)<1:parser.error('positive ensemble sizes required')
    run(Path(__file__).resolve().parents[1],args)


if __name__=='__main__':main()
