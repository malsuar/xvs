"""Independent contamination/recovery pilot with a count-based noise budget."""
import argparse
import hashlib
import json
import subprocess
import sys
import time as clock
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d,percentile_filter
from scipy.linalg import qr
from blind_co_recovery import load_models,generate,write_csv,TRUE_PHASE
from xds.molecular import C,bin_average
from xds.moving_shadow import MovingShadow,REFERENCE_PERIOD
from xds.observing_contamination import (CountBudget,electron_budget,time_whiten,
    spectral_metric,ObservedCOModel,ReferenceOrbitBank)
from xds.orbit_search import fit_null,empirical_tail,binomial_interval


FAMILIES=['atmosphere','stellar','telluric','combined']


def templates(root,old):
    """Physical stellar flux and sky spectra at R=1e5, no invented lines."""
    sky=np.loadtxt(root/'data/observing/skycalc_paranal.csv',delimiter=',',skiprows=1)
    star=np.loadtxt(root/'data/observing/phoenix_5700_4p5.csv',delimiter=',',skiprows=1)
    x=C*np.log(old['wavenumber']/4300);step=np.median(np.diff(x))
    fine=np.arange(x[0]-100,x[-1]+100,.25);wave=1e4/(4300*np.exp(fine/C))
    raw=np.interp(wave,star[:,0],star[:,1])
    continuum=percentile_filter(raw,percentile=95,size=801,mode='nearest')
    local=np.clip(raw/continuum,.01,1.1)
    # Disk rotation with linear limb darkening; PHOENIX surface flux is used
    # as a surrogate local profile, not a mu-resolved stellar atmosphere.
    velocities=np.arange(-2.,2.00001,.25);z=velocities/2
    kernel=2*(1-.6)*np.sqrt(np.maximum(1-z*z,0))+.5*np.pi*.6*(1-z*z)
    kernel/=kernel.sum()
    disk=np.convolve(np.pad(local,(8,8),mode='edge'),kernel,mode='valid')
    sigma=step/np.sqrt(8*np.log(2))/.25
    edges=np.r_[x-step/2,x[-1]+step/2]
    def observe(y):return bin_average(fine,gaussian_filter1d(y,sigma,mode='nearest'),edges)
    transmission=observe(np.interp(wave*1000,sky[:,0],sky[:,1]))
    radiance=observe(np.interp(wave*1000,sky[:,0],sky[:,2]))
    local=observe(local);disk=observe(disk)
    tau=-np.log(np.clip(transmission,1e-6,1))
    nuisance=np.array([tau,np.gradient(tau,x),tau*tau,1-disk,np.gradient(disk,x),np.gradient(np.gradient(disk,x),x)])
    # Approximate monochromatic Ks normalization with a 5700 K Planck slope.
    wavelength=1e4/old['wavenumber'];a=14387.76877/5700
    slope=(2.159/wavelength)**3*np.expm1(a/2.159)/np.expm1(a/wavelength)
    return dict(x=x,wavelength=wavelength,transmission=transmission,radiance=radiance,
        local=local,disk=disk,tau=tau,nuisance=nuisance,stellar_relative=disk*slope)


def contamination(tpl,time,rng,family,offset=0.):
    """Correlated post-correction residuals, signed in absorption units.

    Telluric optical-depth error is 0.1--0.5%, wavelength error 3--15 m/s.
    Stellar RM/phenomenological CLV residual retains 3--15% of a surrogate.
    These ranges are sensitivity assumptions, not measured pipeline residuals.
    """
    x=tpl['x'];t=np.asarray(time)
    orbit_phase=2*np.pi*t/720
    a=(30/365.25)**(2/3)*149597870.7
    chord=a*np.sin(orbit_phase)/695700
    mu=np.sqrt(1-chord**2)
    vp=(2*np.pi*a/(720*3600))*np.sin(orbit_phase)
    def shift(profile,v):return np.interp(x-v,x,profile,left=profile[0],right=profile[-1])
    cov=np.exp(-.5*((t[:,None]-t[None,:])/.7)**2)+1e-9*np.eye(len(t))
    z=np.linalg.cholesky(cov)@rng.normal(size=(len(t),3))
    stellar=np.zeros((len(t),len(x)));telluric=np.zeros_like(stellar)
    retention=rng.uniform(.03,.15)
    optical_error=rng.uniform(.001,.005);drift=rng.uniform(.003,.015)
    intensity=(1-.4*(1-mu)-.25*(1-mu)**2)/(1-.4/3-.25/6)
    depth=(1.1*71492/695700)**2*intensity
    for k in range(len(t)):
        disk=shift(tpl['disk'],-vp[k])
        local=1-(1-shift(tpl['local'],2*chord[k]-vp[k]))*(1+.2*(1-mu[k]))
        stellar[k]=retention*depth[k]/(1-depth[k])*(local/disk-1)*(1+.1*z[k,0])
        telluric[k]=optical_error*z[k,1]*shift(tpl['tau'],offset-vp[k]+.1*t[k])
        telluric[k]+=drift*z[k,2]*shift(np.gradient(tpl['tau'],x),offset-vp[k]+.1*t[k])
    if family not in ['stellar','combined']:stellar[:]=0
    if family not in ['telluric','combined']:telluric[:]=0
    return stellar+telluric,stellar,telluric


def setup_analysis(old,shadow,model,bank,sigma,keep,tpl,protected):
    q,transform=spectral_metric(old['basis'],sigma,keep,tpl['nuisance'] if protected else None)
    v,r=qr(transform.T@shadow.rotation,mode='economic')
    moving=MovingShadow(shadow.axes,shadow.library@r.T,
        np.array([p(shadow.centers) for p in shadow.position]),shadow.centers,v,shadow.rp)
    observed=ObservedCOModel(model,transform)
    return dict(q=q,transform=transform,model=observed,bank=ReferenceOrbitBank(bank,moving))


def summarize(root,records,tag):
    decisions=[];rates=[]
    for level in ['control','39m_K6']:
        for method in ['unprotected','protected']:
            selected=[r for r in records if r['level']==level and r['method']==method]
            calibration={f:np.array([r['statistic'] for r in selected if r['stream']=='calibration' and r['family']==f]) for f in FAMILIES}
            if not all(len(v) for v in calibration.values()):continue
            for row in selected:
                if row['stream']=='calibration':continue
                p={f:empirical_tail(v,row['statistic']) for f,v in calibration.items()}
                row=dict(row,p_atmosphere=p['atmosphere'],p_robust=max(p.values()))
                row.update(flag_atmosphere=row['p_atmosphere']<=.05,flag_robust=row['p_robust']<=.05)
                decisions.append(row)
            for stream,labels in [('validation',FAMILIES),('injection',['clean_fixed','combined_fixed','combined_random'])]:
                for label in labels:
                    group=[r for r in decisions if r['level']==level and r['method']==method and r['stream']==stream and r['family']==label]
                    if not group:continue
                    for rule in ['atmosphere','robust']:
                        k=sum(r['flag_'+rule] for r in group);lo,hi=binomial_interval(k,len(group))
                        rates.append(dict(level=level,method=method,stream=stream,family=label,calibration=rule,
                            n=len(group),flags=k,fraction=k/len(group),ci95_low=lo,ci95_high=hi))
    if decisions:write_csv(root/f'docs/co_observing_{tag}_decisions.csv',decisions)
    if rates:write_csv(root/f'docs/co_observing_{tag}_rates.csv',rates)
    write_csv(root/f'docs/co_observing_{tag}_calibration.csv',[r for r in records if r['stream']=='calibration'])


def run(args):
    root=Path(__file__).resolve().parents[1]
    old,shadow,models,bank=load_models(root)
    model=models[3];tpl=templates(root,old)
    keep=(tpl['transmission']>=.3)&(tpl['disk']>=.2)
    budgets=[];budget_arrays={}
    for label,area,mag in [('8m_K6',49.,6.),('39m_K6',978.,6.),('39m_K3',978.,3.)]:
        # Detector linearity is checked BEFORE adopting the count forecast.
        # The 39m bright case is retained only as an invalid photon upper bound.
        spec=CountBudget(area_m2=area,k_magnitude=mag,live_seconds=240. if area==49 else 84*1.427,
            reads_per_exposure=8 if area==49 else 84,read_noise_e=12.)
        values=electron_budget(tpl['wavelength'],tpl['stellar_relative'],tpl['transmission'],tpl['radiance'],spec)
        budget_arrays[label]=values
        budgets.append(dict(label=label,**asdict(spec),retained_bins=int(keep.sum()),total_bins=len(keep),
            median_sigma_ppm=float(np.median(values['sigma'][keep])*1e6),
            q16_sigma_ppm=float(np.quantile(values['sigma'][keep],.16)*1e6),
            q84_sigma_ppm=float(np.quantile(values['sigma'][keep],.84)*1e6),
            median_source_e=float(np.median(values['source'][keep])),
            median_sky_e=float(np.median(values['sky'][keep])),median_thermal_e=float(np.median(values['thermal'][keep])),
            peak_source_e_per_read_pixel=float(2*np.max(values['source'])/(spec.reads_per_exposure*spec.pixels_per_bin)),
            linearity_limit_e=12900.,valid_linearity=bool(2*np.max(values['source'])/(spec.reads_per_exposure*spec.pixels_per_bin)<12900.)))
    write_csv(root/'docs/co_observing_budgets.csv',budgets)
    noise_profiles={'control':budget_arrays['39m_K6']['sigma']*11.394654178131821e-6/np.median(budget_arrays['39m_K6']['sigma'][keep]),
        '39m_K6':budget_arrays['39m_K6']['sigma']}
    config=dict(seed=63150,n_calibration=args.n_calibration,n_validation=args.n_validation,n_injections=args.n_injections,
        n_reference=12,n_controls=3,n_periods=73,n_phases=1024,alpha=.05,
        kept_bins=int(keep.sum()),families=FAMILIES,levels=list(noise_profiles),
        assumptions='R100000; count table specifies live time/reads; one-hour shared OOT; fixed median airmass1.3')
    h=hashlib.sha256(json.dumps(config,sort_keys=True).encode())
    for name in ['scripts/observing_contamination_pilot.py','src/xds/observing_contamination.py',
                 'src/xds/orbit_search.py','src/xds/nonlinear_atmosphere.py','src/xds/moving_shadow.py',
                 'data/observing/skycalc_paranal.csv','data/observing/phoenix_5700_4p5.csv']:
        h.update((root/name).read_bytes())
    config['fingerprint']=h.hexdigest()
    folder=root/'outputs/co_observing_runs'/args.tag;folder.mkdir(parents=True,exist_ok=True)
    path=folder/f'config_worker{args.worker}.json'
    if path.exists() and json.loads(path.read_text())!=config:raise RuntimeError('configuration changed: choose new tag')
    path.write_text(json.dumps(config,indent=2)+'\n')
    raw=folder/f'trials_worker{args.worker}.jsonl'
    records=[json.loads(r) for r in raw.read_text().splitlines()] if raw.exists() else []
    completed={(r['level'],r['method'],r['stream'],r['family'],r['draw']) for r in records}
    diagnostics=[]
    started=clock.monotonic()
    with raw.open('a') as handle:
        for level,sigma in noise_profiles.items():
            analyses={name:setup_analysis(old,shadow,model,bank,sigma,keep,tpl,name=='protected') for name in ['unprotected','protected']}
            p=np.zeros((2,len(old['time']),3));moon=shadow.signal(p,REFERENCE_PERIOD,TRUE_PHASE,old['time'],order=16)
            for name,analysis in analyses.items():
                observed_signal=time_whiten(moon@analysis['transform'],12)
                null=fit_null(analysis['model'],np.zeros_like(observed_signal),1.)
                direct_null=dict(null,residual=observed_signal)
                cond=analysis['bank'].conditional(direct_null,1.)
                value,f=analysis['bank'].score_point(cond,REFERENCE_PERIOD,TRUE_PHASE)
                diagnostics.append(dict(level=level,method=name,full_signal_norm=float(np.linalg.norm(observed_signal)),
                    known_orbit_kinematic_norm=float(np.sqrt(value)),best_blocking=f))
            for stream,labels,count,sid in [('calibration',FAMILIES,args.n_calibration,0),
                ('validation',FAMILIES,args.n_validation,1),
                ('injection',['clean_fixed','combined_fixed','combined_random'],args.n_injections,2)]:
                for fi,label in enumerate(labels):
                    for draw in range(count):
                        if draw%args.workers!=args.worker:continue
                        if all((level,name,stream,label,draw) in completed for name in analyses):continue
                        rng=np.random.default_rng(np.random.SeedSequence([63150,sid,fi,draw]))
                        phase=float(rng.uniform()) if label=='combined_random' else TRUE_PHASE
                        truth,_,signal_norm=generate(shadow,model,rng,'wind_pulse',phase,stream=='injection')
                        family=('atmosphere' if label=='clean_fixed' else 'combined') if stream=='injection' else label
                        contaminants,stellar,telluric=contamination(tpl,old['time'],rng,family)
                        # Full-bin draws make bare/protected and noise levels paired.
                        white=rng.normal(size=(len(old['time']),int(keep.sum())))
                        reference=rng.normal(size=int(keep.sum()))/np.sqrt(12)
                        for name,analysis in analyses.items():
                            key=(level,name,stream,label,draw)
                            if key in completed:continue
                            q=analysis['q']
                            data=truth@analysis['transform']+(contaminants[:,keep]/sigma[keep])@q
                            data+=(white+reference)@q
                            null=fit_null(analysis['model'],time_whiten(data,12),1.)
                            search,_=analysis['bank'].search(null,1.)
                            row=dict(level=level,method=name,stream=stream,family=label,draw=draw,
                                statistic=search['statistic'],grid_statistic=search['grid_max'],period=search['period'],phase=search['phase'],
                                blocking=search['blocking'],physical_signal_norm=signal_norm,
                                phase_error=float(abs((search['phase']-phase+.5)%1-.5)) if stream=='injection' else None,
                                stellar_peak_ppm=float(np.max(abs(stellar[:,keep]))*1e6),
                                telluric_peak_ppm=float(np.max(abs(telluric[:,keep]))*1e6),
                                null_start_spread=null['start_spread'],refinement_failures=sum(not c['refinement_success'] for c in search['candidates']))
                            handle.write(json.dumps(row)+'\n');handle.flush();records.append(row);completed.add(key)
                        if draw%8==0 or draw==count-1:print(level,stream,label,draw+1,'/',count,'elapsed',round(clock.monotonic()-started,1),'s',flush=True)
            if args.workers==1:summarize(root,records,args.tag)
    if args.worker==0:
        write_csv(root/f'docs/co_observing_{args.tag}_information.csv',diagnostics)
        (root/f'docs/co_observing_{args.tag}_config.json').write_text(json.dumps(config,indent=2)+'\n')
        np.savez_compressed(root/f'outputs/co_observing_{args.tag}_spectra.npz',**tpl,keep=keep,
            sigma_control=noise_profiles['control'],sigma_real=noise_profiles['39m_K6'],time=old['time'])
    print('COMPLETE',len(records),'searches',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot')
    parser.add_argument('--n-calibration',type=int,default=64)
    parser.add_argument('--n-validation',type=int,default=32)
    parser.add_argument('--n-injections',type=int,default=24)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--worker',type=int,default=-1)
    args=parser.parse_args()
    if args.workers<1:raise ValueError('positive worker count required')
    if args.workers==1:
        args.worker=0;run(args)
    elif args.worker>=0:run(args)
    else:
        root=Path(__file__).resolve().parents[1]
        folder=root/'outputs/co_observing_runs'/args.tag;folder.mkdir(parents=True,exist_ok=True)
        processes=[];handles=[]
        try:
            for worker in range(args.workers):
                log=(folder/f'worker{worker}.log').open('a');handles.append(log)
                command=[sys.executable,__file__,*sys.argv[1:],'--worker',str(worker)]
                processes.append(subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT))
            while any(p.poll() is None for p in processes):
                if any(p.poll() not in [None,0] for p in processes):raise RuntimeError('worker failed: inspect logs')
                clock.sleep(15)
                count=sum(len(p.read_text().splitlines()) for p in folder.glob('trials_worker*.jsonl'))
                print('Checkpointed searches:',count,flush=True)
            if any(p.returncode for p in processes):raise RuntimeError('worker failed: inspect logs')
        finally:
            for p in processes:
                if p.poll() is None:p.terminate()
            for p in processes:p.wait()
            for h in handles:h.close()
        records=[]
        configs=[]
        for worker in range(args.workers):
            records.extend(json.loads(r) for r in (folder/f'trials_worker{worker}.jsonl').read_text().splitlines())
            configs.append(json.loads((folder/f'config_worker{worker}.json').read_text()))
        if any(c!=configs[0] for c in configs):raise RuntimeError('worker procedures differ')
        expected=2*2*(4*(args.n_calibration+args.n_validation)+3*args.n_injections)
        if len(records)!=expected:raise RuntimeError(f'incomplete sample: {len(records)} != {expected}')
        summarize(root,records,args.tag)
        print('COMPLETE',len(records),'searches',flush=True)
