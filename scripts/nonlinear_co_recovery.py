"""Bounded, coherent nonlinear CO recovery; a controlled known-orbit experiment."""
import argparse
import csv
import hashlib
import itertools
import inspect
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from molecular_injection_recovery import build_molecular_geometry,synthesize
from xds.nonlinear_atmosphere import CoherentCOModel,fit_coherent,physical_parameters,temporal_basis


AXES=[np.linspace(-1,1,n) for n in (5,7,9)]
RANK=128


def direct(root,geometry,p):
    t,x,v=physical_parameters(np.asarray(p))
    return synthesize(root,geometry=geometry,temperature=t,co_mixing_ratio=x,wind_km_s=v,
        shift_interpolation='cubic',return_limb_components=True)


def spectral_fields(result):
    return np.stack([result['limb_reference'],result['limb_signal']],axis=2)


def emulator_cache_key(root,geometry):
    fingerprint=hashlib.sha256()
    for path in ['src/xds/molecular.py','src/xds/detectability.py','scripts/molecular_injection_recovery.py','data/exomol/CO_Li2015/manifest.json']:
        fingerprint.update((root/path).read_bytes())
    for function in [emulator_cache_key,build_emulator,direct,spectral_fields,physical_parameters]:
        fingerprint.update(inspect.getsource(function).encode())
    for name in ['full','missing','shifts','radii','rp','scale','time']:
        value=np.ascontiguousarray(geometry[name]);fingerprint.update(str(value.shape).encode());fingerprint.update(value.tobytes())
    fingerprint.update(json.dumps(dict(axes=[a.tolist() for a in AXES],rank=RANK,basis_seed=63120)).encode())
    return fingerprint.hexdigest()[:16]


def build_emulator(root,geometry):
    key=emulator_cache_key(root,geometry)
    cache=root/'outputs/co_nonlinear_cache'/key;cache.mkdir(parents=True,exist_ok=True)
    base=direct(root,geometry,np.zeros((2,3)))
    if (cache/'basis.npz').exists():basis=np.load(cache/'basis.npz')['basis']
    else:
        rng=np.random.default_rng(63120)
        points=[np.array(p) for p in itertools.product([-1.,1.],repeat=3)]+[np.zeros(3)]
        points.extend(rng.uniform(-1,1,(9,3)))
        rows=[]
        for i,p in enumerate(points):
            fields=spectral_fields(direct(root,geometry,np.tile(p,(2,1)))).reshape(-1,base['reference'].shape[-1])
            norm=np.linalg.norm(fields,axis=1);keep=norm>1e-22
            rows.append(fields[keep]/norm[keep,None])
            print(f'Spectral basis snapshot {i+1}/{len(points)}',flush=True)
        a=np.concatenate(rows)
        basis=np.linalg.qr(a.T@rng.normal(size=(len(a),RANK)),mode='reduced')[0]
        # A power iteration improves the approximate range; independent direct
        # spectra below set the actual accuracy requirement, not the rank alone.
        basis=np.linalg.qr(a.T@(a@basis),mode='reduced')[0]
        np.savez_compressed(cache/'basis.npz',basis=basis)
    library=np.zeros((*map(len,AXES),2,len(geometry['time']),2,RANK))
    indices=list(itertools.product(*(range(len(a)) for a in AXES)))
    for n,index in enumerate(indices):
        path=cache/('node_'+'_'.join(map(str,index))+'.npz')
        if path.exists():coeff=np.load(path)['coefficients']
        else:
            p=np.array([a[i] for a,i in zip(AXES,index)])
            fields=spectral_fields(direct(root,geometry,np.tile(p,(2,1))))
            coeff=fields@basis
            np.savez_compressed(path,coefficients=coeff)
        library[index]=coeff
        if n%15==0 or n==len(indices)-1:print(f'Nonlinear grid {n+1}/{len(indices)}',flush=True)
    return basis,library,base,key


def direct_trajectory(root,geometry,p):
    """Independent line-by-line injection, frozen within each exposure."""
    fields=[]
    for k in range(len(geometry['time'])):
        one=dict(geometry,full=geometry['full'][:,k:k+1],missing=geometry['missing'][:,k:k+1],time=geometry['time'][k:k+1])
        fields.append(spectral_fields(direct(root,one,p[:,k])))
    return np.concatenate(fields,axis=1)


def validate_emulator(root,geometry,basis,library,base):
    model=CoherentCOModel(AXES,library,geometry['time'],1)
    control_noise=base['noise']/130.8285010091114
    checks=[];rng=np.random.default_rng(63121)
    for i in range(10):
        p=rng.uniform(-.95,.95,(2,3));exact=spectral_fields(direct(root,geometry,p))
        for blocking in [0.,1.]:
            prediction=model.evaluate(p.ravel(),blocking,subtract_reference=False)[0]@basis.T
            truth=np.sum(exact[:,:,0]+blocking*exact[:,:,1],axis=0)
            checks.append(dict(test=i,blocking=blocking,error_control_noise_norm=float(np.linalg.norm(prediction-truth)/control_noise),
                max_error_ppm=float(np.max(abs(prediction-truth))*1e6)))
    return checks


def fit_pair(model,data,noise,rng,n_starts=3):
    starts=[np.zeros(model.n_parameters)]
    starts.extend(rng.uniform(-.65,.65,(n_starts-1,model.n_parameters)))
    null,null_runs=fit_coherent(model,data,noise,starts)
    # Include the best null point in the joint starts so the nested alternative
    # is not accidentally handicapped by unrelated initial conditions.
    joint,joint_runs=fit_coherent(model,data,noise,[null['coefficients'],*starts[1:]],with_moon=True)
    if joint['chi2']>null['chi2']:
        joint=dict(null,blocking=0.)
    return null,joint,null_runs,joint_runs


def serializable_fit(fit):
    return {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in fit.items()}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--build-only',action='store_true')
    parser.add_argument('--n-mc',type=int,default=16);args=parser.parse_args()
    if args.n_mc<1:parser.error('--n-mc must be at least one')
    root=Path(__file__).resolve().parents[1];out=root/'outputs';out.mkdir(exist_ok=True)
    geometry=build_molecular_geometry(nr=52,n_scale_heights=16.,split_limbs=True)
    basis,library,base,key=build_emulator(root,geometry)
    checks=validate_emulator(root,geometry,basis,library,base)
    (root/'docs/co_nonlinear_emulator_checks.json').write_text(json.dumps(dict(cache_key=key,checks=checks),indent=2)+'\n')
    max_error=max(c['error_control_noise_norm'] for c in checks)
    print('Maximum independent emulator error / control noise:',max_error,flush=True)
    if max_error>.1:raise RuntimeError('Emulator error too large for recovery; refine grid before inference')
    np.savez_compressed(out/'co_nonlinear_emulator.npz',basis=basis,library=library,
        time=base['time'],wavenumber=base['wavenumber'],reference=base['reference'],signal=base['signal'],noise=base['noise'])
    if args.build_only:return
    run_experiment(root,geometry,basis,library,base,key,checks,args.n_mc)


def run_experiment(root,geometry,basis,library,base,key,checks,nmc):
    out=root/'outputs';time=base['time'];b0=base['reference'];s=base['signal']
    native_noise=base['noise'];control_noise=native_noise/130.8285010091114
    models={n:CoherentCOModel(AXES,library,time,n) for n in [1,3,5]}
    tangent_controls=[]
    for n,model in models.items():
        _,jac,moon=model.evaluate(np.zeros(model.n_parameters),0.)
        design=jac.reshape(-1,model.n_parameters)
        norms=np.linalg.norm(design,axis=0)
        u,singular,_=np.linalg.svd(design/norms,full_matrices=False)
        u=u[:,singular>singular[0]*1e-8]
        q=moon.ravel()-u@(u.T@moon.ravel())
        tangent_controls.append(dict(n_controls=n,rank=u.shape[1],
            native_local_fisher_snr=float(np.linalg.norm(q)/native_noise),
            control_local_fisher_snr=float(np.linalg.norm(q)/control_noise)))
    rng=np.random.default_rng(63122)
    # Smooth held-out physical atmosphere: each sector has distinct low-order
    # histories, not a time-dependent mixture of two endpoint spectra.
    truth_coeff=np.array([[[-.25,-.10,.10],[-.25,-.10,.15],[-.30,-.10,.25]],
                         [[ .25, .10,-.05],[ .25, .10,-.10],[ .30, .10,-.20]]])
    p=models[3].trajectory(truth_coeff.ravel())
    physical=direct_trajectory(root,geometry,p)
    weather=np.sum(physical[:,:,0],axis=0)-b0
    weather_moon=weather+np.sum(physical[:,:,1],axis=0)
    pulse_p=p.copy();pulse_p[:,:,2]+=.6*np.exp(-.5*((time+1.1)/.35)**2)[None,:]
    pulse_fields=direct_trajectory(root,geometry,pulse_p)
    pulse=np.sum(pulse_fields[:,:,0],axis=0)-b0
    cases={'native_moon':s,'control_moon':s,'control_atmosphere':weather,'control_atmosphere_moon':weather_moon,
           'control_wind_pulse':pulse,'control_wind_pulse_moon':pulse+np.sum(pulse_fields[:,:,1],axis=0)}
    rows=[];fits=[]
    for case,truth in cases.items():
        noise=native_noise if case=='native_moon' else control_noise
        data=truth@basis
        for n,model in models.items():
            null,joint,nr,jr=fit_pair(model,data,noise,rng)
            separation=float(np.sqrt(max(null['chi2']-joint['chi2'],0)))
            rows.append(dict(case=case,n_temporal_controls=n,n_atmospheric_parameters=model.n_parameters,
                noise_ppm=noise*1e6,null_chi2=null['chi2'],joint_chi2=joint['chi2'],
                root_delta_chi2=separation,recovered_blocking=joint['blocking'],
                null_start_chi2_spread=max(r['chi2'] for r in nr)-min(r['chi2'] for r in nr),
                joint_start_chi2_spread=max(r['chi2'] for r in jr)-min(r['chi2'] for r in jr),
                all_starts_success=all(r['success'] for r in nr+jr),
                joint_active_bounds=joint['active_bounds']))
            fits.append(dict(case=case,n_controls=n,null=serializable_fit(null),joint=serializable_fit(joint),
                null_starts=[serializable_fit(r) for r in nr],joint_starts=[serializable_fit(r) for r in jr]))
            print(case,n,'separation',separation,'blocking',joint['blocking'],flush=True)
    mc=[];model=models[3]
    for k in range(nmc):
        noise_draw=rng.normal(size=(len(time),basis.shape[1]))
        for label,truth in [('no_moon',weather),('moon',weather_moon)]:
            data=truth@basis+control_noise*noise_draw
            null,joint,nr,jr=fit_pair(model,data,control_noise,rng,n_starts=2)
            mc.append(dict(draw=k,injection=label,root_delta_chi2=float(np.sqrt(max(null['chi2']-joint['chi2'],0))),
                blocking=joint['blocking'],all_starts_success=all(r['success'] for r in nr+jr),
                start_chi2_spread=max(max(r['chi2'] for r in nr)-min(r['chi2'] for r in nr),
                    max(r['chi2'] for r in jr)-min(r['chi2'] for r in jr))))
        print(f'Nonlinear noise pair {k+1}/{nmc}',flush=True)
    profiles=[]
    for case in ['control_moon','control_atmosphere_moon','control_wind_pulse']:
        previous=np.array(next(f for f in fits if f['case']==case and f['n_controls']==3)['joint']['coefficients'])
        for fraction in np.linspace(0,1,11):
            best,runs=fit_coherent(model,cases[case]@basis,control_noise,
                [previous,np.zeros(model.n_parameters)],fixed_blocking=float(fraction))
            previous=best['coefficients']
            profiles.append(dict(case=case,blocking=float(fraction),chi2=best['chi2'],
                all_starts_success=all(r['success'] for r in runs)))
    # Re-bin a direct 2.5-minute physical simulation to assess the frozen-state
    # five-minute approximation, including its covariance with moon motion.
    fine_geometry=build_molecular_geometry(nr=52,n_scale_heights=16.,split_limbs=True,n_exposures=96)
    fine_p=np.einsum('gck,tk->gtc',truth_coeff,temporal_basis(fine_geometry['time'],3))
    fine=direct_trajectory(root,fine_geometry,fine_p)
    fine_total=np.sum(fine[:,:,0]+fine[:,:,1],axis=0).reshape(48,2,-1).mean(axis=1)
    temporal_error=float(np.linalg.norm(fine_total-(weather_moon+b0))/control_noise)
    direct_fit_checks=[]
    for case,hypothesis in [('control_moon','null'),('control_atmosphere_moon','joint'),('control_wind_pulse','null')]:
        saved=next(f for f in fits if f['case']==case and f['n_controls']==3)[hypothesis]
        coefficients=np.array(saved['coefficients']);fraction=saved['blocking']
        exact=direct_trajectory(root,geometry,model.trajectory(coefficients))
        prediction=np.sum(exact[:,:,0]+fraction*exact[:,:,1],axis=0)-b0
        approximate=model.evaluate(coefficients,fraction)[0]@basis.T
        direct_fit_checks.append(dict(case=case,hypothesis=hypothesis,
            emulator_error_control_noise_norm=float(np.linalg.norm(prediction-approximate)/control_noise),
            direct_chi2=float(np.sum(((prediction-cases[case])/control_noise)**2)),emulated_chi2=saved['chi2'],
            minimum_absolute_absorption=float(np.min(prediction+b0))))
    summary=dict(cache_key=key,seed=63122,n_noise_pairs=nmc,known_orbit=True,
        native_noise_ppm=native_noise*1e6,control_noise_ppm=control_noise*1e6,
        control='physical opaque moon, all Gaussian noise reduced by 130.8285; not instrument forecast',
        bounds=dict(temperature_K=[560,720],CO_mixing_ratio=[1e-4/3,3e-4],wind_km_s=[-1.5,1.5],blocking=[0,1]),
        temporal_controls=[1,3,5],time_window_h=[-2,2],no_independent_depth_scale=True,
        physical_model='same isothermal LTE CO, fixed Rp/P0 and two velocity sectors; not equilibrium/circulation',
        n_grid_nodes=[len(a) for a in AXES],spectral_rank=RANK,
        max_independent_emulator_error_control_noise_norm=max(c['error_control_noise_norm'] for c in checks),
        frozen_atmosphere_exposure_check_control_noise_norm=temporal_error,
        local_fisher_controls=tangent_controls,independent_direct_fit_checks=direct_fit_checks,
        wind_pulse_scope='true nonlinear uniform wind addition 0.9 km/s, Gaussian center -1.1 h, sigma .35 h; outside spline family',
        root_delta_chi2_is_not_calibrated_sigma=True)
    for name,data in [('co_nonlinear_recovery.csv',rows),('co_nonlinear_noise_recovery.csv',mc),('co_nonlinear_profiles.csv',profiles)]:
        with (root/'docs'/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
    (root/'docs/co_nonlinear_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (root/'docs/co_nonlinear_fits.json').write_text(json.dumps(fits,indent=2)+'\n')
    np.savez_compressed(out/'co_nonlinear_truth.npz',time=time,truth_parameters=p,weather=weather,weather_moon=weather_moon)
    plot_results(out,models,basis,s,rows,mc,fits,p,control_noise)
    plot_saved_diagnostics(root)


def plot_results(out,models,basis,signal,rows,mc,fits,truth_parameters,noise):
    fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for case,label,color in [('control_moon','Physical moon, idealized noise','#2575a8'),
                             ('control_atmosphere_moon','Moon + smooth atmosphere','#438343'),
                             ('control_atmosphere','Smooth atmosphere, no moon','#a44d56'),
                             ('control_wind_pulse','Unmodelled wind pulse, no moon','#a8752c')]:
        selected=[r for r in rows if r['case']==case]
        ax[0,0].plot([r['n_temporal_controls'] for r in selected],[r['root_delta_chi2'] for r in selected],
            'o-',label=label,color=color)
    ax[0,0].set(xlabel='Temporal control points per parameter / limb',ylabel='sqrt(delta chi2)',title='Noiseless nested-model separation')
    ax[0,0].legend(fontsize=8)
    chosen=next(f for f in fits if f['case']=='control_moon' and f['n_controls']==3)
    model=models[3];atmosphere=model.evaluate(np.array(chosen['null']['coefficients']),0.)[0]@basis.T
    ax[0,1].plot(model.time,np.linalg.norm(signal,axis=1)/noise,label='Physical lunar residual')
    ax[0,1].plot(model.time,np.linalg.norm(signal-atmosphere,axis=1)/noise,label='Remaining after bounded smooth fit')
    ax[0,1].set(xlabel='Time [h]',ylabel='Spectral norm / control noise',title='Atmosphere-only attempt to imitate the moon')
    ax[0,1].legend(fontsize=8)
    for label,color in [('no_moon','#a44d56'),('moon','#2575a8')]:
        v=[r['root_delta_chi2'] for r in mc if r['injection']==label]
        ax[1,0].hist(v,bins=np.linspace(0,max(7,max(v)+.1),13),histtype='step',label=label,color=color)
    ax[1,0].set(xlabel='sqrt(delta chi2), not calibrated sigma',ylabel='Noise realizations',title='Bounded nonlinear recovery: paired noise')
    ax[1,0].legend(fontsize=8)
    recovered=np.array(next(f for f in fits if f['case']=='control_atmosphere_moon' and f['n_controls']==3)['joint']['coefficients'])
    trajectory=model.trajectory(recovered)
    for g,color in enumerate(['#2575a8','#d27d28']):
        ax[1,1].plot(model.time,physical_parameters(truth_parameters)[0][g],color=color,label=f'Limb {g+1}: injected')
        ax[1,1].plot(model.time,physical_parameters(trajectory)[0][g],ls='--',color=color,label=f'Limb {g+1}: recovered')
    ax[1,1].set(xlabel='Time [h]',ylabel='Temperature [K]',title='Held-out smooth atmospheric injection')
    ax[1,1].legend(fontsize=8)
    for ext in ['png','pdf']:fig.savefig(out/f'co_nonlinear_recovery.{ext}',dpi=180)
    plt.close(fig)


def plot_saved_diagnostics(root):
    out=root/'outputs'
    with (root/'docs/co_nonlinear_profiles.csv').open() as f:profiles=list(csv.DictReader(f))
    with (root/'docs/co_nonlinear_recovery.csv').open() as f:rows=list(csv.DictReader(f))
    fits=json.loads((root/'docs/co_nonlinear_fits.json').read_text())
    saved=np.load(out/'co_nonlinear_emulator.npz');truth=np.load(out/'co_nonlinear_truth.npz')
    model=CoherentCOModel(AXES,saved['library'],saved['time'],3)
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    for case,label,color in [('control_moon','Moon only','#2575a8'),
        ('control_atmosphere_moon','Smooth atmosphere + moon','#438343'),
        ('control_wind_pulse','Wind pulse, no moon','#a8752c')]:
        selected=[p for p in profiles if p['case']==case]
        chi=np.array([float(p['chi2']) for p in selected])
        ax[0].plot([float(p['blocking']) for p in selected],chi-chi.min(),'o-',ms=3,label=label,color=color)
    ax[0].set(xlabel='Lunar blocking fraction f',ylabel='Profile chi2 - minimum chi2',title='Known-orbit noiseless profiles')
    ax[0].legend(fontsize=8)
    cases=['control_moon','control_atmosphere','control_atmosphere_moon','control_wind_pulse','control_wind_pulse_moon']
    selected=[next(r for r in rows if r['case']==c and int(r['n_temporal_controls'])==3) for c in cases]
    x=np.arange(len(cases))
    ax[1].bar(x-.18,[np.sqrt(float(r['null_chi2'])) for r in selected],.36,label='Atmosphere only',color='#777777')
    ax[1].bar(x+.18,[np.sqrt(float(r['joint_chi2'])) for r in selected],.36,label='Atmosphere + moon',color='#2575a8')
    ax[1].set_xticks(x,['Moon\nonly','Atmosphere\nonly','Atmosphere\n+ moon','Wind pulse\nonly','Wind pulse\n+ moon'],fontsize=8)
    ax[1].set(ylabel='Noiseless residual norm / control noise',title='Does either hypothesis actually fit?')
    ax[1].legend(fontsize=8)
    for ext in ['png','pdf']:fig.savefig(out/f'co_nonlinear_profiles.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(3,2,figsize=(11,8),constrained_layout=True)
    labels=['Temperature [K]','CO mixing ratio / 1e-4','LOS wind [km/s]']
    limits=[(560,720),(1/3,3),(-1.5,1.5)]
    for col,(case,hypothesis) in enumerate([('control_moon','null'),('control_atmosphere_moon','joint')]):
        fit=next(f for f in fits if f['case']==case and f['n_controls']==3)[hypothesis]
        values=physical_parameters(model.trajectory(np.array(fit['coefficients'])))
        for i in range(3):
            for g,color in enumerate(['#2575a8','#d27d28']):
                v=values[i][g]/(1e-4 if i==1 else 1.)
                axes[i,col].plot(model.time,v,color=color,label=f'Sector {g+1}: fitted')
                if col==1:
                    true=physical_parameters(truth['truth_parameters'])[i][g]/(1e-4 if i==1 else 1.)
                    axes[i,col].plot(model.time,true,'o',markevery=6,ms=3,color=color,label=f'Sector {g+1}: injected')
            axes[i,col].set(xlabel='Time [h]',ylabel=labels[i],ylim=limits[i])
    axes[0,0].set_title('Atmosphere-only attempt to imitate the moon')
    axes[0,1].set_title('Noiseless atmosphere + moon recovery')
    axes[0,0].legend(fontsize=8);axes[0,1].legend(fontsize=8)
    for ext in ['png','pdf']:fig.savefig(out/f'co_nonlinear_trajectories.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
