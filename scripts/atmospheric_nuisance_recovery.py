"""Known-orbit CO recovery against uniform and two-limb atmospheric alternatives."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from molecular_injection_recovery import build_molecular_geometry,synthesize
from xds.molecular import nuisance_basis,project_nuisance


MODES=['temperature','log_abundance','wind','limb_temperature','limb_log_abundance','limb_wind']
STEPS=np.array([10.,.1,.2,10.,.1,.2])


def atmosphere_parameters(mode,amount):
    t=np.full(2,640.);lx=np.zeros(2);v=np.zeros(2)
    contrast=np.array([-1.,1.]) if mode.startswith('limb_') else np.ones(2)
    if mode.endswith('temperature'):t+=amount*contrast
    elif mode.endswith('log_abundance'):lx+=amount*contrast
    elif mode.endswith('wind'):v+=amount*contrast
    else:raise ValueError(mode)
    return dict(temperature=t,co_mixing_ratio=1e-4*np.exp(lx),wind_km_s=v)


def make_subspaces(b,derivatives,curvatures):
    definitions={
        'depth':[b],
        'depth_T':[b,derivatives[0]],
        'depth_X':[b,derivatives[1]],
        'depth_wind':[b,derivatives[2]],
        'uniform_linear':[b,*derivatives[:3]],
        'two_limb_linear':[b,*derivatives],
        'two_limb_curvature':[b,*derivatives,*curvatures],
    }
    return {name:nuisance_basis(modes) for name,modes in definitions.items()}


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def atmospheric_synthesis(root,**kwargs):
    # Smooth sub-grid shifts are necessary for meaningful wind curvatures.
    return synthesize(root,shift_interpolation='cubic',**kwargs)


def main():
    root=Path(__file__).resolve().parents[1];out=root/'outputs';out.mkdir(exist_ok=True)
    geometry=build_molecular_geometry(nr=40,n_scale_heights=12.,split_limbs=True)
    base=atmospheric_synthesis(root,geometry=geometry)
    linear_interpolation_base=synthesize(root,geometry=geometry)
    b,s,noise=base['reference'],base['signal'],base['noise']
    derivatives=[];curvatures=[];half_derivatives=[];half_curvatures=[];step_checks=[]
    for mode,step in zip(MODES,STEPS):
        approximations=[]
        for factor in [1.,.5]:
            plus=atmospheric_synthesis(root,geometry=geometry,**atmosphere_parameters(mode,step*factor))['reference']
            minus=atmospheric_synthesis(root,geometry=geometry,**atmosphere_parameters(mode,-step*factor))['reference']
            approximations.append(((plus-minus)/(2*step*factor),(plus+minus-2*b)/(step*factor)**2))
        (d,c),(hd,hc)=approximations
        derivatives.append(d);curvatures.append(c);half_derivatives.append(hd);half_curvatures.append(hc)
        step_checks.append(dict(mode=mode,step=float(step),
            derivative_relative_change=float(np.linalg.norm(d-hd)/np.linalg.norm(hd)),
            curvature_relative_change=float(np.linalg.norm(c-hc)/np.linalg.norm(hc))))
        print('Finished atmospheric mode:',mode,flush=True)
    # Use the smaller steps for reported inference; larger steps are an audit.
    spaces=make_subspaces(b,half_derivatives,half_curvatures)
    coarse_spaces=make_subspaces(b,derivatives,curvatures)
    names=list(spaces)
    components={name:project_nuisance(s,basis) for name,(basis,_) in spaces.items()}
    snr={name:float(np.linalg.norm(q)/noise) for name,q in components.items()}
    lab_scale=5/snr['depth']
    units=np.array([(components[name]/np.linalg.norm(components[name])).ravel() for name in names])
    rng=np.random.default_rng(63118);nmc=2000;scores=[]
    for start in range(0,nmc,32):
        n=rng.normal(size=(min(32,nmc-start),s.size))
        scores.append(n@units.T)
    scores=np.concatenate(scores)
    recovery=[]
    example=lab_scale*s+noise*rng.normal(size=s.shape)
    for i,name in enumerate(names):
        coarse=project_nuisance(s,coarse_spaces[name][0]);q=components[name]
        example_projected=project_nuisance(example,spaces[name][0])
        example_z=float(np.sum(example*q)/(noise*np.linalg.norm(q)))
        fitted_amplitude=example_z/snr[name]
        joint_residual=project_nuisance(example-fitted_amplitude*s,spaces[name][0])
        delta_chi2=float((np.sum(example_projected**2)-np.sum(joint_residual**2))/noise**2)
        assert np.isclose(delta_chi2,example_z**2,atol=1e-7)
        recovery.append(dict(model=name,nuisance_rank_min=int(spaces[name][1].min()),
            nuisance_rank_max=int(spaces[name][1].max()),native_snr=snr[name],
            retained_norm_vs_depth=snr[name]/snr['depth'],
            controlled_expected_z=lab_scale*snr[name],
            controlled_mean_z=float((scores[:,i]+lab_scale*snr[name]).mean()),
            controlled_fraction_z_above3=float(np.mean(scores[:,i]+lab_scale*snr[name]>3)),
            null_mean_z=float(scores[:,i].mean()),null_std_z=float(scores[:,i].std(ddof=1)),
            controlled_example_z=example_z,controlled_example_delta_chi2=delta_chi2,
            controlled_example_fitted_multiplier=fitted_amplitude,
            multiplier_for_expected_z5=5/snr[name],
            half_step_snr_relative_change=float(abs(np.linalg.norm(coarse)/np.linalg.norm(q)-1)),
            half_step_component_relative_difference=float(np.linalg.norm(coarse-q)/np.linalg.norm(q))))

    alternatives={
        'uniform_T_plus55K':dict(temperature=695.),
        'uniform_X_times2p3':dict(co_mixing_ratio=2.3e-4),
        'uniform_wind_plus1p3':dict(wind_km_s=1.3),
        'limb_T_contrast35K':atmosphere_parameters('limb_temperature',35.),
        'limb_logX_contrast0p45':atmosphere_parameters('limb_log_abundance',.45),
        'limb_wind_contrast0p7':atmosphere_parameters('limb_wind',.7),
        'uniform_combined':dict(temperature=675.,co_mixing_ratio=1.6e-4,wind_km_s=.6),
        'limb_combined':dict(temperature=np.array([615.,665.]),
            co_mixing_ratio=1e-4*np.exp(np.array([-.3,.3])),wind_km_s=np.array([-.4,.4])),
    }
    # An explicit stress test: variable covering fraction of a second spectrum.
    # This is not identical to a globally uniform T(t), X(t), or wind(t).
    envelope=np.exp(-.5*((base['time']+1.1)/.35)**2)
    finite=[];alternative_maps={}
    for label,parameters in alternatives.items():
        alternative=atmospheric_synthesis(root,geometry=geometry,**parameters)['reference']-b
        alternative_maps[label]=envelope[:,None]*alternative
        for evolution,perturbation in [('static',alternative),('transient_mixture',alternative_maps[label])]:
            for i,name in enumerate(names):
                q=components[name];norm=np.linalg.norm(q)
                projected=project_nuisance(perturbation,spaces[name][0])
                z=float(np.sum(perturbation*q)/(noise*norm))
                finite.append(dict(alternative=label,evolution=evolution,recovery_model=name,
                    atmospheric_residual_snr=float(np.linalg.norm(perturbation)/noise),
                    postfit_residual_snr=float(np.linalg.norm(projected)/noise),
                    expected_false_lunar_z=z,
                    expected_controlled_moon_plus_atmosphere_z=lab_scale*snr[name]+z,
                    signed_overlap_after_nuisance=float(np.sum(projected*q)/(max(np.linalg.norm(projected),1e-30)*norm))))
        print('Finished finite alternative:',label,flush=True)

    # Atmospheric columns have a common fixed geometric boundary; test the
    # hottest alternative and a strong abundance case on a taller domain.
    taller=build_molecular_geometry(nr=48,n_scale_heights=14.,split_limbs=True)
    taller_base=atmospheric_synthesis(root,geometry=taller)
    extent=[]
    for label in ['uniform_T_plus55K','uniform_X_times2p3']:
        high=atmospheric_synthesis(root,geometry=taller,**alternatives[label])['reference']-taller_base['reference']
        low=alternative_maps[label]/envelope[:,None]
        extent.append(dict(alternative=label,
            relative_map_difference=float(np.linalg.norm(low-high)/np.linalg.norm(high))))

    summary=dict(seed=63118,n_noise_realizations=nmc,n_exposures=s.shape[0],n_bins=s.shape[1],
        noise_per_bin_ppm=noise*1e6,native_full_snr=float(np.linalg.norm(s)/noise),
        reference_height_scale_temperature_k=640.,domain_scale_heights=12.,n_radial_cells=40,
        per_exposure_free_nuisance=True,known_lunar_trajectory=True,
        controlled_moon_multiplier=lab_scale,svd_relative_tolerance=1e-8,
        doppler_interpolation='cubic spline; previous pilot retains linear default',
        interpolation_residual_map_relative_change=float(np.linalg.norm(s-linear_interpolation_base['signal'])/np.linalg.norm(s)),
        perturbations=step_checks,outer_boundary_checks=extent,
        conventions='positive wind adds to model velocity; limb order is v_rot<0 then v_rot>=0',
        curvature_scope='separate single-mode second derivatives; no cross derivatives or nonlinear retrieval',
        finite_transient_scope='two-column spectral mixture, Gaussian covering fraction centered -1.1 h, sigma .35 h')
    write_csv(root/'docs/co_atmospheric_recovery.csv',recovery)
    write_csv(root/'docs/co_atmospheric_false_positives.csv',finite)
    (root/'docs/co_atmospheric_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    np.savez_compressed(out/'co_atmospheric_nuisance.npz',time=base['time'],wavenumber=base['wavenumber'],
        signal=s,reference=b,noise=noise,models=np.array(names),
        components=np.array([components[n] for n in names]),null_scores=scores,
        finite_alternatives=np.array(list(alternative_maps)),finite_maps=np.array(list(alternative_maps.values())))
    plot_results(out,base,recovery,finite,components,scores,alternative_maps,lab_scale)
    print(json.dumps(recovery,indent=2),flush=True)


def plot_results(out,base,recovery,finite,components,scores,alternatives,lab_scale):
    labels=['Depth','Depth + T','Depth + CO','Depth + wind','Uniform: all','Two limbs: linear','Two limbs: curvature']
    colors=['#555555','#2575a8','#d27d28','#438343','#9c485c','#348c99','#80579a']
    plt.rcParams.update({'font.size':11,'axes.titlesize':12})
    fig,ax=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
    snr=np.array([r['native_snr'] for r in recovery])
    ax[0,0].barh(labels,snr,color=colors);ax[0,0].invert_yaxis()
    ax[0,0].set(xlabel='Native known-trajectory S/N',title='Lunar information after atmospheric fitting')
    for i,value in enumerate(snr):ax[0,0].text(value+.0003,i,f'{value:.4f}',va='center',fontsize=9)
    ax[0,0].set_xlim(0,snr.max()*1.2)
    for name,color in zip(['depth','uniform_linear','two_limb_linear','two_limb_curvature'],[colors[0],*colors[4:]]):
        q=components[name]
        ax[0,1].plot(base['time'],np.linalg.norm(q,axis=1)/base['noise'],label=name.replace('_',' '),color=color)
    ax[0,1].set(xlabel='Time [h]',ylabel='Residual norm / noise per exposure',title='Native lunar residual')
    ax[0,1].legend(fontsize=8)
    for name,color in zip(['depth','uniform_linear','two_limb_linear','two_limb_curvature'],[colors[0],*colors[4:]]):
        i=list(components).index(name)
        z=scores[:,i]+lab_scale*snr[i]
        ax[1,0].hist(z,bins=np.linspace(-4,9,60),histtype='step',density=True,label=name.replace('_',' '),color=color)
    ax[1,0].axvline(3,color='0.5',ls=':',lw=1)
    ax[1,0].set(xlabel='Recovered signed lunar score Z',ylabel='Density',title='Same amplified moon, increasingly flexible atmosphere')
    ax[1,0].legend(fontsize=8)
    selected=[r for r in finite if r['evolution']=='transient_mixture']
    modes=list(components);cases=list(alternatives)
    matrix=np.array([[next(r['signed_overlap_after_nuisance'] for r in selected if r['alternative']==case and r['recovery_model']==mode)
                      for mode in modes] for case in cases])
    vmax=max(abs(matrix).max(),1e-6)
    im=ax[1,1].imshow(matrix,cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto')
    case_labels=['T +55 K','CO x2.3','Wind +1.3 km/s','Limb T +/-35 K','Limb ln(CO) +/-0.45',
                 'Limb wind +/-0.7 km/s','Uniform combined','Limb combined']
    ax[1,1].set_yticks(np.arange(len(cases)),case_labels,fontsize=8)
    ax[1,1].set_xticks(np.arange(len(modes)),['Depth','+T','+CO','+wind','Uniform','Limbs','Curvature'],rotation=40,ha='right',fontsize=8)
    for (i,j),value in np.ndenumerate(matrix):
        ax[1,1].text(j,i,f'{value:.2f}',ha='center',va='center',fontsize=7,color='white' if abs(value)>.6*vmax else 'black')
    ax[1,1].set_title('No moon: overlap after atmospheric fitting')
    fig.colorbar(im,ax=ax[1,1],label='Signed overlap with lunar template',shrink=.85)
    for ext in ['png','pdf']:fig.savefig(out/f'co_atmospheric_degeneracy.{ext}',dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    k=int(np.argmax(np.linalg.norm(base['signal'],axis=1)))
    band=(base['wavenumber']>4302.8)&(base['wavenumber']<4304.5)
    x=base['wavenumber'][band];signal=base['signal'][k,band]*1e6
    ax[0].plot(x,signal,color='black',lw=2,label='Native lunar residual')
    for name,color in zip(['depth','uniform_linear','two_limb_linear','two_limb_curvature'],[colors[0],*colors[4:]]):
        remaining=components[name][k,band]*1e6
        label=name.replace('_',' ')
        ax[0].plot(x,signal-remaining,color=color,label=label)
        ax[1].plot(x,remaining,color=color,label=label)
    for a in ax:
        a.axhline(0,color='0.7',lw=.7)
        a.set(xlabel='Wavenumber [cm^-1]',ylabel='Residual [ppm]')
    ax[0].set_title(f'Atmospheric approximation at t={base["time"][k]:.2f} h')
    ax[1].set_title('Lunar component remaining after fitting')
    ax[0].legend(fontsize=8);ax[1].legend(fontsize=8)
    for ext in ['png','pdf']:fig.savefig(out/f'co_atmospheric_line_shapes.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
