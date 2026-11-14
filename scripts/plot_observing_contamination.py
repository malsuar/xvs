"""Figures for count budgets, correlated contaminants and held-out search tests."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from observing_contamination_pilot import templates,contamination,FAMILIES
from nonlinear_co_recovery import AXES
from xds.moving_shadow import MovingShadow,REFERENCE_PERIOD
from xds.molecular import load_co_lines,C
from xds.orbit_search import empirical_tail


def read(path):
    with path.open() as handle:return list(csv.DictReader(handle))


def save(fig,root,name):
    for extension in ['png','pdf']:fig.savefig(root/f'outputs/{name}.{extension}',dpi=180)
    plt.close(fig)


def budget_figure(root,tag):
    old=np.load(root/'outputs/co_nonlinear_emulator.npz');tpl=templates(root,old)
    rows=read(root/'docs/co_observing_budgets.csv')
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    wave=tpl['wavelength'];selection=(wave>=2.318)&(wave<=2.326)
    ax=axes[0,0]
    ax.plot(wave[selection],tpl['disk'][selection],color='#336e9a',lw=1,label='PHOENIX 5700 K')
    ax.plot(wave[selection],tpl['transmission'][selection],color='#b96535',lw=1,label='SkyCalc Paranal')
    ax.set(xlabel='Vacuum wavelength [um]',ylabel='Normalized flux / transmission',title='Physical line templates, R=100000')
    ax.legend(fontsize=8)
    ax=axes[0,1]
    centers=np.arange(3);values=np.array([float(r['median_sigma_ppm']) for r in rows])
    lower=np.array([float(r['q16_sigma_ppm']) for r in rows]);upper=np.array([float(r['q84_sigma_ppm']) for r in rows])
    for i,r in enumerate(rows):
        valid=r['valid_linearity']=='True'
        ax.errorbar(i,values[i],yerr=[[values[i]-lower[i]],[upper[i]-values[i]]],fmt='o' if valid else 'x',color='#336e9a' if valid else '#ae495b',capsize=4)
    ax.axhline(11.394654,color='0.4',ls=':',label='Previous laboratory median')
    ax.set_xticks(centers,['8 m, K=6','39 m, K=6','39 m, K=3\nnonlinear: invalid'])
    ax.set(yscale='log',ylabel='Single-exposure uncertainty [ppm]',title='Counts, background, read noise and duty cycle')
    ax.legend(fontsize=8)
    time=old['time'];combined,stellar,telluric=contamination(tpl,time,np.random.default_rng(63151),'combined')
    ax=axes[1,0]
    keep=(tpl['transmission']>=.3)&(tpl['disk']>=.2)
    for data,color,label in [(stellar,'#336e9a','Stellar residual'),(telluric,'#b96535','Telluric residual')]:
        ax.plot(time,np.max(abs(data[:,keep]),axis=1)*1e6,color=color,label=label)
    ax.set(yscale='log',xlabel='Time [h]',ylabel='Maximum retained-bin residual [ppm]',title='One assumed residual realization, not a forecast')
    ax.legend(fontsize=8)
    ax=axes[1,1]
    covariance=np.eye(48)+np.ones((48,48))/12
    correlation=covariance/np.sqrt(np.diag(covariance)[:,None]*np.diag(covariance)[None,:])
    im=ax.imshow(correlation,origin='lower',extent=[-2,2,-2,2],cmap='viridis',vmin=0,vmax=1,aspect='auto')
    fig.colorbar(im,ax=ax,label='Noise correlation')
    ax.set(xlabel='Time [h]',ylabel='Time [h]',title='Finite shared reference: off-diagonal correlation 1/13')
    save(fig,root,f'co_observing_{tag}_budget')


def inference_figures(root,tag):
    rates=read(root/f'docs/co_observing_{tag}_rates.csv')
    cal=read(root/f'docs/co_observing_{tag}_calibration.csv')
    colors={'atmosphere':'#336e9a','stellar':'#82569d','telluric':'#b96535','combined':'#438343'}
    level_names={'control':'Laboratory noise','39m_K6':'39 m, K=6 count budget'}
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for i,level in enumerate(['control','39m_K6']):
        for ax,method in zip(axes[i],['unprotected','protected']):
            maximum=0.
            for family in FAMILIES:
                values=np.array([float(r['statistic']) for r in cal if r['level']==level and r['method']==method and r['family']==family])
                maximum=max(maximum,float(np.sqrt(values.max())))
                grid=np.r_[0.,np.geomspace(max(values.max()*1e-5,1e-9),max(values.max()*1.05,1e-8),300)]
                tail=[empirical_tail(values,t) for t in grid]
                ax.step(np.sqrt(grid),tail,where='post',label=family,color=colors[family])
            ax.axhline(.05,color='0.4',ls=':')
            ax.set(xscale='log' if maximum>40 else 'linear',xlim=(.1 if maximum>40 else 0,maximum*1.06),
                yscale='log',ylim=(.01,1.1),xlabel='sqrt(maximum bounded score), not sigma',
                ylabel='Empirical tail probability',title=f'{level_names[level]}: {method}')
            ax.legend(fontsize=8)
    save(fig,root,f'co_observing_{tag}_null_tails')
    fig,axes=plt.subplots(2,2,figsize=(12,7.5),constrained_layout=True)
    settings=[('unprotected','atmosphere','#b96535','x','Unprotected, atmospheric null'),
        ('unprotected','robust','#82569d','o','Unprotected, all nulls'),
        ('protected','atmosphere','#336e9a','x','Protected, atmospheric null'),
        ('protected','robust','#438343','s','Protected, all nulls')]
    for i,level in enumerate(['control','39m_K6']):
        for ax,stream,labels in [(axes[i,0],'validation',FAMILIES),
            (axes[i,1],'injection',['clean_fixed','combined_fixed','combined_random'])]:
            for j,(method,calibration,color,marker,label) in enumerate(settings):
                selected=[next(r for r in rates if r['level']==level and r['method']==method and r['calibration']==calibration and r['stream']==stream and r['family']==f) for f in labels]
                values=np.array([float(r['fraction']) for r in selected]);lo=np.array([float(r['ci95_low']) for r in selected]);hi=np.array([float(r['ci95_high']) for r in selected])
                ax.errorbar(np.arange(len(labels))+(j-1.5)*.13,values,yerr=[values-lo,hi-values],fmt=marker,color=color,ms=4,capsize=2,label=label)
            names=['Atmosphere','Stellar','Telluric','Combined'] if stream=='validation' else ['No contamination\nfixed phase','Combined\nfixed phase','Combined\nrandom phase']
            ax.set_xticks(np.arange(len(labels)),names,fontsize=8)
            ax.axhline(.05,color='0.5',ls=':')
            ax.set(ylim=(-.025,1.025),ylabel='Candidate fraction',title=f'{level_names[level]}: '+('false alarms' if stream=='validation' else 'injection flags'))
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=8,loc='outside upper center',ncol=2)
    save(fig,root,f'co_observing_{tag}_validation')
    decisions=read(root/f'docs/co_observing_{tag}_decisions.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
    for ax,level in zip(axes,['control','39m_K6']):
        for method,color,marker in [('unprotected','#82569d','o'),('protected','#438343','s')]:
            selected=[r for r in decisions if r['level']==level and r['method']==method and r['stream']=='injection' and r['family']=='combined_fixed']
            faces=[color if r['flag_robust']=='True' else 'none' for r in selected]
            ax.scatter([float(r['period'])/21.238292071543107 for r in selected],
                [float(r['phase_error']) for r in selected],facecolors=faces,edgecolors=color,marker=marker,label=method,s=30)
        ax.axvline(1,color='0.5',ls=':')
        ax.set(xscale='log',xlim=(.34,4.2),ylim=(-.01,.51),xlabel='Maximum-score period / injected period',
            ylabel='Circular phase error [cycles]',title=level_names[level]+': combined contamination')
        ax.set_xticks([.5,1,2,4],['0.5','1','2','4']);ax.legend(fontsize=8)
    fig.suptitle('Filled symbols: all-null candidate flag, not proof of the correct orbit',fontsize=11)
    save(fig,root,f'co_observing_{tag}_orbital_errors')


def line_maps(root,tag):
    old=np.load(root/'outputs/co_nonlinear_emulator.npz');moving=np.load(root/'outputs/co_blind_shadow.npz')
    shadow=MovingShadow(AXES,moving['library'],moving['left'],moving['centers'],moving['rotation'],float(moving['rp']))
    tpl=templates(root,old);time=old['time']
    moon=shadow.signal(np.zeros((2,len(time),3)),REFERENCE_PERIOD,.7473,time,order=16)@old['basis'].T
    combined,_,_=contamination(tpl,time,np.random.default_rng(63151),'combined')
    keep=(tpl['transmission']>=.3)&(tpl['disk']>=.2)
    # This illustrates spectral nuisance removal BEFORE optional compression
    # to the fixed planetary sufficient subspace. No data-selected line is used.
    from xds.observing_contamination import CountBudget,electron_budget
    budget=CountBudget(area_m2=978.,live_seconds=84*1.427,reads_per_exposure=84,read_noise_e=12.)
    sigma=electron_budget(tpl['wavelength'],tpl['stellar_relative'],tpl['transmission'],tpl['radiance'],budget)['sigma']
    u,s,_=np.linalg.svd(tpl['nuisance'][:,keep].T/sigma[keep,None],full_matrices=False)
    u=u[:,s>s[0]*1e-9]
    white=combined[:,keep]/sigma[keep]
    protected=np.full_like(combined,np.nan)
    protected[:,keep]=(white-(white@u)@u.T)*sigma[keep]
    lines=load_co_lines(root/'data/exomol/CO_Li2015',640.)
    eligible=(lines.wavenumber>old['wavenumber'][0]*np.exp(60/C))&(lines.wavenumber<old['wavenumber'][-1]*np.exp(-60/C))
    rest=lines.wavenumber[np.argmax(np.where(eligible,lines.strength,-1))]
    velocity=tpl['x']-C*np.log(rest/4300);selected=abs(velocity)<40
    fig,axes=plt.subplots(1,3,figsize=(12,4),sharey=True,constrained_layout=True)
    for ax,data,title in zip(axes,[moon,combined,protected],
        ['Physical moon only','Stellar + telluric residual','After fixed nuisance removal']):
        patch=data[:,selected]*1e6;maximum=np.nanmax(abs(patch))
        im=ax.imshow(patch,origin='lower',aspect='auto',extent=[velocity[selected][0],velocity[selected][-1],-2,2],
            cmap='RdBu_r',vmin=-maximum,vmax=maximum)
        ax.set(xlabel='Model velocity [km/s]',title=title)
        fig.colorbar(im,ax=ax,label='Residual [ppm]',fraction=.046)
    axes[0].set_ylabel('Time [h]')
    fig.suptitle(f'Noiseless illustration at CO rest wavelength {1e4/rest:.5f} um; independent color scales',fontsize=11)
    save(fig,root,f'co_observing_{tag}_line_maps')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot');parser.add_argument('--budget-only',action='store_true')
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'legend.fontsize':8})
    budget_figure(root,args.tag)
    if not args.budget_only:
        inference_figures(root,args.tag)
        line_maps(root,args.tag)
