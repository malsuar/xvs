"""Paper-style diagnostics for the calibrated finite phase/period search."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from xds.orbit_search import empirical_tail
from xds.moving_shadow import REFERENCE_PERIOD,period_limits


COLORS={'reference':'#2575a8','smooth':'#438343','wind_pulse':'#d27d28'}


def read_csv(path):
    with path.open() as f:return list(csv.DictReader(f))


def calibration_plot(root,tag,records,examples=None):
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'legend.fontsize':8})
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for ax,n in zip(axes,[3,5]):
        groups={family:np.array([r['statistic'] for r in records if r['stream']=='calibration' and r['family']==family and r['n_controls']==n])
                for family in COLORS}
        if not all(len(a) for a in groups.values()):continue
        for family,values in groups.items():
            order=np.sort(values);tail=(1+np.arange(len(order),0,-1))/(len(order)+1)
            ax.step(np.sqrt(order),tail,where='post',color=COLORS[family],label=family.replace('_',' '))
        fixed=np.sort([r['oracle_statistic'] for r in records if r['stream']=='calibration' and r['family']=='reference' and r['n_controls']==n])
        ax.step(np.sqrt(fixed),(1+np.arange(len(fixed),0,-1))/(len(fixed)+1),where='post',ls='--',color=COLORS['reference'],label='reference: fixed orbit')
        grid=np.linspace(0,max(a.max() for a in groups.values())*1.15,400)
        envelope=[max(empirical_tail(a,t) for a in groups.values()) for t in grid]
        ax.step(np.sqrt(grid),envelope,where='post',color='black',lw=1.4,label='worst represented tail')
        ax.axhline(.05,color='0.5',ls=':',lw=1)
        if examples:
            value=next(float(r['statistic']) for r in examples if r['case']=='direct_offgrid' and int(r['n_controls'])==n)
            ax.axvline(np.sqrt(value),color='#9c485c',ls='-.',lw=1,label='direct off-grid injection')
        ax.set(yscale='log',ylim=(.006,1.1),xlabel='sqrt(maximum bounded score), not sigma',
            ylabel='Empirical tail probability',title=f'{n} temporal controls: search penalty and atmospheric nulls')
        ax.legend(loc='upper right',fontsize=7)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/co_blind_{tag}_calibration.{ext}',dpi=180)
    plt.close(fig)


def final_plots(root,tag,records):
    rates=read_csv(root/f'docs/co_blind_{tag}_rates.csv');examples=read_csv(root/f'docs/co_blind_{tag}_examples.csv')
    decisions=read_csv(root/f'docs/co_blind_{tag}_decisions.csv')
    calibration_plot(root,tag,records,examples)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),constrained_layout=True)
    settings=[(3,'gaussian','#d27d28','x','Gaussian only, 3'),(5,'gaussian','#a44d56','x','Gaussian only, 5'),
              (3,'robust','#2575a8','o','Atmospheric calibration, 3'),(5,'robust','#438343','s','Atmospheric calibration, 5'),
              (0,'consensus','#222222','D','Consensus')]
    labels=['reference','smooth','wind_pulse']
    for j,(n,rule,color,marker,title) in enumerate(settings):
        selected=[next(r for r in rates if r['stream']=='validation' and r['family']==f and int(r['n_controls'])==n and r['rule']==rule) for f in labels]
        values=np.array([float(r['fraction']) for r in selected]);low=np.array([float(r['ci95_low']) for r in selected]);high=np.array([float(r['ci95_high']) for r in selected])
        axes[0].errorbar(np.arange(3)+(j-2)*.12,values,yerr=[values-low,high-values],fmt=marker,color=color,ms=4,capsize=2,label=title)
    axes[0].axhline(.05,color='0.5',ls=':',lw=1)
    axes[0].set_xticks(np.arange(3),['Reference','Smooth atmosphere','Wind pulse'])
    axes[0].set(ylabel='Held-out candidate fraction',ylim=(-.025,1.025),title='False alarms: independent no-moon validation')
    axes[0].legend(fontsize=7,loc='upper left')
    labels=['showcase_reference','showcase_smooth','showcase_wind_pulse','uniform_phase']
    for j,(n,rule,color,marker,title) in enumerate(settings[2:]):
        selected=[next(r for r in rates if r['stream']=='injection' and r['family']==f and int(r['n_controls'])==n and r['rule']==rule) for f in labels]
        values=np.array([float(r['fraction']) for r in selected]);low=np.array([float(r['ci95_low']) for r in selected]);high=np.array([float(r['ci95_high']) for r in selected])
        axes[1].errorbar(np.arange(4)+(j-1)*.15,values,yerr=[values-low,high-values],fmt=marker,color=color,ms=4,capsize=2,label=title)
    axes[1].set_xticks(np.arange(4),['Reference\nfixed phase','Smooth\nfixed phase','Wind pulse\nfixed phase','Smooth\nrandom phase'],fontsize=8)
    axes[1].set(ylabel='Candidate fraction',ylim=(-.025,1.025),title='Injected moons: flags are not correct-orbit probabilities')
    axes[1].legend(fontsize=7,loc='upper left')
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/co_blind_{tag}_validation.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    images=[];folder=root/'outputs/co_blind_runs'/tag
    maxstat=max(float(r['statistic']) for r in examples if r['case']=='direct_offgrid')
    for ax,n in zip(axes[0],[3,5]):
        data=np.load(folder/f'grid_direct_offgrid_{n}.npz');values=np.sqrt(data['statistic'])
        image=ax.pcolormesh(data['phases'],data['periods'],values,cmap='magma',shading='nearest',vmin=0,vmax=np.sqrt(maxstat))
        row=next(r for r in examples if r['case']=='direct_offgrid' and int(r['n_controls'])==n)
        ax.scatter([.7473],[REFERENCE_PERIOD],marker='x',color='#62d78c',s=65,label='Injected')
        ax.scatter([float(row['phase'])],[float(row['period'])],facecolors='none',edgecolors='cyan',s=70,label='Score maximum')
        ax.set(yscale='log',xlabel='Orbital phase at mid-transit [cycles]',ylabel='Period [h]',title=f'Direct off-grid injection: {n} controls')
        ax.set_yticks([8,12,20,40,80],['8','12','20','40','80'])
        ax.legend(fontsize=8,loc='upper left');images.append(image)
    fig.colorbar(images[0],ax=list(axes[0]),label='sqrt(score), not posterior probability',shrink=.9)
    bins=np.geomspace(*period_limits(),17)
    ax=axes[1,0]
    for n,color in [(3,'#2575a8'),(5,'#438343')]:
        values=[float(r['period']) for r in decisions if r['stream']=='injection' and r['family']=='showcase_reference' and int(r['n_controls'])==n]
        ax.hist(values,bins=bins,histtype='step',color=color,label=f'{n} controls')
    ax.axvline(REFERENCE_PERIOD,color='black',ls=':',lw=1,label='Injected period')
    ax.set(xscale='log',xlabel='Period of maximum score [h]',ylabel='Realizations',title='Fixed off-grid phase: all noisy injections')
    ax.set_xticks([8,12,20,40,80],['8','12','20','40','80'])
    ax.legend(fontsize=8)
    ax=axes[1,1]
    for n,color,marker in [(3,'#2575a8','o'),(5,'#438343','s')]:
        selected=[r for r in decisions if r['stream']=='injection' and r['family']=='showcase_reference' and int(r['n_controls'])==n]
        ax.scatter([float(r['period']) for r in selected],[float(r['phase']) for r in selected],
            color=color,marker=marker,s=20,label=f'{n} controls',alpha=.65)
    ax.scatter([REFERENCE_PERIOD],[.7473],marker='x',s=60,color='black',label='Injected orbit')
    ax.set(xlabel='Period at maximum score [h]',ylabel='Phase at maximum score [cycles]',title='Noisy joint orbital maxima, not confidence regions')
    ax.legend(fontsize=8)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/co_blind_{tag}_orbits.{ext}',dpi=180)
    plt.close(fig)
    from matplotlib.lines import Line2D
    fig,axes=plt.subplots(2,2,figsize=(10,7.5),constrained_layout=True)
    labels=['showcase_reference','showcase_smooth','showcase_wind_pulse','uniform_phase']
    titles=['Reference atmosphere, fixed phase','Smooth atmosphere, fixed phase','Wind pulse + moon, fixed phase','Random phase, including invisible moons']
    for ax,label,title in zip(axes.ravel(),labels,titles):
        for n,color,marker in [(3,'#2575a8','o'),(5,'#438343','s')]:
            selected=[r for r in decisions if r['stream']=='injection' and r['family']==label and int(r['n_controls'])==n]
            faces=[color if r['flag_consensus']=='True' else 'none' for r in selected]
            ax.scatter([float(r['period'])/REFERENCE_PERIOD for r in selected],
                [float(r['phase_circular_error']) for r in selected],facecolors=faces,edgecolors=color,marker=marker,s=22,alpha=.8)
        ax.axvline(1,color='0.6',ls=':',lw=.8)
        ax.set(xscale='log',xlim=(period_limits()[0]/REFERENCE_PERIOD*.95,period_limits()[1]/REFERENCE_PERIOD*1.05),
            ylim=(-.012,.512),xlabel='Recovered period / injected period',ylabel='Circular phase error [cycles]',title=title)
        ax.set_xticks([.5,1,2,4],['0.5','1','2','4'])
    handles=[Line2D([],[],ls='',marker='o',color='#2575a8',label='3 controls'),
        Line2D([],[],ls='',marker='s',color='#438343',label='5 controls'),
        Line2D([],[],ls='',marker='o',color='black',label='Consensus flag'),
        Line2D([],[],ls='',marker='o',color='black',markerfacecolor='none',label='Not flagged')]
    fig.legend(handles=handles,loc='outside upper center',ncol=4,fontsize=8)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/co_blind_{tag}_orbital_errors.{ext}',dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot');args=parser.parse_args()
    root=Path(__file__).resolve().parents[1];raw=root/'outputs/co_blind_runs'/args.tag/'trials.jsonl'
    records=[json.loads(line) for line in raw.read_text().splitlines() if line]
    if (root/f'docs/co_blind_{args.tag}_rates.csv').exists():final_plots(root,args.tag,records)
    else:calibration_plot(root,args.tag,records)


if __name__=='__main__':main()
