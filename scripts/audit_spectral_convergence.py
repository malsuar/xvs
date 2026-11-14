"""Spectral and representative-cell-velocity convergence stress controls."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
from audit_grazing_convergence import evaluate


def main():
    root=Path(__file__).resolve().parents[1]; (root/'outputs').mkdir(exist_ok=True)
    rows=[]; phases=[0.,.25,.5,.75]
    cases=[(12.,1.5),(25.,1.5),(25.,.5),(25.,.2)]
    fig,axes=plt.subplots(2,4,figsize=(16,7),constrained_layout=True)
    for col,(veq,sigma) in enumerate(cases):
        for b,label in [(0.,'central'),(1.057,'grazing')]:
            refs=[]
            for nv in [101,201,401,801,1601]:
                row,_,_,_=evaluate(480,16,b,moon_radius=.03,fractional=True,order=16,
                    n_velocities=nv,velocity_limit=40.,v_eq=veq,line_sigma=sigma,
                    all_exposure_phases=True)
                row.update(v_eq=veq,line_sigma=sigma,n_velocities=nv,
                           velocity_limit=40.,kind='spectral')
                refs.append(row)
            ref=refs[-1]
            for row in refs:
                row['snr_error']=max(abs(row[f'snr_phase_{p}']/ref[f'snr_phase_{p}']-1) for p in phases)
            fine,_,_,_=evaluate(960,32,b,moon_radius=.03,fractional=True,order=32,
                n_velocities=1601,velocity_limit=40.,v_eq=veq,line_sigma=sigma,
                all_exposure_phases=True)
            fine.update(v_eq=veq,line_sigma=sigma,n_velocities=1601,velocity_limit=40.,kind='spatial',
                        snr_error=max(abs(fine[f'snr_phase_{p}']/ref[f'snr_phase_{p}']-1) for p in phases))
            rows.extend(refs+[fine])
            axes[0,col].loglog([80/(r['n_velocities']-1) for r in refs[:-1]],
                              [max(r['snr_error'],1e-12) for r in refs[:-1]],'o-',label=label)
            axes[1,col].bar(label,100*fine['snr_error'])
            print('CASE',veq,sigma,label,'spectral_errors',[(r['n_velocities'],r['snr_error']) for r in refs],
                  'spatial_change',fine['snr_error'],flush=True)
        axes[0,col].set(title=f'veq={veq:g}, sigma={sigma:g} km/s',xlabel='Input velocity step [km/s]',ylabel='Worst same-phase S/N error')
        axes[0,col].legend(fontsize=8)
        axes[0,col].set_xticks([.1,.2,.4,.8],labels=['0.1','0.2','0.4','0.8'])
        axes[0,col].xaxis.set_minor_formatter(NullFormatter())
        axes[0,col].axhline(.01,color='0.5',ls=':',lw=1)
        axes[1,col].set(ylabel='S/N change doubling spatial grid [%]')
    with (root/'docs/spectral_convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/spectral_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
