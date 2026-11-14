"""Time-grid convergence at fixed exposure cadence in fractional models."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_grazing_convergence import evaluate


def main():
    root=Path(__file__).resolve().parents[1];(root/'outputs').mkdir(exist_ok=True)
    rows=[]; phases=[0.,.25,.5,.75]
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for col,(veq,sigma) in enumerate([(12.,1.5),(25.,.5),(25.,.2)]):
        for b,label in [(0.,'central'),(1.057,'grazing')]:
            case=[]
            for nt in [601,1201,2401,4801]:
                row,_,_,_=evaluate(480,16,b,nt=nt,moon_radius=.03,fractional=True,order=16,
                    n_velocities=801,velocity_limit=40.,v_eq=veq,line_sigma=sigma,all_exposure_phases=True)
                row.update(v_eq=veq,line_sigma=sigma,time_step_seconds=4*3600/(nt-1))
                case.append(row)
            reference=case[-1]
            for row in case:
                row['worst_snr_error']=max(abs(row[f'snr_phase_{p}']/reference[f'snr_phase_{p}']-1) for p in phases)
            rows.extend(case)
            axes[0,col].loglog([r['time_step_seconds'] for r in case[:-1]],
                [max(r['worst_snr_error'],1e-12) for r in case[:-1]],'o-',label=label)
            axes[1,col].plot(phases,[100*(reference[f'snr_phase_{p}']/reference['snr_phase_0.0']-1) for p in phases],'o-',label=label)
            print('CASE',veq,sigma,label,[(r['time_step_seconds'],r['worst_snr_error']) for r in case],flush=True)
        axes[0,col].axhline(.01,color='0.5',ls=':',label='1% tolerance')
        axes[0,col].set(title=f'veq={veq:g}, sigma={sigma:g} km/s',xlabel='Internal time step [s]',ylabel='Worst same-phase S/N error')
        axes[1,col].set(xlabel='Exposure start / duration',ylabel='S/N change from phase zero [%]')
        for ax in axes[:,col]:ax.legend(fontsize=8)
    with (root/'docs/temporal_convergence.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/temporal_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
