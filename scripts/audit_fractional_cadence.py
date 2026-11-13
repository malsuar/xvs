"""Angular origins and exposure phases with fractional-cell integration."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_grazing_convergence import evaluate


def main():
    root=Path(__file__).resolve().parents[1]
    (root/'outputs').mkdir(exist_ok=True)
    rows=[];phases=[0.,.25,.5,.75]
    fig,axes=plt.subplots(3,4,figsize=(16,10),constrained_layout=True)
    cases=[(.06,0.),(.06,1.084),(.03,0.),(.03,1.057)]
    for col,(rm,b) in enumerate(cases):
        reference,_,_,_=evaluate(1920,64,b,moon_radius=rm,fractional=True,order=32,
                                all_exposure_phases=True)
        case=[]
        for na,nr in [(240,8),(480,16),(960,32)]:
            for shift in [0.,.25,.5]:
                row,_,_,_=evaluate(na,nr,b,moon_radius=rm,fractional=True,order=16,
                                   cell_phase=shift,all_exposure_phases=True)
                row.update(moon_radius=rm,angular_phase=shift,kind='spatial',
                           worst_snr_error=max(abs(row[f'snr_phase_{p}']/reference[f'snr_phase_{p}']-1) for p in phases))
                case.append(row)
        for shift in [0.,.25,.5]:
            subset=[r for r in case if r['angular_phase']==shift]
            x=[r['n_azimuths']*r['n_radii'] for r in subset]
            axes[0,col].loglog(x,[max(r['worst_snr_error'],1e-14) for r in subset],'o-',label=f'Angular shift {shift}')
            axes[1,col].loglog(x,[max(r['area_relative_error'],1e-14) for r in subset],'o-',label=f'{shift}')
            fine=subset[-1]
            axes[2,col].plot(phases,[100*(fine[f'snr_phase_{p}']/reference['snr_phase_0.0']-1) for p in phases],'o-',label=f'{shift}')
        reference.update(moon_radius=rm,angular_phase=0.,kind='reference',worst_snr_error=0.)
        temporal,_,_,_=evaluate(960,32,b,nt=2401,moon_radius=rm,fractional=True,order=16,
                                all_exposure_phases=True)
        coarse=next(r for r in case if r['n_azimuths']==960 and r['angular_phase']==0)
        temporal.update(moon_radius=rm,angular_phase=0.,kind='temporal',
                        worst_snr_error=max(abs(temporal[f'snr_phase_{p}']/coarse[f'snr_phase_{p}']-1) for p in phases))
        rows.extend(case+[reference,temporal])
        for k,label in [(0,'S/N error vs same-phase reference'),(1,'Integrated area relative error')]:
            axes[k,col].set(xlabel='Spatial cells',ylabel=label)
            axes[k,col].legend(fontsize=7)
        axes[0,col].set_title(f'Rm={rm:.2f} Rp, b={b:.3f} Rp')
        axes[2,col].set(xlabel='Exposure start / duration',ylabel='S/N change from reference phase zero [%]')
        axes[2,col].legend(fontsize=7)
        spread=(max(reference[f'snr_phase_{p}'] for p in phases)-min(reference[f'snr_phase_{p}'] for p in phases))/np.mean([reference[f'snr_phase_{p}'] for p in phases])
        print('CASE',rm,b,'worst_errors',[(na,max(r['worst_snr_error'] for r in case if r['n_azimuths']==na)) for na in [240,480,960]],
              'cadence_range',spread,'time_change',temporal['worst_snr_error'],flush=True)
    with (root/'docs/fractional_cadence_convergence.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/fractional_cadence_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
