"""Radial sampling origins and smaller-moon convergence controls."""
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
    cases=[(.06,0.),(.06,1.03+.9*.06),(.03,0.),(.03,1.03+.9*.03)]
    fig,axes=plt.subplots(3,4,figsize=(16,10),constrained_layout=True)
    for col,(rm,b) in enumerate(cases):
        case=[]
        for na,nr in [(480,16),(960,32),(1920,64)]:
            for radial in [.25,.5,.75]:
                row,t,f,exact=evaluate(na,nr,b,moon_radius=rm,
                    radial_sample_fraction=radial,all_exposure_phases=True)
                row.update(moon_radius=rm,radial_sample_fraction=radial,kind='spatial')
                case.append(row)
                if na==1920:
                    axes[0,col].plot(t,f,label=f'Radial position {radial}',lw=1)
        axes[0,col].plot(t,exact,'k--',label='Analytic area')
        fine=[r for r in case if r['n_azimuths']==1920]
        refs={p:np.mean([r[f'snr_phase_{p}'] for r in fine]) for p in phases}
        for row in case:
            row['worst_snr_error']=max(abs(row[f'snr_phase_{p}']/refs[p]-1) for p in phases)
        for radial in [.25,.5,.75]:
            subset=[r for r in case if r['radial_sample_fraction']==radial]
            x=[r['n_azimuths']*r['n_radii'] for r in subset]
            axes[1,col].loglog(x,[max(r['area_relative_error'],1e-12) for r in subset],'o-',label=f'{radial}')
            axes[2,col].loglog(x,[max(r['worst_snr_error'],1e-12) for r in subset],'o-',label=f'{radial}')
        temporal,_,_,_=evaluate(1920,64,b,nt=2401,moon_radius=rm,all_exposure_phases=True)
        mid=next(r for r in fine if r['radial_sample_fraction']==.5)
        temporal.update(moon_radius=rm,radial_sample_fraction=.5,kind='temporal',
            worst_snr_error=max(abs(temporal[f'snr_phase_{p}']/mid[f'snr_phase_{p}']-1) for p in phases))
        rows.extend(case+[temporal])
        axes[0,col].set(title=f'Rm={rm:.2f} Rp, b={b:.3f} Rp',xlabel='Time [h]',ylabel='Blocked fraction')
        for k,label in [(1,'Integrated area relative error'),(2,'S/N error vs fine phase ensemble')]:
            axes[k,col].set(xlabel='Spatial cells',ylabel=label);axes[k,col].legend(title='Radial fraction',fontsize=8)
        axes[0,col].legend(fontsize=7)
        print('CASE',rm,b,'max_SNR_error_by_grid',[(na,max(r['worst_snr_error'] for r in case if r['n_azimuths']==na)) for na in [480,960,1920]],
              'max_area_error_by_grid',[(na,max(r['area_relative_error'] for r in case if r['n_azimuths']==na)) for na in [480,960,1920]],
              'time_error',temporal['worst_snr_error'],flush=True)
    with (root/'docs/small_moon_radial_convergence.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/small_moon_radial_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
