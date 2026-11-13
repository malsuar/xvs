"""Combined spatial grids, azimuthal origins and exposure alignment."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_grazing_convergence import evaluate


def main():
    root=Path(__file__).resolve().parents[1]
    rows=[]
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for col,b in enumerate([0.,1.08,1.14]):
        case=[]
        for na,nr in [(480,16),(960,32),(1920,64)]:
            for shift in [0.,.25,.5]:
                row,_,_,_=evaluate(na,nr,b,cell_phase=shift,all_exposure_phases=True)
                row['cell_phase']=shift
                case.append(row)
        # Compare the same observing alignment to a fine-grid phase ensemble.
        phases=[0.,.25,.5,.75]
        fine=[r for r in case if r['n_azimuths']==1920]
        refs={p:np.mean([r[f'snr_phase_{p}'] for r in fine]) for p in phases}
        for row in case:
            row['worst_snr_error']=max(abs(row[f'snr_phase_{p}']/refs[p]-1) for p in phases)
            values=[row[f'snr_phase_{p}'] for p in phases]
            row['exposure_snr_range_fraction']=(max(values)-min(values))/np.mean(values)
        rows.extend(case)
        for shift in [0.,.25,.5]:
            subset=[r for r in case if r['cell_phase']==shift]
            axes[0,col].loglog([r['n_azimuths']*r['n_radii'] for r in subset],
                               [max(r['worst_snr_error'],1e-12) for r in subset],'o-',label=f'Cell shift {shift}')
        for row in fine:
            axes[1,col].plot(phases,[row[f'snr_phase_{p}'] for p in phases],'o-',label=f"Cell shift {row['cell_phase']}")
        axes[0,col].set(title=f'b={b:.2f} Rp',xlabel='Spatial cells',ylabel='Worst S/N error across exposures')
        axes[1,col].set(xlabel='Exposure start / exposure duration',ylabel='Corrected S/N (fine grid)')
        axes[0,col].legend(fontsize=8);axes[1,col].legend(fontsize=8)
        print('CASE',b,[(r['n_azimuths'],r['cell_phase'],r['worst_snr_error'],r['exposure_snr_range_fraction']) for r in case],flush=True)
    with (root/'docs/combined_resolution.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/combined_resolution_audit.{ext}',dpi=180)


if __name__=='__main__':main()
