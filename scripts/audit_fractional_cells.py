"""Compare point and fractional blocking against analytic area and fine S/N."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_grazing_convergence import evaluate


def main():
    root=Path(__file__).resolve().parents[1];(root/'outputs').mkdir(exist_ok=True)
    rows=[];fig,axes=plt.subplots(2,4,figsize=(16,7),constrained_layout=True)
    for col,(rm,b) in enumerate([(.06,0.),(.06,1.084),(.03,0.),(.03,1.057)]):
        ref,_,_,_=evaluate(1920,64,b,moon_radius=rm,fractional=True,order=32)
        case=[]
        for na,nr in [(240,8),(480,16),(960,32)]:
            for method in ['point','fractional']:
                row,_,_,_=evaluate(na,nr,b,moon_radius=rm,fractional=method=='fractional',order=16)
                row.update(moon_radius=rm,method=method,order=16 if method=='fractional' else 0,
                           snr_error=abs(row['snr']/ref['snr']-1))
                case.append(row)
        for method in ['point','fractional']:
            subset=[r for r in case if r['method']==method]
            x=[r['n_azimuths']*r['n_radii'] for r in subset]
            axes[0,col].loglog(x,[max(r['area_relative_error'],1e-14) for r in subset],'o-',label=method)
            axes[1,col].loglog(x,[max(r['snr_error'],1e-14) for r in subset],'o-',label=method)
        # Separate quadrature and representative-point controls on practical grid.
        for radial in [.25,.5,.75]:
            row,_,_,_=evaluate(480,16,b,moon_radius=rm,fractional=True,order=32,radial_sample_fraction=radial)
            row.update(moon_radius=rm,method=f'fractional_radial_{radial}',order=32,
                       snr_error=abs(row['snr']/ref['snr']-1));case.append(row)
        ref.update(moon_radius=rm,method='reference',order=32,snr_error=0.)
        rows.extend(case+[ref])
        axes[0,col].set(title=f'Rm={rm}, b={b}',xlabel='Spatial cells',ylabel='Integrated area relative error')
        axes[1,col].set(xlabel='Spatial cells',ylabel='S/N error vs fine fractional grid')
        for ax in axes[:,col]:ax.legend(fontsize=8)
        print('CASE',rm,b,'reference',ref,'rows',case,flush=True)
    with (root/'docs/fractional_cell_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/fractional_cell_comparison.{ext}',dpi=180)


if __name__=='__main__':main()
