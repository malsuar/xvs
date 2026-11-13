"""Joint numerical refinement of the CO atmospheric nuisance experiment."""
import json
from pathlib import Path
import numpy as np
from atmospheric_nuisance_recovery import (MODES,STEPS,atmosphere_parameters,
    atmospheric_synthesis,build_molecular_geometry,make_subspaces,project_nuisance)


def main():
    root=Path(__file__).resolve().parents[1]
    saved=np.load(root/'outputs/co_atmospheric_nuisance.npz')
    geometry=build_molecular_geometry(na=960,nr=80,nt=2401,shift_step=.125,
        n_scale_heights=12.,split_limbs=True)
    base=atmospheric_synthesis(root,geometry=geometry,dv=.1,wing_sigma=200.)
    b=base['reference'];s=base['signal'];d=[];c=[]
    for mode,step in zip(MODES,STEPS/2):
        plus=atmospheric_synthesis(root,geometry=geometry,dv=.1,wing_sigma=200.,
            **atmosphere_parameters(mode,step))['reference']
        minus=atmospheric_synthesis(root,geometry=geometry,dv=.1,wing_sigma=200.,
            **atmosphere_parameters(mode,-step))['reference']
        d.append((plus-minus)/(2*step));c.append((plus+minus-2*b)/step**2)
        print('Refined nuisance mode:',mode,flush=True)
    spaces=make_subspaces(b,d,c);rows=[]
    for i,name in enumerate(saved['models']):
        q=project_nuisance(s,spaces[name][0]);old=saved['components'][i]
        rows.append(dict(model=str(name),fine_native_snr=float(np.linalg.norm(q)/base['noise']),
            relative_snr_change=float(abs(np.linalg.norm(old)/np.linalg.norm(q)-1)),
            relative_component_difference=float(np.linalg.norm(old-q)/np.linalg.norm(q)),
            fine_rank_min=int(spaces[name][1].min()),fine_rank_max=int(spaces[name][1].max())))
    result=dict(fine_mesh=dict(n_azimuths=960,n_radii=80,n_times=2401,
            spectral_step_km_s=.1,velocity_channel_step_km_s=.125,wing_sigma=200.),
        signal_map_relative_difference=float(np.linalg.norm(saved['signal']-s)/np.linalg.norm(s)),
        models=rows,all_snr_changes_below_1percent=bool(all(r['relative_snr_change']<.01 for r in rows)))
    (root/'docs/co_atmospheric_convergence.json').write_text(json.dumps(result,indent=2)+'\n')
    plot_convergence(root/'outputs',result)
    print(json.dumps(result,indent=2),flush=True)


def plot_convergence(out,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels=['Depth','+ T','+ CO','+ wind','Uniform','Two limbs','Curvatures']
    fig,ax=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for a,key,title in zip(ax,['relative_snr_change','relative_component_difference'],
                           ['Change in retained lunar S/N','Change in retained spectral map']):
        values=[100*r[key] for r in result['models']]
        a.bar(labels,values,color='#348c99')
        a.set(ylabel='Difference from joint refinement [%]',title=title)
        a.tick_params(axis='x',rotation=35)
        for i,value in enumerate(values):a.text(i,value+.002,f'{value:.3f}',ha='center',fontsize=8)
        a.set_ylim(0,max(values)*1.25)
    for ext in ['png','pdf']:fig.savefig(out/f'co_atmospheric_convergence.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
