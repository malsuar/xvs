"""Deterministic physical pilot; no heuristic preselection or detection claims."""
import csv
import json
from pathlib import Path
from dataclasses import replace
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from xds import PlanetConfig,MoonConfig,run_xds_model
from xds.config import AtmosphereConfig,ObservationConfig,InstrumentConfig
from xds.reference import reference_simulation
from xds.detectability import detectability_report,scale_height_km,JUPITER_RADIUS_KM
from xds.dynamics import (moon_period_hours,fluid_roche_limit_rp,planet_density_g_cm3,
                         stable_outer_moon_limit_rp,breakup_rotation_period_hours,
                         transit_duration_hours)


def main():
    root=Path(__file__).resolve().parents[1];(root/'outputs').mkdir(exist_ok=True)
    atmosphere=AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,
        temperature=640.,stellar_radius_solar=1.,mean_molecular_weight=2.3,effective_scale_heights=5.)
    h=5*scale_height_km(atmosphere)/(1.1*JUPITER_RADIUS_KM)
    planet=PlanetConfig(v_eq=12.,line_sigma=1.5,atmosphere_half_width=h)
    observation=ObservationConfig(exposure_time=5/60,continuum_snr_per_resolution_element=300.,
                                  n_lines=200,systematic_floor=10e-6)
    instrument=InstrumentConfig()
    pm=float(moon_period_hours(3.,.7,1.1))
    roche=float(fluid_roche_limit_rp(planet_density_g_cm3(.7,1.1),5.5))
    outer=float(stable_outer_moon_limit_rp(30.,.7,planet_radius_jupiter=1.1))
    prot=2*np.pi*1.1*JUPITER_RADIUS_KM/12/3600
    assert roche<3<outer and prot>breakup_rotation_period_hours(.7,1.1)
    stellar_axis_au=(30/365.25)**(2/3)
    teq=5772*np.sqrt((695700/149597870.7)/(2*stellar_axis_au))
    metadata=dict(planet_mass_mj=.7,planet_radius_rj=1.1,stellar_mass_msun=1.,stellar_radius_rsun=1.,
        stellar_temperature_k=5772.,planet_period_days=30.,approx_zero_albedo_teq_k=float(teq),
        adopted_temperature_k=640.,moon_density_g_cm3=5.5,annulus_thickness_rp=h,
        moon_period_h=pm,roche_rp=roche,outer_stability_rp=outer,rotation_period_h=prot,
        approximate_central_transit_h=float(transit_duration_hours(30.,1.)),
        observing_window_h=4.,lunar_phases=[0.,.25,.5,.75],exposure_phases=[0.,.25,.5,.75])
    (root/'docs/pilot_assumptions.json').write_text(json.dumps(metadata,indent=2)+'\n')
    base=reference_simulation(planet,4.,observation,instrument)
    fine=replace(base,n_azimuths=960,n_radii=32,n_times=2*(base.n_times-1)+1,
                 n_velocities=2*(base.n_velocities-1)+1,blocking_quadrature_order=32)
    rows=[];fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for col,rm in enumerate([.03,.06,.09]):
        for kind,b,color in [('central',0.,'C0'),('grazing',1+h+.9*rm,'C1')]:
            for phase in [0.,.25,.5,.75]:
                moon=MoonConfig(radius=rm,semi_major_axis=3.,period=pm,phase=phase,
                    inclination=np.rad2deg(np.arcsin(b/3)),longitude_node=90.)
                scores={}
                for stage,sim in [('base',base),('refined',fine)]:
                    result=run_xds_model(planet,moon,sim)
                    scores[stage]=[]
                    for cadence in [0.,.25,.5,.75]:
                        obs=replace(observation,exposure_start_offset=cadence*observation.exposure_time)
                        report=detectability_report(result,atmosphere,obs,instrument)
                        scores[stage].append(report.matched_filter_snr)
                errors=[]
                for j,cadence in enumerate([0.,.25,.5,.75]):
                    x,y=scores['base'][j],scores['refined'][j]
                    no_signal=x==0 and y==0
                    err=0. if no_signal else (abs(x/y-1) if y>0 else float('inf'))
                    errors.append(err)
                    rows.append(dict(moon_radius_rp=rm,encounter=kind,minimum_projected_separation_rp=b,
                        lunar_phase=phase,exposure_phase=cadence,base_snr=x,refined_snr=y,
                        relative_change=err,no_signal=no_signal,passes_1percent=err<.01))
                refined=np.array(scores['refined'])
                axes[0,col].errorbar(phase+(0.015 if kind=='grazing' else -.015),np.mean(refined),
                    yerr=[[np.mean(refined)-min(refined)],[max(refined)-np.mean(refined)]],
                    fmt='o',color=color,label=kind if phase==0 else None,capsize=3)
                axes[1,col].plot(phase,max(errors)*100,'o',color=color)
                print('PILOT',rm,kind,phase,'SNR_range',min(refined),max(refined),
                      'worst_change',max(errors),flush=True)
        axes[0,col].set(title=f'Rm={rm:.2f} Rp',xlabel='Lunar phase [cycles]',ylabel='Refined S/N: cadence range')
        axes[0,col].legend(fontsize=8)
        axes[1,col].axhline(1,color='0.5',ls=':',label='1% tolerance')
        axes[1,col].set(xlabel='Lunar phase [cycles]',ylabel='Worst joint refinement change [%]')
        axes[1,col].legend(fontsize=8)
    with (root/'docs/pilot_phase_convergence.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/pilot_phase_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
