"""Joint refinement and production/audit orchestration consistency."""
import csv
from pathlib import Path
from dataclasses import replace
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from audit_grazing_convergence import evaluate
from xds import PlanetConfig, MoonConfig, run_xds_model
from xds.config import AtmosphereConfig, ObservationConfig, InstrumentConfig
from xds.reference import reference_simulation
from xds.detectability import detectability_report, atmospheric_line_depth
from xds.observing import observe_residual


def main():
    root=Path(__file__).resolve().parents[1];(root/'outputs').mkdir(exist_ok=True)
    atmosphere=AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,temperature=500.,stellar_radius_solar=1.)
    observation=ObservationConfig(exposure_time=5/60,continuum_snr_per_resolution_element=300.,n_lines=200,systematic_floor=10e-6)
    instrument=InstrumentConfig();phases=[0.,.25,.5,.75];rows=[]
    cases=[(12.,1.5,0.),(12.,1.5,1.057),(25.,.2,0.),(25.,.2,1.057)]
    fig,axes=plt.subplots(2,4,figsize=(16,7),constrained_layout=True)
    for col,(veq,sigma,b) in enumerate(cases):
        planet=PlanetConfig(v_eq=veq,line_sigma=sigma)
        moon=MoonConfig(radius=.03,semi_major_axis=3.,period=21.238292071543107,
                        inclination=np.rad2deg(np.arcsin(b/3)),longitude_node=90.,phase=.25)
        base=reference_simulation(planet,4.,observation,instrument)
        fine=replace(base,n_azimuths=base.n_azimuths*2,n_radii=base.n_radii*2,
                     n_times=2*(base.n_times-1)+1,n_velocities=2*(base.n_velocities-1)+1,
                     blocking_quadrature_order=32)
        results=[]
        for stage,sim in [('base',base),('refined',fine)]:
            prod=run_xds_model(planet,moon,sim)
            audit,t,frac,exact,amap=evaluate(sim.n_azimuths,sim.n_radii,b,nt=sim.n_times,
                fractional=True,order=sim.blocking_quadrature_order,moon_radius=.03,
                n_velocities=sim.n_velocities,velocity_limit=sim.velocity_max,v_eq=veq,line_sigma=sigma,return_result=True)
            assert np.array_equal(t,prod.time) and np.array_equal(amap.velocity,prod.velocity)
            map_error=float(np.max(abs(prod.delta_absorption-amap.delta_absorption)))
            fraction_error=float(np.max(abs(prod.combined_fraction-frac)))
            assert map_error<1e-12 and fraction_error<1e-12
            assert np.allclose(prod.absorption,prod.reference_absorption+prod.delta_absorption,atol=1e-14)
            scores=[]
            for phase in phases:
                obs=replace(observation,exposure_start_offset=phase*observation.exposure_time)
                report=detectability_report(prod,atmosphere,obs,instrument)
                sampled=observe_residual(amap,atmospheric_line_depth(atmosphere),obs,instrument)
                expected=np.linalg.norm(sampled.signal/sampled.noise)*np.sqrt(obs.n_lines)
                snr_diff=abs(report.matched_filter_snr/expected-1)
                assert snr_diff<1e-10
                scores.append(report.matched_filter_snr)
                rows.append(dict(v_eq=veq,line_sigma=sigma,impact_rp=b,stage=stage,
                    n_azimuths=sim.n_azimuths,n_radii=sim.n_radii,n_times=sim.n_times,
                    n_velocities=sim.n_velocities,order=sim.blocking_quadrature_order,
                    exposure_phase=phase,snr=report.matched_filter_snr,
                    production_audit_map_max_error=map_error,production_audit_fraction_max_error=fraction_error,
                    production_audit_snr_relative_error=snr_diff,analytic_area_error=audit['area_relative_error']))
            results.append(np.array(scores))
            axes[0,col].plot(prod.time,prod.combined_fraction,label=stage,lw=1)
        change=abs(results[0]/results[1]-1)
        axes[1,col].plot(phases,100*change,'o-')
        axes[1,col].axhline(1,color='0.5',ls=':',label='1% tolerance')
        axes[0,col].set(title=f'veq={veq:g}, sigma={sigma:g}, b={b:g}',xlabel='Time [h]',ylabel='Blocked fraction')
        axes[1,col].set(xlabel='Exposure start / duration',ylabel='Joint S/N refinement change [%]')
        for ax in axes[:,col]:ax.legend(fontsize=8)
        print('CASE',veq,sigma,b,'joint_max',max(change),'map_error',map_error,flush=True)
    with (root/'docs/joint_production_convergence.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/joint_production_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
