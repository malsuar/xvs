"""Compare stellar continuum weighting with the archived uniform pilot."""
import csv
from dataclasses import replace
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from xds import PlanetConfig,MoonConfig,run_xds_model
from xds.config import AtmosphereConfig,ObservationConfig,InstrumentConfig,StellarSurfaceConfig
from xds.reference import reference_simulation
from xds.detectability import scale_height_km,JUPITER_RADIUS_KM,detectability_report
from xds.dynamics import moon_period_hours


def main():
    root=Path(__file__).resolve().parents[1]
    with (root/'docs/pilot_phase_convergence.csv').open() as f:pilot=list(csv.DictReader(f))
    atm=AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,temperature=640.)
    h=5*scale_height_km(atm)/(1.1*JUPITER_RADIUS_KM)
    planet=PlanetConfig(v_eq=12.,line_sigma=1.5,atmosphere_half_width=h)
    obs=ObservationConfig(exposure_time=5/60,continuum_snr_per_resolution_element=300.,n_lines=200,systematic_floor=10e-6)
    inst=InstrumentConfig();base=reference_simulation(planet,4.,obs,inst)
    fine=replace(base,n_azimuths=960,n_radii=32,n_times=2*(base.n_times-1)+1,
                 n_velocities=2*(base.n_velocities-1)+1,blocking_quadrature_order=32)
    star=StellarSurfaceConfig(695700/(1.1*JUPITER_RADIUS_KM),
        (30/365.25)**(2/3)*149597870.7/(1.1*JUPITER_RADIUS_KM),720.,u1=.4,u2=.25)
    rows=[];fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for col,rm in enumerate([.03,.06,.09]):
        for kind,b in [('central',0.),('grazing',1+h+.9*rm)]:
            for phase in [.25,.75]:
                moon=MoonConfig(radius=rm,semi_major_axis=3.,period=float(moon_period_hours(3,.7,1.1)),
                    inclination=np.rad2deg(np.arcsin(b/3)),longitude_node=90.,phase=phase)
                results={stage:run_xds_model(planet,moon,sim,star) for stage,sim in [('base',base),('refined',fine)]}
                candidates=[r for r in pilot if float(r['moon_radius_rp'])==rm and r['encounter']==kind and float(r['lunar_phase'])==phase]
                ratios=[]
                for cadence in [0.,.25,.5,.75]:
                    o=replace(obs,exposure_start_offset=cadence*obs.exposure_time)
                    reports={stage:detectability_report(result,atm,o,inst) for stage,result in results.items()}
                    uniform=float(next(r for r in candidates if float(r['exposure_phase'])==cadence)['refined_snr'])
                    new=reports['refined'].matched_filter_snr
                    ratios.append(new/uniform)
                    rows.append(dict(radius_rp=rm,encounter=kind,lunar_phase=phase,exposure_phase=cadence,
                        uniform_snr=uniform,limb_darkened_snr=new,ratio=new/uniform,
                        joint_relative_change=abs(reports['base'].matched_filter_snr/new-1),
                        limb_darkened_peak_ppm=reports['refined'].peak_xds_depth*1e6))
                if kind=='central' and phase==.75:
                    r=results['refined']
                    axes[0,col].plot(r.time,np.trapz(r.reference_absorption,r.velocity,axis=1)/(planet.line_depth*planet.line_sigma*np.sqrt(2*np.pi)),label='Atmospheric reference')
                    axes[0,col].axhline(1,color='0.5',ls=':',label='Uniform-star reference')
                axes[1,col].plot([0.,.25,.5,.75],ratios,'o-' if phase==.25 else 's--',label=f'{kind}, lunar phase {phase}')
                print('CASE',rm,kind,phase,'ratio',min(ratios),max(ratios),flush=True)
                del results
        axes[0,col].set(title=f'Rm={rm:.2f} Rp',xlabel='Time [h]',ylabel='Reference equivalent-width ratio')
        axes[1,col].set(xlabel='Exposure start / duration',ylabel='Limb-darkened / uniform S/N')
        for ax in axes[:,col]:ax.legend(fontsize=7)
    with (root/'docs/stellar_surface_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/stellar_surface_comparison.{ext}',dpi=180)


if __name__=='__main__':main()
