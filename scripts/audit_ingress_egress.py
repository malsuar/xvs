"""Joint stellar/lunar cell integration around selected limb encounters."""
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
from xds.penumbra import circle_overlap
from xds.stellar_cells import stellar_cell_weights


def main():
    root=Path(__file__).resolve().parents[1];(root/'outputs').mkdir(exist_ok=True)
    atm=AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,temperature=640.)
    h=5*scale_height_km(atm)/(1.1*JUPITER_RADIUS_KM)
    planet=PlanetConfig(v_eq=12,line_sigma=1.5,atmosphere_half_width=h)
    obs=ObservationConfig(exposure_time=5/60,continuum_snr_per_resolution_element=300,n_lines=200,systematic_floor=10e-6)
    inst=InstrumentConfig()
    star=StellarSurfaceConfig(695700/(1.1*JUPITER_RADIUS_KM),(30/365.25)**(2/3)*149597870.7/(1.1*JUPITER_RADIUS_KM),720,u1=.4,u2=.25)
    limb_radius=1+h/2
    tc=720/(2*np.pi)*np.arcsin(np.sqrt(star.radius_rp**2-limb_radius**2)/star.orbital_separation_rp)
    pm=float(moon_period_hours(3,.7,1.1));rows=[]
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for col,(name,center) in enumerate([('ingress',-tc),('egress',tc)]):
        # Controlled lunar contact at the stellar boundary, not random phase.
        phase=(np.arcsin(limb_radius/3)/(2*np.pi)-center/pm)%1
        moon=MoonConfig(radius=.03,semi_major_axis=3,period=pm,inclination=90,longitude_node=0,phase=phase)
        base=replace(reference_simulation(planet,.8,obs,inst),time_center=center,stellar_cell_integration=True)
        fine=replace(base,n_azimuths=960,n_radii=32,n_times=2*(base.n_times-1)+1,n_velocities=2*(base.n_velocities-1)+1,blocking_quadrature_order=32)
        models={}
        for mode,sim in [('point_star',replace(base,stellar_cell_integration=False)),('joint',base),('refined',fine)]:
            r=run_xds_model(planet,moon,sim,star);models[mode]=r
            axes[0,col].plot(r.time,np.trapz(r.reference_absorption,r.velocity,axis=1)/(1.5*np.sqrt(2*np.pi)),label=mode)
        for cadence in [0.,.25,.5,.75]:
            o=replace(obs,exposure_start_offset=cadence*obs.exposure_time)
            scores={k:detectability_report(v,atm,o,inst).matched_filter_snr for k,v in models.items()}
            change=abs(scores['joint']/scores['refined']-1)
            axes[1,col].plot(cadence,100*change,'o',color='C0')
            rows.append(dict(encounter=name,center_h=center,lunar_phase=phase,exposure_phase=cadence,
                point_star_snr=scores['point_star'],joint_snr=scores['joint'],refined_snr=scores['refined'],joint_relative_change=change))
        # Independent uniform-star circle-overlap control at matching epochs.
        uniform=replace(star,u1=0,u2=0);errs=[]
        for t in np.linspace(center-.4,center+.4,31):
            vis,_=stellar_cell_weights(480,16,h,t,uniform,np.array([3.,1.,0.]),.03,order=16)
            d=abs(star.orbital_separation_rp*np.sin(2*np.pi*t/720))
            exact=(circle_overlap(star.radius_rp,1+h,np.array([d]))-circle_overlap(star.radius_rp,1.,np.array([d])))/(np.pi*((1+h)**2-1))
            errs.append(abs(vis.mean()-exact[0]))
        print(name,'max_joint_change',max(x['joint_relative_change'] for x in rows if x['encounter']==name),'max_uniform_fraction_error',max(errs),flush=True)
        axes[0,col].set(title=name,xlabel='Time from planetary mid-transit [h]',ylabel='Normalized reference equivalent width')
        axes[0,col].legend(fontsize=8)
        axes[1,col].axhline(1,color='0.5',ls=':')
        axes[1,col].set(xlabel='Exposure start / duration',ylabel='Joint S/N refinement change [%]')
    with (root/'docs/ingress_egress_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/ingress_egress_comparison.{ext}',dpi=180)


if __name__=='__main__':main()
