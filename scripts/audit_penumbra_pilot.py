"""Compare direct-ray pilot with a distinct irradiance-response toy model."""
import csv
import json
from dataclasses import replace
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from xds import PlanetConfig,MoonConfig,run_xds_model
from xds.config import AtmosphereConfig,ObservationConfig,InstrumentConfig
from xds.reference import reference_simulation
from xds.detectability import scale_height_km,JUPITER_RADIUS_KM,detectability_report
from xds.dynamics import moon_period_hours
from xds.penumbra import irradiance_blocked_fractions
from xds.spectra import gaussian_profiles


def main():
    root=Path(__file__).resolve().parents[1];(root/'outputs').mkdir(exist_ok=True)
    atm=AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,temperature=640.)
    h=5*scale_height_km(atm)/(1.1*JUPITER_RADIUS_KM)
    star_radius=695700/(1.1*JUPITER_RADIUS_KM)
    star_distance=(30/365.25)**(2/3)*149597870.7/(1.1*JUPITER_RADIUS_KM)
    planet=PlanetConfig(v_eq=12,line_sigma=1.5,atmosphere_half_width=h)
    # Entire annulus remains backed by the uniform stellar disk in this window.
    maximum_background_radius=star_distance*np.sin(2*np.pi*2/(30*24))+(1+h)
    assert maximum_background_radius<star_radius
    obs=ObservationConfig(exposure_time=5/60,continuum_snr_per_resolution_element=300,n_lines=200,systematic_floor=10e-6)
    inst=InstrumentConfig();sim=reference_simulation(planet,4.,obs,inst)
    rows=[];fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for col,rm in enumerate([.03,.06,.09]):
        for kind,b in [('central',0.),('grazing',1+h+.9*rm)]:
            for phase in [.25,.75]:
                moon=MoonConfig(radius=rm,semi_major_axis=3.,period=float(moon_period_hours(3,.7,1.1)),
                     phase=phase,inclination=np.rad2deg(np.arcsin(b/3)),longitude_node=90.)
                direct=run_xds_model(planet,moon,sim)
                kernels=gaussian_profiles(direct.velocity,direct.v_los,planet.line_sigma)
                models={}
                for order in ([8,16] if phase==.75 else [8]):
                    delta=direct.delta_absorption.copy();fraction=direct.combined_fraction.copy()
                    for k,p in enumerate(direct.moon_positions):
                        if p[0]<=0:continue
                        f=irradiance_blocked_fractions(sim.n_azimuths,sim.n_radii,h,p,rm,
                            star_radius,star_distance,order=order)
                        selected=f>0
                        delta[k]=-(direct.terminator.weights[selected]*f[selected])@kernels[selected]
                        fraction[k]=direct.terminator.weights@f
                    models[order]=replace(direct,delta_absorption=delta,combined_fraction=fraction)
                for cadence in [0.,.25,.5,.75]:
                    o=replace(obs,exposure_start_offset=cadence*obs.exposure_time)
                    a=detectability_report(direct,atm,o,inst)
                    s=detectability_report(models[8],atm,o,inst)
                    fine=detectability_report(models[max(models)],atm,o,inst)
                    rows.append(dict(radius_rp=rm,encounter=kind,lunar_phase=phase,exposure_phase=cadence,
                        direct_ray_snr=a.matched_filter_snr,irradiance_toy_snr=fine.matched_filter_snr,
                        toy_to_direct_ratio=fine.matched_filter_snr/a.matched_filter_snr,
                        quadrature_relative_change=abs(s.matched_filter_snr/fine.matched_filter_snr-1),
                        direct_peak_ppm=1e6*a.peak_xds_depth,toy_peak_ppm=1e6*fine.peak_xds_depth))
                if phase==.75:
                    if kind=='central':
                        axes[0,col].plot(direct.time,direct.combined_fraction,label='Direct ray')
                        axes[0,col].plot(direct.time,models[16].combined_fraction,label='Irradiance toy')
                    subset=[r for r in rows if r['radius_rp']==rm and r['encounter']==kind and r['lunar_phase']==phase]
                    axes[1,col].plot([r['exposure_phase'] for r in subset],[r['toy_to_direct_ratio'] for r in subset],'o-',label=kind)
                print('CASE',rm,kind,phase,'ratio',rows[-1]['toy_to_direct_ratio'],'quad',rows[-1]['quadrature_relative_change'],flush=True)
        axes[0,col].set(title=f'Rm={rm:.2f} Rp (illumination)',xlabel='Time [h]',ylabel='Blocked fraction')
        axes[1,col].set(xlabel='Exposure start / duration',ylabel='Irradiance toy / direct-ray S/N')
        for ax in axes[:,col]:ax.legend(fontsize=8)
    with (root/'docs/penumbra_pilot.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (root/'docs/penumbra_assumptions.json').write_text(json.dumps(dict(stellar_radius_rp=star_radius,
        stellar_distance_rp=star_distance,stellar_angular_radius_approx_rad=star_radius/star_distance,
        maximum_approx_blur_rp=3*star_radius/star_distance,
        maximum_background_radius_rp=float(maximum_background_radius),
        direct_ray_uniform_star_fully_backed=True),indent=2)+'\n')
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/penumbra_pilot.{ext}',dpi=180)


if __name__=='__main__':main()
