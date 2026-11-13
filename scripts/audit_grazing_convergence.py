"""Separate angular/radial convergence of circular grazing lunar contacts."""
import csv
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from audit_alignment_contacts import overlap
from xds import MoonConfig
from xds.config import AtmosphereConfig, InstrumentConfig, ObservationConfig
from xds.detectability import atmospheric_line_depth
from xds.moon_orbit import moon_positions
from xds.observing import observe_residual
from xds.rotation import line_of_sight_velocity
from xds.spectra import gaussian_profiles
from xds.terminator import build_projected_annulus
from xds.cell_blocking import blocked_cell_fractions


def evaluate(na, nr, b, nt=1201, cell_phase=0., exposure_phase=0., all_exposure_phases=False,
             moon_radius=.12, radial_sample_fraction=.5, fractional=False, order=16,
             n_velocities=401, velocity_limit=30., v_eq=12., line_sigma=1.5,
             return_result=False):
    time = np.linspace(-2., 2., nt)
    velocity = np.linspace(-velocity_limit, velocity_limit, n_velocities)
    # Normal=(b/a,0,sqrt(1-(b/a)^2)); conjunction at mean anomaly pi/2.
    moon = MoonConfig(radius=moon_radius, semi_major_axis=3., period=21.238292071543107,
                      inclination=np.rad2deg(np.arcsin(b/3)), longitude_node=90., phase=.25)
    positions = moon_positions(time, moon)
    grid = build_projected_annulus(na, nr, .03,radial_sample_fraction=radial_sample_fraction)
    # Shift azimuths by a fraction of one cell without moving the physical orbit.
    angle = cell_phase*2*np.pi/na
    y,z = grid.positions[:,1].copy(),grid.positions[:,2].copy()
    grid.positions[:,1] = np.cos(angle)*y-np.sin(angle)*z
    grid.positions[:,2] = np.sin(angle)*y+np.cos(angle)*z
    vlos = line_of_sight_velocity(grid.positions,v_eq,90.,0.,(-1,0,0),"observer")
    kernels = gaussian_profiles(velocity,vlos,line_sigma)
    signal = np.zeros((nt,len(velocity)))
    fractions = np.zeros(nt)
    # Parallel star/observer rays and sky-plane elements: union of both masks
    # is the projected disk wherever moon x != 0 (all epochs used here).
    assert np.all(abs(positions[:,0]) > 0)
    for k,p in enumerate(positions):
        if fractional:
            blocked=blocked_cell_fractions(na,nr,.03,p,moon_radius,order=order,angular_phase=cell_phase)
            mask=blocked>0
            fractions[k]=grid.weights @ blocked
            if np.any(mask): signal[k]=-((grid.weights[mask]*blocked[mask]) @ kernels[mask])
            continue
        mask = np.sum((grid.positions[:,1:]-p[1:])**2,axis=1) <= moon_radius**2
        fractions[k] = grid.weights[mask].sum()
        if np.any(mask): signal[k] = -(grid.weights[mask] @ kernels[mask])
    d = np.linalg.norm(positions[:,1:],axis=1)
    exact = (overlap(1.03,moon_radius,d)-overlap(1.,moon_radius,d))/(np.pi*(1.03**2-1))
    result = SimpleNamespace(time=time,velocity=velocity,delta_absorption=signal)
    atmosphere = AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,
                                 stellar_radius_solar=1.,temperature=500.)
    obs = ObservationConfig(continuum_snr_per_resolution_element=300., exposure_time=5/60,
                            n_lines=200,systematic_floor=10e-6,
                            exposure_start_offset=exposure_phase*5/60)
    measured = observe_residual(result,atmospheric_line_depth(atmosphere),obs,InstrumentConfig())
    snr = np.linalg.norm(measured.signal/measured.noise)*np.sqrt(obs.n_lines)
    area_error = abs(np.trapz(fractions,time)/np.trapz(exact,time)-1)
    active = np.flatnonzero(fractions>0)
    duration = 0. if len(active)==0 else time[active[-1]]-time[active[0]]
    # Outer contact span includes any intervening gap for a central encounter.
    theta = np.arcsin(np.sqrt(((1.03+moon_radius)**2-b*b)/(9-b*b)))
    exact_span = 2*moon.period/(2*np.pi)*theta
    record = dict(impact_rp=b,n_azimuths=na,n_radii=nr,n_times=nt,snr=float(snr),
                  area_relative_error=float(area_error),
                  max_fraction_error=float(np.max(abs(fractions-exact))),
                  outer_span_error_minutes=float(60*(duration-exact_span)))
    if all_exposure_phases:
        from dataclasses import replace
        for phase in [0.,.25,.5,.75]:
            observed=observe_residual(result,atmospheric_line_depth(atmosphere),
                                      replace(obs,exposure_start_offset=phase*5/60),InstrumentConfig())
            record[f'snr_phase_{phase}']=float(np.linalg.norm(observed.signal/observed.noise)*np.sqrt(obs.n_lines))
    if return_result:
        return record,time,fractions,exact,result
    return record,time,fractions,exact


def main():
    root=Path(__file__).resolve().parents[1]; out=root/'outputs'; out.mkdir(exist_ok=True)
    rows=[]
    fig,axes=plt.subplots(3,3,figsize=(13,10),constrained_layout=True)
    for column,b in enumerate([0.,1.08,1.14]):
        ref,t,f,e=evaluate(1920,64,b)
        rows.append(dict(ref,series='reference',snr_relative_error=0.))
        axes[0,column].plot(t,e,'k',label='Analytic area')
        axes[0,column].plot(t,f,'--',label='1920 x 64')
        axes[0,column].set(title=f'Projected impact b={b:.2f} Rp',xlabel='Time [h]',ylabel='Blocked fraction')
        axes[0,column].legend(fontsize=8)
        for series,grids,color in [('angular',[(n,64) for n in [120,240,480,960]],'C0'),
                                    ('radial',[(1920,n) for n in [2,4,8,16,32]],'C1')]:
            subset=[]
            for na,nr in grids:
                row,_,_,_=evaluate(na,nr,b)
                row.update(series=series,snr_relative_error=abs(row['snr']/ref['snr']-1))
                rows.append(row); subset.append(row)
            x=[r['n_azimuths'] if series=='angular' else r['n_radii'] for r in subset]
            axes[1,column].loglog(x,[max(r['area_relative_error'],1e-12) for r in subset],'o-',color=color,label=series)
            axes[2,column].loglog(x,[max(r['snr_relative_error'],1e-12) for r in subset],'o-',color=color,label=series)
        # Reference time refinement provides a separate temporal control.
        finer,_,_,_=evaluate(1920,64,b,nt=2401)
        finer.update(series='time_check',snr_relative_error=abs(finer['snr']/ref['snr']-1))
        rows.append(finer)
        for k,label in [(1,'Integrated area relative error'),(2,'Corrected S/N relative error')]:
            axes[k,column].set(xlabel='Cells in refined direction',ylabel=label)
            axes[k,column].legend(fontsize=8)
        print('CASE',b,'reference',ref,'time_check_snr_error',finer['snr_relative_error'],flush=True)
    with (root/'docs/grazing_convergence.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for ext in ['png','pdf']:fig.savefig(out/f'grazing_spatial_snr_convergence.{ext}',dpi=180)


if __name__=='__main__':main()
