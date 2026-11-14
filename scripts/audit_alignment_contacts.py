"""Analytic contact and projected-area convergence for aligned/polar orbits."""
import csv
from pathlib import Path
from dataclasses import replace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from xds import MoonConfig, PlanetConfig, SimulationConfig, run_xds_model
from xds.moon_orbit import moon_positions, spin_orbit_angle
from xds.rotation import spin_axis
from xds.terminator import build_projected_annulus


def overlap(radius, moon_radius, d):
    d = np.asarray(d)
    value = np.zeros_like(d)
    inside = d <= abs(radius-moon_radius)
    value[inside] = np.pi*min(radius, moon_radius)**2
    cross = (d > abs(radius-moon_radius)) & (d < radius+moon_radius)
    x = d[cross]
    a = np.clip((x*x+radius**2-moon_radius**2)/(2*x*radius), -1, 1)
    b = np.clip((x*x+moon_radius**2-radius**2)/(2*x*moon_radius), -1, 1)
    value[cross] = radius**2*np.arccos(a)+moon_radius**2*np.arccos(b)-.5*np.sqrt(
        np.maximum(0, (-x+radius+moon_radius)*(x+radius-moon_radius)*(x-radius+moon_radius)*(x+radius+moon_radius)))
    return value


def main():
    root = Path(__file__).resolve().parents[1]
    out = root/'outputs'; out.mkdir(exist_ok=True)
    h, rm, am, period = .03, .12, 3., 21.238292071543107
    moon = MoonConfig(radius=rm, semi_major_axis=am, period=period, phase=0., inclination=0.)
    # One complete passage through the front projected limb near conjunction.
    time = np.linspace(0, period/8, 1201)
    pos = moon_positions(time, moon)
    d = np.linalg.norm(pos[:,1:], axis=1)
    exact = (overlap(1+h, rm, d)-overlap(1, rm, d))/(np.pi*((1+h)**2-1))
    contacts = period/(2*np.pi)*np.arcsin(np.array([1-rm, 1+h+rm])/am)
    rows, curves = [], []
    for na, nr in [(120,4),(240,8),(480,16),(960,32)]:
        grid = build_projected_annulus(na,nr,h)
        fraction = np.empty(len(time))
        for k,p in enumerate(pos):
            blocked = np.sum((grid.positions[:,1:]-p[1:])**2,axis=1) <= rm**2
            fraction[k] = grid.weights @ blocked
        active = np.flatnonzero(fraction > 0)
        row = dict(n_azimuths=na,n_radii=nr,
                   max_absolute_fraction_error=float(np.max(abs(fraction-exact))),
                   integrated_fraction_relative_error=float(abs(np.trapz(fraction,time)/np.trapz(exact,time)-1)),
                   entry_error_minutes=float(60*(time[active[0]]-contacts[0])),
                   exit_error_minutes=float(60*(time[active[-1]]-contacts[1])))
        rows.append(row); curves.append(fraction)
    with (root/'docs/contact_convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    fig,axes=plt.subplots(1,3,figsize=(15,4),constrained_layout=True)
    axes[0].plot(time,exact,'k',label='Analytic annulus overlap')
    for row,curve in zip(rows,curves):
        axes[0].plot(time,curve,label=f"{row['n_azimuths']} x {row['n_radii']}",alpha=.7)
    axes[0].set(xlabel='Time since conjunction [h]',ylabel='Blocked fraction'); axes[0].legend(fontsize=8)
    axes[1].loglog([r['n_azimuths']*r['n_radii'] for r in rows],
                  [r['max_absolute_fraction_error'] for r in rows],'o-')
    axes[1].set(xlabel='Spatial elements',ylabel='Maximum absolute fraction error')
    theta=np.linspace(0,2*np.pi,500)
    axes[2].plot(np.cos(theta),np.sin(theta),'k',lw=1)
    for inc,label in [(0.,'Aligned, mutual angle 0 deg'),(90.,'Polar, mutual angle 90 deg')]:
        orbit=moon_positions(np.linspace(0,period,500),replace(moon,inclination=inc))
        axes[2].plot(orbit[:,1],orbit[:,2],label=label)
    axes[2].set(xlabel='Sky y/Rp',ylabel='Sky z/Rp',aspect='equal'); axes[2].legend(fontsize=8)
    for ext in ['png','pdf']: fig.savefig(out/f'alignment_contact_audit.{ext}',dpi=180)
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    spin=spin_axis(90,0,convention='observer')
    for ax,inc in zip(axes,[0.,90.]):
        m=replace(moon,inclination=inc,phase=.08)
        sim=SimulationConfig(geometry='projected_annulus',spin_convention='observer',
                             n_azimuths=480,n_radii=8,n_times=481,n_velocities=301,
                             transit_duration=16)
        result=run_xds_model(PlanetConfig(line_sigma=1.5),m,sim)
        im=ax.imshow(result.delta_absorption,origin='lower',aspect='auto',
                     extent=[-25,25,-8,8],cmap='Blues_r',vmin=-.08,vmax=0)
        ax.set(xlabel='Velocity toward observer [km/s]',ylabel='Time [h]',
               title=f'Mutual angle {spin_orbit_angle(m,spin):.0f} deg')
    fig.colorbar(im,ax=axes,label='Dimensionless residual')
    for ext in ['png','pdf']: fig.savefig(out/f'alignment_residual_maps.{ext}',dpi=180)
    print('Analytic contacts [h]',contacts)
    print(rows)


if __name__=='__main__': main()
