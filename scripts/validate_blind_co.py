"""Independent moving-orbit spectra and quadrature checks against production."""
import json
from pathlib import Path
import numpy as np
from nonlinear_co_recovery import AXES,direct
from molecular_injection_recovery import build_molecular_geometry
from xds.moving_shadow import MovingShadow,separation_from_period,period_limits,REFERENCE_PERIOD
from xds import MoonConfig,StellarSurfaceConfig
from xds.moon_orbit import moon_positions
from xds.terminator import build_projected_annulus
from xds.rotation import line_of_sight_velocity
from xds.cell_blocking import blocked_cell_fractions
from xds.stellar_surface import background_intensity
from xds.molecular import bin_average


def independent_geometry(base,period,phase,nt=1201):
    nr=len(base['radii']);na=480;h=np.sqrt(1+2*nr*(base['radii'][0]**2-1))-1
    grid=build_projected_annulus(na,nr,h);rp=base['rp']
    star=StellarSurfaceConfig(695700/rp,(30/365.25)**(2/3)*149597870.7/rp,720.,u1=.4,u2=.25)
    velocity=line_of_sight_velocity(grid.positions,12.,90.,0.,(-1.,0.,0.),'observer')
    shifts=base['shifts'];step=shifts[1]-shifts[0]
    index=np.clip(np.floor((velocity-shifts[0])/step).astype(int),0,len(shifts)-2)
    fraction=(velocity-shifts[index])/step;flat=np.repeat(np.arange(nr),na)*len(shifts)+index
    moon=MoonConfig(radius=.09,semi_major_axis=float(separation_from_period(period)),period=float(period),
        phase=float(phase),inclination=0.,longitude_node=90.)
    time=np.linspace(-2,2,nt);positions=moon_positions(time,moon)
    missing=np.zeros((2,nt,nr,len(shifts)))
    for k,position in enumerate(positions):
        blocked=blocked_cell_fractions(na,nr,h,position,.09,order=16)
        weight=grid.weights*blocked*background_intensity(grid.positions,time[k],star)
        for g,region in enumerate([velocity<0,velocity>=0]):
            w=weight*region
            missing[g,k]=(np.bincount(flat,weights=w*(1-fraction),minlength=nr*len(shifts))+
                np.bincount(flat+1,weights=w*fraction,minlength=nr*len(shifts))).reshape(nr,-1)
    edges=np.linspace(-2,2,49)
    missing=np.array([bin_average(time,m.reshape(nt,-1).T,edges).T.reshape(48,nr,-1) for m in missing])
    return dict(base,missing=missing)


def main():
    root=Path(__file__).resolve().parents[1];saved=np.load(root/'outputs/co_blind_shadow.npz')
    shadow=MovingShadow(AXES,saved['library'],saved['left'],saved['centers'],saved['rotation'],float(saved['rp']))
    basis=saved['basis'];base=build_molecular_geometry(nr=52,n_scale_heights=16.,split_limbs=True)
    noise=base['noise']/130.8285010091114
    lo,hi=period_limits();cases=[(lo,.75),(REFERENCE_PERIOD,.7473),(hi,.75),(13.,.2841),(50.,.21)]
    checks=[]
    for i,(period,phase) in enumerate(cases):
        p=np.zeros((2,3)) if i%2==0 else np.array([[-.23,-.16,.21],[.19,.22,-.27]])
        geometry=independent_geometry(base,period,phase)
        physical=direct(root,geometry,p)['signal']
        trajectory=np.broadcast_to(p[:,None,:],(2,48,3))
        coarse=shadow.signal(trajectory,period,phase,base['time'])@basis.T
        predicted=shadow.signal(trajectory,period,phase,base['time'],order=16)@basis.T
        norm=np.linalg.norm(physical)
        check=dict(period_h=float(period),phase=float(phase),signal_control_noise_norm=float(norm/noise),
            error_control_noise_norm=float(np.linalg.norm(physical-predicted)/noise),
            relative_signal_error=float(np.linalg.norm(physical-predicted)/norm) if norm>1e-20 else None,
            quadrature8_vs16_control_noise_norm=float(np.linalg.norm(coarse-predicted)/noise))
        checks.append(check);print(check,flush=True)
    passed=all(c['error_control_noise_norm']<.05 and (c['relative_signal_error'] is None or
        c['signal_control_noise_norm']<.1 or c['relative_signal_error']<.01) for c in checks)
    result=dict(checks=checks,passed=passed,exposure_order=16,criterion='absolute error <0.05 control-noise units; <1% for signal norm>0.1')
    (root/'docs/co_blind_forward_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    if not passed:raise RuntimeError('Moving shadow failed production comparison')


if __name__=='__main__':main()
