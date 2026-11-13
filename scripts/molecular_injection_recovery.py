"""CO line-by-line first recovery control against per-exposure depth changes."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from xds import PlanetConfig,MoonConfig,StellarSurfaceConfig
from xds.config import AtmosphereConfig,ObservationConfig,InstrumentConfig
from xds.reference import reference_simulation
from xds.terminator import build_projected_annulus
from xds.rotation import line_of_sight_velocity
from xds.moon_orbit import moon_positions
from xds.cell_blocking import blocked_cell_fractions
from xds.stellar_surface import background_intensity
from xds.detectability import atmospheric_line_depth,scale_height_km,JUPITER_RADIUS_KM
from xds.dynamics import moon_period_hours
from xds.molecular import load_co_lines,co_absorption,bin_average,depth_orthogonal_component,C


def synthesize(root,na=480,nr=32,nt=1201,dv=.2,shift_step=.25,wing_sigma=100.,n_scale_heights=10.):
    atm=AtmosphereConfig(planet_mass_jupiter=.7,planet_radius_jupiter=1.1,temperature=640.,effective_scale_heights=n_scale_heights)
    hp=scale_height_km(atm);rp=1.1*JUPITER_RADIUS_KM;h=n_scale_heights*hp/rp
    star=StellarSurfaceConfig(695700/rp,(30/365.25)**(2/3)*149597870.7/rp,720.,u1=.4,u2=.25)
    planet=PlanetConfig(v_eq=12.,atmosphere_half_width=h)
    moon=MoonConfig(radius=.09,semi_major_axis=3.,period=float(moon_period_hours(3,.7,1.1)),
                    inclination=0.,longitude_node=90.,phase=.75)
    obs=ObservationConfig(exposure_time=5/60,continuum_snr_per_resolution_element=300.,n_lines=1,systematic_floor=10e-6)
    grid=build_projected_annulus(na,nr,h)
    velocity=line_of_sight_velocity(grid.positions,12.,90.,0.,(-1,0,0),'observer')
    vlim=np.ceil(max(abs(velocity))/shift_step)*shift_step
    shifts=np.arange(-vlim,vlim+shift_step/2,shift_step)
    index=np.clip(np.floor((velocity+vlim)/shift_step).astype(int),0,len(shifts)-2)
    frac=(velocity-shifts[index])/shift_step
    layer=np.repeat(np.arange(nr),na);flat=layer*len(shifts)+index
    time=np.linspace(-2,2,nt);positions=moon_positions(time,moon)
    full=np.zeros((nt,nr,len(shifts)));missing=np.zeros_like(full)
    for k,p in enumerate(positions):
        intensity=background_intensity(grid.positions,time[k],star)
        weights=grid.weights*intensity
        blocked=blocked_cell_fractions(na,nr,h,p,.09,order=16)
        for dest,w in [(full,weights),(missing,weights*blocked)]:
            dest[k]=(np.bincount(flat,weights=w*(1-frac),minlength=nr*len(shifts))+
                     np.bincount(flat+1,weights=w*frac,minlength=nr*len(shifts))).reshape(nr,-1)
    tedges=np.arange(49)*obs.exposure_time-2
    full=bin_average(time,full.reshape(nt,-1).T,tedges).T.reshape(48,nr,-1)
    missing=bin_average(time,missing.reshape(nt,-1).T,tedges).T.reshape(48,nr,-1)
    xmin=C*np.log(4250/4300);xmax=C*np.log(4400/4300)
    x=np.arange(xmin-50,xmax+50+dv/2,dv);nu=4300*np.exp(x/C)
    dvinst=C/100000
    edges=xmin+np.arange(int((xmax-xmin)/dvinst)+1)*dvinst
    lines=load_co_lines(root/'data/exomol/CO_Li2015',640.)
    reference=np.zeros((48,len(edges)-1));signal=np.zeros_like(reference)
    radii=np.sqrt(1+(np.arange(nr)+.5)*((1+h)**2-1)/nr)
    profiles=[];columns=[]
    for j,r in enumerate(radii):
        pbar=.001*np.exp(-(r-1)*rp/hp)
        # Isothermal tangent-column approximation; constant CO mixing 1e-4.
        column=1e-4*(pbar*1e5/(1.380649e-23*640.)/1e6)*np.sqrt(2*np.pi*(r*rp*1e5)*(hp*1e5))
        local=co_absorption(nu,lines,640.,pbar,column,wing_sigma=wing_sigma)
        profiles.append(local);columns.append(column)
        broadened=gaussian_filter1d(local,dvinst/np.sqrt(8*np.log(2))/dv,mode='constant',truncate=8.)
        library=np.array([np.interp(x-s,x,broadened,left=0,right=0) for s in shifts])
        library=bin_average(x,library,edges)
        reference+=full[:,j]@library
        signal-=missing[:,j]@library
    # Exact annular area, rather than the earlier thin-annulus approximation.
    scale=atmospheric_line_depth(atm)*(1+h/2)
    noise=float(np.hypot(1/300*np.sqrt((1/60)/obs.exposure_time),obs.systematic_floor))
    return dict(reference=scale*reference,signal=scale*signal,noise=noise,
        time=(tedges[:-1]+tedges[1:])/2,wavenumber=4300*np.exp((edges[:-1]+edges[1:])/2/C),
        input_wavenumber=nu,local_profiles=np.array(profiles),column_cm2=np.array(columns),
        n_lines=len(lines.wavenumber),partition=lines.partition,
        retained_strength_fraction=lines.retained_strength_fraction,annulus_area_scale=scale)


def main():
    root=Path(__file__).resolve().parents[1];out=root/'outputs';out.mkdir(exist_ok=True)
    model=synthesize(root)
    s,b,noise=model['signal'],model['reference'],model['noise']
    q=depth_orthogonal_component(s,b)
    snr_full=float(np.linalg.norm(s)/noise);snr_kin=float(np.linalg.norm(q)/noise)
    rng=np.random.default_rng(63117);nmc=2000
    # Actual independent Gaussian samples, projected onto two known templates.
    zfull=[];zkin=[]
    full_unit=(s/np.linalg.norm(s)).ravel();kin_unit=(q/np.linalg.norm(q)).ravel()
    for start in range(0,nmc,32):
        n=rng.normal(size=(min(32,nmc-start),s.size))
        zfull.extend(n@full_unit);zkin.extend(n@kin_unit)
    zfull=np.array(zfull);zkin=np.array(zkin)
    laboratory_scale=5/snr_kin
    # Known-orbit test: signed Gaussian score, no orbital-template search.
    rows=[]
    for name,amplitude in [('native',1.),('controlled_kinematic_snr5',laboratory_scale)]:
        z=zkin+amplitude*snr_kin
        rows.append(dict(injection=name,amplitude_multiplier=amplitude,
            expected_full_snr=amplitude*snr_full,expected_depth_profiled_snr=amplitude*snr_kin,
            mean_recovered_amplitude=float(z.mean()/snr_kin),
            amplitude_std=float(z.std(ddof=1)/snr_kin),fraction_signed_z_above_3=float(np.mean(z>3))))
    summary=dict(n_lines=model['n_lines'],partition=model['partition'],
        retained_strength_fraction=model['retained_strength_fraction'],
        annulus_area_scale_ppm=model['annulus_area_scale']*1e6,noise_per_bin_ppm=noise*1e6,
        n_exposures=s.shape[0],n_spectral_bins=s.shape[1],native_full_snr=snr_full,
        native_depth_profiled_snr=snr_kin,kinematic_norm_fraction=snr_kin/snr_full,
        controlled_amplitude_multiplier=laboratory_scale,null_z_mean=float(zkin.mean()),
        null_z_std=float(zkin.std(ddof=1)),null_fraction_z_above_3=float(np.mean(zkin>3)),
        seed=63117,n_noise_realizations=nmc,
        min_column_cm2=float(min(model['column_cm2'])),max_column_cm2=float(max(model['column_cm2'])))
    # Depth-only injection of the best depth approximation to the lunar signal.
    depth=s-q
    summary['depth_only_expected_kinematic_z']=float(np.sum(depth*q)/(noise*np.linalg.norm(q)))
    depth_scale=5/(np.sum(depth*s)/(noise*np.linalg.norm(s)))
    rows.append(dict(injection='controlled_depth_only_full_score5',amplitude_multiplier=depth_scale,
        expected_full_snr=5.,expected_depth_profiled_snr=0.,
        mean_recovered_amplitude=float(zkin.mean()/snr_kin),
        amplitude_std=float(zkin.std(ddof=1)/snr_kin),fraction_signed_z_above_3=float(np.mean(zkin>3))))
    # Explicit injected noisy sequence and nested least-squares fits.
    example=laboratory_scale*s+noise*rng.normal(size=s.shape)
    z_example=float(np.sum(example*q)/(noise*np.linalg.norm(q)))
    fitted_amplitude=z_example/snr_kin
    depth_fit=example-depth_orthogonal_component(example,b)
    remainder=example-fitted_amplitude*s
    joint_fit=fitted_amplitude*s+remainder-depth_orthogonal_component(remainder,b)
    delta_chi2=float(np.sum(((example-depth_fit)/noise)**2)-np.sum(((example-joint_fit)/noise)**2))
    assert np.isclose(delta_chi2,z_example**2,atol=1e-7)
    summary.update(controlled_example_z=z_example,controlled_example_fitted_amplitude=fitted_amplitude,
        controlled_example_delta_chi2=delta_chi2,depth_only_control_multiplier=float(depth_scale),
        atmospheric_temperature_k=640.,base_pressure_bar=.001,co_mixing_ratio=1e-4,
        n_scale_heights=10.,n_radii=32,amplitude_normalization='exact annular projected area / stellar disk area',
        microturbulence_km_s=1.,isotope='12C-16O',resolving_power=100000,
        exposure_minutes=5.,continuum_snr_reference=300.,reference_exposure_minutes=1.,
        species_line_gain_multiplier=1.,stellar_limb_u1=.4,stellar_limb_u2=.25,
        planet_mass_mj=.7,planet_radius_rj=1.1,moon_radius_rp=.09,moon_separation_rp=3.,
        lunar_phase_at_midtransit=.75,geometry='aligned_circular_annulus',
        retrieval='known orbit, unrestricted signed amplitude, one depth nuisance per exposure')
    with (root/'docs/co_injection_recovery.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (root/'docs/co_injection_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    np.savez_compressed(out/'co_molecular_injection.npz',**model,kinematic_component=q)
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    sel=(model['input_wavenumber']>=4250)&(model['input_wavenumber']<=4400)
    axes[0,0].plot(1e4/model['input_wavenumber'][sel],model['local_profiles'][0,sel],lw=.7,label='Lowest radial cell')
    axes[0,0].plot(1e4/model['input_wavenumber'][sel],model['local_profiles'][-1,sel],lw=.7,label='Highest radial cell')
    axes[0,0].set(xlabel='Wavelength [micron]',ylabel='Local CO absorption');axes[0,0].legend(fontsize=8)
    axes[0,1].plot(model['time'],np.linalg.norm(s,axis=1)/noise,label='All lunar residual')
    axes[0,1].plot(model['time'],np.linalg.norm(q,axis=1)/noise,label='After depth removal')
    axes[0,1].set(xlabel='Time [h]',ylabel='Template norm / noise per exposure');axes[0,1].legend(fontsize=8)
    k=np.argmax(np.linalg.norm(s,axis=1));nu=model['wavenumber']
    band=(nu>4300)&(nu<4310)
    axes[1,0].plot(nu[band],s[k,band]*1e6,label='Kinematic lunar template')
    axes[1,0].plot(nu[band],depth[k,band]*1e6,label='Best depth-only approximation')
    axes[1,0].set(xlabel='Wavenumber [cm^-1]',ylabel='Residual [ppm]');axes[1,0].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,1.18))
    for values,label in [(zkin,'Noise / depth-only null'),(zkin+snr_kin,'Native injection'),(zkin+5,'Controlled injection, expected Z=5')]:
        axes[1,1].hist(values,bins=45,density=True,histtype='step',label=label)
    axes[1,1].set(xlabel='Signed known-trajectory score after depth removal',ylabel='Density');axes[1,1].legend(fontsize=7,loc='upper center',bbox_to_anchor=(.5,1.18))
    for ext in ['png','pdf']:fig.savefig(out/f'co_injection_recovery.{ext}',dpi=180)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
