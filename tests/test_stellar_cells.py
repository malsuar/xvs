import numpy as np
from xds.config import StellarSurfaceConfig
from xds.stellar_cells import stellar_cell_weights
from xds.penumbra import circle_overlap


def test_uniform_stellar_overlap_and_joint_bounds():
    star=StellarSurfaceConfig(8.85,359.44,720.)
    t=720/(2*np.pi)*np.arcsin(8.85/359.44)
    vis,joint=stellar_cell_weights(120,8,.03,t,star,np.array([-3.,1.,0.]),.06,order=24)
    exact=(circle_overlap(8.85,1.03,np.array([8.85]))-circle_overlap(8.85,1.,np.array([8.85])))/(np.pi*(1.03**2-1))
    assert abs(vis.mean()-exact[0])<1e-5
    assert np.all(joint>=0) and np.all(joint<=vis)
    dark,_=stellar_cell_weights(120,8,.03,4.,star,np.array([3.,1.,0.]),.06)
    assert not np.any(dark)


def test_fully_backed_joint_versus_lunar_fraction():
    from xds.cell_blocking import blocked_cell_fractions
    p=np.array([3.,1.01,0.]);star=StellarSurfaceConfig(8.85,359.44,720.)
    vis,joint=stellar_cell_weights(120,8,.03,0.,star,p,.06,order=24)
    lunar=blocked_cell_fractions(120,8,.03,p,.06,order=24)
    assert np.allclose(vis,1)
    assert np.allclose(joint,lunar,atol=1e-10)


def test_joint_intersection_is_not_product_of_cell_fractions():
    # Coincident projected stellar/lunar disks: their intersection is the disk
    # itself, not the square of its partial-cell coverage.
    star=StellarSurfaceConfig(.06,1.5,10.)
    t=10/(2*np.pi)*np.arcsin(1.01/1.5)
    vis,joint=stellar_cell_weights(120,8,.03,t,star,np.array([3.,1.01,0.]),.06,order=24)
    assert np.any((vis>0)&(vis<1))
    assert np.allclose(vis,joint,atol=1e-12)
    assert np.max(abs(vis*vis-joint))>.01


def test_production_limb_window_and_disabled_masks():
    from dataclasses import replace
    from xds.config import PlanetConfig,MoonConfig,SimulationConfig
    from xds.spectra import run_xds_model
    star=StellarSurfaceConfig(8.85,359.44,720.,u1=.4,u2=.25)
    center=-720/(2*np.pi)*np.arcsin(np.sqrt(8.85**2-1.015**2)/359.44)
    phase=(np.arcsin(1.015/3)/(2*np.pi)-center/21.24)%1
    sim=SimulationConfig(geometry='projected_annulus',spin_convention='observer',
        blocking_integration='fractional',stellar_cell_integration=True,
        time_center=center,transit_duration=.8,n_times=9,n_velocities=41,
        n_azimuths=120,n_radii=4)
    moon=MoonConfig(radius=.1,semi_major_axis=3,period=21.24,phase=phase,inclination=90)
    r=run_xds_model(PlanetConfig(),moon,sim,star)
    assert np.min(r.delta_absorption)<0
    assert np.allclose(r.reference_absorption[0],0)
    assert np.allclose(r.absorption,r.reference_absorption+r.delta_absorption)
    off=run_xds_model(PlanetConfig(),moon,replace(sim,use_illumination_shadow=False,use_observer_occultation=False),star)
    assert np.allclose(off.delta_absorption,0)
