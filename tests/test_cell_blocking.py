import numpy as np
from xds.cell_blocking import blocked_cell_fractions
from xds import PlanetConfig,MoonConfig,SimulationConfig,run_xds_model


def test_fractional_area_concentric_and_empty():
    f=blocked_cell_fractions(120,8,.03,[3,0,0],1.015)
    assert np.all((f>=0)&(f<=1))
    assert np.isclose(f.mean(),(1.015**2-1)/(1.03**2-1))
    assert not np.any(blocked_cell_fractions(120,8,.03,[3,3,0],.03))
    assert np.allclose(blocked_cell_fractions(120,8,.03,[3,0,0],2),1)


def test_fractional_forward_reference_and_masks():
    from dataclasses import replace
    sim=SimulationConfig(geometry='projected_annulus',spin_convention='observer',
        blocking_integration='fractional',n_azimuths=120,n_radii=4,n_times=31,n_velocities=41)
    planet=PlanetConfig();moon=MoonConfig(radius=.3,inclination=0,phase=.08)
    result=run_xds_model(planet,moon,sim)
    assert np.max(result.combined_fraction)>0
    assert np.all(result.delta_absorption<=1e-14)
    off=run_xds_model(planet,moon,replace(sim,use_illumination_shadow=False,use_observer_occultation=False))
    assert np.allclose(off.absorption,off.reference_absorption)


def test_partial_overlap_rotation_and_quadrature():
    # Known disk intersection area, including a crossing of the angular seam.
    from math import acos, sqrt, pi
    def area(r,q,d):
        if d>=r+q:return 0.
        if d<=abs(r-q):return pi*min(r,q)**2
        return (r*r*acos((d*d+r*r-q*q)/(2*d*r))+
                q*q*acos((d*d+q*q-r*r)/(2*d*q))-
                .5*sqrt((-d+r+q)*(d+r-q)*(d-r+q)*(d+r+q)))
    d,q=1.057,.03
    exact=(area(1.03,q,d)-area(1.,q,d))/(pi*(1.03**2-1))
    for angle in [0.,.7,3.14]:
        p=[3.,d*np.cos(angle),d*np.sin(angle)]
        f=blocked_cell_fractions(240,8,.03,p,q,order=32)
        assert abs(f.mean()/exact-1)<1e-4


def test_fractional_deficit_equivalent_width_and_mask_partition():
    from dataclasses import replace
    sim=SimulationConfig(geometry='projected_annulus',spin_convention='observer',
        blocking_integration='fractional',n_azimuths=120,n_radii=4,
        n_times=31,n_velocities=501,velocity_min=-40,velocity_max=40)
    planet=PlanetConfig(line_depth=2.5,line_sigma=1.5)
    moon=MoonConfig(radius=.3,inclination=0,phase=.08)
    combined=run_xds_model(planet,moon,sim)
    ill=run_xds_model(planet,moon,replace(sim,use_observer_occultation=False))
    obs=run_xds_model(planet,moon,replace(sim,use_illumination_shadow=False))
    assert np.max(combined.combined_fraction)>0
    assert np.allclose(combined.delta_absorption,ill.delta_absorption+obs.delta_absorption)
    integral=np.trapz(-combined.delta_absorption,combined.velocity,axis=1)
    expected=2.5*1.5*np.sqrt(2*np.pi)*combined.combined_fraction
    assert np.allclose(integral,expected,atol=1e-12)
