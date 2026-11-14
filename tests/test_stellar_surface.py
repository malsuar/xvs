import numpy as np
from xds.config import StellarSurfaceConfig,PlanetConfig,MoonConfig,SimulationConfig
from xds.stellar_surface import background_intensity
from xds.spectra import run_xds_model


def test_stellar_flux_normalization_and_visibility():
    star=StellarSurfaceConfig(10.,360.,720.,u1=.4,u2=.25)
    r2=(np.arange(10000)+.5)/10000
    p=np.column_stack([np.zeros(len(r2)),10*np.sqrt(r2),np.zeros(len(r2))])
    assert abs(background_intensity(p,0,star).mean()-1)<1e-6
    assert background_intensity(np.array([[0,11,0]]),0,star)[0]==0


def test_uniform_star_and_time_dependent_reference():
    sim=SimulationConfig(geometry='projected_annulus',spin_convention='observer',
        blocking_integration='fractional',n_azimuths=120,n_radii=4,n_times=31,n_velocities=101,transit_duration=4)
    p=PlanetConfig();m=MoonConfig(radius=.3,phase=.08,inclination=0)
    star=StellarSurfaceConfig(8.85,359.44,720.)
    old=run_xds_model(p,m,sim)
    uniform=run_xds_model(p,m,sim,star)
    assert np.array_equal(old.delta_absorption,uniform.delta_absorption)
    from dataclasses import replace
    ld=run_xds_model(p,m,sim,replace(star,u1=.4,u2=.25))
    assert not np.allclose(ld.reference_absorption[0],ld.reference_absorption[15])
    assert np.allclose(ld.absorption,ld.reference_absorption+ld.delta_absorption)
