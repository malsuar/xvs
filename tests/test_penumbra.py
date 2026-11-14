import numpy as np
from xds.penumbra import circle_overlap,irradiance_blocked_fractions
from xds.cell_blocking import blocked_cell_fractions


def test_uniform_disk_overlap_limits():
    assert np.allclose(circle_overlap(2,1,np.array([0.,.5,4.])),[np.pi,np.pi,0.])
    # Central irradiance deficit when projected stellar blur exceeds the moon.
    assert np.isclose(circle_overlap(.1,.03,np.array([0.]))[0]/(np.pi*.1**2),.09)


def test_point_source_limit_and_downstream():
    p=np.array([3.,1.,0.])
    f=irradiance_blocked_fractions(120,8,.03,p,.06,0.,1e12)
    ref=blocked_cell_fractions(120,8,.03,p,.06)
    assert np.allclose(f,ref,atol=1e-9)
    assert not np.any(irradiance_blocked_fractions(120,8,.03,-p,.06,9.,360.))


def test_penumbra_cell_quadrature_bounds_and_convergence():
    p=np.array([3.,1.03,0.])
    values=[]
    for order in [4,8,16]:
        f=irradiance_blocked_fractions(240,8,.03,p,.03,8.85,360.,order=order)
        assert np.all((f>=0)&(f<=1))
        values.append(f.mean())
    assert abs(values[-1]/values[-2]-1)<.005
