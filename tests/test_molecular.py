import numpy as np
from xds.molecular import COLines,co_absorption,bin_average,depth_orthogonal_component


def test_exact_linear_bin_integration():
    x=np.linspace(0,1,11);edges=np.array([.025,.225,.725,.975])
    result=bin_average(x,np.array([2*x+1,3*x]),edges)
    center=(edges[1:]+edges[:-1])/2
    assert np.allclose(result,np.array([2*center+1,3*center]))


def test_depth_projection_removes_depth_but_retains_shapes():
    ref=np.array([[1.,2.,3.],[3.,1.,2.]])
    depth=ref*np.array([[.7],[-.2]])
    assert np.allclose(depth_orthogonal_component(depth,ref),0)
    signal=np.array([[3.,1.,0.],[0.,4.,1.]])
    q=depth_orthogonal_component(signal,ref)
    assert np.allclose(np.sum(q*ref,axis=1),0)
    assert np.linalg.norm(q)>0


def test_co_thin_limit_integrated_line_strength():
    lines=COLines(np.array([4300.]),np.array([1e-22]),np.array([.07]),np.array([.65]),1.,1.)
    nu=np.linspace(4299.5,4300.5,10001)
    absorption=co_absorption(nu,lines,640.,0.,1e14)
    assert abs(np.trapz(absorption,nu)/1e14/1e-22-1)<1e-5
    thicker=co_absorption(nu,lines,640.,0.,1e20)
    assert np.all(thicker>=absorption)
    assert np.max(thicker)<=1
