import numpy as np
from xds.molecular import (COLines,co_absorption,bin_average,depth_orthogonal_component,
                          nuisance_basis,project_nuisance,cubic_shifted_bin_average)


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


def test_nuisance_projection_rank_units_and_depth_equivalence():
    b=np.array([[1.,2.,3.,4.],[4.,3.,2.,1.]])
    s=np.array([[3.,1.,0.,2.],[0.,4.,1.,2.]])
    basis,rank=nuisance_basis([b,1000*b,np.zeros_like(b)])
    assert np.all(rank==1)
    q=project_nuisance(s,basis)
    assert np.allclose(q,depth_orthogonal_component(s,b))
    assert np.allclose(project_nuisance(q,basis),q)
    assert np.allclose(project_nuisance(b,basis),0)


def test_nuisance_nested_subspace_and_gaussian_recovery():
    rng=np.random.default_rng(29)
    modes=rng.normal(size=(3,2,30));s=rng.normal(size=(2,30))
    small,_=nuisance_basis(modes[:1]);big,rank=nuisance_basis(modes)
    q=project_nuisance(s,big)
    assert np.all(rank==3)
    assert np.linalg.norm(q)<=np.linalg.norm(project_nuisance(s,small))
    assert np.allclose(np.einsum('tv,mtv->mt',q,modes),0,atol=1e-12)
    noise=rng.normal(size=(2000,2,30))
    z=np.einsum('ntv,tv->n',noise,q)/np.linalg.norm(q)
    assert abs(np.mean(z))<.08 and abs(np.std(z)-1)<.08
    y=2*s+4*modes[0]-3*modes[1]
    amplitude=np.sum(y*q)/np.sum(q*q)
    assert np.isclose(amplitude,2)


def test_cubic_shifted_bins_exact_polynomial_and_wind_derivative():
    x=np.linspace(-3,3,61);edges=np.array([-1.,-.3,.2,1.]);shifts=np.array([-.1,0.,.2])
    values=cubic_shifted_bin_average(x,x*x+2*x+3,edges,shifts)
    lo=edges[:-1][None,:]-shifts[:,None];hi=edges[1:][None,:]-shifts[:,None]
    expected=(hi**3-lo**3)/(3*(hi-lo))+hi+lo+3
    assert np.allclose(values,expected,atol=1e-12)
    for step in [.1,.05]:
        plus=cubic_shifted_bin_average(x,x*x+2*x+3,edges,np.array([step]))
        minus=cubic_shifted_bin_average(x,x*x+2*x+3,edges,np.array([-step]))
        zero=cubic_shifted_bin_average(x,x*x+2*x+3,edges,np.array([0.]))
        assert np.allclose((plus+minus-2*zero)/step**2,2,atol=1e-10)
