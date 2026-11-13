import numpy as np
from xds.nonlinear_atmosphere import CoherentCOModel,fit_coherent,physical_parameters,temporal_basis


def toy_model(n_controls=3):
    axes=[np.linspace(-1,1,5)]*3;time=np.linspace(-1.8,1.8,5)
    a,b,c=np.meshgrid(*axes,indexing='ij')
    library=np.zeros((5,5,5,2,5,2,9))
    for g in range(2):
        for k in range(5):
            library[:,:,:,g,k,0,3*g]=1+a+.1*a*a
            library[:,:,:,g,k,0,3*g+1]=1+b+.1*a*b
            library[:,:,:,g,k,0,3*g+2]=1+c+.1*c*c
            library[:,:,:,g,k,1,6+g]=-.1*(1+.1*a+.1*b)*(1+.02*k)
    return CoherentCOModel(axes,library,time,n_controls)


def test_coherent_trajectories_preserve_bounds():
    t=np.linspace(-2,2,101)
    for n in [1,3,5]:
        b=temporal_basis(t,n)
        assert np.all(b>=-1e-15) and np.allclose(b.sum(axis=1),1)
        controls=np.linspace(-1,1,n)
        assert np.max(abs(b@controls))<=1+1e-15
    temp,abundance,wind=physical_parameters(np.array([[-1,-1,-1],[1,1,1]]))
    assert np.allclose(temp,[560,720])
    assert np.allclose(abundance,[1e-4/3,3e-4])
    assert np.allclose(wind,[-1.5,1.5])


def test_nonlinear_analytic_jacobian_and_blocking_dependence():
    m=toy_model();x=np.linspace(-.3,.3,m.n_parameters)
    prediction,jac,ds=m.evaluate(x,.4)
    h=1e-5
    for i in range(len(x)):
        e=np.zeros_like(x);e[i]=h
        numerical=(m.evaluate(x+e,.4)[0]-m.evaluate(x-e,.4)[0])/(2*h)
        assert np.allclose(jac[...,i],numerical,atol=1e-9)
    assert np.allclose(m.evaluate(x,1.)[0]-m.evaluate(x,0.)[0],ds)
    assert np.allclose(prediction,m.evaluate(x,0.)[0]+.4*ds)
    assert not np.allclose(ds,m.evaluate(np.zeros_like(x),.4)[2])


def test_bounded_nonlinear_recovery_and_null_nestedness():
    m=toy_model(1);x=np.array([.1,-.2,.15,-.1,.2,-.15])
    truth=m.evaluate(x,.6)[0]
    null,_=fit_coherent(m,truth,.01,[np.zeros(6)])
    joint,runs=fit_coherent(m,truth,.01,[np.zeros(6),np.full(6,.4)],with_moon=True)
    assert all(r['success'] for r in runs)
    assert joint['chi2']<1e-12 and joint['chi2']<=null['chi2']
    assert abs(joint['blocking']-.6)<1e-6
    fixed,_=fit_coherent(m,truth,.01,[np.zeros(6)],fixed_blocking=.6)
    assert fixed['chi2']<1e-12 and fixed['blocking']==.6
