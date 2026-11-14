import numpy as np
from xds.observing_contamination import CountBudget,electron_budget,time_whiten,spectral_metric


def test_shared_reference_whitening():
    rng=np.random.default_rng(12);n=7;nr=3
    covariance=np.eye(n)+np.ones((n,n))/nr
    w=time_whiten(np.eye(n),nr)
    assert np.allclose(w@covariance@w.T,np.eye(n),atol=1e-14)
    x=rng.normal(size=(n,4))
    assert np.allclose(np.sum(time_whiten(x,nr)**2),np.sum(x*np.linalg.solve(covariance,x)))


def test_count_budget_flux_and_area_scaling():
    a=np.ones(8);w=np.linspace(2.3,2.35,8)
    first=electron_budget(w,a,a,0*a,CountBudget(area_m2=49))
    second=electron_budget(w,a,a,0*a,CountBudget(area_m2=98))
    faint=electron_budget(w,a,a,0*a,CountBudget(k_magnitude=8.5))
    assert np.allclose(second['source'],2*first['source'])
    assert np.allclose(faint['source'],first['source']/10)
    assert np.all(first['variance']>first['source'])
    assert np.all(first['sigma']>0)


def test_spectral_metric_preserves_likelihood_and_removes_nuisance():
    rng=np.random.default_rng(4);b=np.linalg.qr(rng.normal(size=(50,5)))[0]
    sigma=np.linspace(1,3,50);keep=np.arange(50)%5!=0
    nuisance=rng.normal(size=(2,50))
    q,transform=spectral_metric(b,sigma,keep,nuisance)
    assert np.allclose(q.T@q,np.eye(5))
    assert np.max(abs((nuisance[:,keep]/sigma[keep])@q))<1e-13
    assert np.allclose((b[keep]/sigma[keep,None]).T@q,transform)
    coeff=rng.normal(size=5)
    assert np.allclose((b[keep]@coeff/sigma[keep])@q,coeff@transform)


def test_reference_metric_sparse_scores_match_full_space():
    from xds.observing_contamination import ReferenceOrbitBank
    from xds.orbit_search import bounded_improvement
    class Shadow:
        library=np.zeros((1,1,1,2,3))
        rotation=np.eye(3)
        def orbit_weights(self,p,phase,t,order):
            return np.array([0,0]),np.array([0,2]),np.array([0,0]),np.array([[1.,.2],[.3,.8]])
    bank=object.__new__(ReferenceOrbitBank)
    bank.shadow=Shadow();bank.time=np.arange(4);bank.order=16;bank.n_reference=2;bank.n_templates=1
    bank.oi=np.array([0,0]);bank.ti=np.array([0,2]);bank.g=np.array([0,0])
    bank.weights=np.array([[1.,.2],[.3,.8]])
    bank.groups=[(0,t,np.flatnonzero(bank.ti==t)) for t in range(4)]
    rng=np.random.default_rng(5);states=rng.normal(size=(1,4,2,3))
    original=rng.normal(size=(4,3));whitened=time_whiten(original,2)
    u=np.linalg.qr(rng.normal(size=(12,2)))[0]
    condition=dict(states=states,ur=time_whiten(u.reshape(4,3,2),2),
        residual=time_whiten(whitened,2),utr=u.T@whitened.ravel(),noise=1.)
    values=np.zeros((4,3));values[0]=bank.weights[0]@states[0,0];values[2]=bank.weights[1]@states[0,2]
    v=time_whiten(values,2).ravel();q=v-u@(u.T@v)
    expected=bounded_improvement(q@whitened.ravel(),q@q)
    point=bank.score_point(condition,1.,0.)
    grid=bank.score_grid(condition)
    assert np.allclose(point,expected)
    assert np.allclose([grid[0][0],grid[1][0]],expected)


def test_contaminant_calibration_is_separate_and_conservative():
    import csv
    import sys
    import tempfile
    from pathlib import Path
    scripts=Path(__file__).resolve().parents[1]/'scripts'
    sys.path.insert(0,str(scripts))
    try:
        from observing_contamination_pilot import summarize,FAMILIES
    finally:sys.path.remove(str(scripts))
    records=[]
    for method in ['unprotected','protected']:
        for fi,family in enumerate(FAMILIES):
            for draw in range(64):
                records.append(dict(level='control',method=method,stream='calibration',family=family,
                    draw=draw,statistic=float(draw+fi)))
        records.append(dict(level='control',method=method,stream='validation',family='combined',draw=0,statistic=64.))
    with tempfile.TemporaryDirectory() as folder:
        root=Path(folder);(root/'docs').mkdir()
        summarize(root,records,'test')
        with (root/'docs/co_observing_test_decisions.csv').open() as handle:rows=list(csv.DictReader(handle))
        assert len(rows)==2 and all(r['stream']=='validation' for r in rows)
        assert all(r['flag_atmosphere']=='True' and r['flag_robust']=='False' for r in rows)
        assert all(float(r['p_robust'])>=float(r['p_atmosphere']) for r in rows)
