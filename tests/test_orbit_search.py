import numpy as np
from xds.orbit_search import bounded_improvement,empirical_tail,binomial_interval
from xds.moving_shadow import separation_from_period,REFERENCE_PERIOD,period_limits,stellar_moments
from xds.dynamics import moon_period_hours


def test_bounded_score_equals_direct_scalar_least_squares():
    rng=np.random.default_rng(15)
    for _ in range(20):
        y=rng.normal(size=50);s=rng.normal(size=50)
        c=y@s;d=s@s;stat,f=bounded_improvement(c,d)
        assert 0<=f<=1
        assert np.isclose(stat,np.sum(y*y)-np.sum((y-f*s)**2))
    assert bounded_improvement(-2.,3.)[0]==0
    assert bounded_improvement(2.,0.)[0]==0


def test_empirical_p_is_nonzero_and_binomial_limits():
    calibration=np.arange(128.)
    assert empirical_tail(calibration,1000.)==1/129
    assert empirical_tail(calibration,-1.)==1.
    assert binomial_interval(0,64)[0]==0 and binomial_interval(64,64)[1]==1
    lo,hi=binomial_interval(3,64)
    assert lo<3/64<hi


def test_blind_period_changes_separation_consistently():
    p=np.geomspace(*period_limits(),20)
    a=separation_from_period(p)
    assert np.allclose(moon_period_hours(a,.7,1.1),p)
    assert np.isclose(separation_from_period(REFERENCE_PERIOD),3.)


def test_stellar_moments_against_direct_quadratic_limb_darkening():
    rp=1.1*71492.;rstar=695700/rp;axis=(30/365.25)**(2/3)*149597870.7/rp
    for t in [-2.,0.,2.]:
        for center in [-1.1,1.1]:
            coeff=stellar_moments(t,center,rp)
            for du,z in [(.07,.04),(-.06,.02)]:
                mu=np.sqrt(1-((center+du+axis*np.sin(2*np.pi*t/720))**2+z*z)/rstar**2)
                exact=(1-.4*(1-mu)-.25*(1-mu)**2)/(1-.4/3-.25/6)
                approximate=coeff@np.array([1,du/.09,(du/.09)**2,(z/.09)**2])
                assert abs(exact-approximate)<1e-5


def test_summary_separates_calibration_and_uses_conservative_consensus():
    import csv
    import sys
    import tempfile
    from pathlib import Path
    scripts=Path(__file__).resolve().parents[1]/'scripts'
    sys.path.insert(0,str(scripts))
    try:
        from blind_co_recovery import summarize,FAMILIES,INJECTIONS
    finally:sys.path.remove(str(scripts))
    records=[]
    for stream,labels,count in [('calibration',FAMILIES,8),('validation',FAMILIES,1),('injection',INJECTIONS,1)]:
        for label in labels:
            for draw in range(count):
                for n in [3,5]:
                    value=draw*(2 if n==5 else 1) if stream=='calibration' else 8.
                    records.append(dict(stream=stream,family=label,draw=draw,n_controls=n,statistic=value,
                        oracle_statistic=value,period=22.,phase=.746,blocking=.8,true_period=21.,true_phase=.7473,
                        null_start_spread=0.,refinement_failures=0))
    observed=[dict(case='test',n_controls=n,statistic=8.,oracle_statistic=8.) for n in [3,5]]
    with tempfile.TemporaryDirectory() as name:
        root=Path(name);(root/'docs').mkdir();folder=root/'raw';folder.mkdir()
        (folder/'candidate_fits.json').write_text('[]')
        summarize(root,folder,records,observed,dict(alpha=.25,n_calibration=8,n_periods=17,n_phases=128),'test')
        with (root/'docs/co_blind_test_decisions.csv').open() as f:rows=list(csv.DictReader(f))
        assert all(r['stream']!='calibration' for r in rows)
        assert all(r['flag_consensus']=='False' for r in rows)
        assert any(r['flag_robust']=='True' for r in rows)
        assert all(float(r['p_consensus'])>=float(r['p_robust']) for r in rows)


def test_sparse_bank_score_matches_full_space_projection_and_phase_wrap():
    from xds.orbit_search import OrbitBank
    rng=np.random.default_rng(71);time=np.linspace(-1.5,1.5,4)
    rotation=np.linalg.qr(rng.normal(size=(3,2)))[0]
    states=rng.normal(size=(2,4,2,2))
    class ToyShadow:
        library=np.zeros((1,2,2))
        def atmospheric_states(self,p):return states,None
        def orbit_weights(self,periods,phases,time,order=16):
            periods=np.atleast_1d(periods);phases=np.atleast_1d(phases)
            oi=np.repeat(np.arange(len(periods)),len(time));ti=np.tile(np.arange(len(time)),len(periods))
            a=2*np.pi*(time[ti]/periods[oi]+phases[oi])
            return oi,ti,np.zeros(len(oi),dtype=int),np.column_stack([np.cos(a),np.sin(a)])
    shadow=ToyShadow();shadow.rotation=rotation
    bank=OrbitBank(shadow,time,3,4)
    u=np.linalg.qr(rng.normal(size=(12,2)))[0];residual=rng.normal(size=(4,3))
    null=dict(parameters=np.zeros((2,4,3)),u=u,residual=residual)
    condition=bank.conditional(null,.7);scores,f=bank.score_grid(condition)
    for i,(p,phase) in enumerate(zip(bank.periods,bank.phases)):
        value,fraction=bank.score_point(condition,p,phase)
        assert np.isclose(value,scores[i]) and np.isclose(fraction,f[i])
        assert np.isclose(value,bank.score_point(condition,p,phase+1)[0])
        _,ti,_,w=shadow.orbit_weights([p],[phase],time)
        full=(np.einsum('tk,tkr->tr',w,states[0])@rotation.T/.7).ravel()
        q=full-u@(u.T@full)
        exact,_=bounded_improvement(residual.ravel()@q,q@q)
        assert np.isclose(value,exact)
    found,_=bank.search(null,.7)
    assert found['statistic']>=scores.max()-1e-10
