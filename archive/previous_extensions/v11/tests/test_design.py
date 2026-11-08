from pathlib import Path
import sys,re,json
import numpy as np
from scipy.special import expit,logit
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from attrition_power_and_observables import paired_probabilities,required_power_normal,sim_experiment
from archival_lag_audit import audit

def test_exact_equivalence():
    f0=.990359;c=.1456747;m0=.85;tau=10.
    nu=c*tau;q=(1-f0)/f0;K=m0*q**(1/nu)
    tt=np.linspace(0,60,1703)
    f=expit(logit(f0)-c*tt)
    m=m0*np.exp(-tt/tau)
    h=1/(1+(K/m)**nu)
    assert np.max(np.abs(h-f))<2e-15

def test_positive_exponential_constraint():
    d0=17.56-4.68;d14=17.56-5.46;d37=17.56-13.41
    assert d37/d0 < (d14/d0)**(37/14)
    rng=np.random.default_rng(31)
    for _ in range(100):
        lam=rng.uniform(0,.5,5);weights=rng.dirichlet(np.ones(5))
        a=lambda t: np.sum(weights*np.exp(-lam*t))
        assert a(37)+1e-12>=a(14)**(37/14)

def test_paired_contact_marginals_and_covariance():
    for h in [0,14,37]:
        f,p,sel,hom=paired_probabilities(h)
        assert np.isclose(sel.sum(),1.) and np.isclose(hom.sum(),1.)
        assert np.isclose(sel[2]+sel[3],p)
        assert np.isclose(hom[2]+hom[3],p)
        cov=sel[3]-p*p
        assert np.isclose(cov,f*(1-f)*(.75-.15)**2)
        assert cov>0

def test_information_increases_away_from_pure_phenotype():
    vals=[]
    for h in [0,14,37]:
        _,_,a,b=paired_probabilities(h)
        _,d,_,_,_=required_power_normal(a,b)
        vals.append(d)
    assert vals[0]<vals[1]<vals[2]

def test_sharp_missing_assignment_bound():
    for N in range(2,10):
        for n in range(0,N+1):
            for S in range(n+1):
                low=S/N;high=(S+N-n)/N
                vals=[(S+m)/N for m in range(N-n+1)]
                assert min(vals)==low and max(vals)==high

def test_power_drops_when_many_cells_missing():
    a=sim_experiment(np.random.default_rng(12),Nstudies=350,effect_logodds=.85,loss=0,mode='uniform')
    b=sim_experiment(np.random.default_rng(12),Nstudies=350,effect_logodds=.85,loss=.4,mode='uniform')
    assert a['manski_joint']>.9 and b['manski_joint']<.15

def test_bounds_reduce_spurious_recovery_declarations():
    a=sim_experiment(np.random.default_rng(2026),Nstudies=800,effect_logodds=0,loss=.4,mode='phenotype_biased')
    assert a['raw_joint']>a['manski_joint']
    assert a['manski_joint']<.025

def test_independent_kinetic_robustness():
    a=sim_experiment(np.random.default_rng(8),kinetic='logistic_exponential',Nstudies=350,effect_logodds=.85,loss=.1)
    b=sim_experiment(np.random.default_rng(8),kinetic='cloglog_weibull',Nstudies=350,effect_logodds=.85,loss=.1)
    assert a['manski_joint']>.85 and b['manski_joint']>.85

def test_archive_unit_count():
    r=audit(ROOT/'results/archival_occupancy_validation_input.csv')
    assert r['n_independent_pairs']==3
    assert r['hold_independent_biological_replicates']=={'Tr14':1,'Tr37':1}
    assert np.isclose(r['mean_paired_difference_h'],-10.32666666666)
    assert not r['identifiable_between_replicate_washout_CI']

def test_figures_and_source():
    if not (ROOT/'overleaf/main.tex').exists():
        # The public code archive intentionally excludes the private manuscript.
        assert (ROOT/'analysis/attrition_power_and_observables.py').exists()
        return
    main=(ROOT/'overleaf/main.tex').read_text()
    sup=(ROOT/'overleaf/supplement.tex').read_text()
    a=re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',main,re.S).group(1)
    a=re.sub(r'\\[a-zA-Z]+(?:\{.*?\})?','',a)
    assert len(re.findall(r"\b[\w'-]+\b",a))<=200
    for word in ['v7','v8','v9','v10','v11','overleaf','manuscript source','reviewer','hostile','extensions/','frozen protocol','tranche']:
        assert word not in main.lower(),word
        assert word not in sup.lower(),word
    for f in re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}',main+'\n'+sup):
        assert (ROOT/'overleaf/figures'/f).exists(),f
    assert '\\renewcommand*{\\bibfont}' in main

def test_simulation_csv_summary():
    import pandas as pd
    d=pd.read_csv(ROOT/'results/attrition_power_frontier.csv')
    assert len(d)==96
    row=d[(d.kinetic=='logistic_exponential')&(d['mode']=='uniform')&(d.effect_logodds==.85)&(d.requested_loss==.4)].iloc[0]
    assert np.isclose(row.manski_joint,.030)
    r=json.loads((ROOT/'results/design_results_summary.json').read_text())
    assert r['attrition_designs']==96
