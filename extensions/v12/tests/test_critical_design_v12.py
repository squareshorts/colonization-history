"""Regression and integrity checks for the analytic manuscript inputs."""
import csv, json, re, sys
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'analysis'))
import critical_design_v12 as A
import archival_within_estimator_v12 as B

def test_missing_outcome_bounds_with_two_arm_dropout():
    x=A.simulate(np.random.default_rng(10),'both_uniform_equal',.20,.85,'logistic_exponential',studies=300,blocks=12,founders=48)
    assert x['exact_bounds_cover_generated_full_data']
    assert .185 < x['actual_pulse_loss'] < .215
    assert .185 < x['actual_sham_loss'] < .215

def test_sham_half_dropout():
    x=A.simulate(np.random.default_rng(11),'sham_half_uniform',.40,.85,'cloglog_weibull',studies=240,blocks=10,founders=48)
    assert x['exact_bounds_cover_generated_full_data']
    assert .385<x['actual_pulse_loss']<.415
    assert .185<x['actual_sham_loss']<.215

def test_phenotype_recovery_nominal_not_realized():
    x=A.simulate(np.random.default_rng(12),'pulse_only_phenotype',.4,0.,'logistic_exponential',studies=250,blocks=20,founders=48)
    assert x['exact_bounds_cover_generated_full_data']
    assert .25<x['actual_pulse_loss']<.32
    assert abs(x['actual_sham_loss'])<1e-12

def test_sham_more_loss_case_not_trimming():
    x=A.simulate(np.random.default_rng(13),'sham_more_uniform',.4,.85,'logistic_exponential',studies=150,blocks=12,founders=48)
    assert x['actual_sham_loss'] > x['actual_pulse_loss']
    assert x['exact_bounds_cover_generated_full_data']

def test_paired_contact_same_single_contact_mean():
    f,p,sel,hom=A.paired_probabilities(14.,.1,.03)
    assert np.all(sel>0) and np.all(hom>0)
    assert abs(sum(sel)-1)<1e-12 and abs(sum(hom)-1)<1e-12
    assert abs(sel[2]+sel[3]-p)<1e-12
    assert abs(hom[2]+hom[3]-p)<1e-12
    assert sel[3]-p*p>0

def test_repeated_contact_realistic_requirements():
    for hold,target in ((0,22633507),(14,459080),(37,14901)):
        _,_,sel,hom=A.paired_probabilities(hold,.1,.03)
        n,*_=A.normal_power_sample_size(sel,hom)
        assert abs(n-target)/target<.001

def test_empirical_nonmonotone_within_estimator():
    x=B.run();r=x['same_estimator_normalized_advantage']
    assert 1.13<r['14h']<1.15
    assert .27<r['37h']<.30
    assert x['biological_interval_for_holding_time_response'].startswith('not identifiable')
    assert all(i<0 for i in x['technical_window_adversarial_paired_difference_h'])

def test_coincident_technical_minima_are_different_windows():
    x=json.loads((BASE/'results/within_estimator_lag_audit.json').read_text())['technical_minima_coincidence']
    assert round(x['FC1-1']['lag_h'],2)==round(x['FC2-2']['lag_h'],2)==4.97
    assert x['FC1-1']['fit_cell_eq'] != x['FC2-2']['fit_cell_eq']

def test_published_results_grid_and_bounds():
    r=json.loads((BASE/'results/v12_results_summary.json').read_text())
    assert r['all_generated_true_effects_inside_observed_missing_bounds']
    assert r['designs']==280
    assert r['studies_per_design']==1200
    assert r['max_null_bounded_joint'] < .01
    dat=list(csv.DictReader((BASE/'results/v12_two_arm_loss_grid.csv').open()))
    def find(s,k,e,l):
        return next(x for x in dat if x['scenario']==s and x['kinetic']==k and float(x['effect'])==e and float(x['nominal_loss'])==l)
    assert float(find('both_uniform_equal','logistic_exponential',.85,.1)['bounded_joint']) < .1
    assert float(find('pulse_only_uniform','logistic_exponential',.85,.1)['bounded_joint']) > .9

def test_manuscript_no_internal_workflow_language():
    src='\n'.join((BASE/'overleaf'/p).read_text() for p in ['main.tex','supplement.tex'])
    for pattern in [r'v[0-9]+\s*[-–]\s*v[0-9]+',r'\bOverleaf\b',r'\bhostile[- ]review\b',r'\b(reviewer|referee|submission package)\b',r'\bTODO\b',r'\bPLACEHOLDER\b',r'\btranche\b',r'\bfrozen protocol\b']:
        assert not re.search(pattern,src,re.I),pattern

def test_bibliography_encoding_fixed():
    b=(BASE/'overleaf/references.bib').read_text()
    assert 'Delbr{\\"u}ck' in b
    assert 'str{\\"o}m' in b
    assert 'Klingm{\\"u}ller' in b

def test_abstract_and_figures():
    tex=(BASE/'overleaf/main.tex').read_text()
    abstract=tex.split('\\begin{abstract}',1)[1].split('\\end{abstract}',1)[0]
    cleaned=re.sub(r'\\[a-zA-Z]+','',abstract).replace('\\%','%')
    assert len(cleaned.split())<=200
    for name in ('main.tex','supplement.tex'):
        doc=(BASE/'overleaf'/name).read_text()
        graphics=re.findall(r'\\includegraphics\[[^]]*\]\{([^}]+)\}',doc)
        assert graphics
        for g in graphics:assert (BASE/'overleaf/figures'/g).is_file(),g
    assert '\\label{fig:contacts}' in tex
    assert '\\ref{fig:contacts}' in tex

