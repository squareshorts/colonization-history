import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from discriminating_test_v9 import (analytical_equivalence,design_rank,five_carrier_rank,
    calibration_geometry,observable_probability,dataset,summarize)

class CriticalPredictionTest(unittest.TestCase):
    def test_macro_equivalence_across_all_hold_times(self):
        s,c=analytical_equivalence()
        self.assertLess(s['lag_curve_max_abs_difference_h_0_to_60h'],1e-12)
        np.testing.assert_allclose(c['selection'],c['internal'],atol=1e-12)
    def test_lee_published_points_matched_by_both(self):
        s,_=analytical_equivalence()
        for i,v in [(0,4.68),(2,5.46),(5,13.41)]:
            self.assertAlmostEqual(s['lags_selection_h'][i],v,10)
            self.assertAlmostEqual(s['lags_intracellular_h'][i],v,10)
    def test_global_observational_equivalence_without_data_fitting(self):
        rng=np.random.default_rng(202611083)
        f0=rng.uniform(.05,.95,100)
        c=rng.uniform(.04,.25,100)
        tau=rng.uniform(6.,25.,100)
        m0=rng.uniform(.4,.95,100)
        nu=c*tau
        q=(1-f0)/f0
        K=m0*np.power(q,1/nu)
        times=np.linspace(0,60,60)
        sel=1/(1+q[:,None]*np.exp(c[:,None]*times))
        m=m0[:,None]*np.exp(-times/tau[:,None])
        hh=1/(1+np.power(K[:,None]/m,nu[:,None]))
        np.testing.assert_allclose(sel,hh,rtol=1e-12,atol=1e-12)
    def test_hill_rate_tautology_explicit(self):
        s,_=analytical_equivalence()
        self.assertAlmostEqual(s['hill_exponent']/s['intracellular_state_tau_h'],
                               s['relative_growth_disadvantage_per_h'],12)
    def test_two_extra_contrasts_are_minimal_in_ideal_design(self):
        d=design_rank()
        self.assertEqual(d['rank_pre'],2)
        self.assertEqual(d['rank_full'],4)
        self.assertEqual(d['rank_without_pulse'],3)
        self.assertEqual(d['rank_without_washed_surface'],3)
    def test_fifth_carrier_defeats_identification(self):
        d=five_carrier_rank()
        self.assertEqual(d['rank_with_sister_and_surface_controls'],4)
        self.assertEqual(d['rank_with_one_selective_localization_assay'],5)
        x=np.array([[1,1,0,0,1],[0,0,1,1,0],[1,0,0,0,1],[0,0,1,0,0],[1,0,0,0,1]])
        np.testing.assert_allclose(x @ np.array(d['null_direction']),np.zeros(5))
    def test_calibration_ill_conditioning(self):
        r=calibration_geometry()
        nominal=next(x for x in r if x['selection_leakage_kappa']==0 and x['trace_retention_q']==1)
        badly=next(x for x in r if x['selection_leakage_kappa']==.99 and x['trace_retention_q']==.6)
        self.assertLess(nominal['matrix_condition_number'],3)
        self.assertGreater(badly['matrix_condition_number'],100)
    def test_probability_simulation_matches_kinetics(self):
        # Independently simulate attachment times and exponential detach waits.
        probs=observable_probability(np.array([[np.log(.04/.96)]]),np.array([[.22]]),
               np.array([[1.25]]),np.array([[1.1]]),np.array([[2.6]]),
               np.array([[1.2]]),np.array([[.6]]))[0,0]
        rng=np.random.default_rng(70677)
        N=200000
        z=rng.random(N)<.22
        pa=1/(1+np.exp(-(np.log(.04/.96)+1.25*z+1.1)))
        arrivals=rng.random(N)<pa
        hazard=2.6*np.exp(-1.2*1.1-.6*z)
        off=rng.exponential(1/hazard)
        empirical=np.mean(arrivals & (off>5/60))
        self.assertLess(abs(empirical-probs),.003)
    def test_selection_null_rate_and_mother_daughter_inheritance(self):
        n=1200; k=32
        a=summarize(dataset('intracellular_inherited',n,k,48,912))
        b=summarize(dataset('intracellular_no_inheritance',n,k,48,913))
        c=summarize(dataset('immutable_selection',n,k,48,914))
        self.assertGreater(a['joint_positive_fraction'],.93)
        self.assertGreater(b['mother_positive_fraction'],.95)
        self.assertLess(b['daughter_positive_fraction'],.055)
        self.assertLess(c['joint_positive_fraction'],.02)
        self.assertLess(c['mother_positive_fraction'],.055)
        self.assertLess(c['daughter_positive_fraction'],.055)
    def test_contaminating_selection_creates_false_inheritance_signature(self):
        a=summarize(dataset('selection_with_loss',1200,32,48,1002,selection_survival_odds=1))
        b=summarize(dataset('selection_with_loss',1200,32,48,1003,selection_survival_odds=8))
        self.assertLess(a['joint_positive_fraction'],.03)
        self.assertGreater(b['joint_positive_fraction'],.65)
    def test_cell_coating_can_mimic_inheritance(self):
        a=summarize(dataset('intracellular_inherited',1200,32,48,315))
        b=summarize(dataset('cell_coating',1200,32,48,316))
        self.assertLess(abs(a['joint_positive_fraction']-b['joint_positive_fraction']),.04)
    def test_washout_of_microscopic_predicted_effect(self):
        a=summarize(dataset('intracellular_inherited',1000,32,48,1234,hold_h=0))
        b=summarize(dataset('intracellular_inherited',1000,32,48,1235,hold_h=37))
        self.assertGreater(a['joint_positive_fraction'],.9)
        self.assertLess(b['joint_positive_fraction'],.15)
    def test_seed_reproducibility(self):
        a=dataset('intracellular_inherited',10,8,48,9123)
        b=dataset('intracellular_inherited',10,8,48,9123)
        np.testing.assert_array_equal(a['mother_effect'],b['mother_effect'])
        np.testing.assert_array_equal(a['daughter_detected'],b['daughter_detected'])

if __name__=='__main__':unittest.main()
