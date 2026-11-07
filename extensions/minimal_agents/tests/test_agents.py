"""Deterministic and falsifiability checks for exploratory agent model."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'src'))
import unittest, numpy as np
from dataclasses import replace
from minimal_agents import (Parameters,CAPABILITIES,run_field,initial_memory,
                            daughter_transfer_probability,sigmoid_hill)
from oscillatory_signatures import simulate,corr_vs_lag

class ModelChecks(unittest.TestCase):
    def setUp(self):
        self.p=Parameters(horizon_h=12.)
    def test_state_zero_without_prior_exposure(self):
        self.assertEqual(initial_memory(0,self.p,CAPABILITIES['passive']),0.0)
    def test_exposed_state_fades_in_suspension(self):
        c=CAPABILITIES['cell_memory']
        self.assertGreater(initial_memory(0,self.p,c),initial_memory(14,self.p,c))
        self.assertGreater(initial_memory(14,self.p,c),initial_memory(37,self.p,c))
    def test_no_decay_maintains_transfer_state(self):
        c=CAPABILITIES['cell_no_decay']
        self.assertEqual(initial_memory(0,self.p,c),initial_memory(37,self.p,c))
    def test_passive_no_history_signature(self):
        c=CAPABILITIES['passive']
        a=run_field(self.p,c,0.,110,transfer=False)
        b=run_field(self.p,c,37.,110,transfer=True)
        self.assertAlmostEqual(a['restricted_threshold_h'],b['restricted_threshold_h'])
        self.assertEqual(a['n_at_8h'],b['n_at_8h'])
        self.assertEqual(a['short_fraction_all_contacts'],b['short_fraction_all_contacts'])
    def test_environment_only_no_transfer_on_fresh(self):
        c=CAPABILITIES['environment']
        a=run_field(self.p,c,0.,117,transfer=False)
        b=run_field(self.p,c,37.,117,transfer=True)
        self.assertAlmostEqual(a['restricted_threshold_h'],b['restricted_threshold_h'])
        self.assertEqual(a['n_at_8h'],b['n_at_8h'])
        self.assertEqual(a['short_fraction_all_contacts'],b['short_fraction_all_contacts'])
    def test_results_reproducible(self):
        c=CAPABILITIES['cell_memory']
        self.assertEqual(run_field(self.p,c,14.,122),run_field(self.p,c,14.,122))
    def test_capacity_and_range(self):
        c=CAPABILITIES['cell_memory']
        x=run_field(self.p,c,0.,99)
        self.assertLessEqual(x['final_count'],self.p.carrying_capacity)
        self.assertGreaterEqual(x['final_count'],0)
        self.assertTrue(0<=x['short_fraction_all_contacts']<=1)
    def test_early_stopping_does_not_change_lag(self):
        p=Parameters(horizon_h=18.)
        c=CAPABILITIES['cell_memory']
        a=run_field(p,c,0.,124)
        b=run_field(replace(p,stop_at_lag_and_8h=False),c,0.,124)
        self.assertEqual(a['lag_threshold_h'],b['lag_threshold_h'])
        self.assertEqual(a['n_at_8h'],b['n_at_8h'])
    def test_fresh_daughter_transfer_distinguishes_inheritance(self):
        p=self.p
        a=daughter_transfer_probability(p,CAPABILITIES['cell_memory'],.8)
        b=daughter_transfer_probability(p,CAPABILITIES['cell_no_inheritance'],.8)
        self.assertGreater(a,b*10)
        self.assertAlmostEqual(b,p.stable_attach_naive)
    def test_bounded_hill(self):
        xs=np.linspace(0,1,100)
        y=sigmoid_hill(xs,self.p)
        self.assertTrue(np.all((y>=0)&(y<=1)))
        self.assertTrue(np.all(np.diff(y)>=0))
    def test_oscillatory_signature_requires_additional_state_in_tested_family(self):
        x,y,dt=simulate('scalar',n_agents=100,seed=99213)
        corr=corr_vs_lag(x,y,dt)
        a=float(corr[5.0].mean())
        x,y,dt=simulate('oscillator',n_agents=100,seed=99213)
        corr=corr_vs_lag(x,y,dt)
        b=float(corr[5.0].mean())
        self.assertLess(a,0)
        self.assertGreater(b,0)

if __name__=='__main__':unittest.main()
