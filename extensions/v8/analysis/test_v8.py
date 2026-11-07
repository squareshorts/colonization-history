import sys,unittest,json,csv
from pathlib import Path
from math import isclose
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from hostile_checks_v8 import model_audit,lhs_values
from selection_counterexample_v8 import f,lag,model
from lag_robustness_v8 import margin

class HostileAuditChecks(unittest.TestCase):
    def test_source_values_and_holdout(self):
        a=model_audit(); self.assertAlmostEqual(a['linear_exponential']['prediction_37_h'],6.640221060438822,places=9)
        self.assertAlmostEqual(a['observations']['37h'],13.41)
    def test_logconvexity_counterexample(self):
        a=model_audit();r=a['nonnegative_exponential_mixture']
        self.assertLess(r['observed_normalized_37'],r['minimum_possible_normalized_37_from_normalized_14'])
    def test_selection_population_interpolation(self):
        for t,expected in [(0,4.68),(14,5.46),(37,13.41)]:self.assertAlmostEqual(lag(t),expected,places=9)
        self.assertGreater(f(0),f(14));self.assertGreater(f(14),f(37)); self.assertGreater(model['relative_suspension_fitness_cost_per_hour'],0)
    def test_worstcase_one_hour_perturbations(self):
        a=json.loads((ROOT/'results'/'v8'/'lag_uncertainty_robustness.json').read_text())
        by={x['each_lag_perturbed_within_plusminus_h']:x for x in a['perturbations']}
        self.assertLess(by[1.]['maximum_attainable_margin'],0)
        self.assertGreater(by[2.]['maximum_attainable_margin'],0)
    def test_priors_independent_and_reproducible(self):
        p=lhs_values(n=5,seed=42);q=lhs_values(n=5,seed=42)
        self.assertEqual(p,q); self.assertTrue(all(4<=x.free_memory_tau_h<=35 for x in p))
    def test_broad_sweep_has_each_mechanism(self):
        a=json.loads((ROOT/'results'/'v8'/'agent_global_prior_sweep.json').read_text())
        self.assertEqual(a['n_parameter_draws'],160)
        for v in a['model'].values():self.assertEqual(v['parameter_draws'],160)
        self.assertEqual(a['model']['selection']['all_three_criteria_count'],0)
        self.assertEqual(a['model']['cell_memory']['all_three_criteria_count'],39)
    def test_endpoint_mismatch_is_not_numeric_fit(self):
        with (ROOT/'results'/'v8'/'endpoint_robustness_v8_summary.csv').open() as fh:
            a=list(csv.DictReader(fh))
        r=next(x for x in a if x['model']=='cell_memory' and x['arm']=='0h')
        self.assertGreater(float(r['mean_restricted_first_passage_h']),0)
        self.assertLess(float(r['log_linear_onset_h']),0)
if __name__=='__main__':unittest.main()
