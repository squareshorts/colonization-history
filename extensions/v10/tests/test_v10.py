"""Deterministic v10 algebra, finite-population and simulation tests."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.special import logit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
import generality_attrition_v10 as v


class V10Tests(unittest.TestCase):
    def test_general_lag_reparam_is_exact_along_path(self):
        t=np.linspace(0,60,1001)
        f0,c,m0,tau=.4,.17,.7,13.
        m=v.latent_state(t,m0,tau)
        f=v.latent_selection(t,f0,c)
        back_m=-tau*np.log(m/m0)
        back_f=(logit(f0)-logit(f))/c
        for name in ['weibull','stretched','biexponential','logistic','delayed']:
            actual=v.lag_family(t,name)
            self.assertLess(np.max(np.abs(v.lag_family(back_m,name)-actual)),1e-12)
            self.assertLess(np.max(np.abs(v.lag_family(back_f,name)-actual)),1e-12)

    def test_disjoint_holdout_uses_no_interpolation_test_values(self):
        x=v.reparameterization_stress(n_trials=20,n_test=150)
        self.assertEqual(len(x),20)
        self.assertLess(np.mean([r['rmse_state_h'] for r in x]),.05)
        self.assertLess(np.mean([r['rmse_selection_h'] for r in x]),.05)

    def test_attrition_sharp_bounds_and_extremes(self):
        lo,hi=v.missing_bounds(12,34,40)
        self.assertAlmostEqual(lo,12/40)
        self.assertAlmostEqual(hi,18/40)
        self.assertAlmostEqual(hi-lo,6/40)

    def test_missing_outcome_bounds_cover_arbitrary_censoring(self):
        rng=np.random.default_rng(12)
        for _ in range(2000):
            y=rng.integers(0,2,size=80)
            observed=rng.random(80)<rng.uniform(.01,1)
            lo,hi=v.missing_bounds(int(y[observed].sum()),int(observed.sum()),80)
            self.assertLessEqual(lo-1e-12,y.mean())
            self.assertGreaterEqual(hi+1e-12,y.mean())

    def test_difference_bounds_cover_full_group_difference(self):
        rng=np.random.default_rng(13)
        for _ in range(700):
            a=rng.integers(0,2,size=47);b=rng.integers(0,2,size=53)
            ra=rng.random(47)>.2;rb=rng.random(53)>.4
            lo,hi=v.difference_bounds(a[ra].sum(),ra.sum(),len(a),b[rb].sum(),rb.sum(),len(b))
            self.assertLessEqual(lo-1e-12,a.mean()-b.mean())
            self.assertGreaterEqual(hi+1e-12,a.mean()-b.mean())

    def test_lossless_equivalence_of_statistical_tests(self):
        rng=np.random.default_rng(14)
        a=rng.binomial(48,.25,(20,32));b=rng.binomial(48,.20,(20,32))
        k=np.full((20,32),48)
        res=v.clone_inference(a,k,b,k,N=48)
        self.assertTrue(np.array_equal(res['raw_detect'],res['bound_detect']))
        self.assertTrue(np.allclose(res['raw_effect'],res['bound_mean_lower']))

    def test_selection_loss_cannot_be_confident_if_robust_bound_negative(self):
        for model in ('exp','weibull'):
            r=v.simulate_clones('phenotype_biased_recovery',model,n_studies=200,
                                recovery_high=.98,recovery_low=.45,seed=118)
            self.assertGreater(np.mean(r['mother']['raw_detect']),.1)
            self.assertEqual(int(np.sum(r['mother']['bound_detect'])),0)
            self.assertEqual(int(np.sum(r['daughter']['bound_detect'])),0)

    def test_control_daughter_versus_inherited_effect(self):
        for model in ('exp','weibull'):
            inherit=v.simulate_clones('inherited_state',model,n_studies=180,seed=25)
            no=v.simulate_clones('no_daughter_inheritance',model,n_studies=180,seed=26)
            self.assertGreater(np.mean(inherit['daughter']['bound_detect']),.8)
            self.assertLess(np.mean(no['daughter']['bound_detect']),.1)

    def test_inherited_coating_is_undistinguished(self):
        for model in ('exp','weibull'):
            coat=v.simulate_clones('inherited_coating',model,n_studies=140,seed=123)
            self.assertGreater(np.mean(coat['daughter']['bound_detect']),.8)

    def test_fifth_assay_rank_and_selectivity(self):
        self.assertEqual(np.linalg.matrix_rank(v.carrier_matrix(1,1)),4)
        self.assertEqual(np.linalg.matrix_rank(v.carrier_matrix(1,0)),5)
        self.assertGreater(np.linalg.cond(v.carrier_matrix(1,.98)),
                           40*np.linalg.cond(v.carrier_matrix(1,0)))

    def test_bounded_A_C_contains_true_with_valid_calibration_and_intervals(self):
        rng=np.random.default_rng(26)
        for _ in range(2500):
            A,C=rng.uniform(0,.6,2)
            gamma=rng.uniform(.95,1)
            kappa=rng.uniform(0,.2)
            U=A+C
            V=gamma*A+kappa*C
            w=rng.uniform(.0001,.06)
            b=v.bounded_AC((U-w,U+w),(V-w,V+w),(.95,1),(0,.2))
            self.assertLessEqual(b['A_interval_conservative'][0]-1e-11,A)
            self.assertGreaterEqual(b['A_interval_conservative'][1]+1e-11,A)
            self.assertLessEqual(b['C_interval_conservative'][0]-1e-11,C)
            self.assertGreaterEqual(b['C_interval_conservative'][1]+1e-11,C)

    def test_calibration_ranges_overlap_means_no_guarantee(self):
        r=v.bounded_AC((.3,.4),(.2,.3),(.3,.5),(.4,.6))
        self.assertFalse(r['identified'])

    def test_invalid_calibration_can_destroy_coverage(self):
        r=v.location_interval_coverage(n_sim=350,seed=20261114)
        self.assertGreater(r[0]['coverage_both_A_C'],.93)
        self.assertLess(r[1]['coverage_both_A_C'],.5)
        self.assertLess(r[2]['coverage_both_A_C'],.5)


if __name__=='__main__':unittest.main()
