import csv,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from observed_establishment import moving_average,time_to_sustained_fold,paired_summary

class TestObservedEstablishment(unittest.TestCase):
    def test_constant_no_crossing(self):
        t=np.arange(2400)/120
        y=np.ones_like(t)*3
        val,base,target=time_to_sustained_fold(t,y)
        self.assertTrue(np.isnan(val))
        self.assertAlmostEqual(base,3)
        self.assertAlmostEqual(target,15)

    def test_boundary_smoothing_preserves_constant(self):
        np.testing.assert_allclose(moving_average(np.full(100,4.),20),4.,atol=1e-13)

    def test_stable_crossing_within_record(self):
        t=np.arange(3600)/120
        y=np.where(t<11,2.,12.)
        val,base,target=time_to_sustained_fold(t,y)
        self.assertGreater(val,10.8)
        self.assertLess(val,11.1)
        self.assertEqual(base,2)
        self.assertEqual(target,10)
        self.assertLess(val, t[-1])

    def test_spike_rejected_by_15_min_rule(self):
        t=np.arange(2000)/120
        y=np.full_like(t,2.)
        y[(t>=6)&(t<6.08)]=24.
        val,_,_=time_to_sustained_fold(t,y)
        self.assertTrue(np.isnan(val))

    def test_floor_affects_nearly_empty_baseline(self):
        t=np.arange(2400)/120
        y=np.where(t<8,.4,3.)
        val,base,target=time_to_sustained_fold(t,y)
        self.assertTrue(np.isnan(val))
        self.assertAlmostEqual(base,1.)
        self.assertAlmostEqual(target,5.)
        val2,base2,target2=time_to_sustained_fold(t,y,minimum_baseline=0.)
        self.assertTrue(np.isfinite(val2))
        self.assertLess(val2,8.2)

    def test_archival_pair_summary_matches_report(self):
        with (ROOT/'results/observed_establishment.csv').open(newline='') as f:
            rows=list(csv.DictReader(f))
        rows=[{'record':r['record'],'time_h':float(r['time_h'])} for r in rows]
        v=paired_summary(rows)
        self.assertAlmostEqual(v['paired_mean_h'],-10.3638888667,places=6)
        self.assertTrue(all(d<0 for d in v['paired_differences_h']))
        self.assertLess(v['descriptive_95pct_interval_h'][1],0.)
        self.assertTrue(any(r['record']=='FC1-2' for r in rows))

    def test_all_archived_crossings_observed(self):
        with (ROOT/'results/observed_establishment.csv').open(newline='') as f:
            rows=list(csv.DictReader(f))
        for r in rows:
            self.assertGreater(float(r['time_h']),0)
            self.assertLess(float(r['time_h'])+.25,float(r['last_record_h']))

    def test_sensitivity_grid_shows_consistency(self):
        with (ROOT/'results/observed_establishment_sensitivity.csv').open(newline='') as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),6)
        for r in rows:
            self.assertTrue(all(float(r[f'pair_{i}_h'])<0 for i in (1,2,3)))

if __name__=='__main__':unittest.main()
