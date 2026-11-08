import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis'))
from endpoint_observation_comparison import absolute_sustained_time,normalized_advantage,paired_differences,summarize_from_records

class TestEndpoints(unittest.TestCase):
    def test_crossing_when_initially_above(self):
        t=np.arange(600)/120
        x=np.full_like(t,40.)
        value,status=absolute_sustained_time(t,x,reference_px=1.,threshold=30)
        self.assertEqual(value,0.)
        self.assertEqual(status,'at_start')
    def test_record_right_censoring(self):
        t=np.arange(600)/120
        x=np.full_like(t,25.)
        value,status=absolute_sustained_time(t,x,reference_px=1.,threshold=50)
        self.assertTrue(np.isnan(value))
        self.assertEqual(status,'after_end')
    def test_transient_not_maintained(self):
        t=np.arange(1200)/120
        x=np.full_like(t,3.)
        x[(t>3)&(t<3.07)]=120.
        value,status=absolute_sustained_time(t,x,reference_px=1.,threshold=20)
        self.assertTrue(np.isnan(value))
        self.assertEqual(status,'after_end')
    def test_observed_crossing(self):
        t=np.arange(1200)/120
        x=np.where(t<4,3.,60.)
        value,status=absolute_sustained_time(t,x,reference_px=1.,threshold=30)
        self.assertTrue(3.9<value<4.2)
        self.assertEqual(status,'observed')
    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            absolute_sustained_time(np.arange(4),np.arange(4),0,20)
    def test_censoring_gives_bound(self):
        r=[]
        for j in (1,2,3):
            r.append({'record':f'FC1-{j}','absolute50_h':None if j==2 else 25.,'record_end_h':26.74})
            r.append({'record':f'FC2-{j}','absolute50_h':15.1,'record_end_h':19.})
        d=paired_differences(r,'absolute50_h')
        self.assertEqual(d[1]['status'],'earlier_censored')
        self.assertAlmostEqual(d[1]['upper_bound_h'],-11.64)
    def test_normalized_advantage(self):
        self.assertAlmostEqual(normalized_advantage(25.2444444444,6.4222222222,7.7166666667),.931227863,places=7)
    def test_golden_absolute30(self):
        known={
           'FC1-1':24.933333333333334,'FC2-1':.075,
           'FC1-2':26.425,'FC2-2':14.141666666666667,
           'FC1-3':24.375,'FC2-3':5.05,
           'Tr14':7.716666666666667,'Tr37':16.141666666666666,
           'FC3-3':0.0}
        rs=[{'record':k,'absolute30_h':v,'record_end_h':30.} for k,v in known.items()]
        by=summarize_from_records([{**r, 'relative5_h':r['absolute30_h'],
           'absolute20_h':r['absolute30_h'],'absolute50_h':r['absolute30_h']} for r in rs])
        self.assertAlmostEqual(by['absolute30_h']['r14'],0.9312278630460449,places=7)
        self.assertAlmostEqual(by['absolute30_h']['r37'],0.48361865407319954,places=7)

if __name__=='__main__':unittest.main()

class TestArchivedCrossChecks(unittest.TestCase):
    def test_saved_censoring_statuses(self):
        import csv
        with (ROOT/'results/endpoint_observation_comparison.csv').open(newline='') as f:
            rows={r['record']:r for r in csv.DictReader(f)}
        self.assertEqual(rows['FC2-1']['absolute20_status'],'at_start')
        self.assertEqual(rows['FC1-2']['absolute50_status'],'after_end')
        self.assertEqual(rows['FC1-2']['absolute50_h'],'')
        self.assertAlmostEqual(float(rows['FC2-2']['absolute50_h']),15.1083333333,places=6)
    def test_readout_ratios_in_saved_results(self):
        import json
        summary=json.loads((ROOT/'results/endpoint_observation_summary.json').read_text())
        self.assertAlmostEqual(summary['absolute30_h']['r14'],.9312278630460449,places=7)
        self.assertAlmostEqual(summary['absolute30_h']['r37'],.48361865407319954,places=7)
        self.assertAlmostEqual(summary['relative5_h']['r37'],.7000804073974807,places=7)
    def test_absolute_sensitivity_sign_and_censoring(self):
        import csv
        with (ROOT/'results/absolute_threshold_sensitivity.csv').open(newline='') as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),27)
        self.assertEqual({int(r['threshold_cell_eq']) for r in rows},{20,30,50})
        for r in rows:
            for col in ('series1_fc2_minus_fc1_h_or_upper_bound',
                        'series2_fc2_minus_fc1_h_or_upper_bound',
                        'series3_fc2_minus_fc1_h_or_upper_bound'):
                self.assertLess(float(r[col]),0)
