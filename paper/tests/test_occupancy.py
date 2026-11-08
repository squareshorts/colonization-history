import json,sys,unittest
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'analysis'))
from occupancy_estimator import fit_lag,DATASETS

class EdgeCorrectedOccupancy(unittest.TestCase):
    def test_constant_trace_stays_constant_to_boundaries(self):
        counts=np.full(120,40.,dtype=float)
        kernel=np.ones(20)
        smoothed=np.convolve(counts,kernel,mode="same")/np.convolve(np.ones_like(counts),kernel,mode="same")
        np.testing.assert_allclose(smoothed,40.,rtol=0,atol=1e-13)
    def test_exponential_lag_remains_finite(self):
        t=np.arange(600)*.008333333333333333
        y=np.exp(.7*(t-1))+25
        lag,n=fit_lag(t,y,smooth=20,lo=20,hi=300)
        self.assertTrue(np.isfinite(lag));self.assertGreater(n,20)
    def test_third_flowcell_included(self):
        self.assertEqual(sum(label=="FC3-3" for label,_ in DATASETS),1)
        summary=json.loads((R/"results/shared_calibration_summary.json").read_text())
        self.assertAlmostEqual(summary["FC3_single_series_3_lag_h"],.9359044992,places=8)
    def test_paired_effects_persist_edge_correction(self):
        s=json.loads((R/"results/edge_corrected_summary.json").read_text())
        self.assertAlmostEqual(s["corrected_per_record"]["paired_mean_h"],-10.69076006456,places=8)
        self.assertAlmostEqual(s["corrected_shared"]["paired_mean_h"],-11.66481609207,places=8)
        self.assertTrue(all(x<0 for x in s["corrected_shared"]["paired_FC2_minus_FC1_h"]))
    def test_large_series_2_lag_extrapolated(self):
        import csv
        with (R/"results/shared_calibration.csv").open() as f: rows={r["dataset"]:r for r in csv.DictReader(f)}
        self.assertEqual(int(rows["FC1-2"]["n_standard_fit"]),152)
        self.assertGreater(float(rows["FC1-2"]["original_lag_h"]),20.)

if __name__=="__main__": unittest.main()
