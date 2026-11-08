import csv, json, math, re, sys, unittest
from pathlib import Path
import numpy as np
from scipy.stats import binom
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'analysis'))
import fine_recovery as F

def csv_rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

class V13Tests(unittest.TestCase):
    def test_expected_attrition_identity(self):
        e1,e0,p1,p0=.1,.04,.16,.095
        lhs=(1-e1)*p1-(1-e0)*p0-e0
        rhs=(1-e1)*p1-p0-e0*(1-p0)
        self.assertAlmostEqual(lhs,rhs)
    def test_daughter_threshold(self):
        r=F.expected_baselines(seed=119,draws=250000)[0]
        self.assertTrue(.052<r['daughter_sham_loss_limit_at_pulse_10pct']<.056)
        self.assertTrue(r['mother_sham_loss_limit_at_pulse_10pct']>r['daughter_sham_loss_limit_at_pulse_10pct'])
    def test_exact_sharp_bounds(self):
        # True finite assigned fraction must lie in sharp binary bounds.
        N,n,S=48,43,7
        for M in range(N-n+1):
            true=(S+M)/N
            self.assertLessEqual(S/N,true)
            self.assertLessEqual(true,(S+N-n)/N)
    def test_binomial_tail_validity(self):
        m=150; q=binom.ppf(.995,m,.10)
        self.assertGreaterEqual(binom.cdf(q,m,.10),.995)
        self.assertLess(q,m)
    def test_calibration_is_additional_assumption(self):
        x=csv_rows(R/'results/assisted_calibration.csv')
        y=next(a for a in x if a['loss_mechanism']=='uniform' and a['epsilon']=='0.1' and a['effect']=='0.85' and a['external_sham_upper_p']=='0.1')
        self.assertEqual(y['assumption_valid'],'True')
        self.assertGreater(float(y['cap_joint']),float(y['worst_joint']))
    def test_invalid_upper_can_lose_coverage(self):
        rows=csv_rows(R/'results/assisted_calibration.csv')
        x=next(a for a in rows if a['loss_mechanism']=='uniform' and a['epsilon']=='0.1' and a['effect']=='0.85' and a['external_sham_upper_p']=='0.03')
        self.assertEqual(x['assumption_valid'],'False')
        self.assertLess(float(x['mean_missing_sham_tail_coverage']),.99)
    def test_fine_grid_has_both_laws_and_arm_controls(self):
        rows=csv_rows(R/'results/fine_loss_detection.csv')
        self.assertEqual(len(rows),96)
        self.assertEqual({x['kinetic'] for x in rows},set(F.KINETICS))
        self.assertEqual({x['scenario'] for x in rows},set(('pulse_only_uniform','both_uniform_equal','sham_half_uniform')))
    def test_small_loss_cliff(self):
        rows=csv_rows(R/'results/fine_loss_detection.csv')
        def find(x):return float(next(z['bounded_joint'] for z in rows if z['kinetic']=='logistic_exponential' and z['scenario']=='both_uniform_equal' and z['nominal_loss']==str(x) and z['effect']=='0.85'))
        self.assertGreater(find(.01),.90)
        self.assertLess(find(.05),.3)
        self.assertLess(find(.1),.01)
    def test_shared_area_controls(self):
        a=json.loads((R/'results/shared_calibration_summary.json').read_text())
        self.assertEqual(a['shared_area_reference_px'],262)
        self.assertTrue(all(z<0 for z in a['paired_FC2_minus_FC1_h']))
        self.assertTrue(1 < a['advantage_ratio_14h'] < 1.1)
    def test_common_scale_variation_preserves_difference(self):
        x=csv_rows(R/'results/shared_calibration.csv')
        d={r['dataset']:r for r in x}
        for i in (1,2,3):
            self.assertLess(float(d[f'FC2-{i}']['shared_max_h']),float(d[f'FC1-{i}']['shared_min_h']))
    @unittest.skipUnless((R/'overleaf/main.tex').is_file(),'Manuscript sources distributed separately as the Overleaf package')
    def test_full_reference_authors(self):
        bib=(R/'overleaf/references.bib').read_text()
        self.assertNotIn('and others',bib)
        self.assertIn('and Wong, Gerard C. L.',bib)
    @unittest.skipUnless((R/'overleaf/main.tex').is_file(),'Manuscript sources distributed separately as the Overleaf package')
    def test_funding_not_invented_and_roles(self):
        t=(R/'overleaf/main.tex').read_text()
        self.assertNotIn('supervision;',t)
        self.assertIn('The research received no specific grant funding',t)
    @unittest.skipUnless((R/'overleaf/main.tex').is_file(),'Manuscript sources distributed separately as the Overleaf package')
    def test_no_internal_language_in_scientific_text(self):
        for file in ['main.tex','supplement.tex']:
            text=(R/'overleaf'/file).read_text()
            for pattern in (r'\bv(?:7|8|9|10|11|12|13)\b',r'Overleaf package',r'hostile.review',r'editorial workflow',r'to be inserted',r'placeholder'):
                self.assertIsNone(re.search(pattern,text,flags=re.I),(file,pattern))
    @unittest.skipUnless((R/'overleaf/main.tex').is_file(),'Manuscript sources distributed separately as the Overleaf package')
    def test_abstract_below_200(self):
        tex=(R/'overleaf/main.tex').read_text()
        t=tex.split('\\begin{abstract}',1)[1].split('\\end{abstract}',1)[0]
        words=re.findall(r"\b[\w]+(?:[-'][\w]+)*\b",re.sub(r'\\[%&]',' ',t))
        self.assertLessEqual(len(words),200)
    @unittest.skipUnless((R/'overleaf/main.tex').is_file(),'Manuscript sources distributed separately as the Overleaf package')
    def test_crossref_figure_paths(self):
        for name in ['main.tex','supplement.tex']:
            tex=(R/'overleaf'/name).read_text()
            for p in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',tex):
                self.assertTrue((R/'overleaf'/'figures'/p).exists(),(name,p))
    def test_counts_previous_repeats(self):
        rows=csv_rows(R/'results/repeated_contact_information.csv')
        k=[x for x in rows if x['p_high']=='0.1' and x['p_low']=='0.03']
        self.assertEqual(len(k),3)
        self.assertGreater(float(k[0]['normal_oracle_n80']),10000000)

if __name__=='__main__':unittest.main()
