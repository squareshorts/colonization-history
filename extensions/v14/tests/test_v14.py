import csv,json,re,sys,unittest
from pathlib import Path
from math import sqrt
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'analysis'))
import calibration_stress as CS

def records(filename):
 with (R/'results'/filename).open(newline='') as h:return list(csv.DictReader(h))
class AuditV14(unittest.TestCase):
 def test_original_34_percent_not_duplicated(self):
  a=next(x for x in records('v12_two_arm_loss_grid.csv') if x['kinetic']=='cloglog_weibull' and x['scenario']=='sham_half_uniform' and x['effect']=='0.85' and x['nominal_loss']=='0.1')
  b=next(x for x in records('fine_loss_detection.csv') if x['kinetic']=='cloglog_weibull' and x['scenario']=='sham_half_uniform' and x['effect']=='0.85' and x['nominal_loss']=='0.1')
  c=next(x for x in records('assisted_calibration.csv') if x['loss_mechanism']=='uniform' and x['external_sham_upper_p']=='0.1' and x['effect']=='0.85' and x['epsilon']=='0.1')
  self.assertEqual(int(a['studies']),1200);self.assertEqual(int(b['studies']),1600);self.assertEqual(int(c['studies']),2800)
  self.assertAlmostEqual(float(a['bounded_joint']),.3433333333333333)
  self.assertAlmostEqual(float(b['bounded_joint']),.304375)
  self.assertAlmostEqual(float(c['cap_joint']),.34285714285714286)
 def test_strong_monte_carlo_replication(self):
  x=json.loads((R/'results/v14_cloglog_repeat_check.json').read_text())
  self.assertEqual(x['per_replicate']*len(x['samples']),32000)
  self.assertAlmostEqual(x['combined'],.328125)
 def test_exact_expected_bound_and_equal_loss(self):
  p0=.09443261036368705;d=.15927438472747357;m=.19514708121837185
  for p1 in [d,m]:
   threshold=(p1-p0)/(1+p1-p0)
   for epsilon in (threshold-.001,threshold,threshold+.001):
    lower=(1-epsilon)*p1-p0-epsilon*(1-p0)
    self.assertAlmostEqual(lower,(p1-p0)-epsilon*(1+p1-p0))
    if epsilon<threshold:self.assertGreater(lower,0)
    if epsilon>threshold:self.assertLess(lower,0)
  self.assertAlmostEqual((d-p0)/(1+d-p0),.06089334202025243)
  self.assertAlmostEqual((m-p0)/(1+m-p0),.09149917941614943)
 def test_invalid_cap_adversarial_false_attribution(self):
  a=records('null_asymmetric_calibration_stress.csv')
  h=next(x for x in a if x['phenotype_low']=='0.8' and x['sham_loss']=='0.4' and x['pulse_loss']=='0.0')
  self.assertEqual(int(h['studies']),6000)
  self.assertEqual(int(h['invalid_u003_joint_count']),929)
  self.assertAlmostEqual(float(h['invalid_u003_joint']),.15483333333333332)
  self.assertEqual(float(h['valid_u010_joint']),0.)
  self.assertLess(float(h['invalid_u003_coverage_mother']),.01)
 def test_equal_loss_does_not_misstate_false_attribution(self):
  x=next(y for y in records('assisted_calibration.csv') if y['epsilon']=='0.1' and y['effect']=='0.0' and y['loss_mechanism']=='uniform' and y['external_sham_upper_p']=='0.03')
  self.assertEqual(float(x['cap_joint']),0.0)
  self.assertLess(float(x['mean_missing_sham_tail_coverage']),.95)
 def test_calibration_bound_is_outside_kinetic_parameter_distribution(self):
  # A specific admissible main-simulator parameter tuple exceeds the fixed cap.
  from v12_model_core import attachment_probability
  p=attachment_probability(-1.8,1.,1.6,0.,.6,.1,'logistic_exponential')
  self.assertGreater(float(p),.10)
 def test_fine_grid_factorial_has_96(self):
  a=records('fine_loss_detection.csv')
  self.assertEqual(len(a),8*3*2*2)
 @unittest.skipUnless((R/'overleaf/main.tex').is_file(),'LaTeX comes in separate Overleaf ZIP')
 def test_manuscript_clean(self):
  s=(R/'overleaf/main.tex').read_text();t=(R/'overleaf/supplement.tex').read_text()
  abstract=s.split('\\begin{abstract}',1)[1].split('\\end{abstract}',1)[0]
  self.assertLessEqual(len(abstract.strip().split()),200)
  self.assertNotIn('sharply reduced joint',abstract)
  self.assertIn('\\label{eq:equalattrition}',s)
  self.assertIn('15.5\\%',s)
  self.assertIn('15.48\\%',t)
  self.assertIn('32,000-study',t)
  self.assertIn('\\cite{ImbensManski2004}',s)
  self.assertNotIn('\\cite{LeeTrimming2009,ImbensManski2004}',s)
  self.assertTrue(all(not re.search(p,s+'\n'+t,flags=re.I) for p in [r'\\bv(?:7|8|9|10|11|12|13|14)\\b','hostile.review','previous version','overleaf packages']))
 def test_small_stress_replay(self):
  a=CS.simulate_single(seed=119,studies=150,sham_loss=.30,pulse_loss=0.,effect=0.,phenotype_low=.80,phenotype_high=.99)
  self.assertIn('invalid_u003_joint',a)
  self.assertLess(a['invalid_u003_coverage_mother'],a['valid_u010_coverage_mother'])
if __name__=='__main__':unittest.main()
