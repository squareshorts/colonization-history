"""Calculate additional identifiability boundaries; append to v9 outputs."""
import csv,json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'analysis'))
from discriminating_test_v9 import five_carrier_rank,calibration_geometry
path=root/'results'/'critical_test_v9.json'
obj=json.loads(path.read_text())
obj['five_carrier_state_vs_coating_nonidentifiability']=five_carrier_rank()
obj['calibration_geometry_scope']='Gaussian additive four-carrier contrast errors, field SD=0.12, n=32; illustrative scale only'
r=calibration_geometry()
with (root/'results'/'assay_design_condition_numbers.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=r[0]);w.writeheader();w.writerows(r)
obj['design_ill_conditioning_examples']=[x for x in r if (x['selection_leakage_kappa'],x['trace_retention_q']) in [(0.,1.),(.8,.75),(.99,.6)]]
path.write_text(json.dumps(obj,indent=2))
print('five carriers rank',obj['five_carrier_state_vs_coating_nonidentifiability'])
print('geometry examples',obj['design_ill_conditioning_examples'])
