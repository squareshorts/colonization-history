"""Minimal no-intracellular-memory counterexample for lag washout.

One selected immutable adhesive type suffers a relative growth-rate cost during
suspension. Population composition changes by differential amplification; no
individual cell acquires or loses a molecular memory state.

Uses the same phenomenological linear lag mapping as the exponential benchmark.
This is exact fit with two parameters and TWO remaining observations, not a
validated mechanistic explanation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'v8';OUT.mkdir(parents=True,exist_ok=True)
Lnaive,L0,L14,L37=17.56,4.68,5.46,13.41
A0=Lnaive-L0;r14=(Lnaive-L14)/A0;r37=(Lnaive-L37)/A0
z14=1/r14-1;z37=1/r37-1
# (exp(c*37)-1)/(exp(c*14)-1)=z37/z14
c=brentq(lambda c:np.expm1(c*37)/np.expm1(c*14)-z37/z14,.00001,0.5,xtol=1e-13)
f0=1-z14/np.expm1(c*14)
def f(t):return 1/(1+(1/f0-1)*np.exp(c*t))
def lag(t):return Lnaive-A0*f(t)/f0
model={'relative_suspension_fitness_cost_per_hour':float(c),'fraction_of_high_adhesiveness_type_immediately_after_selection':float(f0),
       'high_adhesive_type_fraction_14h':float(f(14)), 'high_adhesive_type_fraction_37h':float(f(37)),
       'predicted_lags_h':{str(t):float(lag(t)) for t in [0,7,14,21,24,30,37,48]},
       'model':'two immutable adhesive phenotypes; relative planktonic growth cost changes population composition during hold; observable lag advantage proportional to high-adhesive frequency',
       'warnings':['Both free parameters tuned to two nonzero hold-time points. This is exact interpolation, not biological validation or prospective prediction.','Requires experimentally unmeasured initial composition and relative growth cost.','The lag mapping is phenomenological; not derived from a flow-cell agent simulator.']}
(OUT/'selection_composition_counterexample.json').write_text(json.dumps(model,indent=2))
print(json.dumps(model,indent=2))
