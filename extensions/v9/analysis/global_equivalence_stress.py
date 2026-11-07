"""Parameter-global identity test, independent of Lee et al. point fits."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
rng=np.random.default_rng(202611082)
n=1000
f0=rng.uniform(.05,.95,n)
c=rng.uniform(.04,.25,n)
tau=rng.uniform(6.,25.,n)
m0=rng.uniform(.4,.95,n)
nu=c*tau
q=(1-f0)/f0
K=m0*np.power(q,1/nu)
t=np.linspace(0,60,121)
f=1/(1+q[:,None]*np.exp(c[:,None]*t[None,:]))
m=m0[:,None]*np.exp(-t[None,:]/tau[:,None])
h=1/(1+np.power(K[:,None]/m,nu[:,None]))
maxgap=np.max(np.abs(f-h),axis=1)
out={'parameter_draws':n,'sampling':'independent parameter stress for f0 in [0.05,0.95], c in [0.04,0.25]/h, tau in [6,25]h, m0 in [0.4,0.95]; no fitted bacterial observations',
     'holding_times_count':len(t),'max_abs_fraction_difference':float(maxgap.max()),
     'median_max_abs_fraction_difference':float(np.median(maxgap)),
     'all_parameter_draws_satisfy_1e-12':bool(np.all(maxgap<1e-12))}
(ROOT/'results'/'global_equivalence_stress.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
