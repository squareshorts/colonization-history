"""Worst-case deterministic perturbation audit; never a statistical CI."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution
ROOT=Path(__file__).resolve().parents[1]
L=np.array([17.56,4.68,5.46,13.41],float)

def margin(z):
    n, l0, l14, l37=z
    if n<=max(l0,l14,l37)+1e-9:return -1e5
    d0=n-l0; d14=n-l14; d37=n-l37
    if min(d0,d14,d37)<=0:return -1e5
    r14=d14/d0;r37=d37/d0
    return float(r37-r14**(37/14))

report={'origin':'Lee et al 2018 reported fitted lag point estimates; uniform bounded perturbations are not uncertainty distributions',
        'baseline_margin':margin(L),'positive_margin_means_positive_exponential_mixture_might_fit':True,'perturbations':[]}
for delta in [0.0,0.25,0.5,1.,2.,3.,4.]:
    if delta==0:z=L;m=margin(L)
    else:
        x=differential_evolution(lambda z:-margin(z),[(a-delta,a+delta) for a in L],seed=9821,popsize=16,maxiter=200,polish=True,tol=1e-10)
        z=x.x;m=margin(z)
    # Certified loose interval envelope: independently maximize r37 and
    # minimize r14. If its margin is negative, the incompatibility holds for
    # EVERY possible combination inside the supplied lag intervals.
    if delta == 0:
        upper37 = (L[0]-L[3])/(L[0]-L[1]);lower14 = (L[0]-L[2])/(L[0]-L[1])
    else:
        upper37 = ((L[0]+delta)-(L[3]-delta))/((L[0]+delta)-(L[1]+delta))
        lower14 = ((L[0]-delta)-(L[2]+delta))/((L[0]-delta)-(L[1]-delta))
    certified_margin=float(upper37-lower14**(37/14))
    report['perturbations'].append({'certified_upper_margin_via_interval_envelopes':certified_margin,
            'incompatibility_guaranteed_by_interval_envelope':bool(certified_margin<0),
            'each_lag_perturbed_within_plusminus_h':delta,
            'maximum_attainable_margin':m,'any_compatible_positive_exponential_mixture':bool(m>=-1e-9),
            'maximizing_lags':list(map(float,z))})
print(json.dumps(report,indent=2))
path=ROOT/'results'/'v8'/'lag_uncertainty_robustness.json'
path.write_text(json.dumps(report,indent=2))
