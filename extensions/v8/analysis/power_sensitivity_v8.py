"""Hypothetical power sensitivity from three heterogeneous paired field contrasts.

NOT a validated sample-size recommendation: n=3 variance is unstable, units
may be nonexchangeable and original exposures were not randomized.
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import nct,t
ROOT=Path(__file__).resolve().parents[1]
a=np.array([.136187,.038469,-.126263])
sd=float(a.std(ddof=1))
rows=[]
for effect in [0.05,0.10,0.15,0.20]:
    for inflation in [1.,1.5,2.]:
        s=sd*inflation
        for n in range(3,501):
            df=n-1;crit=t.ppf(.975,df);noncentral=effect*np.sqrt(n)/s
            power=float(nct.sf(crit,df,noncentral)+nct.cdf(-crit,df,noncentral))
            if power>=.8:break
        rows.append({'hypothesized_true_mean_difference':effect,'sd_inflation':inflation,
             'assumed_pair_sd':s,'required_independent_pairs_for_nominal_80percent_power':n,
             'achieved_power':power})
result={'observed_pair_differences':a.tolist(),'estimated_sd_n3':sd,
        'alpha_two_sided':.05,'target_power':.8,'method':'Student-t noncentral distribution, hypothetical independent paired observations',
        'warning':'n=3 variance is unstable; original series differ in timing, and no randomized exposure assignment supports exact t assumptions. These are sensitivity numbers, NOT validated prospective sample sizes.',
        'results':rows}
(ROOT/'results'/'v8'/'power_sensitivity_v8.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
