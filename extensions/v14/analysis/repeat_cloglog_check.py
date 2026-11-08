"""Reproduce a 32,000-study independent check of half-rate sham-loss detection."""
import json
from pathlib import Path
import numpy as np
from v12_model_core import simulate
R=Path(__file__).resolve().parents[1]

def run(seed=20261125):
    rng=np.random.default_rng(seed)
    rates=[]
    for _ in range(4):
        r=simulate(rng,'sham_half_uniform',.10,.85,'cloglog_weibull',studies=8000)
        rates.append(r['bounded_joint'])
    summary={'seed':seed,'per_replicate':8000,'samples':rates,'combined':float(np.mean(rates)),
             'original_1200':.3433333333333333,'fine_1600':.304375}
    (R/'results/v14_cloglog_repeat_check.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=='__main__':run()
