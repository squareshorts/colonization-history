"""Precommitted parameter sensitivity to distinguish capability from tuning."""
import csv,json,itertools
from pathlib import Path
from dataclasses import replace
import numpy as np
from minimal_agents import *

OUT=Path(__file__).resolve().parent.parent/'results'

def run(reps=20):
    base=Parameters()
    grid=list(itertools.product([7.,10.,14.],[.075,.115,.160],[4.,7.,10.]))
    results=[]
    for i,(tau,mid,hill) in enumerate(grid):
        p=replace(base,free_memory_tau_h=tau,memory_midpoint=mid,memory_hill=hill)
        x=summarize(replicated(p,CAPABILITIES['cell_memory'],reps=reps,seed=80151+i*53))
        a={k:float(v['restricted_threshold_h']['mean']) for k,v in x.items()}
        # Qualitative signatures based on independent endpoint, not numerical fitting:
        # faster zero-gap transfer, 14-h persistence, and delayed 37-h washout.
        z=dict(tau=tau,mid=mid,hill=hill,**a)
        z['transfer_advantage_h']=a['fresh']-a['transfer_0h']
        z['washout_delay_h']=a['transfer_37h']-a['transfer_0h']
        z['short_hold_change_h']=a['transfer_14h']-a['transfer_0h']
        z['all_three_qualitative']=(z['transfer_advantage_h']>=4 and z['washout_delay_h']>=4 and z['short_hold_change_h']<=3 and z['short_hold_change_h']>=-1)
        results.append(z)
        print(i+1,len(grid),z['all_three_qualitative'],flush=True)
    OUT.mkdir(exist_ok=True)
    with (OUT/'robustness_grid.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    (OUT/'robustness_grid.json').write_text(json.dumps({'replicates':reps,'results':results},indent=2))
    print('pass',sum(x['all_three_qualitative'] for x in results),'of',len(results))
    return results
if __name__=='__main__':run()
