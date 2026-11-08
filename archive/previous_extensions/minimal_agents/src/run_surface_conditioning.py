"""Environment-only positive control: naive cells on previously visited surface."""
import csv,json
from pathlib import Path
import numpy as np
from minimal_agents import Parameters,CAPABILITIES,run_field
OUT=Path(__file__).resolve().parent.parent/'results'

def run(n=60):
    p=Parameters();c=CAPABILITIES['environment'];rows=[]
    for k in range(n):
        seed=40101+19*k
        for label,fresh in [('fresh_surface',True),('conditioned_surface',False)]:
            r=run_field(p,c,0.,seed,transfer=False,fresh_surface=fresh)
            rows.append(dict(field=k,surface=label,**r))
    with (OUT/'environment_positive_control.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    res={}
    for label in ('fresh_surface','conditioned_surface'):
        a=[x for x in rows if x['surface']==label]
        res[label]={'n_fields':len(a),'mean_n8':float(np.mean([x['n_at_8h'] for x in a])),
                    'mean_restricted_lag_h':float(np.mean([x['restricted_threshold_h'] for x in a])),
                    'cross_fraction':float(np.mean([x['threshold_reached'] for x in a]))}
    (OUT/'environment_positive_control.json').write_text(json.dumps(res,indent=2))
    print(json.dumps(res,indent=2));return res
if __name__=='__main__':run()
