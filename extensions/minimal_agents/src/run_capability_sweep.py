"""Fixed capability set with independent stochastic experimental fields.

Pre-specified comparisons: 0h transferred vs fresh, 14h persistence, 37h
washout, and early lineage expansion. No model fit to published observations.
"""
from pathlib import Path
import argparse,csv,json,sys
import numpy as np
from dataclasses import asdict
from minimal_agents import Parameters,CAPABILITIES,replicated,summarize

OUT=Path(__file__).resolve().parent.parent/'results'
MODELS=['passive','selection','environment','cell_no_decay','cell_linear_feedback',
        'cell_no_inheritance','cell_no_ongoing_sensing','cell_one_bit','cell_memory','cell_plus_environment']

def run(reps=70,seed=7007,models=MODELS):
    out={};raw_rows=[]; summary_rows=[]
    for j,name in enumerate(models):
        print('SIM',name,flush=True)
        runs=replicated(Parameters(),CAPABILITIES[name],reps=reps,seed=seed,model=name)
        out[name]=summarize(runs)
        for arm,records in runs.items():
            for idx,rec in enumerate(records):
                raw_rows.append({'model':name,'arm':arm,'field':idx,**rec})
            s=out[name][arm]
            summary_rows.append(dict(model=name,arm=arm,n_fields=reps,
                restricted_lag_mean=s['restricted_threshold_h']['mean'],
                crossing_fraction=s['threshold_reached']['mean'],
                median_lag_conditional=s['lag_threshold_h']['median'],
                mean_cells_at_8h=s['n_at_8h']['mean'],
                median_cells_at_8h=s['n_at_8h']['median'],
                deep_fraction_at_8h=s['deep_fraction_8h']['mean'],
                short_contact_fraction=s['short_fraction_all_contacts']['mean']))
    OUT.mkdir(exist_ok=True)
    with (OUT/'agent_field_observations.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(raw_rows[0]));w.writeheader();w.writerows(raw_rows)
    with (OUT/'agent_summaries.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary_rows[0]));w.writeheader();w.writerows(summary_rows)
    (OUT/'agent_summaries.json').write_text(json.dumps({'parameters':asdict(Parameters()),'fields_per_arm':reps,'seed':seed,'models':out},indent=2))
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--reps',type=int,default=70);ap.add_argument('--seed',type=int,default=7007)
    ap.add_argument('--models',nargs='*',default=MODELS)
    a=ap.parse_args();res=run(a.reps,a.seed,a.models)
    for name,d in res.items():
        print(name,[round(d[k]['restricted_threshold_h']['mean'],2) for k in d],flush=True)
