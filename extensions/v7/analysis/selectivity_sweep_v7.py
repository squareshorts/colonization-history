"""Independent synthetic sweep of assay cross-reactivity under absent acquired state.

This guards against reporting only a hand-picked incorrect sensitivity. Results
remain within a dimensionless additive observation model, not bacterial data.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from calibration_bounds_v7 import Scenario,aggregate,BOUNDS,NOMINAL


def run(outdir=Path('results/selectivity_sweep_v7'),replicates=500,ns=(32,128),seed=20261107):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    qvals=np.round(np.linspace(0.,.30,13),3)
    rows=[]
    for n in ns:
        for q in qvals:
            s=Scenario('cross_reactivity_sweep',A=0.,S=.8,E=0.,R=.7,I=-.12,
                       qA=.78,qS=float(q),qE=.14,qR=.67)
            v=aggregate(s,n,replicates,rng)
            rows.append({'n_fields':n,'q_selection_true':float(q),
                         'inside_stipulated_bounds':bool(BOUNDS['S'][0] <= q <= BOUNDS['S'][1]),
                         'point_false_A':v['point_A_positive_fraction'],
                         'bounded_false_A':v['set_A_positive_fraction'],
                         'feasible_rate':v['feasible_rate'],
                         'bounded_A_coverage':v['set_A_coverage_conditional_on_feasibility']})
    with (outdir/'selectivity_sweep_v7.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (outdir/'selectivity_sweep_v7.json').open('w') as f:
        json.dump({'n_independent_fields':list(ns),'replicates':replicates,
                   'q_nominal':NOMINAL['S'],'q_bounds':BOUNDS['S'],
                   'seed':seed,'rows':rows},f,indent=2)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.7,4.8))
    for n in ns:
        rs=[r for r in rows if r['n_fields']==n]
        ax.plot([r['q_selection_true'] for r in rs],
                [r['point_false_A'] for r in rs],marker='o',label=f'Point calibration, n={n}')
        ax.plot([r['q_selection_true'] for r in rs],
                [r['bounded_false_A'] for r in rs],marker='s',linestyle='--',label=f'Bounded calibration, n={n}')
    ax.axvline(NOMINAL['S'],color='0.4',linewidth=1,linestyle=':')
    ax.axvline(BOUNDS['S'][1],color='0.5',linewidth=1,linestyle='-.')
    ax.set(xlabel='True cross-reactivity of selected cells in cellular assay',
           ylabel='Fraction falsely declaring acquired-cell contribution',
           ylim=(-.025,1.025),title='A calibration boundary, not sample size, controls false attribution')
    ax.legend(fontsize=9,loc='upper left')
    fig.tight_layout();fig.savefig(outdir/'false_attribution_sweep.png',dpi=200)
    print('swept',len(rows),'cells')
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--replicates',type=int,default=500)
    p.add_argument('--n',type=int,nargs='+',default=[32,128])
    p.add_argument('--outdir',type=Path,default=Path('results/selectivity_sweep_v7'))
    p.add_argument('--seed',type=int,default=20261107)
    a=p.parse_args()
    run(a.outdir,a.replicates,a.n,a.seed)
