"""Purely illustrative oscillator vs scalar-OU capability test.

Lee et al. 2018 report cAMP / type IV pilus cross-correlation sign reversal
near +/-5 h. A linear scalar Markov OU state with an instantaneous negative
readout cannot produce that sign change. A coupled 2-state OU oscillator can.
No fluorescence data is fitted. No inference of actual biochemical necessity.
"""
from pathlib import Path
import numpy as np,json,csv

def simulate(kind='oscillator',n_agents=160,seed=99213,dt=.05,horizon=55.,burn=10.):
    rng=np.random.default_rng(seed)
    n=int(round(horizon/dt))
    x=np.zeros((n_agents,n),np.float32);y=np.zeros_like(x)
    ga=.22;om=.60;noise=1.0
    for j in range(1,n):
        a=x[:,j-1]
        b=y[:,j-1]
        x[:,j]=a+dt*(-ga*a-(om*b if kind=='oscillator' else 0.))+np.sqrt(dt)*noise*rng.standard_normal(n_agents)
        y[:,j]=b+dt*((om*a if kind=='oscillator' else 0.)-.6*b)+np.sqrt(dt)*noise*rng.standard_normal(n_agents)
    xx=x[:,int(round(burn/dt)):]
    yy=y[:,int(round(burn/dt)):]
    tfp=(-.8*xx+.35*yy).astype(float)
    return xx.astype(float),tfp,dt

def corr_vs_lag(x,y,dt,lags=np.arange(0,8.01,.5)):
    # Unit-of-inference is a lineage: time-normalized within lineage, then
    # average across independently simulated agents/lineages.
    x=(x-x.mean())/max(x.std(),1e-8)
    y=(y-y.mean())/max(y.std(),1e-8)
    d={}
    for lag in lags:
        k=int(round(lag/dt))
        z=np.mean(x[:,:-k or None]*y[:,k:] if k else x*y,axis=1)
        d[float(lag)]=z
    return d

def run():
    out=Path(__file__).resolve().parent.parent/'results'
    out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(101)
    rows=[]
    for kind in ('scalar','oscillator'):
        x,y,dt=simulate(kind)
        for lag,values in corr_vs_lag(x,y,dt).items():
            bootstrap=np.mean(rng.choice(values,(1200,len(values)),replace=True),axis=1)
            rows.append(dict(model=kind,lag_h=lag,mean=float(values.mean()),
                ci_lo=float(np.quantile(bootstrap,.025)),ci_hi=float(np.quantile(bootstrap,.975)),n_lineages=len(values)))
    with (out/'oscillator_correlations.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'oscillator_correlations.json').write_text(json.dumps(rows,indent=2))
    print([(z['model'],z['lag_h'],round(z['mean'],2),round(z['ci_lo'],2),round(z['ci_hi'],2)) for z in rows if z['lag_h'] in (0.,5.,8.)])
    return rows
if __name__=='__main__':run()
