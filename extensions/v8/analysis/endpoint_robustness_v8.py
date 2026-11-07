"""Endpoint sensitivity in the exploratory nonspatial-agent model.

This analysis does not claim to reproduce Lee et al.'s exact exponential-lag
fitting procedure, which used actual image-derived trajectories. A transparent
log-linear extrapolation is a second synthetic readout, with diagnostics.
"""
import sys, json, csv
from dataclasses import replace
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
_agent_source=ROOT/'agents'/'src'
if not _agent_source.is_dir():
    _agent_source=ROOT.parent/'minimal_agents'/'src'
sys.path.insert(0,str(_agent_source))
from minimal_agents import Parameters,CAPABILITIES,run_field
OUT=ROOT/'results'/'v8';OUT.mkdir(parents=True,exist_ok=True)

def log_onset(curve,dt):
    t=np.arange(len(curve))*dt
    y=np.asarray(curve,float)
    # Predefined operational mid-growth region; only when counts genuinely grow.
    ok=(t>=1)&(y>=8)&(y<=120)
    # Restrict to rising limb: earliest contiguous region before first 120 crossing.
    crossings=np.flatnonzero(y>=120)
    if len(crossings):ok &= np.arange(len(y))<=crossings[0]
    if ok.sum()<8:return np.nan,np.nan,int(ok.sum())
    x=t[ok];lny=np.log(y[ok]);
    slope,intercept=np.polyfit(x,lny,1)
    res=lny-(slope*x+intercept)
    ss_tot=np.sum((lny-np.mean(lny))**2)
    r2=1-np.sum(res**2)/ss_tot if ss_tot>0 else np.nan
    if slope<=0 or not np.isfinite(r2):return np.nan,r2,int(ok.sum())
    onset=-intercept/slope
    if not (-50<onset<40):return np.nan,r2,int(ok.sum())
    return float(onset),float(r2),int(ok.sum())

def run(n=48,seed=87000):
    p=replace(Parameters(),stop_at_lag_and_8h=False)
    rows=[]; summary=[]
    for model in ('passive','cell_memory','cell_one_bit','cell_no_inheritance','selection'):
        for j,(arm,hold,tr) in enumerate((('fresh',0.,False),('0h',0.,True),('14h',14.,True),('37h',37.,True))):
            curves=[]; first=[]
            for r in range(n):
                run=run_field(p,CAPABILITIES[model],hold,seed+j*100000+r*73,transfer=tr,store_curve=True,
                              fresh_surface=True)
                curves.append(np.asarray(run['curve'],float))
                first.append(float(run['restricted_threshold_h']))
            av=np.mean(np.asarray(curves),axis=0)
            logon,r2,points=log_onset(av,p.dt_h)
            # Ensemble first-crossing readout, rather than mean individual lag.
            cross=np.flatnonzero(av>=100)
            mean_cross=float(cross[0]*p.dt_h) if len(cross) else None
            row={'model':model,'arm':arm,'n_fields':n,'mean_restricted_first_passage_h':float(np.mean(first)),
                'ensemble_threshold_h':mean_cross,'log_linear_onset_h':logon if np.isfinite(logon) else None,
                'log_linear_fit_r2':r2 if np.isfinite(r2) else None,'fit_points':points,
                'max_ensemble_count':float(np.max(av))}
            summary.append(row)
            print(model,arm,row['mean_restricted_first_passage_h'],row['log_linear_onset_h'],row['log_linear_fit_r2'],flush=True)
            for t,y in zip(np.arange(len(av))*p.dt_h,av):rows.append({'model':model,'arm':arm,'time_h':float(t),'mean_count':float(y)})
    for name,values in [('endpoint_robustness_v8_summary.csv',summary),('endpoint_robustness_v8_curves.csv',rows)]:
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
    (OUT/'endpoint_robustness_v8.json').write_text(json.dumps({'scope':'Uncalibrated stochastic simulations; log-linear onset is an explicit alternative readout, not published Lee fit','n_fields':n,'seed':seed,'summary':summary},indent=2))
if __name__=='__main__':
    import argparse
    z=argparse.ArgumentParser();z.add_argument('--n',type=int,default=16);z.add_argument('--seed',type=int,default=87000);a=z.parse_args();run(a.n,a.seed)
