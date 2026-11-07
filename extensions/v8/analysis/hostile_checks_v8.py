"""Adversarial checks of published lag summaries and synthetic-agent sufficiency.

Independent audit 2026-11-07. Published fitted lag point values are not
independent biological replicates. Every synthetic rate is an illustrative prior.
"""
import json, csv, sys
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.stats import qmc
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
_agent_source=ROOT/'agents'/'src'
if not _agent_source.is_dir():
    _agent_source=ROOT.parent/'minimal_agents'/'src'
sys.path.insert(0,str(_agent_source))
from minimal_agents import Parameters, CAPABILITIES, run_field

OUT=ROOT/'results'/'v8'
OUT.mkdir(parents=True,exist_ok=True)
OBS = {'naive':17.56, '0h':4.68, '14h':5.46, '37h':13.41}


def model_audit():
    n=OBS['naive']; d0=n-OBS['0h']; d14=n-OBS['14h']; d37=n-OBS['37h']
    r14=d14/d0; r37=d37/d0
    tau=-14/np.log(r14)
    d37_hat=d0*np.exp(-37/tau)
    # Any positive mixture of simple exponential components with an additive
    # observation mapping must be log convex, hence r37 >= r14**(37/14).
    lower_r37_for_mixture=r14**(37/14)
    eff_short=-np.log(r14)/14
    eff_late=-np.log(r37/r14)/(37-14)
    # Logistic-in-time effective readout. Calibrated to all three transfer lags.
    k=np.log((1/r37-1)/(1/r14-1))/np.log(37/14)
    t50=14/((1/r14-1)**(1/k))
    def lag_logistic(t):return n-d0/(1+(t/t50)**k)
    def lag_exponential(t):return n-d0*np.exp(-t/tau)
    # A priori fixed intracellular state tau=10 h; nonlinear Hill readout.
    # Normalized effect: (1+x)/(1+x exp(k*t/10)); x>0.
    exponents=[0.05,0.1,0.2,0.5,1.,2.,3.,4.,5.,7.,10.]
    hill_rows=[]
    for h in exponents:
        den=r14*np.exp(h*14/10)-1
        if den<=0:continue
        x=(1-r14)/den
        if x<=0:continue
        def ratio(t):return (1+x)/(1+x*np.exp(h*t/10))
        hill_rows.append({'hill_exponent':h,'x':x,'pred_37h':n-d0*ratio(37),'pred_24h':n-d0*ratio(24)})
    full={'source':'Lee et al. PNAS 2018 Fig 4; published fit point values, not independent replicate estimates',
        'observations':OBS, 'advantages':{'0':d0,'14':d14,'37':d37},
        'normalized_advantages':{'14':r14,'37':r37},
        'linear_exponential':{'tau_fitted_to_0_14_h':tau,'prediction_37_h':lag_exponential(37),
            'observed_37_h':OBS['37h'],'prediction_error_observed_minus_predicted_h':OBS['37h']-lag_exponential(37),
            'effective_rate_0_to_14_hinv':eff_short,'effective_rate_14_to_37_hinv':eff_late,
            'late_to_early_rate_ratio':eff_late/eff_short},
        'nonnegative_exponential_mixture':{'minimum_possible_normalized_37_from_normalized_14':lower_r37_for_mixture,
            'observed_normalized_37':r37,'violates_log_convexity_at_point_estimates':bool(r37+1e-12<lower_r37_for_mixture)},
        'two_parameter_effective_hill_in_time':{'h':k,'half_decay_h':t50,'prediction_24h':lag_logistic(24),
                'note':'Two freely fitted parameters and two nonzero-hold observations: exact interpolation is not validation'},
        'single_state_tau10_h_hill_family_fitted_only_at_14h':hill_rows,
        'note':'All calculations are descriptive algebra on published fit summaries; do not use their within-trajectory CIs as biological replicate errors.'}
    (OUT/'published_lag_model_audit.json').write_text(json.dumps(full,indent=2))
    return full


def lhs_values(n=160,seed=2026110709):
    from scipy.special import expit
    u=qmc.LatinHypercube(d=8,seed=seed).random(n)
    defaults=Parameters()
    ps=[]
    for v in u:
        # A broad, explicitly constructed design distribution, without fitting
        # or screening on the held-out 37-h outcome.
        ps.append(replace(defaults,
            free_memory_tau_h=float(np.exp(np.log(4)+v[0]*np.log(35/4))),
            memory_midpoint=float(np.exp(np.log(0.035)+v[1]*np.log(.32/.035))),
            memory_hill=float(1+11*v[2]),
            contact_acquisition=float(.25+.55*v[3]),
            on_surface_acquisition_h=float(.3+1.4*v[4]),
            stable_attach_naive=float(.001+.008*v[5]),
            inflow_per_h=float(15+30*v[6]),
            daughter_ret_naive=float(.45+.36*v[7])
        ))
    return ps


def batch_sweep(n=160,reps=5,seed=2026110709):
    ps=lhs_values(n,seed)
    rows=[]
    for i,p in enumerate(ps):
        for model in ['cell_memory','cell_one_bit','cell_no_inheritance','selection']:
            c=CAPABILITIES[model]
            by_arm={}
            for j,(name,hold,tr) in enumerate([('fresh',0,False),('0h',0,True),('14h',14,True),('37h',37,True)]):
                # Challenge randomness independent across arms. Same parameter
                # and fixed seed formula for every model & param draw.
                ys=[]
                for k in range(reps):
                    idseed=seed+10000000*i+100000*j+1297*k
                    y=run_field(p,c,hold,idseed,transfer=tr)
                    ys.append(y['restricted_threshold_h'])
                by_arm[name]=float(np.mean(ys))
            da=by_arm['fresh']-by_arm['0h']
            db=by_arm['fresh']-by_arm['14h']
            dc=by_arm['fresh']-by_arm['37h']
            # The 0/14 conditions are treated as preliminary calibration targets;
            # 37h is a RETROSPECTIVE PSEUDO-HOLDOUT known from the literature,
            # not a prospectively unseen outcome.
            fit_early=(da>6 and db>6 and db/da>=0.65) if da>0 else False
            pass_late=(dc/da<=0.5 and by_arm['37h']>=by_arm['14h']+4) if da>0 else False
            pass_all=fit_early and pass_late
            row=dict(draw=i,model=model,seed=seed,reps=reps,**{f't_{k}':v for k,v in by_arm.items()},
                     advantage_0=da,advantage_14=db,advantage_37=dc,fit_early=fit_early,
                     pass_late=pass_late,pass_all=pass_all,
                     **{k:getattr(p,k) for k in ['free_memory_tau_h','memory_midpoint','memory_hill','contact_acquisition','on_surface_acquisition_h','stable_attach_naive','inflow_per_h','daughter_ret_naive']})
            rows.append(row)
        if (i+1)%40==0: print('Finished candidate parameter sets:',i+1,flush=True)
    with (OUT/'agent_global_prior_sweep.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    summary={}
    for model in ['cell_memory','cell_one_bit','cell_no_inheritance','selection']:
        mr=[r for r in rows if r['model']==model]
        early=[r for r in mr if r['fit_early']]
        summary[model]={'parameter_draws':len(mr),'fields_per_arm':reps,
             'early_qualifying_draws':len(early),'early_fraction':len(early)/n,
             'late_washout_given_early':sum(x['pass_late'] for x in early)/len(early) if early else None,
             'late_washout_count_given_early':sum(x['pass_late'] for x in early),
             'all_three_criteria_count':sum(x['pass_all'] for x in mr),
             'all_three_criteria_fraction':sum(x['pass_all'] for x in mr)/n,
             'median_early_advantage_0_h':float(np.median([r['advantage_0'] for r in mr])),
             'median_advantage_37_h':float(np.median([r['advantage_37'] for r in mr]))}
    full={'scope':'Latin-hypercube synthetic parameter perturbations; not an empirical prior, not a biological fit',
          'selection_comparator':'Immutable adhesive selection; no growth cost or phenotype change in holding tube',
          'pseudo_holdout_warning':'37h exposure was known to model developer, so this is not blinded or prospective validation',
          'criteria':{'early0':'fresh minus 0h > 6 synthetic h',
                      'early14':'fresh minus 14h > 6 synthetic h and >=0.65 of immediate advantage',
                      'late37':'37h advantage <=0.5 of immediate advantage and 37h threshold at least 4 synthetic h later than 14h'},
          'n_parameter_draws':n,'reps_per_arm':reps,'seed':seed,'model':summary}
    (OUT/'agent_global_prior_sweep.json').write_text(json.dumps(full,indent=2))
    return full


def main():
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--draws',type=int,default=160);a.add_argument('--reps',type=int,default=5);a.add_argument('--seed',type=int,default=2026110709);z=a.parse_args()
    m=model_audit()
    print(json.dumps(m['linear_exponential'],indent=2))
    print(json.dumps(m['nonnegative_exponential_mixture'],indent=2))
    if z.draws>0:
        g=batch_sweep(z.draws,z.reps,z.seed)
        print(json.dumps(g['model'],indent=2))
if __name__=='__main__':main()
