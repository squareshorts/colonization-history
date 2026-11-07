"""Identifiability, missing-outcome bounds and held-out experimental discrimination.

An explicitly synthetic and deliberately adversarial theoretical experiment.
No biological rates are inferred or calibrated from Lee et al. 2018.
The shape-conjugacy result is a standard mathematical reparameterization.

Main experiments:
(1) Monotone latent-path invariance with arbitrary observation functions.
(2) Sharp missing-data bounds for binary, tracked mother/daughter units.
(3) Finite-sample synthetic paired-clone challenge in two independent kinetic laws.
(4) Location-specific fifth-carrier assay and interval-specificity uncertainty.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import expit,logit
from scipy.stats import t as t_dist, norm

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'; FIGURES=ROOT/'figures'


def latent_selection(t, f0=.31, cost=.12):
    return expit(logit(f0)-cost*np.asarray(t))


def latent_state(t, m0=.83, tau=9.3):
    return m0*np.exp(-np.asarray(t)/tau)


def lag_family(t, family:str, params=None):
    t=np.asarray(t,dtype=float)
    lo,hi=4.4,17.8
    if family=='weibull': response=1-np.exp(-(t/25.)**2.2)
    elif family=='stretched':response=1-np.exp(-(t/30.)**.65)
    elif family=='biexponential':response=1-(.53*np.exp(-t/7.)+.47*np.exp(-t/44.))
    elif family=='logistic':
        b=expit(-4.1); response=(expit((t-26.)/6.)-b)/(1-b)
    elif family=='delayed':response=1-np.exp(-(np.maximum(t-9.,0)/20.)**1.5)
    else:raise ValueError(family)
    return lo+(hi-lo)*response


def reparameterization_stress(n_trials=300, seed=20261110, n_train=25, n_test=600):
    """Disjoint pseudo-times; transformations re-learn arbitrary smooth readouts.
    The target family is a synthetic known function, not external biological data.
    Reports interpolation error; exact equivalence is an analytic identity.
    """
    rng=np.random.default_rng(seed)
    families=['weibull','stretched','biexponential','logistic','delayed']
    records=[]
    for i in range(n_trials):
        family=families[i%len(families)]
        t_train=np.linspace(0,60,n_train)
        t_test=rng.uniform(0,60,n_test)
        f0=rng.uniform(.06,.9); c=rng.uniform(.04,.32)
        m0=rng.uniform(.25,.95); tau=rng.uniform(5.,30.)
        truth=lag_family(t_test,family)
        x_m=latent_state(t_train,m0,tau); x_f=latent_selection(t_train,f0,c)
        y=lag_family(t_train,family)
        # Sorting is required for a spline defined in the latent coordinate.
        r_m=PchipInterpolator(x_m[::-1],y[::-1])(latent_state(t_test,m0,tau))
        r_f=PchipInterpolator(x_f[::-1],y[::-1])(latent_selection(t_test,f0,c))
        err=float(np.sqrt(np.mean((r_m-truth)**2)))
        errf=float(np.sqrt(np.mean((r_f-truth)**2)))
        # If a common fitted latent-time response had been the only observation,
        # the two models' noise-free predictons differ through interpolation only.
        records.append({'replicate':i,'family':family,'rmse_state_h':err,
                        'rmse_selection_h':errf,
                        'max_between_readouts_h':float(np.max(np.abs(r_m-r_f)))})
    return records


def missing_bounds(success_recovered, n_recovered, n_assigned):
    """Exact finite-population identification bounds for an outcome in [0,1].

    Known: assigned N, recovered n, sum S of observed binary successes.
    Unknown: potential binary outcomes of N-n lost, independently of why lost.
    The identified set for the full-assigned mean is [S/N,(S+N-n)/N].
    These are sharp: assign all missing outcomes zero / one.
    """
    N=np.asarray(n_assigned,dtype=float)
    n=np.asarray(n_recovered,dtype=float)
    s=np.asarray(success_recovered,dtype=float)
    if np.any(N<=0) or np.any(n<0) or np.any(n>N) or np.any(s<0) or np.any(s>n):
        raise ValueError('invalid assigned/recovered/success counts')
    return s/N,(s+N-n)/N


def difference_bounds(success_pulse,recovered_pulse,assigned_pulse,
                      success_sham,recovered_sham,assigned_sham):
    pa,pb=missing_bounds(success_pulse,recovered_pulse,assigned_pulse)
    ca,cb=missing_bounds(success_sham,recovered_sham,assigned_sham)
    return pa-cb,pb-ca


def clone_inference(success_pulse,recovered_pulse,success_sham,recovered_sham,
                    N=48,alpha=.025):
    """Cluster-level t inference: approximate sampling error, not distribution-free.

    Returns unadjusted recovered-only contrast and attrition-robust lower-limit
    of full randomized-group contrast. The latter is valid as an identified
    bound, but its t sampling interval needs independent blocks / approximation.
    """
    p=np.divide(success_pulse,recovered_pulse,
                out=np.zeros_like(success_pulse,dtype=float),where=np.asarray(recovered_pulse)>0)
    c=np.divide(success_sham,recovered_sham,
                out=np.zeros_like(success_sham,dtype=float),where=np.asarray(recovered_sham)>0)
    d=p-c
    lb,ub=difference_bounds(success_pulse,recovered_pulse,N,success_sham,recovered_sham,N)
    def low95(x):
        k=x.shape[-1]
        m=x.mean(axis=-1)
        se=x.std(axis=-1,ddof=1)/np.sqrt(k)
        return m-t_dist.ppf(1-alpha,k-1)*se
    return {'raw_effect':d.mean(axis=-1),'raw_detect':low95(d)>0,
            'bound_mean_lower':lb.mean(axis=-1),'bound_mean_upper':ub.mean(axis=-1),
            'bound_ci_lower':low95(lb),'bound_detect':low95(lb)>0,
            'mean_recovery_pulse':np.mean(recovered_pulse/N,axis=-1),
            'mean_recovery_sham':np.mean(recovered_sham/N,axis=-1)}


def kinetic_rates(alpha,z,beta_z,mem,model,det_rate,eta,phenotype_det):
    """Two independent contact/dwell observation laws, neither fitted to data.

    Outcome is stable attachment at five minutes after an initial contact.
    'exp' uses logistic contact and exponential detachment; 'weibull' uses
    complementary log-log contact and Weibull detachment. Both have monotone
    effects but quantitatively different nuisance structure.
    """
    time=5/60
    if model=='exp':
        attach=expit(alpha+beta_z*z+mem)
        det=det_rate*np.exp(-eta*mem-phenotype_det*z)
        surv=np.exp(-det*time)
    elif model=='weibull':
        attach=-np.expm1(-np.exp(alpha+beta_z*z+mem))
        shape=1.55
        det=det_rate*np.exp(-eta*mem-phenotype_det*z)
        surv=np.exp(-(det*time)**shape)
    else:raise ValueError(model)
    return np.clip(attach*surv,0,1)


def simulate_clones(carrier,kinetic,n_studies=1500,blocks=32,assigned=48,seed=20261111,
                    recovery_high=1.,recovery_low=1.,mem_scale=1.,
                    hold_h=0.):
    """Randomized within-clone pulse/sham comparison of traceable founder units.

    One prespecified daughter per surviving founder; no population proliferation.
    Complete recovery implies no compositional sorting by construction.
    Differential survival is allowed and is *recorded* per assigned arm.
    """
    allowed={'inherited_state','no_daughter_inheritance','pure_selection',
             'phenotype_biased_recovery','inherited_coating','induced_switch'}
    if carrier not in allowed:raise ValueError(carrier)
    rng=np.random.default_rng(seed)
    shape=(n_studies,blocks)
    zfrac=expit(logit(rng.uniform(.06,.50,(n_studies,1)))+rng.normal(0,.45,shape))
    base=rng.uniform(-3.5,-1.7,(n_studies,1))+rng.normal(0,.35,shape)
    beta_z=rng.uniform(.5,1.65,(n_studies,1))
    det_rate=rng.uniform(.6,4.,(n_studies,1))
    eta=rng.uniform(.25,1.2,(n_studies,1))
    phenotype_det=rng.uniform(.10,.90,(n_studies,1))
    m0=rng.uniform(.6,.95,(n_studies,1))
    tau=rng.uniform(6,20,(n_studies,1))
    med=rng.uniform(.15,.50,(n_studies,1))
    expn=rng.uniform(1.5,4.5,(n_studies,1))
    mem=m0*np.exp(-hold_h/tau)
    response=mem**expn/(mem**expn+med**expn)
    effect=rng.uniform(.7,1.7,(n_studies,1))*response*mem_scale
    mom=effect if carrier in {'inherited_state','no_daughter_inheritance','inherited_coating','induced_switch'} else np.zeros_like(effect)
    daughter=effect*rng.uniform(.4,.9,(n_studies,1)) if carrier in {'inherited_state','inherited_coating','induced_switch'} else np.zeros_like(effect)
    # 'induced_switch' deliberately represents acquisition of an adhesive
    # phenotype without identifying an intracellular molecular store.
    # Its kinetic effect can be indistinguishable from an inherited state.
    nhA=rng.binomial(assigned,zfrac)
    nhB=rng.binomial(assigned,zfrac)
    if carrier=='phenotype_biased_recovery':
        rhi=float(recovery_high); rlo=float(recovery_low)
    else:rhi,rlo=1.,1.
    # Sham is a complete-recovery control; experimental pulse may lose cells.
    pulse_hi=rng.binomial(nhA,rhi)
    pulse_lo=rng.binomial(assigned-nhA,rlo)
    sham_hi=nhB;sham_lo=assigned-nhB
    pulse_n=pulse_hi+pulse_lo
    sham_n=sham_hi+sham_lo
    
    def score(hi,lo,boost):
        p_hi=kinetic_rates(base,1,beta_z,boost,kinetic,det_rate,eta,phenotype_det)
        p_lo=kinetic_rates(base,0,beta_z,boost,kinetic,det_rate,eta,phenotype_det)
        return rng.binomial(hi,p_hi)+rng.binomial(lo,p_lo)
    results={}
    for who,boost in [('mother',mom),('daughter',daughter)]:
        # Daughter endpoints tracked one-to-one to assigned mothers. Units lost
        # during preparation have unobserved daughter outcome; no replacement.
        sA=score(pulse_hi,pulse_lo,boost)
        sB=score(sham_hi,sham_lo,0.)
        results[who]=clone_inference(sA,pulse_n,sB,sham_n,N=assigned,alpha=.025)
    data={
       'carrier':carrier,'kinetic':kinetic,'blocks':blocks,'cells_per_arm_clone':assigned,
       'n_studies':n_studies,'recovery_high':rhi,'recovery_low':rlo,
       'mother':results['mother'],'daughter':results['daughter'],
       'mean_true_exposure_effect_amplitude':float(np.mean(effect)) if carrier in {'inherited_state','no_daughter_inheritance','inherited_coating','induced_switch'} else 0.,
    }
    return data


def summarize_sim(data):
    d={'carrier':data['carrier'],'kinetic':data['kinetic'],
       'n_studies':data['n_studies'],'blocks':data['blocks'],
       'cells_per_arm_clone':data['cells_per_arm_clone'],
       'recovery_high':data['recovery_high'],'recovery_low':data['recovery_low']}
    for who in ('mother','daughter'):
        for key in ['raw_effect','raw_detect','bound_mean_lower','bound_detect','mean_recovery_pulse','mean_recovery_sham']:
            d[f'{who}_{key}']=float(np.mean(data[who][key]))
    for method in ('raw','bound'):
        d[f'{method}_joint_detect']=float(np.mean(data['mother'][f'{method}_detect']&data['daughter'][f'{method}_detect']))
    return d


def carrier_matrix(gamma=1.,kappa=0.):
    """Five carrier columns (A,S,C,E,R); exterior readout affected by gamma*A+kappa*C."""
    return np.array([[1,1,1,0,0],[0,0,0,1,1],[1,0,1,0,0],[0,0,0,1,0],
                     [gamma,0,kappa,0,0]],float)


def bounded_AC(U_interval,V_interval,gamma_interval,kappa_interval):
    """Sharp extrema of A=(V-kU)/(gamma-k), C=U-A on rectangular inputs,
    restricted to physically admissible nonnegative A,C.
    
    For gamma_min>kappa_max, coordinate extrema are attained at corners
    of feasible input boxes; to handle clipped feasibility explicitly, enumerate
    vertices by optimizing bounded inequalities with SLSQP as fallback. The
    closed-form unconstrained corner envelope is guaranteed conservative.
    """
    ul,uh=U_interval;vl,vh=V_interval;gl,gh=gamma_interval;kl,kh=kappa_interval
    if not(gl>kh>=0 and vl<=vh and ul<=uh):return {'identified':False,'reason':'overlapping or invalid calibration intervals'}
    vals=[]
    for u in (ul,uh):
        for v in (vl,vh):
            for g in (gl,gh):
                for k in (kl,kh):
                    a=(v-k*u)/(g-k);c=u-a
                    vals.append((a,c))
    a_min,a_max=min(x[0] for x in vals),max(x[0] for x in vals)
    c_min,c_max=min(x[1] for x in vals),max(x[1] for x in vals)
    return {'identified':True,'A_interval_conservative':[max(0.,a_min),min(1.,a_max)],
            'C_interval_conservative':[max(0.,c_min),min(1.,c_max)],
            'gamma_range':[gl,gh],'kappa_range':[kl,kh]}


def location_sweep(seed=20261112,n=3000,n_fields=32,sigma=.12):
    rng=np.random.default_rng(seed)
    true=np.array([.25,.4,.18,.14,.32]) # (A,S,C,E,R)
    scenarios=[(1.,0.),(1.,.2),(1.,.6),(1.,.9),(1.,.98),(1.,1.),(.98,.2)]
    records=[]
    z=norm.ppf(1-.05/4)
    for gamma,kappa in scenarios:
        mat=carrier_matrix(gamma,kappa)
        y=mat@true
        Y=y+rng.normal(0,sigma/np.sqrt(n_fields),(n,5))
        if np.linalg.matrix_rank(mat)==5:
            xhat=np.linalg.solve(mat,Y.T).T
            rmse=float(np.sqrt(np.mean((xhat[:,0]-true[0])**2)))
            fp=float(np.mean(xhat[:,0]<0))
        else:rmse=None;fp=None
        records.append({'gamma':gamma,'kappa':kappa,'separation':gamma-kappa,
                        'rank':int(np.linalg.matrix_rank(mat)),
                        'condition_number':float(np.linalg.cond(mat)),
                        'rmse_A':rmse,'negative_estimate_rate':fp,
                        'fields':n_fields,'field_sd':sigma})
    # Interval-calibrated assay example with full uncertainty in gamma and kappa.
    gamma=.98;kappa=.12
    y=carrier_matrix(gamma,kappa)@true
    half=z*sigma/np.sqrt(n_fields)
    U,V=y[2],y[4]
    ex=bounded_AC((U-half,U+half),(V-half,V+half),(.95,1.),(0.,.2))
    ex.update({'true_A':float(true[0]),'true_C':float(true[2]),
               'observed_nominal_U':float(U),'observed_nominal_V':float(V),
               'simultaneous_normal_interval_halfwidth':float(half),
               'field_sd':sigma,'n_fields':n_fields,
               'warning':'Only valid if specificity intervals, noise law and additive response are validated; no biological calibration supplied.'})
    return records,ex



def location_interval_coverage(n_sim=4000, seed=20261114, fields=32, sigma=.12):
    """Sensitivity-interval coverage under valid and deliberately invalid assay calibration.
    Report *joint* nominal Gaussian coverage for two independent mean readouts;
    percentages are conditional on exactly the stated Gaussian model.
    """
    rng=np.random.default_rng(seed)
    A,C=.25,.18
    U=A+C
    half=norm.ppf(1-.05/4)*sigma/np.sqrt(fields)
    scenarios=[('valid_calibration',.98,.12),
               ('invalid_exterior_activity',1.,.65),
               ('invalid_intracellular_activity',.60,.10)]
    out=[]
    for label,g,k in scenarios:
        V=g*A+k*C
        Yu=rng.normal(U,sigma/np.sqrt(fields),n_sim)
        Yv=rng.normal(V,sigma/np.sqrt(fields),n_sim)
        indicators=[]
        positive=[]
        for u,v in zip(Yu,Yv):
            b=bounded_AC((u-half,u+half),(v-half,v+half),(.95,1.),(0.,.2))
            a0,a1=b['A_interval_conservative']
            c0,c1=b['C_interval_conservative']
            indicators.append(a0-1e-10<=A<=a1+1e-10 and c0-1e-10<=C<=c1+1e-10)
            positive.append(a0>0)
        out.append({'scenario':label,'true_gamma':g,'true_kappa':k,
                    'calibration_gamma':[.95,1.],'calibration_kappa':[0.,.2],
                    'n_sim':n_sim,'fields':fields,'coverage_both_A_C':float(np.mean(indicators)),
                    'declares_positive_A':float(np.mean(positive)),
                    'CI_bonferroni_halfwidth':half})
    return out


def run(n_trials=300, n_studies=1500):
    RESULTS.mkdir(parents=True,exist_ok=True); FIGURES.mkdir(parents=True,exist_ok=True)
    equivalences=reparameterization_stress(n_trials=n_trials)
    with (RESULTS/'v10_latent_holdout.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(equivalences[0]));w.writeheader();w.writerows(equivalences)
    allrows=[]
    pairs=[('inherited_state',1.,1.),('no_daughter_inheritance',1.,1.),
           ('pure_selection',1.,1.),('phenotype_biased_recovery',.98,.45),
           ('phenotype_biased_recovery',.998,.98),('inherited_coating',1.,1.),
           ('induced_switch',1.,1.)]
    for kinetic in ('exp','weibull'):
        for i,(carrier,rhi,rlo) in enumerate(pairs):
            seed=20261113+104*i+(0 if kinetic=='exp' else 10000)
            data=simulate_clones(carrier,kinetic,n_studies=n_studies,seed=seed,
                                 recovery_high=rhi,recovery_low=rlo)
            row=summarize_sim(data);allrows.append(row)
    with (RESULTS/'v10_kinetic_stress.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    location,example=location_sweep()
    location_coverage=location_interval_coverage()
    with (RESULTS/'v10_location_selectivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(location[0]));w.writeheader();w.writerows(location)
    summary={'analytic_theorem':'For an injective latent trajectory x(t), any scalar target y(t) has unique along-path readout g(x)=y(x_inverse(x)). This is a standard reparameterization, not new mathematics.',
       'scope':'Observational non-identifiability ONLY for the selected scalar endpoint and unconstrained readout; additional signals may identify carriers.',
       'monotone_output_families':sorted(set(r['family'] for r in equivalences)),
       'holdout_n':len(equivalences),'holdout_rmse_state_mean_h':float(np.mean([r['rmse_state_h'] for r in equivalences])),
       'holdout_rmse_selection_mean_h':float(np.mean([r['rmse_selection_h'] for r in equivalences])),
       'holdout_max_between_readouts_h':float(max(r['max_between_readouts_h'] for r in equivalences)),
       'kinetic_rows':allrows,'location_rows':location,'location_example':example,
       'location_interval_coverage':location_coverage,
       'attrition_identification':'P(R)=1-epsilon; observed S successes in n recovered of N assigned. Full-group success probability in [S/N,(S+N-n)/N]. Difference bounds subtract arm extremes. The maximal difference between recovered and full-arm means is epsilon, not a bound on randomization imbalance. Mother/daughter bounds need one-to-one tracking and no differential proliferation.',
       'sources':{'HorowitzManski2000':'10.1080/01621459.2000.10473902','Manski2005':'10.1016/j.ijar.2004.10.006','Meshkat2011':'10.1016/j.mbs.2011.06.001'},
       'not_biological_data':True}
    (RESULTS/'v10_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    make_figures(allrows,equivalences,location)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('kinetic_rows','location_rows','location_example','location_interval_coverage')},indent=2))
    print('LOCATION COVERAGE',json.dumps(location_coverage,indent=2))
    for v in allrows:print(v['kinetic'],v['carrier'],v['recovery_high'],v['recovery_low'],
                           'raw mother',round(v['mother_raw_detect'],3),
                           'raw d',round(v['daughter_raw_detect'],3),
                           'bound mother',round(v['mother_bound_detect'],3),
                           'bound d',round(v['daughter_bound_detect'],3))
    return summary


def make_figures(rows,eq,location):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5.3))
    families=['weibull','stretched','biexponential','logistic','delayed']
    t=np.linspace(0,60,300)
    for f in families:ax.plot(t,lag_family(t,f),label=f)
    ax.set(xlabel='Holding interval (arbitrary hours)',ylabel='Synthetic lag (arbitrary hours)',
           title='Different monotone lag shapes admit equivalent latent representations')
    ax.legend(fontsize=8,ncol=2);fig.tight_layout()
    fig.savefig(FIGURES/'v10_general_lag_classes.png',dpi=180);plt.close(fig)
    # raw vs attrition robust conditional detection across kinetics and carrier.
    fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True)
    labels=['Inherited state','No daughter inheritance','Selection, complete recovery','Selection, high attrition',
            'Selection, <2% loss','Inherited coating','Induced phenotype switch']
    for ax,k in zip(axs,['exp','weibull']):
        vals=[r for r in rows if r['kinetic']==k]
        xpos=np.arange(len(vals))
        ax.bar(xpos-.2,[100*r['raw_joint_detect'] for r in vals],width=.4,label='Recovered-only test')
        ax.bar(xpos+.2,[100*r['bound_joint_detect'] for r in vals],width=.4,label='Attrition-bounded test')
        ax.set_title('Exponential dwell' if k=='exp' else 'Weibull dwell')
        ax.set_xticks(xpos,labels,rotation=60,ha='right',fontsize=8)
        ax.set_ylim(0,105);ax.grid(axis='y',alpha=.2)
    axs[0].set_ylabel('Joint mother and daughter detections (%)')
    axs[0].legend(fontsize=8)
    fig.suptitle('Synthetic cross-kinetic discrimination: failures remain explicit',fontsize=12)
    fig.tight_layout();fig.savefig(FIGURES/'v10_cross_kinetic_detections.png',dpi=190,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.2,4.8))
    ordered=sorted((z for z in location if z['rank']==5),key=lambda z:z['separation'],reverse=True)
    x=[z['separation'] for z in ordered]
    y=[z['rmse_A'] for z in ordered]
    ax.plot(x,y,'o-',label='RMSE of estimated intracellular component')
    ax.invert_xaxis();ax.set(xlabel='Difference in assay response to intracellular and exterior carriers',
       ylabel='Synthetic estimate RMSE',title='Localization contrast becomes unstable near equal sensitivities')
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(FIGURES/'v10_localization_conditioning.png',dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=300);p.add_argument('--studies',type=int,default=1500)
    a=p.parse_args();run(a.trials,a.studies)
