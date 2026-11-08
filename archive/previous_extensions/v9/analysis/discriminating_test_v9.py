"""Prospective-design stress test of within-clone bacterial history interventions.

Computational validation ONLY. Analytic model: exact observational equivalence
of selection turnover and Hill-transformed exponentially decaying internal state.
Independent simulator: stochastic single-cell attachment and detachment under
paired randomized surface-contact pulses with clone-level heterogeneity.

No simulated rates or treatment effects are measured bacterial quantities.
No raw experimental cell lineages are fitted. The published lags were used to
construct the illustrative equivalence, not a blind external validation.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.special import expit, logit
from scipy.stats import t as student_t
from scipy.stats import beta as beta_dist
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
L_NAIVE, L0, L14, L37 = 17.56, 4.68, 5.46, 13.41


def analytical_equivalence():
    """Exactly equal lag curves under different internal-state interpretations."""
    d0 = L_NAIVE - L0
    r14 = (L_NAIVE-L14)/d0
    r37 = (L_NAIVE-L37)/d0
    z14, z37 = 1/r14-1, 1/r37-1
    cost = brentq(lambda x: np.expm1(x*37)/np.expm1(x*14)-z37/z14,1e-6,.5,xtol=1e-14)
    f0 = 1-z14/np.expm1(cost*14)
    q=(1-f0)/f0
    tau = 10.0  # chosen coordinate; only the product n/tau is identified
    nu = cost*tau
    m0 = 0.85
    k=m0*q**(1/nu)

    def share(t):
        return 1/(1+q*np.exp(cost*np.asarray(t)))
    def state(t):
        return m0*np.exp(-np.asarray(t)/tau)
    def hill(m):
        return 1/(1+(k/np.asarray(m))**nu)
    def lag_s(t):
        return L_NAIVE-d0*share(t)/f0
    def lag_m(t):
        return L_NAIVE-d0*hill(state(t))/hill(state(0.))
    t = np.linspace(0,60,1201)
    maximum=float(np.max(np.abs(lag_s(t)-lag_m(t))))
    data={
      'source':'Published Lee et al 2018 reported lag point values used for algebra; not independent biological replicates',
      'selection_frequency_initial':float(f0),
      'relative_growth_disadvantage_per_h':float(cost),
      'intracellular_state_initial':float(m0),
      'intracellular_state_tau_h':tau,
      'hill_exponent':float(nu),
      'hill_midpoint':float(k),
      'lag_curve_max_abs_difference_h_0_to_60h':maximum,
      'holds_h':[0,7,14,21,30,37,48,60],
      'lags_selection_h':[float(lag_s(v)) for v in [0,7,14,21,30,37,48,60]],
      'lags_intracellular_h':[float(lag_m(v)) for v in [0,7,14,21,30,37,48,60]],
      'proof':'f(t)=[1+q exp(c t)]^-1 and h(m0 exp(-t/tau))=[1+(K/m0)^nu exp(nu t/tau)]^-1; set (K/m0)^nu=q and nu/tau=c',
      'scope':'Exact equivalence of phenomenological lag curves ONLY. Not equivalence of bacterial cAMP or lineage data.'
    }
    return data, {'t':t,'selection':lag_s(t),'internal':lag_m(t)}


def design_rank():
    """Expected contrasts of A,S,E,R under additive calibrated intervention idealization."""
    pre=np.array([[1,1,0,0],[0,0,1,1],[1,1,1,1]],float)
    full=np.vstack([pre,[1,0,0,0],[0,0,1,0]])
    no_pulse=np.vstack([pre,[0,0,1,0]])
    no_wash=np.vstack([pre,[1,0,0,0]])
    out={'column_order':['acquired_cell','selection','environment_trace','viable_carryover'],
         'pre_rows':['exposed_cells_on_virgin','fresh_cells_on_conditioned','both_history'],
         'additional_rows':['randomized_no_selection_contact_pulse','washed_and_viable_free_conditioned_surface'],
         'pre_matrix':pre.astype(int).tolist(),'full_matrix':full.astype(int).tolist(),
         'rank_pre':int(np.linalg.matrix_rank(pre)),'rank_full':int(np.linalg.matrix_rank(full)),
         'rank_without_pulse':int(np.linalg.matrix_rank(no_pulse)),
         'rank_without_washed_surface':int(np.linalg.matrix_rank(no_wash)),
         'assumption':'All measurements are on a common additive scale; pulse eliminates sorting/selection; washing eliminates all viable carryover while preserving deposited trace; no interactions.'}
    return out


def observable_probability(base_logit, phenotype_fraction, phenotype_logodds, response, det_h,
                           det_suppression, phenotype_det_suppression, assay_h=5/60):
    """Independent kinetic Bernoulli simulator marginalizing heterogeneous siblings.

    Early stable attachment comprises contact success times survival of detachment
    hazard through the observation interval. Pre-existing adhesive phenotype z is
    fixed per individual but variably represented in the cell suspension.
    """
    out = 0
    for z in (0,1):
        pa=expit(base_logit+phenotype_logodds*z+response)
        hazard=det_h*np.exp(-det_suppression*response-phenotype_det_suppression*z)
        survival=np.exp(-hazard*assay_h)
        out=out+(phenotype_fraction if z else 1-phenotype_fraction)*pa*survival
    return np.clip(out,0,1)


def dataset(model: str, n_experiments: int, blocks: int, n_cells: int, seed: int,
            selection_survival_odds: float=1., amplitude_scale:float=1., hold_h:float=0.):
    """Independent clone-block assay with genuinely randomized pulse vs sham.

    Assumes no lineage composition change when survival_odds=1. The selection
    stress test violates this with phenotypically differential survival/retention.
    The 'cell_coating' counterexample has EXACTLY the same measured kinetic law
    as an inherited intracellular state by construction (unmeasured carrier).
    """
    allowed={'intracellular_inherited','intracellular_no_inheritance','immutable_selection',
             'selection_with_loss','cell_coating','environment_only','residual_only'}
    if model not in allowed:raise ValueError(model)
    rng=np.random.default_rng(seed)
    N,K=n_experiments,blocks
    shape=(N,K)
    pbase=rng.uniform(.025,.065,(N,1))
    alpha=logit(pbase)+rng.normal(0,.65,shape)
    z0=expit(logit(rng.uniform(.08,.4,(N,1)))+rng.normal(0,.45,shape))
    beta_z=rng.uniform(.65,1.7,(N,1))
    det=rng.uniform(1.,4.,(N,1))
    phenotype_det=rng.uniform(.1,1.,(N,1))
    beta_mem=rng.uniform(.75,1.85,(N,1))*amplitude_scale
    eta=rng.uniform(.3,1.7,(N,1))
    m0=rng.uniform(.6,.95,(N,1))
    tau=rng.uniform(6.,20.,(N,1))
    midpoint=rng.uniform(.12,.5,(N,1))
    hill_exp=rng.uniform(1.3,5.,(N,1))
    inheritance_fraction=rng.uniform(.45,.9,(N,1))
    m=m0*np.exp(-hold_h/tau)
    h=m**hill_exp/(m**hill_exp+midpoint**hill_exp)
    mom_r=beta_mem*h
    daughter_state=m*inheritance_fraction
    daughter_h=daughter_state**hill_exp/(daughter_state**hill_exp+midpoint**hill_exp)
    daughter_r=beta_mem*daughter_h
    acquired=model in {'intracellular_inherited','intracellular_no_inheritance','cell_coating'}
    if not acquired:
        mom_r=np.zeros((N,1)); daughter_r=np.zeros((N,1))
    if model=='intracellular_no_inheritance':daughter_r=np.zeros((N,1))
    # Perturbed condition enriches pre-existing phenotypes if preparation is
    # selectively lossy. In true complete-recovery control the two fractions match.
    if model=='selection_with_loss':
        pulsed_fraction=expit(logit(z0)+np.log(selection_survival_odds))
    else:
        pulsed_fraction=z0
    p0=observable_probability(alpha,z0,beta_z,0,det,eta,phenotype_det)
    pm=observable_probability(alpha,pulsed_fraction,beta_z,mom_r,det,eta,phenotype_det)
    pd=observable_probability(alpha,pulsed_fraction,beta_z,daughter_r,det,eta,phenotype_det)
    # Sham sisters are independent Bernoulli draws conditional on shared clone.
    sham_m=rng.binomial(n_cells,p0)
    pulse_m=rng.binomial(n_cells,pm)
    sham_d=rng.binomial(n_cells,p0)
    pulse_d=rng.binomial(n_cells,pd)
    d_m=(pulse_m-sham_m)/n_cells
    d_d=(pulse_d-sham_d)/n_cells
    # One-sided matched t on biological clones; Bonferroni two endpoints.
    def evaluate(d):
        mean=d.mean(axis=1)
        se=d.std(axis=1,ddof=1)/np.sqrt(K)
        tstat=np.divide(mean,se,out=np.where(mean>0,np.inf,np.where(mean<0,-np.inf,0)),where=se>0)
        p=student_t.sf(tstat,df=K-1)
        return mean,p
    mean_m, pval_m=evaluate(d_m)
    mean_d, pval_d=evaluate(d_d)
    return {'mother_effect':mean_m,'daughter_effect':mean_d,
            'mother_detected':pval_m<.025,'daughter_detected':pval_d<.025,
            'joint_detected':(pval_m<.025)&(pval_d<.025),
            'true_mother_risk_effect':(pm-p0).mean(axis=1),
            'true_daughter_risk_effect':(pd-p0).mean(axis=1),
            'n_experiments':N,'n_blocks':K,'cells_per_arm_per_clone':n_cells,
            'model':model,'selection_survival_odds':selection_survival_odds,
            'hold_h':hold_h}


def summarize(r):
    n=r['n_experiments'];
    def interval(count):
        # Exact binomial interval over independently simulated datasets, not over cells.
        return [float(beta_dist.ppf(.025,count+1,n-count+1)),float(beta_dist.ppf(.975,count+1,n-count+1))]
    out={k:r[k] for k in ['model','n_experiments','n_blocks','cells_per_arm_per_clone','selection_survival_odds','hold_h']}
    for name in ['mother','daughter','joint']:
        count=int(r[f'{name}_detected'].sum())
        out[f'{name}_positive_count']=count
        out[f'{name}_positive_fraction']=count/n
        out[f'{name}_positive_beta_interval']=interval(count)
    for name in ['mother','daughter']:
        out[f'{name}_avg_observed_risk_difference_pp']=float(np.mean(r[f'{name}_effect'])*100)
        out[f'{name}_avg_true_risk_difference_pp']=float(np.mean(r[f'true_{name}_risk_effect'])*100)
    return out


def noisy_lag_classification(eq, n=2500, seed=202611081, noise_sd=.75):
    """A deliberately independent lag-only classifier negative control."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split
    rng=np.random.default_rng(seed)
    holds=np.array([0,5,11,14,21,28,37,48],float)
    pred=eq['intracellular_state_initial']
    q=(1-eq['selection_frequency_initial'])/eq['selection_frequency_initial']
    c=eq['relative_growth_disadvantage_per_h']
    lags=L_NAIVE-(L_NAIVE-L0)/(1+q*np.exp(c*holds))* (1+q)
    y=np.repeat([0,1],n)
    X=np.tile(lags,(2*n,1))+rng.normal(0,noise_sd,(2*n,len(holds)))
    train,test=train_test_split(np.arange(2*n),test_size=.35,random_state=seed,stratify=y)
    fit=LogisticRegression(C=1.,max_iter=2000,random_state=seed).fit(X[train],y[train])
    prob=fit.predict_proba(X[test])[:,1]
    return {'n_per_class':n,'holds_h':holds.tolist(),'measurement_noise_sd_assumed_h':noise_sd,
            'lag_only_holdout_roc_auc':float(roc_auc_score(y[test],prob)),
            'chance_auc_under_exact_equivalence':.5,
            'scope':'Monte Carlo sanity check of analytic equivalence, not a biological classifier benchmark'}


def sweep(n_sim=2000,seed=202611081,cell_n=48):
    eq,curves=analytical_equivalence()
    RESULTS.mkdir(parents=True,exist_ok=True)
    FIGURES.mkdir(parents=True,exist_ok=True)
    rank=design_rank()
    rows=[]
    models=['intracellular_inherited','intracellular_no_inheritance','immutable_selection',
            'cell_coating','selection_with_loss']
    for nblocks in [8,16,32,64]:
        for imod,model in enumerate(models):
            odds=8. if model=='selection_with_loss' else 1.
            rr=dataset(model,n_sim,nblocks,cell_n,seed+100000*nblocks+509*imod,selection_survival_odds=odds)
            rows.append(summarize(rr))
    # Explicit protocol violation severity and low-effect alternatives, same
    # prespecified test: no optimizing the significance threshold on outcomes.
    for odds in [1.,1.5,2.,4.,8.,16.]:
        rr=dataset('selection_with_loss',n_sim,32,cell_n,seed+4431+int(odds*100),selection_survival_odds=odds)
        row=summarize(rr);row['scenario']='selective_loss_sweep';rows.append(row)
    for amp in [.25,.5,.75,1.,1.5]:
        rr=dataset('intracellular_inherited',n_sim,32,cell_n,seed+7770+int(amp*1000),amplitude_scale=amp)
        row=summarize(rr);row['scenario']='strength_sweep';row['amplitude_scale']=amp;rows.append(row)
    for t in [0.,7.,14.,37.]:
        rr=dataset('intracellular_inherited',n_sim,32,cell_n,seed+8770+int(t),hold_h=t)
        row=summarize(rr);row['scenario']='washout_sweep';rows.append(row)
    fields=list(dict.fromkeys(k for d in rows for k in d))
    with (RESULTS/'intervention_power.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    summary={'analytical_equivalence':eq,'minimal_rank_design':rank,
             'lag_only_negative_control':noisy_lag_classification(eq,seed=seed+82),
             'field_level_simulations':{'n_datasets_each':n_sim,'cells_each_arm_per_clone':cell_n,
                 'alpha_per_test_bonferroni_one_sided':.025,
                 'independent_test_family':'Kinetic attachment-plus-survival Bernoulli agents with persistent phenotype heterogeneity and clone random intercepts',
                 'randomization':'Pulse/sham assigned within matched clonal blocks; full recovery required for selection-null contrast',
                 'result_file':'intervention_power.csv'},
             'selected_32block':{r['model']:r for r in rows if r['n_blocks']==32 and 'scenario' not in r},
             'leakage_32block':{str(r['selection_survival_odds']):r for r in rows if r.get('scenario')=='selective_loss_sweep'},
             'weak_effect_32block':{str(r['amplitude_scale']):r for r in rows if r.get('scenario')=='strength_sweep'},
             'hold_32block':{str(r['hold_h']):r for r in rows if r.get('scenario')=='washout_sweep'},
             'failure_boundaries':['The same kinetic signature can be generated by an inherited extracellular coating carried on cells; the randomized clone experiment does not identify intracellular location.',
                 'Differential survival/recovery during the surface pulse can enrich pre-existing phenotypes and yield a false apparent within-clone acquired effect.',
                 'Functional and cross-reactive bacterial assay behavior is uncalibrated; true biological effect sizes are unknown.',
                 'The four-carrier rank criterion is exact only in an additive readout and with genuinely carrier-selective manipulations.']}
    (RESULTS/'critical_test_v9.json').write_text(json.dumps(summary,indent=2))
    return summary,rows,curves


def figures(summary,rows,curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=curves['t']
    fig,ax=plt.subplots(figsize=(7.4,4.1))
    ax.plot(t,curves['selection'],lw=2.5,label='Selection with differential growth',color='#B35F22')
    ax.plot(t,curves['internal'],ls='--',lw=1.6,label='Decaying intracellular state + Hill response',color='#175C82')
    ax.scatter([0,14,37],[L0,L14,L37],s=44,facecolor='white',edgecolor='#1C1C1C',zorder=5,label='Lee et al. reported lag fits')
    ax.axhline(L_NAIVE,c='gray',ls=':',alpha=.65)
    ax.set(xlabel='Suspension holding time (h)',ylabel='Phenomenological colonization lag (h)',xlim=(0,60),ylim=(3.5,18.3))
    ax.legend(fontsize=8,frameon=False); fig.tight_layout()
    fig.savefig(FIGURES/'exact_observational_equivalence.png',dpi=220);plt.close(fig)
    models=['intracellular_inherited','intracellular_no_inheritance','immutable_selection','selection_with_loss','cell_coating']
    labels=['Inherited state','State, no inheritance','Fixed selection','Selection + biased recovery','Inherited coating']
    rr=[summary['selected_32block'][m] for m in models]
    x=np.arange(len(models));w=.36
    fig,ax=plt.subplots(figsize=(9,4.7))
    ax.bar(x-w/2,[z['mother_positive_fraction'] for z in rr],width=w,label='Pulsed mother vs sham',color='#3E7A92')
    ax.bar(x+w/2,[z['daughter_positive_fraction'] for z in rr],width=w,label='Immediate daughter vs sham',color='#C27C44')
    ax.axhline(.025,ls=':',color='black',lw=1,label='Per-endpoint nominal false positive (one-sided)')
    ax.set(ylabel='Fraction of synthetic experiments declaring positive contrast',ylim=(0,1.07),xticks=x,xticklabels=labels)
    ax.tick_params(axis='x',labelrotation=16,labelsize=9)
    ax.legend(fontsize=8,frameon=False,loc='upper right'); fig.tight_layout()
    fig.savefig(FIGURES/'paired_clone_discrimination.png',dpi=220);plt.close(fig)
    odds=[1.,1.5,2.,4.,8.,16.]
    r=[summary['leakage_32block'][str(o)] for o in odds]
    fig,ax=plt.subplots(figsize=(7.3,4.4))
    ax.plot(odds,[x['mother_positive_fraction'] for x in r],marker='o',label='Mother contrast')
    ax.plot(odds,[x['daughter_positive_fraction'] for x in r],marker='s',label='Daughter contrast')
    ax.plot(odds,[x['joint_positive_fraction'] for x in r],marker='^',label='Both contrasts')
    ax.axhline(.025,ls=':',color='black',lw=1)
    ax.set(xscale='log',xlabel='High-adhesion phenotype survival odds ratio during preparation',ylabel='False acquired-state declaration fraction',ylim=(0,1.03))
    ax.set_xticks(odds,labels=[str(a) for a in odds]);ax.legend(frameon=False);fig.tight_layout()
    fig.savefig(FIGURES/'selection_leakage_failure.png',dpi=220);plt.close(fig)
    figure_rows=[z for z in rows if z['n_blocks'] in [8,16,32,64] and 'scenario' not in z]
    fig,ax=plt.subplots(figsize=(7.3,4.4))
    for m,label,mark in [('intracellular_inherited','Inherited scalar state','o'),('intracellular_no_inheritance','Noninherited scalar state','s'),('immutable_selection','Selection negative control','^')]:
        r=sorted([z for z in figure_rows if z['model']==m],key=lambda z:z['n_blocks'])
        ax.plot([z['n_blocks'] for z in r],[z['daughter_positive_fraction'] for z in r],marker=mark,label=label)
    ax.axhline(.025,color='black',ls=':',lw=1)
    ax.set(xlabel='Independent clonal blocks',ylabel='Daughter-pulse positive fraction',ylim=(0,1.03),xticks=[8,16,32,64]);ax.legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(FIGURES/'lineage_block_power.png',dpi=220);plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--datasets',type=int,default=2000)
    parser.add_argument('--cells-per-arm-clone',type=int,default=48)
    parser.add_argument('--seed',type=int,default=202611081)
    args=parser.parse_args()
    s,rows,curves=sweep(args.datasets,args.seed,args.cells_per_arm_clone)
    figures(s,rows,curves)
    print('max lag-curve difference',s['analytical_equivalence']['lag_curve_max_abs_difference_h_0_to_60h'])
    print('contrasts rank',s['minimal_rank_design']['rank_pre'],s['minimal_rank_design']['rank_full'])
    print('lag-only classifier AUC',s['lag_only_negative_control']['lag_only_holdout_roc_auc'])
    for model,v in s['selected_32block'].items():
        print(model,'mother',round(v['mother_positive_fraction'],4),'daughter',round(v['daughter_positive_fraction'],4),'joint',round(v['joint_positive_fraction'],4),'avg_true_pp',(round(v['mother_avg_true_risk_difference_pp'],2),round(v['daughter_avg_true_risk_difference_pp'],2)))
    print('loss contamination',[ (x,round(z['joint_positive_fraction'],3)) for x,z in s['leakage_32block'].items()])
    print('low power',[(x,round(z['joint_positive_fraction'],3)) for x,z in s['weak_effect_32block'].items()])
    print('holding',[ (x,round(z['joint_positive_fraction'],3)) for x,z in s['hold_32block'].items()])

def five_carrier_rank():
    """Prove that inherited surface coating defeats the sister assay as a carrier locator.

    Columns: intracellular A, selected S, environmental E, residual R,
    cell-associated coating C. Sister and mother effects carry A+C alike.
    """
    rows=np.array([[1,1,0,0,1], # conventional cell transfer
                   [0,0,1,1,0], # conditioned recipient
                   [1,0,0,0,1], # randomized contact pulse
                   [0,0,1,0,0], # cleared conditioned recipient
                   [1,0,0,0,1]],float) # daughter assay under full inherited coating
    with_localization=np.vstack([rows,[1,0,0,0,0]])
    return {'columns':['intracellular','selected_type','environmental_trace','residual_viable_cell','cell_coating'],
            'rank_with_sister_and_surface_controls':int(np.linalg.matrix_rank(rows)),
            'rank_with_one_selective_localization_assay':int(np.linalg.matrix_rank(with_localization)),
            'null_direction':[1,0,0,0,-1],
            'caution':'Selective intracellular-versus-cell-coating perturbation has not been validated; molecular colocalization is suggestive but may not equal selective intervention.'}


def calibration_geometry():
    """Noisy linear four-carrier design ill-conditioning with imperfect resets."""
    records=[]
    for kappa in [0.,.2,.5,.8,.95,.99]:
        for q,rho in [(1.,0.),(.9,.1),(.75,.25),(.6,.5)]:
            M=np.array([[1,1,0,0],[0,0,1,1],[1,kappa,0,0],[0,0,q,rho]],float)
            cond=float(np.linalg.cond(M))
            sigma=.12/np.sqrt(32.)
            uncertainty=np.sqrt(np.sum(np.linalg.inv(M)**2,axis=1))*sigma
            records.append({'selection_leakage_kappa':kappa,'trace_retention_q':q,'viable_cell_retention_rho':rho,
                            'matrix_condition_number':cond,
                            'acquired_95pct_halfwidth_at_32_independent_fields':float(1.96*uncertainty[0]),
                            'environment_95pct_halfwidth_at_32_independent_fields':float(1.96*uncertainty[2]),
                            'determinant':float(np.linalg.det(M))})
    return records


if __name__=='__main__':
    main()
