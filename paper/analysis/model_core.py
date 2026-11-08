"""Paired-clone power under two-arm recovery and single-cell information limits.

All parameter distributions are illustrative. No raw bacteria were simulated
or measured in a laboratory. Simulated original-assignment potential outcomes
are used only to verify the sharp missing-outcome bounds, never for detection.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
from scipy.special import expit, logit
from scipy.stats import t, norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
FIG = ROOT / 'figures'
KINETICS = ('logistic_exponential','cloglog_weibull')
SCENARIOS = {
 'pulse_only_uniform': ('uniform',1.,0.),
 'both_uniform_equal': ('uniform',1.,1.),
 'sham_half_uniform': ('uniform',1.,.5),
 'sham_more_uniform': ('uniform',.5,1.),
 'pulse_only_phenotype': ('phenotype',1.,0.),
 'both_phenotype_equal': ('phenotype',1.,1.),
 'sham_half_phenotype': ('phenotype',1.,.5),
}


def attachment_probability(alpha, phenotype, beta, effect, rate, phenotype_detach, kinetic):
    predictor = alpha + beta * phenotype + effect
    hazard = rate * np.exp(-.70*effect - phenotype_detach*phenotype) * (5./60.)
    if kinetic == 'logistic_exponential':
        probability = expit(predictor)
        survival = np.exp(-hazard)
    elif kinetic == 'cloglog_weibull':
        probability = -np.expm1(-np.exp(predictor))
        survival = np.exp(-hazard**1.55)
    else:
        raise ValueError(kinetic)
    return np.clip(probability*survival,0.,1.)


def positive_cluster_t(diffs, alpha=.025):
    # Approximate one-sided independent-clone test; blocks are the units.
    k = diffs.shape[1]
    se = diffs.std(axis=1,ddof=1)/np.sqrt(k)
    return diffs.mean(axis=1) - t.ppf(1-alpha,k-1)*se > 0


def simulate(rng, scenario, nominal_loss, effect, kinetic, studies=1200, blocks=32, founders=48):
    kind, pulse_mult, sham_mult = SCENARIOS[scenario]
    shape=(studies,blocks)
    f=expit(logit(rng.uniform(.07,.48,(studies,1)))+rng.normal(0,.45,shape))
    alpha=rng.uniform(-3.4,-1.8,(studies,1))+rng.normal(0,.3,shape)
    beta=rng.uniform(.5,1.6,(studies,1))
    rate=rng.uniform(.6,4,(studies,1))
    detach=rng.uniform(.1,.9,(studies,1))
    h1=rng.binomial(founders,f);l1=founders-h1
    h0=rng.binomial(founders,f);l0=founders-h0

    lp=nominal_loss*pulse_mult; ls=nominal_loss*sham_mult
    if kind == 'uniform':
        h1r=rng.binomial(h1,1-lp);l1r=rng.binomial(l1,1-lp)
        h0r=rng.binomial(h0,1-ls);l0r=rng.binomial(l0,1-ls)
    else:
        h1r=h1.copy();l1r=rng.binomial(l1,1-lp)
        h0r=h0.copy();l0r=rng.binomial(l0,1-ls)
    n1=h1r+l1r;n0=h0r+l0r
    def draw(h,l,shift):
        ph=attachment_probability(alpha,1.,beta,shift,rate,detach,kinetic)
        pl=attachment_probability(alpha,0.,beta,shift,rate,detach,kinetic)
        return rng.binomial(h,ph)+rng.binomial(l,pl)

    ans={'kinetic':kinetic,'scenario':scenario,'nominal_loss':nominal_loss,
         'requested_pulse_loss':lp,'requested_sham_loss':ls,
         'actual_pulse_loss':float(np.mean(1-n1/founders)),
         'actual_sham_loss':float(np.mean(1-n0/founders)),
         'effect':effect,'studies':studies,'blocks':blocks,'founders_per_arm_block':founders}
    tests={};coverage=[]
    for endp,shift in [('mother',effect),('daughter',.70*effect)]:
        s1=draw(h1r,l1r,shift);s0=draw(h0r,l0r,0.)
        # Only for verifying that the original-assignment contrast is in the bounds.
        missing1=draw(h1-h1r,l1-l1r,shift)
        missing0=draw(h0-h0r,l0-l0r,0.)
        true_diff=(s1+missing1-s0-missing0)/founders
        low=s1/founders-(s0+founders-n0)/founders
        upper=(s1+founders-n1)/founders-s0/founders
        coverage.append(bool(np.all((true_diff>=low-1e-12)&(true_diff<=upper+1e-12))))
        raw=s1/np.maximum(n1,1)-s0/np.maximum(n0,1)
        tests[endp]={'bounded':positive_cluster_t(low),'recovered':positive_cluster_t(raw)}
        ans[endp+'_bound_mean']=float(np.mean(low))
        ans[endp+'_observed_mean']=float(np.mean(raw))
    ans['exact_bounds_cover_generated_full_data']=bool(all(coverage))
    for method in ('bounded','recovered'):
        ans[method+'_joint']=float(np.mean(tests['mother'][method]&tests['daughter'][method]))
        ans[method+'_mother']=float(np.mean(tests['mother'][method]))
        ans[method+'_daughter']=float(np.mean(tests['daughter'][method]))
    return ans


def paired_probabilities(hold,p_hi,p_lo,f0=.990359,c=.145675):
    f=float(expit(logit(f0)-c*hold))
    p=f*p_hi+(1-f)*p_lo
    sel=np.array([f*(1-p_hi)**2+(1-f)*(1-p_lo)**2,
                  f*p_hi*(1-p_hi)+(1-f)*p_lo*(1-p_lo),
                  f*p_hi*(1-p_hi)+(1-f)*p_lo*(1-p_lo),
                  f*p_hi*p_hi+(1-f)*p_lo*p_lo],dtype=float)
    homog=np.array([(1-p)**2,p*(1-p),p*(1-p),p*p],dtype=float)
    return f,p,sel,homog


def normal_power_sample_size(ps,p0,alpha=.05,power=.8):
    ll=np.log(ps/p0)
    a0=float(p0@ll);a1=float(ps@ll)
    v0=float(p0@((ll-a0)**2));v1=float(ps@((ll-a1)**2))
    n=((norm.ppf(1-alpha)*np.sqrt(v0)+norm.ppf(power)*np.sqrt(v1))/(a1-a0))**2
    return float(n),float(a1),float(v0),float(v1)


def oracle_mc(rng,ps,p0,n,studies=5000):
    ll=np.log(ps/p0)
    dnull=rng.multinomial(n,p0,size=studies)@ll
    dalt=rng.multinomial(n,ps,size=studies)@ll
    threshold=np.quantile(dnull,.95)
    return float((dalt>threshold).mean()),float((dnull>threshold).mean())


def run(seed=20261107+12, studies=1200):
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    rec=[]
    for kinetic in KINETICS:
        for scenario in SCENARIOS:
            for loss in (0.,.10,.20,.30,.40):
                for effect in (0.,.45,.85,1.30):
                    rec.append(simulate(rng,scenario,loss,effect,kinetic,studies=studies))
    with (OUT/'two_arm_loss_grid.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rec[0].keys());writer.writeheader();writer.writerows(rec)
    # repeated cell outcomes at simulator-scale rates, with high-signal sensitivity comparator.
    rates=((.10,.03),(.075,.02),(.15,.04),(.75,.15))
    pdat=[]
    for high,low in rates:
        for hold in (0.,14.,37.):
            f,p,ps,p0=paired_probabilities(hold,high,low)
            n,kl,v0,v1=normal_power_sample_size(ps,p0)
            pdat.append({'p_high':high,'p_low':low,'hold_h':hold,'f_high':f,
                         'marginal_p':p,'within_cell_cov':f*(1-f)*(high-low)**2,
                         'kl_nats':kl,'normal_oracle_n80':n})
    with (OUT/'repeated_contact_information.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=pdat[0].keys());w.writeheader();w.writerows(pdat)
    # full Monte Carlo for primary rates, feasible even for the 0-h n>10 million case.
    pmc=[]
    for hold in (0.,14.,37.):
        f,p,ps,p0=paired_probabilities(hold,.10,.03)
        n80,_,_,_=normal_power_sample_size(ps,p0)
        for factor in (.25,.5,1.,2.,4.):
            n=max(1,int(round(n80*factor)))
            power,fpr=oracle_mc(rng,ps,p0,n,studies=6000)
            pmc.append({'hold_h':hold,'n_cells':n,'factor_of_approx_n80':factor,'power':power,'typeI':fpr})
    with (OUT/'v12_repeated_contact_oracle_mc.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=pmc[0].keys());w.writeheader();w.writerows(pmc)
    # Reviewer-accessible summary with explicit denominators and precision.
    subset=[z for z in rec if z['effect']==.85 and z['kinetic']=='logistic_exponential' and z['nominal_loss']==.4]
    max_null_bounded=max(z['bounded_joint'] for z in rec if z['effect']==0.)
    max_null_raw=max(z['recovered_joint'] for z in rec if z['effect']==0.)
    result={'seed':seed,'studies_per_design':studies,'designs':len(rec),
            'max_null_bounded_joint':max_null_bounded,
            'max_null_recovered_joint':max_null_raw,
            'all_generated_true_effects_inside_observed_missing_bounds':all(z['exact_bounds_cover_generated_full_data'] for z in rec),
            'loss_at_nominal_40pct':subset,
            'repeated_contact_primary':[z for z in pdat if z['p_high']==.10],
            'monotone_recovery_trimming':'not applied to the new both-arm loss designs; cross-world monotonicity was neither constructed nor verified'}
    (OUT/'v12_results_summary.json').write_text(json.dumps(result,indent=2))
    print('ALL RESULTS',json.dumps({'designs':len(rec),'null_max_bounded':max_null_bounded,
                                   'n80_primary':[(x['hold_h'],round(x['normal_oracle_n80'])) for x in pdat if x['p_high']==.10]},indent=2))
    # Main recovery plot. All lines are lower-bound tests. Colors denote loss scenarios, panels kinetic laws.
    fig,axes=plt.subplots(1,2,figsize=(10.4,4.1),sharey=True)
    selected=('pulse_only_uniform','both_uniform_equal','sham_half_uniform','both_phenotype_equal')
    names={'pulse_only_uniform':'Pulse loss only, uniform',
           'both_uniform_equal':'Equal loss, uniform',
           'sham_half_uniform':'Sham loss half of pulse',
           'both_phenotype_equal':'Equal selective loss'}
    for ax,k in zip(axes,KINETICS):
        for scenario in selected:
            rr=sorted((z for z in rec if z['scenario']==scenario and z['kinetic']==k and z['effect']==.85), key=lambda z:z['nominal_loss'])
            ax.plot([100*z['nominal_loss'] for z in rr],[100*z['bounded_joint'] for z in rr],marker='o',linewidth=1.7,label=names[scenario])
        ax.set_xlabel('Nominal loss in pulse arm (%)')
        ax.set_title('Logistic / exponential' if k==KINETICS[0] else 'Cloglog / Weibull')
        ax.set_ylim(0,103);ax.set_xticks([0,10,20,30,40]);ax.grid(alpha=.17)
    axes[0].set_ylabel('Joint attribution supported (%)')
    axes[0].legend(fontsize=8,loc='lower center',bbox_to_anchor=(1.05,1.05),ncol=2,frameon=False)
    fig.tight_layout(rect=(0,0,1,.90));fig.savefig(FIG/'v12_recovery_power.pdf');fig.savefig(FIG/'v12_recovery_power.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10.1,4),sharey=True)
    for ax,k in zip(axes,KINETICS):
        for scenario,caption in [('pulse_only_phenotype','Pulse only'),('sham_half_phenotype','Sham half')]:
            rr=sorted((z for z in rec if z['kinetic']==k and z['scenario']==scenario and z['effect']==0.),key=lambda z:z['nominal_loss'])
            ax.plot([100*z['nominal_loss'] for z in rr],[100*z['recovered_joint'] for z in rr],marker='x',linestyle='--',label=caption+' recovered only')
            ax.plot([100*z['nominal_loss'] for z in rr],[100*z['bounded_joint'] for z in rr],marker='o',linestyle='-',label=caption+' bounded')
        ax.set_xlabel('Nominal loss in pulse arm (%)');ax.set_xticks([0,10,20,30,40]);ax.grid(alpha=.15)
        ax.set_title('Logistic / exponential' if k==KINETICS[0] else 'Cloglog / Weibull')
    axes[0].set_ylabel('False joint attribution (%)')
    axes[0].legend(fontsize=8,loc='lower center',bbox_to_anchor=(1.05,1.05),ncol=2,frameon=False)
    fig.tight_layout(rect=(0,0,1,.89));fig.savefig(FIG/'v12_recovery_null.pdf');fig.savefig(FIG/'v12_recovery_null.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.4,4.4))
    for hold in (0.,14.,37.):
        rr=[z for z in pmc if z['hold_h']==hold]
        ax.semilogx([z['n_cells'] for z in rr],[100*z['power'] for z in rr],marker='o',label=f'{hold:g}-h hold')
    ax.axhline(80,color='gray',linestyle=':',linewidth=1,label='80% target')
    ax.set_ylabel('Oracle detection power (%)');ax.set_xlabel('Cells with two repeated attachment outcomes')
    ax.set_ylim(0,105);ax.legend(frameon=False);ax.grid(alpha=.15)
    fig.tight_layout();fig.savefig(FIG/'v12_repeated_contact_power.pdf');fig.savefig(FIG/'v12_repeated_contact_power.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.3,4.3))
    labels=[]; values=[]
    for h in (0.,14.,37.):
        rr=next(z for z in pdat if z['p_high']==.10 and z['hold_h']==h)
        labels.append(f'{h:g} h');values.append(rr['normal_oracle_n80'])
    ax.bar(labels,values)
    ax.set_yscale('log');ax.set_ylabel('Approximate cells for 80% oracle power (log scale)')
    for i,n in enumerate(values):ax.text(i,n*1.08,f'{n:,.0f}',ha='center',va='bottom',fontsize=8)
    fig.tight_layout();fig.savefig(FIG/'v12_sample_size_by_hold.pdf');fig.savefig(FIG/'v12_sample_size_by_hold.png',dpi=180);plt.close(fig)
    return rec,pdat,pmc,result

if __name__=='__main__':
    run()
