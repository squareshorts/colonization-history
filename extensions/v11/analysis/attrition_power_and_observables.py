"""Prospective experimental-design analyses for bacterial history attribution.

All probabilities are illustrative; no assay parameters are fitted to bacteria.
A two-arm randomized clone experiment assigns one founder per endpoint, with no
replacement and one prospectively indexed daughter per founder.  The original
assigned-population treatment effect and the always-observed principal-stratum
contrast are different targets.  The monotone-survival trimming comparison is
included as a conditional sensitivity check, not as identification of the
full assigned-population effect.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.special import expit,logit
from scipy.stats import t,norm
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; FIG=ROOT/'figures'


def kinetic_prob(alpha,z,betaz,eff,rate,phenotype_det,kinetic):
    if kinetic=='logistic_exponential':
        p=expit(alpha+betaz*z+eff)
        survival=np.exp(-rate*np.exp(-.70*eff-phenotype_det*z)*(5/60))
    elif kinetic=='cloglog_weibull':
        p=-np.expm1(-np.exp(alpha+betaz*z+eff))
        survival=np.exp(-(rate*np.exp(-.70*eff-phenotype_det*z)*(5/60))**1.55)
    else: raise ValueError(kinetic)
    return np.clip(p*survival,0,1)


def one_sided_positive(x,alpha=.025):
    se=x.std(axis=1,ddof=1)/np.sqrt(x.shape[1])
    return x.mean(axis=1)-t.ppf(1-alpha,x.shape[1]-1)*se>0


def sim_experiment(rng,kinetic='logistic_exponential',Nstudies=800,B=32,N=48,
                   effect_logodds=0.8,loss=0.1,mode='uniform'):
    shape=(Nstudies,B)
    freq=expit(logit(rng.uniform(.07,.48,(Nstudies,1)))+rng.normal(0,.45,shape))
    alpha=rng.uniform(-3.4,-1.8,(Nstudies,1))+rng.normal(0,.3,shape)
    betaz=rng.uniform(.5,1.6,(Nstudies,1))
    rate=rng.uniform(.6,4,(Nstudies,1))
    phdet=rng.uniform(.1,.9,(Nstudies,1))
    hi1=rng.binomial(N,freq); lo1=N-hi1
    hi0=rng.binomial(N,freq); lo0=N-hi0
    if mode=='uniform': rhi=rlo=1-loss
    elif mode=='phenotype_biased': rhi=1.;rlo=1-loss
    else:raise ValueError(mode)
    hr=rng.binomial(hi1,rhi);lr=rng.binomial(lo1,rlo)
    nr=hr+lr
    def draw(high,low,eff):
        ph=kinetic_prob(alpha,1,betaz,eff,rate,phdet,kinetic)
        pl=kinetic_prob(alpha,0,betaz,eff,rate,phdet,kinetic)
        return rng.binomial(high,ph)+rng.binomial(low,pl)
    # Mothers and indexed daughters; no selection or reproduction after the loss step.
    results={}
    for endp,eff in [('mother',effect_logodds),('daughter',.70*effect_logodds)]:
        s1=draw(hr,lr,eff);s0=draw(hi0,lo0,0)
        raw=s1/np.maximum(nr,1)-s0/N
        # Assigned-population Manski lower bound, sham 100% observed.
        low=s1/N-s0/N
        # Lee bound for the distinct always-recovered principal stratum.
        # Its interpretation relies on potential observation S(pulse)<=S(sham).
        f=np.maximum(nr/N,1e-9)
        lee=s1/np.maximum(nr,1)-np.minimum(1,(s0/N)/f)
        results[endp]={
            'raw':one_sided_positive(raw),'manski':one_sided_positive(low),
            'lee':one_sided_positive(lee),'raw_mean':raw.mean(axis=1),
            'lb_mean':low.mean(axis=1)}
    out={
        'kinetic':kinetic,'effect_logodds':effect_logodds,'mode':mode,
        'requested_loss':loss,'actual_loss_mean':float(np.mean(1-nr/N)),
        'n_studies':Nstudies,'blocks':B,'assigned_per_arm_block':N,
    }
    for method in ['raw','manski','lee']:
        out[method+'_joint']=float(np.mean(results['mother'][method]&results['daughter'][method]))
        out[method+'_mother']=float(np.mean(results['mother'][method]))
        out[method+'_daughter']=float(np.mean(results['daughter'][method]))
    out['true_effect_proxy_complete_recovery']=None
    return out


def paired_probabilities(t_h,p_high=.75,p_low=.15,f0=.990359,c=.1456747):
    f=expit(logit(f0)-c*t_h)
    p=f*p_high+(1-f)*p_low
    sel=np.array([f*(1-p_high)**2+(1-f)*(1-p_low)**2,
                  f*p_high*(1-p_high)+(1-f)*p_low*(1-p_low),
                  f*p_high*(1-p_high)+(1-f)*p_low*(1-p_low),
                  f*p_high**2+(1-f)*p_low**2])
    hom=np.array([(1-p)**2,p*(1-p),p*(1-p),p**2])
    return f,p,sel,hom


def required_power_normal(ps,p0,alpha=.05,power=.80):
    # Simple-vs-simple oracle log-likelihood ratio, known nuisance parameters.
    ll=np.log(ps/p0)
    m0=float(np.sum(p0*ll));m1=float(np.sum(ps*ll))
    v0=float(np.sum(p0*(ll-m0)**2));v1=float(np.sum(ps*(ll-m1)**2))
    n=((norm.ppf(1-alpha)*np.sqrt(v0)+norm.ppf(power)*np.sqrt(v1))/(m1-m0))**2
    return float(n),float(m1),float(-m0),float(v0),float(v1)


def pair_power_mc(rng,sel,hom,n=500,studies=5000):
    ll=np.log(sel/hom)
    # Use Monte Carlo calibrated threshold under the homogeneous model, same n.
    nullcounts=rng.multinomial(n,hom,size=studies)
    altcounts=rng.multinomial(n,sel,size=studies)
    tn=np.quantile(nullcounts@ll,.95)
    return float(np.mean(altcounts@ll>tn)),float(np.mean(nullcounts@ll>tn))


def run(seed=20261107,nstudies=850):
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    records=[]
    for kinetic in ['logistic_exponential','cloglog_weibull']:
        for effect in [0.,.45,.85,1.30]:
            for mode in ['uniform','phenotype_biased']:
                for loss in [0.,.05,.10,.20,.30,.40]:
                    records.append(sim_experiment(rng,kinetic=kinetic,Nstudies=nstudies,
                                   effect_logodds=effect,loss=loss,mode=mode))
    with (OUT/'attrition_power_frontier.csv').open('w',newline='') as file:
        w=csv.DictWriter(file,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    pairs=[]
    for h in [0,14,37]:
        f,p,sel,hom=paired_probabilities(h)
        nreq,d1,d0,v0,v1=required_power_normal(sel,hom)
        for n in [100,250,500,1000,2500,5000,10000]:
            power,fpr=pair_power_mc(rng,sel,hom,n=n,studies=3000)
            pairs.append({'hold_h':h,'n_cells':n,'f_high':float(f),'p_marginal':float(p),
                'pair_covariance':float(f*(1-f)*(.75-.15)**2),
                'one_contact_KL':0.,'two_contact_KL_selection_vs_hom':d1,
                'two_contact_reverse_KL':d0,'approx_n80':nreq,
                'oracle_lrt_power':power,'oracle_lrt_typeI':fpr})
    with (OUT/'paired_contact_information.csv').open('w',newline='') as file:
        w=csv.DictWriter(file,fieldnames=list(pairs[0]));w.writeheader();w.writerows(pairs)
    # Diagnostic figures use matplotlib default colors.
    fig,axes=plt.subplots(1,2,figsize=(11.6,4.35))
    for ax,mode in zip(axes,['uniform','phenotype_biased']):
        for eff in [0.,.45,.85,1.30]:
            rows=[r for r in records if r['kinetic']=='logistic_exponential' and r['mode']==mode and r['effect_logodds']==eff]
            ax.plot([r['actual_loss_mean']*100 for r in rows], [r['manski_joint']*100 for r in rows],marker='o',label=f'effect={eff}')
            if eff==.85:
                ax.plot([r['actual_loss_mean']*100 for r in rows],[r['raw_joint']*100 for r in rows],linestyle='--',marker='x',label='raw, effect=0.85')
                ax.plot([r['actual_loss_mean']*100 for r in rows],[r['lee_joint']*100 for r in rows],linestyle=':',marker='s',label='trim, effect=0.85')
        ax.set_xlabel('Observed founder attrition (%)');ax.set_ylabel('Joint detection frequency (%)')
        ax.set_title('Uniform loss' if mode=='uniform' else 'Adhesion-dependent loss')
        ax.set_ylim(0,104);ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(FIG/'attrition_power_frontier.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.4))
    for h in [0,14,37]:
        rr=[r for r in pairs if r['hold_h']==h]
        ax.semilogx([r['n_cells'] for r in rr],[100*r['oracle_lrt_power'] for r in rr],marker='o',label=f'Hold {h} h')
    ax.axhline(80,linestyle='--',color='gray',linewidth=1)
    ax.set_xlabel('Cells with two conditional-independent contacts')
    ax.set_ylabel('Oracle selection-model detection (%)')
    ax.set_ylim(0,105);ax.legend();fig.tight_layout()
    fig.savefig(FIG/'paired_contact_information.png',dpi=180);plt.close(fig)
    result={'seed':seed,'n_studies_per_cell':nstudies,'attrition_designs':len(records),
            'kinetics':['logistic_exponential','cloglog_weibull'],
            'null_max_bounded_joint_false':max(r['manski_joint'] for r in records if r['effect_logodds']==0),
            'positive_under_loss':[],
            'paired_contact':[]}
    for mode in ['uniform','phenotype_biased']:
        for kinetic in ['logistic_exponential','cloglog_weibull']:
            for e in [.45,.85,1.30]:
                rr=[r for r in records if r['mode']==mode and r['kinetic']==kinetic and r['effect_logodds']==e]
                result['positive_under_loss'].append({'mode':mode,'kinetic':kinetic,'effect':e,
                    'records':[{key:round(float(z[key]),5) for key in ['requested_loss','actual_loss_mean','raw_joint','manski_joint','lee_joint']} for z in rr]})
    for h in [0,14,37]:
        rr=[r for r in pairs if r['hold_h']==h]
        result['paired_contact'].append({'hold_h':h,'high_phenotype_frequency':rr[0]['f_high'],
             'two_contact_KL_nats':rr[0]['two_contact_KL_selection_vs_hom'],
             'approx_n80':rr[0]['approx_n80'],
             'n1000_oracle_power':next(r['oracle_lrt_power'] for r in rr if r['n_cells']==1000)})
    (OUT/'design_results_summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'max_null_bound_false':result['null_max_bounded_joint_false'], 'pair_summary':result['paired_contact'],'sample_effect':result['positive_under_loss'][2]},indent=2))
    return records,pairs,result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--studies',type=int,default=850);p.add_argument('--seed',type=int,default=20261107)
    a=p.parse_args();run(a.seed,a.studies)
