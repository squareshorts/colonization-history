"""Fine-grid recovery frontier and conservative versus assumption-assisted inference.

Two inference tiers:
 (i) sharp realized-outcome [0,1] bounds requiring no recovery assumptions;
 (ii) an exploratory binomial-tail calibration, valid only when independently
      established conditional missing-control success probabilities are <= u0
      and missing outcomes are conditionally independent given cell type.
The second method is *not* an assumption-free population-level bound.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.special import expit, logit
from scipy.stats import t, binom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from model_core import simulate, attachment_probability, positive_cluster_t
OUT=HERE/'results'; FIG=HERE/'figures'
KINETICS=('logistic_exponential','cloglog_weibull')


def expected_baselines(seed=119,draws=250000):
    rng=np.random.default_rng(seed)
    f=expit(logit(rng.uniform(.07,.48,draws))+rng.normal(0,.45,draws))
    alpha=rng.uniform(-3.4,-1.8,draws)+rng.normal(0,.3,draws)
    beta=rng.uniform(.5,1.6,draws)
    rate=rng.uniform(.6,4,draws)
    detach=rng.uniform(.1,.9,draws)
    out=[]
    for law in KINETICS:
        def mean(effect):
            high=attachment_probability(alpha,1.,beta,effect,rate,detach,law)
            low=attachment_probability(alpha,0.,beta,effect,rate,detach,law)
            return float(np.mean(f*high+(1-f)*low))
        p0=mean(0);pm=mean(.85);pd=mean(.595)
        out.append({'law':law,'mean_sham_p':p0,'mean_pulse_mother_p':pm,'mean_pulse_daughter_p':pd,
                    'daughter_sham_loss_limit_at_pulse_10pct':(.9*pd-p0)/(1-p0),
                    'mother_sham_loss_limit_at_pulse_10pct':(.9*pm-p0)/(1-p0),
                    'daughter_sham_loss_limit_at_pulse_0pct':(pd-p0)/(1-p0)})
    return out


def fine_grid(seed=20261122,studies=1600):
    rng=np.random.default_rng(seed)
    loss=(0.,.01,.02,.03,.04,.05,.075,.10)
    scenarios=('pulse_only_uniform','both_uniform_equal','sham_half_uniform')
    rows=[]
    for law in KINETICS:
        for scenario in scenarios:
            for epsilon in loss:
                for effect in (0.,.85):
                    rows.append(simulate(rng,scenario,epsilon,effect,law,studies=studies))
    with (OUT/'fine_loss_detection.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    return rows


def calibrated_assay(rng,effect=.85,epsilon=0.1,studies=2800,blocks=32,founders=48,
                     mechanism='uniform',calibration_u=.10):
    """Only here are sham single-cell success probabilities bounded by 0.10.

    Known pL=0.03, pH=0.10; odds update for acquired pulse. Missing sham
    successes stochastically dominated by Binomial(M0,u) if u>=0.10 and
    conditional Bernoulli independence is defensible. Finite-sample alpha
    is split between cluster-mean uncertainty (0.020) and missing-data
    upper tail (0.005); these are illustrative approximate tests.
    """
    f=rng.uniform(.07,.48,(studies,blocks))
    a1h=rng.binomial(founders,f); a1l=founders-a1h
    a0h=rng.binomial(founders,f); a0l=founders-a0h
    if mechanism=='uniform':
        a1hr=rng.binomial(a1h,1-epsilon);a1lr=rng.binomial(a1l,1-epsilon)
        a0hr=rng.binomial(a0h,1-epsilon);a0lr=rng.binomial(a0l,1-epsilon)
    elif mechanism=='select_low':
        a1hr=a1h.copy();a0hr=a0h.copy()
        a1lr=rng.binomial(a1l,1-epsilon);a0lr=rng.binomial(a0l,1-epsilon)
    else: raise ValueError(mechanism)
    n1=a1hr+a1lr;n0=a0hr+a0lr
    m1=founders-n1;m0=founders-n0
    # Simulated counterfactual missing outcomes are generated solely to audit coverage.
    # They never enter any test statistic.
    N=blocks*founders
    M0=m0.sum(axis=1)
    q=np.asarray(binom.ppf(1-.005,M0,calibration_u),float)
    if np.any(q>M0):raise ValueError('invalid upper quantile')
    p0h=.10;p0l=.03
    null=[];out={}
    for label,shift in [('mother',effect),('daughter',.70*effect)]:
        ph=expit(logit(p0h)+shift);pl=expit(logit(p0l)+shift)
        s1=rng.binomial(a1hr,ph)+rng.binomial(a1lr,pl)
        s0=rng.binomial(a0hr,p0h)+rng.binomial(a0lr,p0l)
        xmiss=rng.binomial(a0h-a0hr,p0h)+rng.binomial(a0l-a0lr,p0l)
        frac_miss=(xmiss.sum(axis=1)<=q)
        raw=s1/np.maximum(n1,1)-s0/np.maximum(n0,1)
        base=(s1-s0)/founders
        worst=base - m0/founders
        # Statistical bound on realized sham missing successes, not a cellwise [0,u] bound.
        # Subtract a study-wise upper quantile for the *total* of missing sham outcomes.
        # Cluster SE covers the 32 independent block terms of the observed contrast.
        tcrit=t.ppf(.98,blocks-1)
        se=base.std(axis=1,ddof=1)/np.sqrt(blocks)
        upper_quantile=q/N
        cap_test=(base.mean(axis=1)-tcrit*se-upper_quantile)>0
        worst_test=positive_cluster_t(worst,alpha=.025)
        recovered_test=positive_cluster_t(raw,alpha=.025)
        out[label]={'cap':cap_test,'worst':worst_test,'recovered':recovered_test}
        out[label+'_missing_sham_coverage']=float(frac_miss.mean())
        out[label+'_missing_sham_upper_fraction']=float(np.mean(upper_quantile))
    record={'effect':effect,'epsilon':epsilon,'studies':studies,'blocks':blocks,
            'founders':founders,'loss_mechanism':mechanism,'external_sham_upper_p':calibration_u,
            'assumption_valid':calibration_u>=.10,
            'realized_pulse_loss':float(np.mean(m1/founders)),'realized_sham_loss':float(np.mean(m0/founders)),
            'mean_missing_sham_tail_coverage':float(np.mean([out[e+'_missing_sham_coverage'] for e in ('mother','daughter')])),
            'mean_assisted_missing_sham_upper':float(np.mean([out[e+'_missing_sham_upper_fraction'] for e in ('mother','daughter')]))}
    for tier in ('cap','worst','recovered'):
        record[tier+'_joint']=float(np.mean(out['mother'][tier]&out['daughter'][tier]))
    return record


def make_figures(fine,calibrated):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    names={'pulse_only_uniform':'Pulse-only loss', 'both_uniform_equal':'Equal loss in both arms',
           'sham_half_uniform':'Sham loss at half pulse rate'}
    fig,ax=plt.subplots(1,2,figsize=(9.1,3.5),sharey=True)
    for a,k in zip(ax,KINETICS):
        for scenario in names:
            rr=sorted((r for r in fine if r['kinetic']==k and r['scenario']==scenario and r['effect']==.85),key=lambda r:r['nominal_loss'])
            a.plot([100*r['nominal_loss'] for r in rr],[100*r['bounded_joint'] for r in rr],marker='o',markersize=4,lw=1.9,label=names[scenario])
        a.set(title='Logistic / exponential' if k==KINETICS[0] else 'Cloglog / Weibull',xlabel='Nominal pulse-arm loss (%)',ylim=(-3,103),xticks=[0,2,4,6,8,10]);a.grid(alpha=.17)
    ax[0].set_ylabel('Joint detection (%)')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',ncol=3,bbox_to_anchor=(.5,1.03),frameon=False,fontsize=8.7)
    fig.tight_layout(rect=[0,0,1,.87]);fig.savefig(FIG/'fine_loss_cliff.pdf',bbox_inches='tight');fig.savefig(FIG/'fine_loss_cliff.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.5,3.6))
    for tier,label in [('worst_joint','Assumption-free bound'),('cap_joint','Calibrated sham-probability tail'),('recovered_joint','Recovered-only comparison')]:
        rr=sorted((r for r in calibrated if r['effect']==.85 and r['loss_mechanism']=='uniform' and r['assumption_valid']),key=lambda r:r['epsilon'])
        ax.plot([100*r['epsilon'] for r in rr],[100*r[tier] for r in rr],marker='o',lw=1.8,label=label)
    ax.set(xlabel='Uniform loss in both arms (%)',ylabel='Joint detection (%)',ylim=(-3,103));ax.grid(alpha=.2)
    ax.legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(FIG/'calibrated_missingness.pdf',bbox_inches='tight');fig.savefig(FIG/'calibrated_missingness.png',dpi=180,bbox_inches='tight');plt.close(fig)
    # A log axis should show marker positions, not filled bars with a false zero.
    import csv
    records=list(csv.DictReader((OUT/'repeated_contact_information.csv').open()))
    arr=[r for r in records if abs(float(r['p_high'])-.10)<1e-6]
    fig,ax=plt.subplots(figsize=(6.3,3.35))
    xx=np.array([float(x['hold_h']) for x in arr]);yy=np.array([float(x['normal_oracle_n80']) for x in arr]);
    ax.plot(xx,yy,ls='none',marker='o',markersize=7)
    for x,y in zip(xx,yy):ax.annotate(f'{y:,.0f}',(x,y),xytext=(8,6),textcoords='offset points',fontsize=9)
    ax.set(xlabel='Suspension hold (h)',ylabel='Approximate paired cells for 80% power',yscale='log',xticks=[0,14,37],xlim=(-3,45),ylim=(5000,1e8));ax.grid(axis='y',alpha=.20)
    fig.tight_layout();fig.savefig(FIG/'repeated_oracle_points.pdf',bbox_inches='tight');fig.savefig(FIG/'repeated_oracle_points.png',dpi=180,bbox_inches='tight');plt.close(fig)


def run(seed=20261122,studies=1600):
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    means=expected_baselines()
    saved=OUT/'fine_loss_detection.csv'
    if saved.exists():
        with saved.open() as fh:
            fine=list(csv.DictReader(fh))
        for r in fine:
            for key in ('effect','nominal_loss','bounded_joint','recovered_joint'):
                r[key]=float(r[key])
    else:
        fine=fine_grid(seed,studies)
    rng=np.random.default_rng(seed+1000)
    calibrated=[]
    for mechanism in ('uniform','select_low'):
        for epsilon in (0.,.02,.05,.10,.20):
            for effect in (0.,.85):
                calibrated.append(calibrated_assay(rng,effect=effect,epsilon=epsilon,mechanism=mechanism,calibration_u=.10))
                if epsilon in (.05,.10):
                    calibrated.append(calibrated_assay(rng,effect=effect,epsilon=epsilon,mechanism=mechanism,calibration_u=.03))
    with (OUT/'assisted_calibration.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=calibrated[0].keys());w.writeheader();w.writerows(calibrated)
    make_figures(fine,calibrated)
    summary={'baseline_mean_probabilities':means,'studies_per_fine_design':studies,'fine_designs':len(fine),
             'calibrated_designs':len(calibrated),'calibrated_studies_per_design':2800,
             'scope':'Synthetic illustrative parameters, not biological power predictions',
             'limitation':'An external p_max on Bernoulli success probabilities is a stochastic constraint, not an individual realized-outcome bound; required conditional independence and valid calibration'}
    (OUT/'recovery_design_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for r in fine:
        if r['kinetic']=='logistic_exponential' and r['scenario']=='both_uniform_equal' and r['effect']==.85:
            print('FINE EQ',r['nominal_loss'],r['bounded_joint'],r['recovered_joint'])
    for r in calibrated:
        if r['effect']==.85 and r['loss_mechanism']=='uniform' and r['assumption_valid']:
            print('CAP',r['epsilon'],r['worst_joint'],r['cap_joint'],r['recovered_joint'])
    return summary
if __name__=='__main__':run()
