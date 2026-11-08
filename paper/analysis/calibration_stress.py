"""Reassess calibration-tail false attribution under asymmetric control loss.

Population is an independent two-phenotype Bernoulli model, with pH=0.10,
pL=0.03. Each assigned founder stays in the denominator, although response
may be missing. Tests approximate a cluster Student-t rule, as in the earlier
analysis, and are not advertised as universal finite-sample error bounds.
Valid calibration u=0.10, invalid u=0.03 share exactly the same generated data.
"""
from pathlib import Path
import csv, json
import numpy as np
from scipy.special import expit,logit
from scipy.stats import t, binom

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'

def simulate_single(seed=20261124, studies=6000, sham_loss=.1, pulse_loss=0.,
                    effect=0., founders=48, blocks=32, invalid=.03, valid=.10,
                    phenotype_low=.07,phenotype_high=.48):
    rng=np.random.default_rng(seed)
    f=rng.uniform(phenotype_low,phenotype_high,(studies,blocks))
    high1=rng.binomial(founders,f);low1=founders-high1
    high0=rng.binomial(founders,f);low0=founders-high0
    # Independently randomized loss in each treatment arm.
    hr1=rng.binomial(high1,1-pulse_loss);lr1=rng.binomial(low1,1-pulse_loss)
    hr0=rng.binomial(high0,1-sham_loss);lr0=rng.binomial(low0,1-sham_loss)
    n1=hr1+lr1;n0=hr0+lr0
    m0=founders-n0;m1=founders-n1
    m0_total=m0.sum(axis=1)
    alpha=.005; tcrit=t.ppf(.98,blocks-1)
    tails={u:np.asarray(binom.ppf(1-alpha,m0_total,u)) for u in (valid,invalid)}
    results={u:[] for u in (valid,invalid)}; worst=[];recovered=[]; coverage={u:[] for u in (valid,invalid)}
    for label,shift in [('mother',effect),('daughter',.7*effect)]:
        p1h=expit(logit(.10)+shift);p1l=expit(logit(.03)+shift)
        x1=rng.binomial(hr1,p1h)+rng.binomial(lr1,p1l)
        x0=rng.binomial(hr0,.10)+rng.binomial(lr0,.03)
        missing_successes=rng.binomial(high0-hr0,.10)+rng.binomial(low0-lr0,.03)
        contrast=(x1-x0)/founders
        se=contrast.std(axis=1,ddof=1)/np.sqrt(blocks)
        score=contrast.mean(axis=1)-tcrit*se
        for u in (valid,invalid):
            results[u].append(score-tails[u]/(blocks*founders)>0)
            coverage[u].append(missing_successes.sum(axis=1)<=tails[u])
        worst_scores=contrast-m0/founders
        se_w=worst_scores.std(axis=1,ddof=1)/np.sqrt(blocks)
        worst.append(worst_scores.mean(axis=1)-t.ppf(.975,blocks-1)*se_w>0)
        r=x1/np.maximum(n1,1)-x0/np.maximum(n0,1)
        se_r=r.std(axis=1,ddof=1)/np.sqrt(blocks)
        recovered.append(r.mean(axis=1)-t.ppf(.975,blocks-1)*se_r>0)
    row={'seed':seed,'studies':studies,'blocks':blocks,'founders_per_arm_block':founders,
         'effect':effect,'sham_loss':sham_loss,'pulse_loss':pulse_loss,
         'phenotype_low':phenotype_low,'phenotype_high':phenotype_high,
         'mean_realized_sham_loss':float(np.mean(m0/founders)),
         'mean_realized_pulse_loss':float(np.mean(m1/founders)),
         'null_recovered_only_joint':float(np.mean(recovered[0]&recovered[1])),
         'null_worst_case_joint':float(np.mean(worst[0]&worst[1]))}
    for u in (valid,invalid):
        tag='valid_u010' if u==valid else 'invalid_u003'
        row[tag+'_joint']=float(np.mean(results[u][0]&results[u][1]))
        row[tag+'_joint_count']=int(np.sum(results[u][0]&results[u][1]))
        row[tag+'_coverage_mother']=float(np.mean(coverage[u][0]))
        row[tag+'_coverage_daughter']=float(np.mean(coverage[u][1]))
    return row

def run():
    data=[]
    for sham_loss,pulse_loss in [(0.,0.),(.05,.05),(.10,.10),(.10,0.),(.20,0.),(.30,0.),(.40,0.)]:
        data.append(simulate_single(seed=20261124+int(sham_loss*100)*100+int(pulse_loss*100),
                                    sham_loss=sham_loss,pulse_loss=pulse_loss))
    for sham_loss,pulse_loss in [(.10,.10),(.20,0.),(.30,0.),(.40,0.)]:
        data.append(simulate_single(seed=20281124+int(sham_loss*100)*100+int(pulse_loss*100),
                                    sham_loss=sham_loss,pulse_loss=pulse_loss,phenotype_low=.80,phenotype_high=.99))
    OUT.mkdir(exist_ok=True,parents=True)
    fn=OUT/'null_asymmetric_calibration_stress.csv'
    with fn.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
    for d in data:
        print(f"enrichment={d['phenotype_low']:.0%}-{d['phenotype_high']:.0%} sham={d['sham_loss']:.0%} pulse={d['pulse_loss']:.0%} invalid={d['invalid_u003_joint']:.3%} valid={d['valid_u010_joint']:.3%} cov_invalid={d['invalid_u003_coverage_mother']:.2%} recovered={d['null_recovered_only_joint']:.2%}")
    return data
if __name__=='__main__':run()
