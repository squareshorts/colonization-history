"""Recalculate summaries, field bootstrap CIs, and publication-ready diagnostics."""
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent.parent
R=ROOT/'results'; F=ROOT/'figures'; F.mkdir(exist_ok=True)
d=pd.read_csv(R/'agent_field_observations.csv')
short=['passive','selection','cell_one_bit','cell_no_inheritance','cell_memory','cell_no_decay']
names={'passive':'Stateless','selection':'Immutable selection','cell_one_bit':'One-bit memory','cell_no_inheritance':'Scalar, no inheritance','cell_memory':'Scalar memory','cell_no_decay':'Scalar, no decay','environment':'Environmental state'}
arms=['fresh','transfer_0h','transfer_14h','transfer_37h']
labels=['Fresh','0 h','14 h','37 h']
means=np.array([[d[(d.model==m)&(d.arm==a)].restricted_threshold_h.mean() for a in arms] for m in short])
fig,ax=plt.subplots(figsize=(9.1,5.8))
for i,m in enumerate(short):
    ax.plot(np.arange(4),means[i],marker='o',linewidth=2,label=names[m])
ax.set_xticks(np.arange(4),labels)
ax.set_ylabel('Restricted time to 100 cells (h)')
ax.set_xlabel('Fresh cells or prior-exposure holding time')
ax.set_title('Capability ablation: transferred agents challenged on fresh surfaces')
ax.set_ylim(0,40)
ax.legend(fontsize=8,frameon=False,ncol=2,loc='upper right')
ax.grid(axis='y',alpha=.18)
fig.text(.10,.02,'Synthetic, hand-set parameters. Censored trajectories assigned 40 h. Not the published lag estimator.',fontsize=8)
fig.tight_layout(rect=(0,.04,1,1))
fig.savefig(F/'agent_capability_ablation.png',dpi=210,bbox_inches='tight');plt.close(fig)

# Separate plots from direct measurement to prevent appearance of exact fitting.
fig,ax=plt.subplots(figsize=(7.6,5.2))
for m in short:
    z=d[d.model==m]
    for a in arms:
        assert len(z[z.arm==a])>0
    fresh=z[z.arm=='fresh'].restricted_threshold_h.mean()
    tr=z[z.arm=='transfer_0h'].restricted_threshold_h.mean()
    t37=z[z.arm=='transfer_37h'].restricted_threshold_h.mean()
    ax.scatter(fresh-tr,t37-tr,s=50,label=names[m])
    ax.annotate(names[m],(fresh-tr,t37-tr),xytext=(5,4),textcoords='offset points',fontsize=8)
ax.axhline(4,linestyle='--',linewidth=1)
ax.axvline(4,linestyle='--',linewidth=1)
ax.set_xlabel('Immediate-transfer advantage over fresh (h)')
ax.set_ylabel('37-h washout delay relative to immediate transfer (h)')
ax.set_title('Minimal behavior constraints: benefit and subsequent washout')
ax.grid(alpha=.15)
fig.tight_layout();fig.savefig(F/'agent_capability_tradeoff.png',dpi=210,bbox_inches='tight');plt.close(fig)

osc=pd.read_csv(R/'oscillator_correlations.csv')
fig,ax=plt.subplots(figsize=(8.1,4.8))
for m in ['scalar','oscillator']:
    x=osc[osc.model==m].sort_values('lag_h')
    ax.plot(x.lag_h,x['mean'],linewidth=2,label='One-state OU' if m=='scalar' else 'Coupled two-state OU')
    ax.fill_between(x.lag_h,x.ci_lo,x.ci_hi,alpha=.15)
ax.axhline(0,linewidth=.7)
ax.set_xlabel('Lag (h)')
ax.set_ylabel('Simulated cAMP-proxy / motility-proxy cross-correlation')
ax.set_title('Additional dynamical state supports oscillatory cross-correlations')
ax.grid(alpha=.17);ax.legend(frameon=False)
fig.tight_layout();fig.savefig(F/'oscillator_vs_scalar.png',dpi=210,bbox_inches='tight');plt.close(fig)

rob=pd.read_csv(R/'robustness_grid.csv')
fig,ax=plt.subplots(figsize=(7.9,4.6))
pass_rates=rob.groupby(['tau','mid'])['all_three_qualitative'].mean().unstack('mid')
im=ax.imshow(pass_rates.values,aspect='auto',vmin=0,vmax=1)
ax.set_xticks(range(len(pass_rates.columns)),[str(x) for x in pass_rates.columns])
ax.set_yticks(range(len(pass_rates.index)),[str(x) for x in pass_rates.index])
ax.set_xlabel('Adhesion threshold in memory state')
ax.set_ylabel('Free-state decay constant (h)')
ax.set_title('Fraction satisfying three qualitative signatures across Hill exponents')
for i in range(len(pass_rates.index)):
    for j in range(len(pass_rates.columns)):
        ax.text(j,i,f'{pass_rates.values[i,j]:.2f}',ha='center',va='center',fontsize=9)
fig.colorbar(im,ax=ax,label='Fraction (3 Hill exponents)')
fig.tight_layout();fig.savefig(F/'local_sufficiency_sensitivity.png',dpi=210,bbox_inches='tight');plt.close(fig)

rng=np.random.default_rng(719)
rows=[]
for m in short:
    v={a:d[(d.model==m)&(d.arm==a)].restricted_threshold_h.to_numpy(float) for a in arms}
    for label,a,b in [('immediate_advantage','fresh','transfer_0h'),('washout_37h','transfer_37h','transfer_0h'),('hold14_penalty','transfer_14h','transfer_0h')]:
        obs=v[a].mean()-v[b].mean()
        boot=(rng.choice(v[a],(3000,len(v[a])),replace=True).mean(axis=1) -
              rng.choice(v[b],(3000,len(v[b])),replace=True).mean(axis=1))
        rows.append(dict(model=m,contrast=label,estimate_h=obs,ci_lower_h=np.quantile(boot,.025),ci_upper_h=np.quantile(boot,.975),n_field_per_arm=len(v[a])))
pd.DataFrame(rows).to_csv(R/'field_bootstrap_contrasts.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda x:f'{x:.2f}'))
