"""Scientific figures for the v8 adversarial analyses."""
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'/'v8'
FIGURES=ROOT/'overleaf'/'figures'
FIGURES.mkdir(parents=True,exist_ok=True)
a=json.loads((RESULTS/'published_lag_model_audit.json').read_text())
b=json.loads((RESULTS/'selection_composition_counterexample.json').read_text())
c=json.loads((RESULTS/'agent_global_prior_sweep.json').read_text())

t=np.linspace(0,40,240)
naive=17.56;d0=12.88
exp_tau=a['linear_exponential']['tau_fitted_to_0_14_h']
f0=b['fraction_of_high_adhesiveness_type_immediately_after_selection'];rate=b['relative_suspension_fitness_cost_per_hour']
selection_lag=lambda x:naive-d0/(f0+(1-f0)*np.exp(rate*x))
kh=a['two_parameter_effective_hill_in_time']['h'];t50=a['two_parameter_effective_hill_in_time']['half_decay_h']
hill_lag=lambda x:naive-d0/(1+(x/t50)**kh)
fig,ax=plt.subplots(figsize=(7.6,4.2))
ax.plot(t,naive-d0*np.exp(-t/exp_tau),label='One-exponential mapping, fitted at 0 and 14 h',linestyle='--')
ax.plot(t,selection_lag(t),label='Selection-composition model, fitted at 0, 14 and 37 h',linestyle='-')
ax.plot(t,hill_lag(t),label='Nonlinear effective hold-time response, interpolated',linestyle=':')
ax.scatter([0,14,37],[4.68,5.46,13.41],color='black',zorder=5,label='Published transfer fit points')
ax.axhline(naive,linestyle=':',color='0.5',alpha=0.7,label='Published naive fit point')
ax.set(xlabel='Suspension interval (h)',ylabel='Colonization lag from published fit (h)',xlim=(0,40),ylim=(0,20))
ax.legend(fontsize=8,loc='upper left',frameon=False)
fig.tight_layout();fig.savefig(FIGURES/'lag_mechanism_counterexamples_v8.png',dpi=220);plt.close(fig)

models=['cell_memory','cell_one_bit','cell_no_inheritance','selection']
labels=['Scalar state','One-bit state','No daughter inheritance','Immutable selection']
fractions=[c['model'][m]['all_three_criteria_fraction'] for m in models]
early=[c['model'][m]['early_fraction'] for m in models]
fig,ax=plt.subplots(figsize=(7,3.5))
y=np.arange(len(models)); ax.barh(y-.16,early,height=.32,label='Immediate and 14-h criteria')
ax.barh(y+.16,fractions,height=.32,label='All criteria, including 37-h washout')
ax.set_yticks(y);ax.set_yticklabels(labels);ax.invert_yaxis();ax.set_xlim(0,1)
ax.set_xlabel('Fraction of 160 synthetic parameter draws')
ax.legend(loc='lower right',fontsize=8,frameon=False)
ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
fig.tight_layout();fig.savefig(FIGURES/'broad_agent_ablation_v8.png',dpi=220);plt.close(fig)

endpoint=list(csv.DictReader((RESULTS/'endpoint_robustness_v8_summary.csv').open()))
models2=['cell_memory','cell_one_bit','cell_no_inheritance','selection']
labels2=['Scalar','One-bit','No inheritance','Selection']
fig,axs=plt.subplots(1,2,figsize=(9.6,3.2),sharex=True)
for m,l in zip(models2,labels2):
    a0=[r for r in endpoint if r['model']==m]
    axs[0].plot([0,14,37],[float(next(r['mean_restricted_first_passage_h'] for r in a0 if r['arm']==arm)) for arm in ['0h','14h','37h']],marker='o',label=l)
    axs[1].plot([0,14,37],[float(next(r['log_linear_onset_h'] for r in a0 if r['arm']==arm)) for arm in ['0h','14h','37h']],marker='o',label=l)
axs[0].set_title('Restricted time to 100 residents')
axs[1].set_title('Log-linear onset, 8–120 residents')
axs[0].set_ylabel('Time (synthetic hours)');axs[1].set_ylabel('Extrapolated time (synthetic hours)')
for ax in axs:ax.set_xlabel('Suspension interval (h)');ax.set_xticks([0,14,37]);ax.grid(alpha=.2)
axs[1].legend(fontsize=8,frameon=False)
fig.tight_layout();fig.savefig(FIGURES/'endpoint_sensitivity_v8.png',dpi=220);plt.close(fig)
print('figures ready')
