from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import expit,logit
from numpy.linalg import cond

D=Path(__file__).resolve().parents[1]; F=D/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.titlesize':10,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

f0=.990359;c=.1456747; LN=17.56;L0=4.68
t=np.linspace(0,60,600);freq=expit(logit(f0)-c*t); L=LN-(LN-L0)*freq/f0
fig,ax=plt.subplots(figsize=(6.45,3.75))
ax.plot(t,L,lw=2.2,color='#1c6481',label='Both constructed carriers (exactly equal)')
ax.scatter([0,14,37],[4.68,5.46,13.41],marker='o',s=55,edgecolor='#c05a28',facecolor='white',linewidth=1.5,zorder=3,label='Published fitted point estimates')
ax.axhline(LN,ls='--',lw=1,color='#555555',label='Published naive lag')
ax.set(xlabel='Holding time in suspension (h)',ylabel='Colonization lag (h)',xlim=(0,60),ylim=(0,19))
ax.legend(loc='upper left',frameon=False)
fig.tight_layout();fig.savefig(F/'lag_carrier_equivalence.pdf');fig.savefig(F/'lag_carrier_equivalence.png',dpi=180);plt.close(fig)

D0=pd.read_csv(D/'results/archival_occupancy_validation_input.csv')
a=D0[(D0.group=='FC1') & D0.lag_h.notna()].lag_h.to_numpy();b=D0[(D0.group=='FC2') & D0.lag_h.notna()].lag_h.to_numpy()
fig,ax=plt.subplots(figsize=(5.6,3.55))
for i,(x,y) in enumerate(zip(a,b)):
 ax.plot([0,1],[x,y],lw=1.35,color='#56636c',alpha=.8,marker='o',ms=7)
 ax.text(1.035,y,f'Pair {i+1}',fontsize=8,va='center')
ax.set_xticks([0,1],['First flow cell','Second flow cell'])
ax.set(xlim=(-.2,1.4),ylim=(0,22),ylabel='Derived occupancy lag (h)',title='Three matched archival series')
fig.tight_layout();fig.savefig(F/'matched_flowcell_lags.pdf');fig.savefig(F/'matched_flowcell_lags.png',dpi=180);plt.close(fig)

p=pd.read_csv(D/'results/attrition_power_frontier.csv')
fig,ax=plt.subplots(figsize=(6.45,4.0))
for model,sty,color in [('logistic_exponential','-','#1c6481'),('cloglog_weibull','--','#ad642c')]:
 d=p[(p.kinetic==model)&(p.effect_logodds==.85)&(p['mode']=='uniform')]
 ax.plot(d.actual_loss_mean*100,d.manski_joint*100,ls=sty,marker='o',color=color,lw=2,label={'logistic_exponential':'Logistic / exponential','cloglog_weibull':'Cloglog / Weibull'}[model])
ax.axhline(80,ls=':',lw=1,color='#666666',label='80% reference')
ax.set(xlabel='Fraction of initially assigned cells lost (%)',ylabel='Joint mother--daughter detection (%)',ylim=(0,106),xlim=(-1,42))
ax.legend(frameon=False,loc='lower left');fig.tight_layout()
fig.savefig(F/'recovery_power_tradeoff.pdf');fig.savefig(F/'recovery_power_tradeoff.png',dpi=180);plt.close(fig)

q=pd.read_csv(D/'results/paired_contact_information.csv')
fig,ax=plt.subplots(figsize=(6.45,4.0))
for hold,color in [(0,'#666666'),(14,'#1c6481'),(37,'#ad642c')]:
 d=q[q.hold_h==hold]
 ax.plot(d.n_cells,d.oracle_lrt_power*100,marker='o',lw=1.8,color=color,label=f'Hold {hold} h')
ax.set_xscale('log');ax.axhline(80,ls=':',color='#666666',lw=1)
ax.set(xlabel='Tracked cells measured twice',ylabel='Power with known model parameters (%)',ylim=(0,106))
ax.legend(frameon=False,loc='lower right');fig.tight_layout()
fig.savefig(F/'repeated_contact_discrimination.pdf');fig.savefig(F/'repeated_contact_discrimination.png',dpi=180);plt.close(fig)

fig,ax=plt.subplots(figsize=(6.45,3.7))
k=np.linspace(0,.975,400)
def matrix(k):return np.array([[1,1,1,0,0],[0,0,0,1,1],[1,0,1,0,0],[0,0,0,1,0],[1,0,k,0,0]])
ax.semilogy(k,[cond(matrix(float(x))) for x in k],lw=2,color='#1c6481')
ax.scatter([0,.9,.98],[4.0489,48.58,250.45],color='#ad642c',s=20,zorder=2)
ax.set(xlabel='Exterior sensitivity of localization readout, $\\kappa$',ylabel='Condition number (log scale)',xlim=(0,1),ylim=(3,1000))
fig.tight_layout();fig.savefig(F/'localization_selectivity.pdf');fig.savefig(F/'localization_selectivity.png',dpi=180);plt.close(fig)
print('saved',list(p.name for p in F.glob('*.pdf')))
