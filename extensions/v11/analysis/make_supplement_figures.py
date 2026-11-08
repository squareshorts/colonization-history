from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path(__file__).resolve().parents[1]
P=pd.read_csv(D/'results/attrition_power_frontier.csv')
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(6.45,3.7))
for k,c in [('logistic_exponential','#1c6481'),('cloglog_weibull','#ad642c')]:
 rows=P[(P['kinetic']==k)&(P['mode']=='phenotype_biased')&(P.effect_logodds==0)]
 ax.plot(rows.actual_loss_mean*100,rows.raw_joint*100,marker='o',label=k.replace('_',' / '),lw=1.8,color=c)
 ax.plot(rows.actual_loss_mean*100,rows.manski_joint*100,linestyle='--',linewidth=1.2,color=c)
ax.set(xlabel='Actual founder attrition (%)',ylabel='False joint mother--daughter declarations (%)',ylim=(-.5,7))
ax.text(1,6.6,'Solid: recovered-only; dashed: attrition-bounded',size=8,va='top')
ax.legend(loc='center left',frameon=False);fig.tight_layout()
fig.savefig(D/'figures/recovery_selection_bias.pdf');fig.savefig(D/'figures/recovery_selection_bias.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(6.45,3.7))
for eff,c in [(.45,'#ad642c'),(.85,'#1c6481'),(1.3,'#47775c')]:
 rows=P[(P['kinetic']=='logistic_exponential')&(P['mode']=='uniform')&(P.effect_logodds==eff)]
 ax.plot(rows.actual_loss_mean*100,rows.manski_joint*100,marker='o',lw=1.7,color=c,label=f'Predictor shift {eff}')
ax.axhline(80,ls=':',lw=1,color='gray');ax.set(xlabel='Founder attrition (%)',ylabel='Bounded joint detection (%)',ylim=(0,105))
ax.legend(frameon=False,loc='lower left');fig.tight_layout()
fig.savefig(D/'figures/recovery_effect_strength.pdf');fig.savefig(D/'figures/recovery_effect_strength.png',dpi=180);plt.close(fig)
