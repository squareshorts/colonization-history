"""Plot the corrected matched flow-cell results; FC3 is one dependent sequence."""
import argparse,csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def run(csvfile,out):
    with open(csvfile,newline='') as f:r={z['dataset']:z for z in csv.DictReader(f)}
    fig,ax=plt.subplots(figsize=(6.2,4.0))
    for j in (1,2,3):
        ids=[f'FC1-{j}',f'FC2-{j}']+(['FC3-3'] if j==3 else [])
        ax.plot(list(range(len(ids))),[float(r[id]['original_lag_h']) for id in ids],'-o',lw=1.9,ms=6,label=f'Pair {j}')
    ax.set_xticks([0,1,2],['FC1','FC2','FC3']);ax.set_xlim(-.15,2.25);ax.set_ylim(0,22.5)
    ax.set_ylabel('Area-based occupancy lag (h)')
    ax.legend(loc='upper right',frameon=False,fontsize=9)
    ax.grid(axis='y',alpha=.18);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--summary',default='results/shared_calibration.csv');p.add_argument('--output',default='figures/matched_flowcell_lags.pdf')
    x=p.parse_args();Path(x.output).parent.mkdir(parents=True,exist_ok=True);run(x.summary,x.output)
