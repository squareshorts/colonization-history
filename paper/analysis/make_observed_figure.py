"""Plot reported fivefold crossings, all based on observed occupied area."""
from pathlib import Path
import csv
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--csv',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    with a.csv.open(newline='') as f:
        rows={r['record']:r for r in csv.DictReader(f)}
    fig,ax=plt.subplots(figsize=(6.35,3.9))
    for i in (1,2,3):
        ids=[f'FC1-{i}',f'FC2-{i}']+(['FC3-3'] if i==3 else [])
        xs=list(range(len(ids)))
        ys=[float(rows[x]['time_h']) for x in ids]
        ax.plot(xs,ys,'-o',lw=1.7,ms=6.5,label=f'Series {i}')
    ax.set(xlabel='Consecutive flow cell',ylabel='Time to sustained fivefold occupancy (h)',ylim=(0,27),xlim=(-.18,2.23))
    ax.set_xticks([0,1,2],['FC1','FC2','FC3'])
    ax.grid(axis='y',alpha=.22)
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.legend(frameon=False,fontsize=9,loc='upper right',ncol=1)
    fig.tight_layout()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.out,bbox_inches='tight')
    fig.savefig(a.out.with_suffix('.png'),dpi=200,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__': main()
