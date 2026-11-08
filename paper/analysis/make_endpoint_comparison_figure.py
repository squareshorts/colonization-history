"""Plot the retrospective relative and common-absolute flow-cell endpoints."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--csv',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    with a.csv.open(newline='') as handle:
        rows={r['record']:r for r in csv.DictReader(handle)}
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.65),sharey=True)
    titles=['Within-record fivefold increase','Shared 30-cell-area threshold']
    keys=['relative5_h','absolute30_h']
    line_handles=[]
    for ax,title,key in zip(axs,titles,keys):
        for s in (1,2,3):
            rr=[rows[f'FC1-{s}'],rows[f'FC2-{s}']]
            yy=[float(z[key]) for z in rr]
            h,=ax.plot([0,1],yy,marker='o',ms=6,lw=1.8,label=f'Series {s}')
            if len(line_handles)<3:line_handles.append(h)
        ax.set_title(title,fontsize=10.5,pad=10)
        ax.set_xticks([0,1],['First exposure','After transfer'])
        ax.set_xlim(-.13,1.13)
        ax.set_ylim(-.6,29.2)
        ax.grid(axis='y',alpha=.19)
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(axis='x',labelsize=9)
    axs[0].set_ylabel('Time to sustained occupancy threshold (h)')
    fig.legend(line_handles,[h.get_label() for h in line_handles],
               loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.52,1.02),fontsize=9)
    fig.tight_layout(rect=(0,0,1,.89),w_pad=1.9)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.out,bbox_inches='tight')
    fig.savefig(a.out.with_suffix('.png'),dpi=220,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':main()
