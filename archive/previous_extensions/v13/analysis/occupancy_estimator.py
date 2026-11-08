"""Recompute occupied-area lag sensitivity from a checked-out public repository.

Source: sq uareshorts/colonization-history, empirical/observables/*/{counts.csv,meta.json}
The occupancy-fit estimator follows the repository's occupancy.py, not the
original Lee et al. published group-level exponential-fit procedure.

Example:
 python recompute_occupancy_lags.py --repo /path/to/colonization-history \
        --out /path/to/technical_fit_sensitivity.csv
"""
import argparse,csv,json
from pathlib import Path
import numpy as np

DATASETS=[
 ('FC1-1','WT YFP-CFP (T=29.5h;FC1)'),
 ('FC2-1','WT YFP-CFP (T=29.5h;FC2)'),
 ('FC1-2','WT YFP-CFP (T=26.75h;FC1)'),
 ('FC2-2','WT YFP-CFP (T=26.75h;FC2)'),
 ('FC1-3','WT YFP-CFP (T1=32.5h,T2=14.75h;FC1)'),
 ('FC2-3','WT YFP-CFP (T1=32.5h,T2=14.75h;FC2)'),
 ('Tr14','WT YFP-CFP (Tr=14h)'),('Tr37','WT YFP-CFP (Tr=37h)')]


def fit_lag(times, counts, smooth=20, lo=20, hi=300):
    cs=np.convolve(counts,np.ones(smooth)/smooth,mode='same')
    base=max(float(np.median(cs[:240])),1.0)
    sel=(cs>=lo)&(cs<=hi)
    if sel.sum()<20 or times[-1]<3:
        return float('nan'),int(sel.sum())
    sl,intercept=np.polyfit(times[sel],np.log(cs[sel]),1)
    lag=(np.log(base)-intercept)/sl if sl>0 else float('nan')
    return float(lag),int(sel.sum())


def run(repo,out):
    obs=Path(repo)/'empirical'/'observables';results=[]
    for short,name in DATASETS:
        d=obs/name
        with (d/'counts.csv').open() as f:rows=list(csv.DictReader(f))
        meta=json.loads((d/'meta.json').read_text())
        t=np.array([float(x['time_h']) for x in rows])
        y=np.array([float(x['area_px']) for x in rows])/float(meta['median_component_area_first2h'])
        baseline,_=fit_lag(t,y)
        valid=[]
        for smooth in [10,20,40]:
            for lo,hi in [(10,250),(20,300),(40,400)]:
                lag,n=fit_lag(t,y,smooth,lo,hi)
                results.append({'dataset':short,'archive_name':name,'smooth':smooth,'lower':lo,'upper':hi,'n_fit':n,'lag_h':lag})
                if np.isfinite(lag):valid.append(lag)
        print(short,'baseline',round(baseline,3),'nfits',len(valid),'range',(round(min(valid),3),round(max(valid),3)))
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,default=Path('technical_fit_sensitivity.csv'))
    a=p.parse_args();run(a.repo,a.out)
