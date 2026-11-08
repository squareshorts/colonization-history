"""Recompute occupancy lags with edge-corrected smoothing from the public data.

Input is empirical/observables/*/{counts.csv,meta.json} in a checkout of
squareshorts/colonization-history. The derived data originate from earlier
segmentation of the Lee et al. Figshare binary-image stacks. This script
performs no new segmentation and cannot estimate biological between-run
uncertainty for the single deposited 14-h and 37-h holding runs.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t
from occupancy_estimator import DATASETS,fit_lag


def run(repo,out):
    root=Path(repo)/'empirical'/'observables';dest=Path(out);dest.mkdir(parents=True,exist_ok=True)
    raw=[]
    for label,name in DATASETS:
        d=root/name;meta=json.loads((d/'meta.json').read_text())
        with (d/'counts.csv').open() as stream:r=list(csv.DictReader(stream))
        raw.append((label,name,float(meta['median_component_area_first2h']),
                    np.array([float(z['time_h']) for z in r]),
                    np.array([float(z['area_px']) for z in r])))
    shared=float(np.median([unit for name,_,unit,_,_ in raw if name.startswith('FC1')]))
    rows=[];grid=[]
    for label,name,unit,t,a in raw:
        own,n=fit_lag(t,a/unit)
        common,_=fit_lag(t,a/shared)
        ownalt=[];commonalt=[]
        for smooth in (10,20,40):
            for low,high in ((10,250),(20,300),(40,400)):
                p,pn=fit_lag(t,a/unit,smooth,low,high)
                q,qn=fit_lag(t,a/shared,smooth,low,high)
                grid.append({'dataset':label,'smooth':smooth,'low':low,'high':high,
                             'per_record_lag_h':p,'shared_lag_h':q,'n_per_record_fit':pn,'n_shared_fit':qn})
                if np.isfinite(p):ownalt.append(p)
                if np.isfinite(q):commonalt.append(q)
        rows.append({'dataset':label,'per_record_area_px':unit,'n_occupancy_samples':len(t),
                     'original_lag_h':own,'shared_lag_h':common,'n_standard_fit':n,
                     'original_min_h':min(ownalt),'original_max_h':max(ownalt),
                     'shared_min_h':min(commonalt),'shared_max_h':max(commonalt),
                     'shared_grid_valid_n':len(commonalt)})
    def writecsv(path,records):
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    writecsv(dest/'shared_calibration.csv',rows)
    writecsv(dest/'technical_fit_grid.csv',grid)
    by={r['dataset']:r for r in rows}
    def summary(key):
        d=np.array([by[f'FC2-{j}'][key]-by[f'FC1-{j}'][key] for j in (1,2,3)])
        mean=float(np.mean(d));ci=student_t.interval(.95,2,loc=mean,scale=np.std(d,ddof=1)/np.sqrt(3))
        naive=float(np.mean([by[f'FC1-{j}'][key] for j in (1,2,3)]))
        immediate=float(np.mean([by[f'FC2-{j}'][key] for j in (1,2,3)]))
        lag14=by['Tr14'][key];lag37=by['Tr37'][key]
        return {'paired_FC2_minus_FC1_h':list(map(float,d)),'paired_mean_h':mean,'descriptive_t95ci_h':list(map(float,ci)),
                'naive_mean_h':naive,'immediate_mean_h':immediate,'lag_14h_single_record_h':lag14,
                'lag_37h_single_record_h':lag37,'advantage_ratio_14h':(naive-lag14)/(naive-immediate),
                'advantage_ratio_37h':(naive-lag37)/(naive-immediate),
                'FC3_single_series_3_lag_h':by['FC3-3'][key]}
    p=summary('original_lag_h');s=summary('shared_lag_h')
    s.update({'method':'Boundary-weight-normalized 20-sample moving average; common cell-area reference 262 pixels',
              'shared_area_reference_px':shared,'FC1_reference_px':[r[2] for r in raw if r[0].startswith('FC1')],
              'shared_technical_worst_FC2_minus_FC1_h':[
                   by[f'FC2-{j}']['shared_max_h']-by[f'FC1-{j}']['shared_min_h'] for j in (1,2,3)],
              'per_record_summary':p,'no_independent_hold_replicates':True,
              'FC1_series2_growth_fit_limited_to_roughly_20_to_33_cell_equivalents':True})
    (dest/'shared_calibration_summary.json').write_text(json.dumps(s,indent=2)+'\n')
    (dest/'edge_corrected_summary.json').write_text(json.dumps({
        'corrected_per_record':p,'corrected_shared':s,
        'per_record_worst_technical_differences':[
            by[f'FC2-{j}']['original_max_h']-by[f'FC1-{j}']['original_min_h'] for j in (1,2,3)],
        'FC3_archival_record':'WT YFP-CFP (T1=32.5h,T2=14.75h;FC3)'},indent=2)+'\n')
    return s

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--repo',required=True);a.add_argument('--out',default='results')
    x=a.parse_args();print(json.dumps(run(x.repo,x.out),indent=2))
