"""Evaluate per-record versus shared FC1-area occupancy calibration.

Reads checked-out derived files from public repo (no raw segmentation in this run).
A common reference is median of FC1 per-record 2-h component-area medians.
CSV acquisition: repository empirical/observables/<record>/{counts.csv,meta.json}.
One FC1 calibration source from a record with early small colonies can bias
absolute cell counts, so this compares 262 px with original record-specific units.
"""
import csv,json,argparse
from pathlib import Path
import numpy as np
from scipy.stats import t as t_dist
from occupancy_estimator import DATASETS, fit_lag

def run(repo,out):
    repo=Path(repo);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    records=[]
    for label,name in DATASETS:
        folder=repo/'empirical'/'observables'/name
        meta=json.loads((folder/'meta.json').read_text())
        rows=list(csv.DictReader((folder/'counts.csv').open()))
        records.append({'label':label,'archive_name':name,'unit':float(meta['median_component_area_first2h']),
                        'times':np.array([float(x['time_h']) for x in rows]),
                        'area':np.array([float(x['area_px']) for x in rows])})
    fc1_refs=[x['unit'] for x in records if x['label'].startswith('FC1')]
    unit=float(np.median(fc1_refs))
    results=[]
    for r in records:
        old,n=fit_lag(r['times'],r['area']/r['unit'])
        shared,n_s=fit_lag(r['times'],r['area']/unit)
        lags=[]
        for smooth in (10,20,40):
            for low,high in ((10,250),(20,300),(40,400)):
                v,_=fit_lag(r['times'],r['area']/unit,smooth,low,high)
                if np.isfinite(v):lags.append(v)
        results.append({'dataset':r['label'],'area_reference_per_record_px':r['unit'],
                        'area_reference_shared_px':unit,'lag_original_h':old,'lag_shared_h':shared,
                        'lag_shared_grid_min_h':min(lags),'lag_shared_grid_max_h':max(lags),
                        'shared_window_valid_n':len(lags)})
    with (out/'shared_calibration.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    by={r['dataset']:r for r in results}
    delta=np.array([by[f'FC2-{j}']['lag_shared_h']-by[f'FC1-{j}']['lag_shared_h'] for j in (1,2,3)])
    mean=float(delta.mean()); ci=t_dist.interval(.95,2,loc=mean,scale=delta.std(ddof=1)/np.sqrt(3))
    nmean=np.mean([by[f'FC1-{j}']['lag_shared_h'] for j in (1,2,3)])
    imean=np.mean([by[f'FC2-{j}']['lag_shared_h'] for j in (1,2,3)])
    summary={'method':'single-cell area calibrated from the median of the three FC1 references',
             'shared_area_reference_px':unit,'FC1_units_px':fc1_refs,
             'paired_differences_h':delta.tolist(),'paired_mean_h':mean,'paired_descriptive_t95ci_h':list(ci),
             'immediate_mean_h':imean,'naive_mean_h':nmean,
             'lag14_h':by['Tr14']['lag_shared_h'],'lag37_h':by['Tr37']['lag_shared_h'],
             'ratio14':(nmean-by['Tr14']['lag_shared_h'])/(nmean-imean),
             'ratio37':(nmean-by['Tr37']['lag_shared_h'])/(nmean-imean),
             'within_record_fit_window_extremes_not_a_biological_confidence_interval':True}
    (out/'shared_calibration_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',default='results')
    a=p.parse_args();run(a.repo,a.out)
