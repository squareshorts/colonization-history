"""Compare within-record fold and common-absolute surface-occupancy endpoints.

Input: area_px and metadata derived from the public Lee et al. image deposit.
Absolute thresholds refer to the same 262-pixel reference (median FC1
calibration) across flow cells. The endpoint is first sustained attainment,
which can be present at time zero or censored at the end of recording.

Usage from a repository checkout:
  python paper/analysis/endpoint_observation_comparison.py --repo . \
      --records paper/results/endpoint_observation_comparison.csv \
      --summary paper/results/endpoint_observation_summary.json \
      --sensitivity paper/results/absolute_threshold_sensitivity.csv
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from observed_establishment import DATASETS, moving_average, read_record, time_to_sustained_fold


def absolute_sustained_time(times, area_px, reference_px, threshold,
                            *, width=20, persistence_minutes=15):
    """First sustained absolute-occupancy threshold; identify censoring.

    time=0 with status `at_start` means that the threshold was already
    attained in the first observation for the required persistence period;
    it is left-censored as a new crossing time. `after_end` means that
    the threshold was not maintained within the recorded interval.
    """
    t, area = np.asarray(times, float), np.asarray(area_px, float)
    if len(t) < 2 or len(t) != len(area) or not np.all(np.diff(t) > 0):
        raise ValueError("time array must be strictly increasing")
    if not np.isfinite(area).all() or (area < 0).any():
        raise ValueError("segmented area must be finite and nonnegative")
    if reference_px <= 0 or threshold <= 0 or width < 1 or persistence_minutes < 0:
        raise ValueError("invalid threshold parameters")
    sm = moving_average(area / reference_px, width)
    period_h = persistence_minutes / 60.
    for i in np.flatnonzero(sm >= threshold):
        end_time = t[i] + period_h
        if t[-1] < end_time:
            continue
        j = np.searchsorted(t, end_time, side="left")
        if np.all(sm[i:j+1] >= threshold):
            return (float(t[i]), "at_start" if i == 0 else "observed")
    return (float("nan"), "after_end")


def normalized_advantage(t_naive, t_transfer, t_hold):
    """Descriptive advantage fraction based on a specified observation law."""
    if not np.isfinite([t_naive, t_transfer, t_hold]).all() or t_naive <= t_transfer:
        raise ValueError("finite times and positive immediate-transfer advantage required")
    return float((t_naive - t_hold) / (t_naive - t_transfer))


def paired_differences(records, key):
    by = {r["record"]: r for r in records}
    d = []
    for s in (1, 2, 3):
        first = by[f"FC1-{s}"]; second = by[f"FC2-{s}"]
        f, t = first[key], second[key]
        if f is None or t is None:
            # A right-censored FC1 still implies an earlier transferred
            # threshold if FC2 attained it before FC1 observation ended.
            if f is None and t is not None and t < first["record_end_h"]:
                d.append({"series":s, "difference_h":None, "upper_bound_h":
                          float(t - first["record_end_h"]),"status":"earlier_censored"})
            else:
                d.append({"series":s,"difference_h":None,"upper_bound_h":None,
                          "status":"indeterminate"})
        else:
            d.append({"series":s,"difference_h":float(t-f),"upper_bound_h":
                      float(t-f),"status":"observed_pair"})
    return d


def summarize_from_records(records):
    by = {r["record"]: r for r in records}
    out = {}
    for key in ("relative5_h", "absolute20_h", "absolute30_h", "absolute50_h"):
        paired = paired_differences(records,key)
        z={"paired":paired}
        if all(row["difference_h"] is not None for row in paired):
            naive=np.mean([by[f"FC1-{i}"][key] for i in (1,2,3)])
            transfer=np.mean([by[f"FC2-{i}"][key] for i in (1,2,3)])
            d=np.array([row["difference_h"] for row in paired])
            z.update({"naive_mean_h":float(naive),"transfer_mean_h":float(transfer),
                      "paired_mean_h":float(d.mean()),
                      "r14":normalized_advantage(naive,transfer,by["Tr14"][key]),
                      "r37":normalized_advantage(naive,transfer,by["Tr37"][key])})
        else:
            z.update({"naive_mean_h":None,"transfer_mean_h":None,"paired_mean_h":None,
                      "r14":None,"r37":None})
        out[key]=z
    return out


def calculate(repo, thresholds=(20,30,50), width=20, persistence_minutes=15):
    calib=[read_record(repo,name)[2] for id,name in DATASETS if id.startswith('FC1')]
    shared=float(np.median(calib))
    records=[]
    for id,name in DATASETS:
        t,area,local=read_record(repo,name)
        sm=moving_average(area/local,width)
        early=max(float(np.median(sm[t<2.])),1.)
        five,_,target=time_to_sustained_fold(
            t,area/local,fold=5,width=width,persistence_minutes=persistence_minutes,
            minimum_baseline=1.)
        item={"record":id,"source_record":name,"record_end_h":round(float(t[-1]),7),
              "local_cell_area_px":local,"shared_cell_area_px":shared,
              "early_baseline_cell_eq":round(early,6),
              "relative_fivefold_target":round(float(target),6),
              "relative5_h":round(float(five),7) if np.isfinite(five) else None}
        for q in thresholds:
            time,status=absolute_sustained_time(t,area,shared,q,
                        width=width,persistence_minutes=persistence_minutes)
            item[f"absolute{q}_h"]=round(float(time),7) if np.isfinite(time) else None
            item[f"absolute{q}_status"]=status
        records.append(item)
    return records,summarize_from_records(records)


def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--repo',type=Path,required=True)
    a.add_argument('--records',type=Path,required=True)
    a.add_argument('--summary',type=Path,required=True)
    a.add_argument('--sensitivity',type=Path,required=True)
    x=a.parse_args()
    records,summ=calculate(x.repo)
    write_csv(x.records,records)
    x.summary.parent.mkdir(parents=True,exist_ok=True)
    x.summary.write_text(json.dumps(summ,indent=2)+'\n')
    rows=[]
    for threshold in (20,30,50):
        for width in (10,20,40):
            for dwell in (0,15,30):
                rec,_=calculate(x.repo,thresholds=(threshold,),width=width,persistence_minutes=dwell)
                diffs=paired_differences(rec,f'absolute{threshold}_h')
                row={'threshold_cell_eq':threshold,'width_samples':width,
                     'persistence_min':dwell}
                for i,d in enumerate(diffs,1):
                    row[f'series{i}_fc2_minus_fc1_h_or_upper_bound']=round(
                         d['upper_bound_h'],4) if d['upper_bound_h'] is not None else None
                rows.append(row)
    write_csv(x.sensitivity,rows)
    print(json.dumps({'shared_area_px':records[0]['shared_cell_area_px'],
          'comparison':summ},indent=2))

if __name__=='__main__':main()
