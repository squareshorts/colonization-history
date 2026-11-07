"""Transcribe/verify aggregate field summaries from connector-fetched GitHub blobs.

This file includes ONLY aggregate derived values, not individual cell records.
For full reproduction, rerun empirical/field_level_v7.py in the original repo.
The blob SHA values make the inspected source snapshots traceable.
"""
from pathlib import Path
import csv,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent / "empirical"))
from field_level_v7 import summarize

BLOBS={
 "29.5h":{"FC1":"588e391436c0c4564bfb01ebea7e39aeb6b7f5be","FC2":"3e3a6a03e0033281fa6832669a0e84c4d3925d2c"},
 "26.75h":{"FC1":"41233ceeaa21036c8453002c1fd28abdf3048b52","FC2":"b8eab5f06ebf18ad33050ee8b7c58c1e75ade779"},
 "32.5h_14.75h":{"FC1":"6241b5c6172beecc754d3a6c17bd5ea8b5818c26","FC2":"6cbad056e481d420cbbeb5a433ada767261fcf10"}
}
# Each tuple (FC1 attached tracks, FC1 five-minute KM, FC2 tracks, FC2 KM).
# Window 0-3 h, varying the minimum observed attachment time (3 s frames).
P={
 "29.5h": {
   2:(34,.08823529411764706,606,.22442244224422428),
   3:(27,.11111111111111113,539,.25231910946196645),
   4:(24,.125,500,.272),
   6:(21,.14285714285714285,426,.3192488262910796),
   10:(14,.21428571428571427,361,.3767313019390581),
 },
 "26.75h":{
   2:(57,.15789473684210528,275,.1963636363636365),
   3:(50,.18,239,.22594142259414216),
   4:(44,.20454545454545453,216,.25),
   6:(37,.24324324324324323,189,.28571428571428603),
   10:(27,.33333333333333337,157,.34394904458598735),
 },
 "32.5h_14.75h":{
   2:(9,.4444444444444446,110,.31818181818181834),
   3:(9,.4444444444444446,102,.3431372549019608),
   4:(7,.5714285714285715,100,.34999999999999987),
   6:(6,.6666666666666667,95,.36842105263157887),
   10:(5,.8,86,.40697674418604646),
 }}
# Additional time window sensitivity at the original >=2-frame definition.
W={
 "29.5h":{"0-1h":(13,.07692307692307696,145,.35172413793103435),
           "1-3h":(21,.09523809523809525,461,.18438177874186543)},
 "26.75h":{"0-1h":(23,.26086956521739135,68,.25),
           "1-3h":(34,.08823529411764702,207,.178743961352657)},
 "32.5h_14.75h":{"0-1h":(6,.33333333333333337,24,.3750000000000001),
           "1-3h":(3,.6666666666666667,86,.30232558139534876)}
}

def run(outdir=Path("results/field_level_v7")):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for pair in P:
        for n,(na,sa,nb,sb) in P[pair].items():
            rows.append(dict(pair=pair,window="0-3h",threshold_frames=n,
                             n_FC1=na,n_FC2=nb,km300_FC1=sa,km300_FC2=sb,
                             delta_FC2_minus_FC1=sb-sa))
        for window,(na,sa,nb,sb) in W[pair].items():
            rows.append(dict(pair=pair,window=window,threshold_frames=2,
                             n_FC1=na,n_FC2=nb,km300_FC1=sa,km300_FC2=sb,
                             delta_FC2_minus_FC1=sb-sa))
    summary=[]
    for window,n in [("0-3h",n) for n in P["29.5h"]]+[("0-1h",2),("1-3h",2)]:
        subset=[r for r in rows if r["window"]==window and r["threshold_frames"]==n]
        summary.append(dict(window=window,threshold_frames=n,
                            **summarize([r["delta_FC2_minus_FC1"] for r in subset])))
    with (outdir/"archival_field_summary_v7.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (outdir/"archival_field_summary_v7.json").open("w") as f:
        json.dump({"source_repo":"https://github.com/squareshorts/colonization-history",
                   "observables":"derived tracks_full.csv, matching tandem flow-cell records",
                   "blob_shas":BLOBS,"frame_seconds":3,"residence_seconds":300,
                   "exploratory":"All analyses condition on >=k consecutive visible frames; no biological mechanism is identified",
                   "field_results":rows,"summary":summary},f,indent=2)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.5,4.7))
    for pair in P:
        subset=[r for r in rows if r["pair"]==pair and r["window"]=="0-3h"]
        ax.plot([r["threshold_frames"] for r in subset],
                [r["delta_FC2_minus_FC1"] for r in subset],marker="o",label=pair)
    ax.axhline(0,color="0.4",linewidth=1)
    ax.set(xlabel="Minimum successive frames classified as attachment",
           ylabel="FC2 − FC1 five-minute survival probability",
           title="Field-level sensitivity of apparent memory to attachment definition")
    ax.legend(title="Matched experiment")
    fig.tight_layout();fig.savefig(outdir/"archival_field_level_v7.png",dpi=200)
    print(json.dumps(summary,indent=2))
    return rows,summary

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser();p.add_argument("--outdir",type=Path,default=Path("results/field_level_v7"))
    a=p.parse_args();run(a.outdir)
