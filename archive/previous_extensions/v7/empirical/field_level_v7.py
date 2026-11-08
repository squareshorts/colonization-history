"""Matched field-level audit of deposited Lee et al. (2018) derived tracks.

Input: original repository's empirical/observables/*/tracks_full.csv files.
Each matched tandem flow-cell series is one pair, not each cell track.
No inference of cellular versus selection history is attempted from residence alone.

Example from repository root:
    python empirical/field_level_v7.py --outdir results/field_level_v7
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

PAIRS = (
    ("29.5h", "WT YFP-CFP (T=29.5h;FC1)", "WT YFP-CFP (T=29.5h;FC2)"),
    ("26.75h", "WT YFP-CFP (T=26.75h;FC1)", "WT YFP-CFP (T=26.75h;FC2)"),
    ("32.5h_14.75h", "WT YFP-CFP (T1=32.5h,T2=14.75h;FC1)",
     "WT YFP-CFP (T1=32.5h,T2=14.75h;FC2)"),
)
WINDOWS = (("0-3h", 0., 3.), ("0-1h", 0., 1.), ("1-3h", 1., 3.))
THRESHOLDS = (2, 3, 4, 6, 10)


def load_arrivals(path: Path):
    """No thinning or reconstruction: read original deposited derived tracks."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing derived tracks: {path}")
    out = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"start", "end", "origin", "right_censored"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Invalid track table: {path}")
        for r in reader:
            if r["origin"] != "arrival":
                continue
            s, e = int(r["start"]), int(r["end"])
            if e < s:
                raise ValueError(f"Negative residence duration: {path}")
            if r["right_censored"] not in ("True", "False"):
                raise ValueError(f"Invalid censor flag: {path}")
            out.append((s, e-s+1, r["right_censored"] == "True"))
    return out


def km_at(arrivals, threshold=2, start_h=0., stop_h=3., residence_s=300.):
    """KM survival after >= threshold frames, at t seconds from first observation.

    This is a conditional-on-observed-attachment estimate. With the current
    deposited data, all relevant tracks end before 300 s or survive beyond it;
    the generic estimator nevertheless handles censoring at intermediate times.
    """
    x = sorted(((length*3., censored) for first, length, censored in arrivals
                if length >= threshold and start_h <= first*3./3600. < stop_h))
    n = len(x)
    if not n:
        return {"n": 0, "survival": None}
    at_risk, survival = n, 1.
    i = 0
    while i < n and x[i][0] <= residence_s:
        tm = x[i][0]
        events, censored = 0, 0
        while i < n and x[i][0] == tm:
            if x[i][1]:
                censored += 1
            else:
                events += 1
            i += 1
        survival *= 1. - events / at_risk
        at_risk -= events + censored
    return {"n": n, "survival": float(survival)}


def exact_sign_flip(values):
    """Exploratory paired sign-flip test; assumes sign exchangeability."""
    v = np.asarray(values, dtype=float)
    obs = abs(v.mean())
    signed = [abs(np.mean(v * s)) for s in itertools.product((-1, 1), repeat=len(v))]
    return float(np.mean(np.array(signed) >= obs-1e-12))


def summarize(diffs):
    d = np.asarray(diffs, dtype=float)
    n = len(d)
    mean = float(d.mean())
    se = float(d.std(ddof=1) / np.sqrt(n)) if n > 1 else None
    ci = [mean-float(student_t.ppf(.975, n-1))*se,
          mean+float(student_t.ppf(.975, n-1))*se] if n > 1 else [None, None]
    return {"n_pairs": n, "mean": mean, "sd_between_pairs": float(d.std(ddof=1)) if n>1 else None,
            "ci95_t": ci, "exact_signflip_p": exact_sign_flip(d),
            "leave_one_pair_out_means": [float(np.delete(d, j).mean()) for j in range(n)]}


def run(repo_root=Path("."), outdir=Path("results/field_level_v7"), make_figure=True):
    base = Path(repo_root) / "empirical" / "observables"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = {}
    for label, name_a, name_b in PAIRS:
        data[label] = (load_arrivals(base / name_a / "tracks_full.csv"),
                       load_arrivals(base / name_b / "tracks_full.csv"))

    records, summaries = [], []
    for threshold in THRESHOLDS:
        for window, t0, t1 in WINDOWS:
            diffs = []
            for label, _, _ in PAIRS:
                fc1, fc2 = data[label]
                a = km_at(fc1, threshold, t0, t1)
                b = km_at(fc2, threshold, t0, t1)
                if a["survival"] is None or b["survival"] is None:
                    raise ValueError(f"Empty matching set in {label} {window} threshold={threshold}")
                d = b["survival"]-a["survival"]
                diffs.append(d)
                records.append({"pair": label, "window": window, "threshold_frames": threshold,
                                "n_FC1": a["n"], "n_FC2": b["n"],
                                "km300_FC1": a["survival"], "km300_FC2": b["survival"],
                                "delta_FC2_minus_FC1": d})
            summaries.append({"window": window, "threshold_frames": threshold,
                              **summarize(diffs)})

    with (outdir / "field_level_v7.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    with (outdir / "field_level_v7_summary.json").open("w") as f:
        json.dump({"source": "GitHub colonization-history, deposited derived tracks_full.csv",
                   "unit": "matched flow-cell series", "residence_endpoint_s": 300,
                   "frame_duration_s": 3, "field_results": records, "summary": summaries}, f, indent=2)
    if make_figure:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5))
        for label, _, _ in PAIRS:
            x = [r for r in records if r["pair"]==label and r["window"]=="0-3h"]
            ax.plot([r["threshold_frames"] for r in x],
                    [r["delta_FC2_minus_FC1"] for r in x], marker="o", label=label)
        ax.axhline(0, lw=1, color="0.35")
        ax.set(xlabel="Minimum consecutive 3-s frames for an attached track",
               ylabel="Difference in five-minute KM survival (FC2 minus FC1)",
               title="Independent flow-cell pairs: operational-definition sensitivity")
        ax.legend(title="Tandem experiment")
        fig.tight_layout()
        fig.savefig(outdir / "field_threshold_sensitivity.png", dpi=200)
        plt.close(fig)
    return records, summaries


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--repo-root", type=Path, default=Path("."))
    p.add_argument("--outdir", type=Path, default=Path("results/field_level_v7"))
    a=p.parse_args()
    r,s=run(a.repo_root,a.outdir)
    print("fields",len(PAIRS),"rows",len(r))
    print(json.dumps(next(z for z in s if z["window"]=="0-3h" and z["threshold_frames"]==2),indent=2))
