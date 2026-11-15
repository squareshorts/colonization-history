"""Occupancy progression between common thresholds, and the decay-mixture bound by readout.

Transferred flow cells begin with roughly ten times the early occupancy of
first-exposure flow cells, so the time to reach a common threshold combines
initial cell delivery with later surface expansion. The time spent between two
common thresholds removes the delivery component: both populations are timed
from the moment they first sustain the same occupied area.

The early occupancy ratio compares occupied area at 4-8 h with the first 2 h;
it measures persistence of surface occupancy and does not follow individual
cells. The expansion rate is the net growth of occupied area, which combines
division, arrival, departure and changes in projected area.

The script also evaluates the positive-decay-mixture constraint
r(37) >= r(14)^(37/14) (electronic supplementary material, section 1) under
every observation rule reported in the article.

Usage from a repository checkout:
  python paper/analysis/matched_progression.py --repo . --out paper/results
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from endpoint_observation_comparison import absolute_sustained_time
from observed_establishment import DATASETS, moving_average, read_record

PUBLISHED_LAGS_H = {"naive": 17.56, "transfer": 4.68, "hold14": 5.46, "hold37": 13.41}
LEVELS = (20, 30, 40, 50)
INTERVALS = ((20, 30), (30, 50), (20, 50), (30, 40), (40, 50))


def shared_reference(repo):
    return float(np.median([read_record(repo, name)[2]
                            for rid, name in DATASETS if rid.startswith("FC1")]))


def crossing_table(repo, width=20, persistence_minutes=15, levels=LEVELS):
    """First sustained attainment of each common threshold, with censoring status."""
    ref = shared_reference(repo)
    out = []
    for rid, name in DATASETS:
        t, area, _ = read_record(repo, name)
        row = {"record": rid, "record_end_h": float(t[-1])}
        for q in levels:
            time, status = absolute_sustained_time(t, area, ref, q, width=width,
                                                   persistence_minutes=persistence_minutes)
            row[f"t{q}_h"] = time if np.isfinite(time) else None
            row[f"t{q}_status"] = status
        out.append(row)
    return out, ref


def progression(row, lo, hi):
    """Time from first sustained `lo` to first sustained `hi` occupancy.

    Returns (value, kind): `observed` for two observed crossings; `lower_bound`
    when `hi` was never sustained (value is the observed remainder of the
    record); `upper_bound` when `lo` was already sustained at the first frame
    (the true crossing of `lo` is at or before time zero, so the interval is at
    least the reported value only if `lo` was crossed exactly at zero; it is
    reported as unavailable instead).
    """
    a, sa = row[f"t{lo}_h"], row[f"t{lo}_status"]
    b, sb = row[f"t{hi}_h"], row[f"t{hi}_status"]
    if sa == "at_start":
        return None, "lower_threshold_at_start"
    if a is None:
        return None, "lower_threshold_not_reached"
    if sb == "after_end":
        return row["record_end_h"] - a, "lower_bound"
    return b - a, "observed"


def progression_table(rows, intervals=INTERVALS):
    out = []
    for row in rows:
        for lo, hi in intervals:
            v, kind = progression(row, lo, hi)
            rate = (np.log(hi / lo) / v) if (v is not None and kind == "observed" and v > 0) else None
            out.append({"record": row["record"], "from_cell_eq": lo, "to_cell_eq": hi,
                        "interval_h": None if v is None else round(float(v), 6),
                        "status": kind,
                        "expansion_rate_per_h": None if rate is None else round(float(rate), 6)})
    return out


def paired_progression(prog, lo, hi):
    by = {(r["record"], r["from_cell_eq"], r["to_cell_eq"]): r for r in prog}
    pairs = []
    for s in (1, 2, 3):
        f = by[(f"FC1-{s}", lo, hi)]
        g = by[(f"FC2-{s}", lo, hi)]
        if f["status"] == "observed" and g["status"] == "observed":
            pairs.append({"series": s, "fc1_h": f["interval_h"], "fc2_h": g["interval_h"],
                          "fc2_minus_fc1_h": round(g["interval_h"] - f["interval_h"], 6),
                          "status": "observed_pair"})
        elif f["status"] == "lower_bound" and g["status"] == "observed":
            # FC1 interval is at least f; the difference is at most g - f.
            pairs.append({"series": s, "fc1_h": f["interval_h"], "fc2_h": g["interval_h"],
                          "fc2_minus_fc1_h": None,
                          "fc2_minus_fc1_upper_bound_h": round(g["interval_h"] - f["interval_h"], 6),
                          "status": "fc1_right_censored"})
        else:
            pairs.append({"series": s, "fc1_h": f["interval_h"], "fc2_h": g["interval_h"],
                          "fc2_minus_fc1_h": None, "status": f"{f['status']}|{g['status']}"})
    return pairs


def ratio(naive, transfer, hold):
    return (naive - hold) / (naive - transfer)


def decay_bound_rows(repo, results_dir):
    """r14, r37 and the positive-decay-mixture lower bound for every readout."""
    results_dir = Path(results_dir)
    summ = json.loads((results_dir / "endpoint_observation_summary.json").read_text())
    lag = json.loads((results_dir / "shared_calibration_summary.json").read_text())
    per = lag["per_record_summary"]
    p = PUBLISHED_LAGS_H
    readouts = [
        ("Published fitted colonization lag",
         ratio(p["naive"], p["transfer"], p["hold14"]), ratio(p["naive"], p["transfer"], p["hold37"])),
        ("Image-derived occupancy lag (record-specific area)",
         per["advantage_ratio_14h"], per["advantage_ratio_37h"]),
        ("Image-derived occupancy lag (common area)",
         lag["advantage_ratio_14h"], lag["advantage_ratio_37h"]),
        ("Recording-specific sustained fivefold occupancy",
         summ["relative5_h"]["r14"], summ["relative5_h"]["r37"]),
        ("Common sustained 20-equivalent occupancy",
         summ["absolute20_h"]["r14"], summ["absolute20_h"]["r37"]),
        ("Common sustained 30-equivalent occupancy",
         summ["absolute30_h"]["r14"], summ["absolute30_h"]["r37"]),
    ]
    rows = []
    for name, r14, r37 in readouts:
        if r14 > 1:
            bound, verdict = None, "violated: advantage increases between 0 and 14 h"
        else:
            bound = r14 ** (37 / 14)
            verdict = "violated" if r37 < bound else "satisfied"
        rows.append({"readout": name, "r14": round(float(r14), 4), "r37": round(float(r37), 4),
                     "mixture_lower_bound_r37": None if bound is None else round(float(bound), 4),
                     "verdict": verdict})
    return rows


def early_occupancy_ratio(t, sm, early=(0., 2.), window=(4., 8.)):
    """Median smoothed occupancy in `window` divided by its median in `early`."""
    e = sm[(t >= early[0]) & (t < early[1])]
    lo, hi = window[0], min(window[1], float(t[-1]))
    w = sm[(t >= lo) & (t < hi)]
    if len(e) == 0 or len(w) == 0 or np.median(e) <= 0:
        return None
    return float(np.median(w) / np.median(e))


def max_expansion_rate(t, sm, span_h=2., floor=5., start_h=2.):
    """Largest least-squares slope of log occupancy over any `span_h` window.

    Only windows starting after `start_h` and with every smoothed value at or
    above `floor` cell equivalents are used, so near-empty fields with
    single-cell noise do not contribute.
    """
    best, at = None, None
    dt = float(np.median(np.diff(t)))
    n = int(round(span_h / dt))
    logy = np.log(np.maximum(sm, 1e-9))
    for i in range(0, len(t) - n, max(1, n // 12)):
        if t[i] < start_h:
            continue
        seg = slice(i, i + n + 1)
        if np.min(sm[seg]) < floor:
            continue
        slope = float(np.polyfit(t[seg], logy[seg], 1)[0])
        if best is None or slope > best:
            best, at = slope, float(t[i])
    return best, at


def dynamics_table(repo, width=20, windows=((4., 8.), (3., 6.), (6., 10.))):
    """Early occupancy ratio and maximal net occupancy expansion rate per record."""
    ref = shared_reference(repo)
    out = []
    for rid, name in DATASETS:
        t, area, _ = read_record(repo, name)
        sm = moving_average(area / ref, width)
        row = {"record": rid, "record_end_h": round(float(t[-1]), 4),
               "early_median_cell_eq": round(float(np.median(sm[t < 2.])), 4)}
        for w in windows:
            r = early_occupancy_ratio(t, sm, window=w)
            row[f"occupancy_ratio_{w[0]:g}_{w[1]:g}h"] = None if r is None else round(r, 4)
        for span, floor in ((2., 5.), (3., 5.), (2., 10.)):
            rate, at = max_expansion_rate(t, sm, span_h=span, floor=floor)
            tag = f"{span:g}h_floor{floor:g}"
            row[f"max_rate_{tag}_per_h"] = None if rate is None else round(rate, 4)
            row[f"max_rate_{tag}_start_h"] = None if at is None else round(at, 3)
        out.append(row)
    return out


def trajectories(repo, width=20):
    """Smoothed occupancy on the common area scale, subsampled every 5 minutes."""
    ref = shared_reference(repo)
    out = []
    for rid, name in DATASETS:
        t, area, _ = read_record(repo, name)
        sm = moving_average(area / ref, width)
        keep = np.arange(0, len(t), 10)
        for i in keep:
            out.append({"record": rid, "time_h": round(float(t[i]), 5),
                        "cell_eq": round(float(sm[i]), 4)})
    return out


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    rows, ref = crossing_table(a.repo)
    prog = progression_table(rows)
    summary = {"shared_area_reference_px": ref,
               "paired": {f"{lo}_to_{hi}": paired_progression(prog, lo, hi) for lo, hi in INTERVALS}}
    sens = []
    for width in (10, 20, 40):
        for dwell in (0, 15, 30):
            r2, _ = crossing_table(a.repo, width=width, persistence_minutes=dwell)
            p2 = progression_table(r2, intervals=((20, 30), (30, 50), (20, 50)))
            for lo, hi in ((20, 30), (30, 50), (20, 50)):
                for pr in paired_progression(p2, lo, hi):
                    sens.append({"width_samples": width, "persistence_min": dwell,
                                 "from_cell_eq": lo, "to_cell_eq": hi, **pr})
    write_csv(a.out / "matched_progression_crossings.csv", rows)
    write_csv(a.out / "matched_progression.csv", prog)
    write_csv(a.out / "matched_progression_sensitivity.csv", sens)
    (a.out / "matched_progression_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_csv(a.out / "decay_bound_by_readout.csv", decay_bound_rows(a.repo, a.out))
    write_csv(a.out / "occupancy_trajectories.csv", trajectories(a.repo))
    write_csv(a.out / "occupancy_dynamics.csv", dynamics_table(a.repo))
    print(json.dumps(summary["paired"], indent=1))


if __name__ == "__main__":
    main()
