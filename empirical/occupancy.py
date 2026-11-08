"""Occupancy lag times per dataset, for validation against Lee et al. (2018).

Attached-cell equivalents are foreground area divided by the median component
area in the first two hours, smoothed over 10 min. An exponential is fitted to
the phase with 20 to 300 cell equivalents; the lag is the time at which that
exponential crosses the median occupancy of the first two hours.
"""
import csv, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as A

PUBLISHED_LAG_H = {"FC1": 17.56, "FC2": 4.68, "Tr14": 5.46, "Tr37": 13.41}


def lag(dataset):
    t, c, m = A.occupancy(dataset)
    k = 20
    cs = np.convolve(c, np.ones(k), mode="same") / np.convolve(np.ones_like(c), np.ones(k), mode="same")
    base = max(float(np.median(cs[:240])), 1.0)
    sel = (cs >= 20) & (cs <= 300)
    if sel.sum() < 20 or t[-1] < 3:
        return base, np.nan, np.nan, float(t[-1])
    sl, ic = np.polyfit(t[sel], np.log(cs[sel]), 1)
    return base, float((np.log(base) - ic) / sl), float(np.log(2) / sl), float(t[-1])


def run():
    rows = []
    groups = dict(A.GROUPS)
    groups.update({k: v for k, v in A.MUTANTS.items()})
    for g, ds in groups.items():
        for n in ds:
            if not (A.OBS / n / "counts.csv").exists():
                continue
            base, tl, dbl, end = lag(n)
            rows.append({"group": g, "dataset": n, "initial_cells": round(base, 2),
                         "lag_h": round(tl, 2) if np.isfinite(tl) else "",
                         "doubling_h": round(dbl, 2) if np.isfinite(dbl) else "",
                         "record_h": round(end, 2), "arrivals_per_h_first3h": round(A.arrival_rate(n), 1),
                         "published_group_lag_h": PUBLISHED_LAG_H.get(g, "")})
    out = A.ROOT.parent / "results" / "empirical" / "occupancy_validation.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    ok = [r for r in rows if r["lag_h"] != "" and r["group"] in ("FC1", "FC2", "FC3", "Tr14", "Tr37")]
    from scipy.stats import spearmanr
    rho = spearmanr([r["initial_cells"] for r in ok], [r["lag_h"] for r in ok])
    print("Spearman initial occupancy vs lag:", rho)
    print("Spearman arrival rate vs lag:", spearmanr([r["arrivals_per_h_first3h"] for r in ok], [r["lag_h"] for r in ok]))
    ex = [r for r in ok if r["group"] != "FC1"]
    print("Spearman arrival rate vs lag, exposed runs:",
          spearmanr([r["arrivals_per_h_first3h"] for r in ex], [r["lag_h"] for r in ex]))
    return rows


if __name__ == "__main__":
    for r in run():
        print(r)
