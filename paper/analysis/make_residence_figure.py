"""Main figure: residence after attachment and five-minute retention per arrival.

  python paper/analysis/make_residence_figure.py --repo . --results paper/results \
      --out paper/figures/residence_composition.pdf
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from residence_composition import EXTRA_FC1, FRAME_S, GROUPS, load_group

COL = {"FC1": "#52514e", "FC2": "#2a78d6", "FC3": "#1baf7a", "Tr14": "#eb6834", "Tr37": "#c98500"}
LAB = {"FC1": "First exposure", "FC2": "After transfer", "FC3": "Second transfer",
       "Tr14": "14-h hold", "Tr37": "37-h hold"}


def km_curve(rows, grid_frames):
    f = np.array([r["frames"] for r in rows if r["frames"] >= 2])
    c = np.array([r["cens"] for r in rows if r["frames"] >= 2])
    out, s, j = [], 1.0, 1
    for k in grid_frames:
        while j < k:
            at = np.sum(f >= j)
            ev = np.sum((f == j) & ~c)
            if at:
                s *= 1 - ev / at
            j += 1
        out.append(s)
    return np.array(out)


def model_curve(weights, rates, ages):
    w, r = np.asarray(weights), np.asarray(rates)
    S = np.exp(-np.outer(ages, r)) @ w
    return S / (np.exp(-r * FRAME_S) @ w)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    d = json.loads((a.results / "residence_composition.json").read_text())
    fits = d["models"]["matched"]["3_types"]["fits"]

    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.1), gridspec_kw={"width_ratios": [1.25, 1], "wspace": .3})
    ax = axs[0]
    frames = np.unique(np.geomspace(2, 1200, 160).astype(int))
    ages = (frames - 1) * FRAME_S
    for g in ("FC1", "FC2"):
        rows = load_group(a.repo, GROUPS[g])
        ax.step(ages, km_curve(rows, frames), where="post", color=COL[g], lw=1.6, label=f"{LAB[g]} (tracks)")
    rw, sh = fits["composition"], fits["rate_shift"]
    fine = np.geomspace(3, 3600, 300)
    for g, key in (("FC1", "weights_FC1"), ("FC2", "weights_FC2")):
        ax.plot(fine, model_curve(rw[key], rw["rates_per_s"], fine), color=COL[g], lw=1.0, ls=(0, (4, 2)))
        rates = np.array(sh["rates_per_s"]) * (np.exp(-sh["delta"]) if g == "FC2" else 1.)
        ax.plot(fine, model_curve(sh[key], rates, fine), color=COL[g], lw=1.0, ls=(0, (1, 1.5)))
    ax.plot([], [], color="#6b6a64", lw=1.0, ls=(0, (4, 2)), label="Composition change or age advance (fit)")
    ax.plot([], [], color="#6b6a64", lw=1.0, ls=(0, (1, 1.5)), label="Rate shift (fit)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3, 3600)
    ax.set_ylim(.02, 1.05)
    ax.axvline(300, color="#9a9890", lw=.6, ls=(0, (2, 2)), zorder=0)
    ax.set_xlabel("Time since attachment (s)", fontsize=8)
    ax.set_ylabel("Fraction of tracks persisting", fontsize=8)
    ax.legend(frameon=False, fontsize=6.8, loc="lower left")
    ax.set_title("Persistence after attachment, matched fields", fontsize=9, pad=4)
    ax.text(-.06, 1.05, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")

    ax = axs[1]
    order = ["FC1", "FC2", "FC3", "Tr14", "Tr37"]
    for i, g in enumerate(order):
        r = d["retention"][g]
        est = r["retained_300s_per_arrival"]
        lo, hi = r["retained_300s_per_arrival_ci95_within_fields"]
        ax.errorbar([i], [est], yerr=[[est - lo], [hi - est]], fmt="o", ms=5.5, color=COL[g],
                    ecolor=COL[g], elinewidth=1.2, capsize=0, zorder=3)
        if len(r["per_field_retained_300s"]) > 1:
            pf = list(r["per_field_retained_300s"].values())
            ax.scatter(np.full(len(pf), i + .22), pf, s=12, facecolor="white", edgecolor=COL[g],
                       linewidth=.9, zorder=2)
        if g == "FC1":
            extra = d["retention"]["FC1_with_first_hour_recording"]["per_field_retained_300s"][EXTRA_FC1]
            ax.scatter([i + .22], [extra], s=16, marker="x", color=COL[g], linewidth=.9, zorder=2)
    ax.set_xticks(range(len(order)), ["First\nexposure", "After\ntransfer", "Second\ntransfer",
                                      "14-h\nhold", "37-h\nhold"], fontsize=6.6)
    ax.set_xlim(-.5, len(order) - .4)
    ax.set_ylim(0, .17)
    ax.set_ylabel("Persisting 5 min after arrival", fontsize=8)
    ax.set_title("Five-minute persistence per arrival", fontsize=9, pad=4)
    ax.text(-.06, 1.05, "b", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    for x in axs:
        x.spines[["top", "right"]].set_visible(False)
        x.tick_params(labelsize=8)
        x.grid(axis="y", alpha=.18, lw=.6)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
