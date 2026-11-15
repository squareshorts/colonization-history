"""Main figure: persistence of linked image components in the matched sequences.

  python paper/analysis/make_residence_figure.py --repo . --results paper/results \
      --out paper/figures/residence_composition.pdf

(a) Kaplan-Meier persistence after attachment for the pooled matched fields,
    with three-component composition (equivalently age-advance) and rate-shift
    fits for each condition. Tracks still present at the end of a one-hour
    window are right-censored; curves stop at the longest follow-up.
(b) Five-minute persistence per arrival in each matched sequence, connected
    within sequence, with 95% intervals from resampling arrivals within that
    field.
(c) Sequence-level differences with their mean and a Student-t interval
    across the three sequences, the between-experiment uncertainty.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from residence_composition import EXTRA_FC1, FRAME_S, GROUPS, load_group

C_FIRST, C_TRANSFER, INK, MUTED = "#52514e", "#2a78d6", "#2b2a27", "#8a8984"
SERIES_MARK = {1: "o", 2: "s", 3: "^"}


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
    return np.array(out), int(f.max())


def model_curve(weights, rates, ages):
    w, r = np.asarray(weights), np.asarray(rates)
    S = np.exp(-np.outer(ages, r)) @ w
    return S / (np.exp(-r * FRAME_S) @ w)


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=.18, lw=.6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    d = json.loads((a.results / "residence_composition.json").read_text())
    fits = d["models"]["matched"]["3_types"]["fits"]
    comp, shift = fits["composition"], fits["rate_shift"]

    fig = plt.figure(figsize=(7.6, 3.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.4, .95, .7], wspace=.5)
    ax = fig.add_subplot(gs[0])
    frames = np.unique(np.geomspace(2, 1200, 160).astype(int))
    fine = np.geomspace(3, 3600, 300)
    for g, col, key, lab in (("FC1", C_FIRST, "weights_FC1", "First exposure"),
                             ("FC2", C_TRANSFER, "weights_FC2", "After transfer")):
        rows = load_group(a.repo, GROUPS[g])
        km, fmax = km_curve(rows, frames)
        keep = frames <= fmax
        ages = (frames[keep] - 1) * FRAME_S
        ax.step(ages, km[keep], where="post", color=col, lw=1.7)
        ax.plot(fine, model_curve(comp[key], comp["rates_per_s"], fine), color=col, lw=1.0, ls=(0, (4, 2)))
        rates = np.array(shift["rates_per_s"]) * (np.exp(-shift["delta"]) if g == "FC2" else 1.)
        ax.plot(fine, model_curve(shift[key], rates, fine), color=col, lw=1.0, ls=(0, (1, 1.5)))
    handles = [plt.Line2D([], [], color=C_FIRST, lw=2.4, label="First exposure"),
               plt.Line2D([], [], color=C_TRANSFER, lw=2.4, label="After transfer"),
               plt.Line2D([], [], color=MUTED, lw=1.7, label="Kaplan-Meier (tracks)"),
               plt.Line2D([], [], color=MUTED, lw=1.0, ls=(0, (4, 2)), label="Composition / age-advance fit"),
               plt.Line2D([], [], color=MUTED, lw=1.0, ls=(0, (1, 1.5)), label="Rate-shift fit")]
    ax.legend(handles=handles, frameon=False, fontsize=6.6, loc="lower left")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(3, 3600); ax.set_ylim(.02, 1.05)
    ax.axvline(300, color="#9a9890", lw=.6, ls=(0, (2, 2)), zorder=0)
    ax.set_xlabel("Time since attachment (s)", fontsize=8)
    ax.set_ylabel("Fraction of linked components persisting", fontsize=8)
    ax.set_title("Matched sequences, pooled", fontsize=9, pad=4)
    ax.text(-.08, 1.04, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    style(ax)

    ax = fig.add_subplot(gs[1])
    for p in d["matched_pairs"]["pairs"]:
        s = p["series"]
        y = [p["first_exposure"], p["after_transfer"]]
        lo = [p["first_exposure_ci95_within_field"][0], p["after_transfer_ci95_within_field"][0]]
        hi = [p["first_exposure_ci95_within_field"][1], p["after_transfer_ci95_within_field"][1]]
        dx = {1: -.08, 2: 0., 3: .08}[s]
        xs = np.array([0, 1]) + dx
        ax.plot(xs, y, color=MUTED, lw=1.0, zorder=1)
        for x, yy, l, h, col in zip(xs, y, lo, hi, (C_FIRST, C_TRANSFER)):
            ax.errorbar([x], [yy], yerr=[[yy - l], [h - yy]], fmt=SERIES_MARK[s], ms=5, color=col,
                        ecolor=col, elinewidth=1.0, capsize=0, zorder=3, mec="white", mew=.6)
    extra = d["retention"]["FC1_with_first_hour_recording"]["per_field_retained_300s"][EXTRA_FC1]
    ax.scatter([-.36], [extra], marker="x", s=22, color=C_FIRST, lw=1.0, zorder=3)
    ax.annotate("unmatched\ncontrol", (-.36, extra), xytext=(0, -6), textcoords="offset points",
                ha="center", va="top", fontsize=6, color=MUTED)
    hs = [plt.Line2D([], [], marker=SERIES_MARK[s], ls="", color=MUTED, ms=5, label=f"Sequence {s}")
          for s in (1, 2, 3)]
    ax.legend(handles=hs, frameon=False, fontsize=6.6, loc="upper left", handletextpad=.2)
    ax.set_xticks([0, 1], ["First\nexposure", "After\ntransfer"], fontsize=7.5)
    ax.set_xlim(-.6, 1.35); ax.set_ylim(0, .17)
    ax.set_ylabel("Persisting 5 min after arrival", fontsize=8)
    ax.set_title("Paired fields", fontsize=9, pad=4)
    ax.text(-.15, 1.04, "b", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    style(ax)

    ax = fig.add_subplot(gs[2])
    pr = d["matched_pairs"]
    for p in pr["pairs"]:
        ax.scatter([0], [p["difference"]], marker=SERIES_MARK[p["series"]], s=26, color=MUTED,
                   edgecolor="white", linewidth=.6, zorder=3)
    lo, hi = pr["t95_interval_between_sequences"]
    m = pr["mean_difference"]
    ax.errorbar([.55], [m], yerr=[[m - lo], [hi - m]], fmt="D", ms=5, color=INK, elinewidth=1.2, capsize=3)
    ax.axhline(0, color="#9a9890", lw=.7, ls=(0, (2, 2)), zorder=0)
    ax.set_xticks([0, .7], ["Pairs", "Mean"], fontsize=7.5)
    ax.set_xlim(-.4, 1.1)
    ax.set_ylabel("Difference, transfer minus first", fontsize=8)
    ax.set_title("Between experiments", fontsize=9, pad=4)
    ax.text(-.3, 1.04, "c", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    style(ax)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
