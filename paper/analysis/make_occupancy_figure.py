"""Main Figure 1: occupancy trajectories, threshold times and early occupancy ratio.

Inputs are produced by matched_progression.py and
endpoint_observation_comparison.py.

  python paper/analysis/make_occupancy_figure.py --results paper/results \
      --out paper/figures/occupancy_endpoints.pdf
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COND = {  # condition -> (label, colour); colours validated for CVD separation
    "FC1": ("First exposure", "#52514e"),
    "FC2": ("After transfer", "#2a78d6"),
    "Tr14": ("14-h hold", "#eb6834"),
    "FC3": ("Second transfer", "#1baf7a"),
    "Tr37": ("37-h hold", "#c98500"),
}
SERIES_MARK = {1: "o", 2: "s", 3: "^"}
FLOOR = 0.3


def cond_of(record):
    return record.split("-")[0]


def read_csv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=.18, lw=.6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    traj = defaultdict(lambda: ([], []))
    for r in read_csv(a.results / "occupancy_trajectories.csv"):
        traj[r["record"]][0].append(float(r["time_h"]))
        traj[r["record"]][1].append(float(r["cell_eq"]))
    ends = {r["record"]: r for r in read_csv(a.results / "endpoint_observation_comparison.csv")}
    dyn = {r["record"]: r for r in read_csv(a.results / "occupancy_dynamics.csv")}

    fig = plt.figure(figsize=(7.4, 6.3))
    gs = fig.add_gridspec(2, 12, height_ratios=[1, 1.02], hspace=.62, wspace=1.6)
    top = [fig.add_subplot(gs[0, 3 * i:3 * i + 3]) for i in range(4)]
    for ax in top[1:]:
        ax.sharey(top[0])
    panels = [("Series 1", ["FC1-1", "FC2-1"]), ("Series 2", ["FC1-2", "FC2-2"]),
              ("Series 3", ["FC1-3", "FC2-3", "FC3-3"]), ("Held before transfer", ["Tr14", "Tr37"])]
    for k, (ax, (title, recs)) in enumerate(zip(top, panels)):
        for q in (20, 30, 50):
            ax.axhline(q, color="#9a9890", lw=.6, ls=(0, (2, 2)), zorder=0)
        ax.axhspan(FLOOR * .9, FLOOR * 1.12, color="#d9d8d3", alpha=.6, lw=0, zorder=0)
        for rec in recs:
            t, y = map(np.asarray, traj[rec])
            label, col = COND[cond_of(rec)]
            ax.plot(t, np.maximum(y, FLOOR), color=col, lw=1.4, solid_capstyle="round")
        ax.set_yscale("log")
        ax.set_ylim(FLOOR * .9, 700)
        ax.set_xlim(0, 30)
        ax.set_xticks([0, 10, 20, 30])
        ax.set_title(title, fontsize=9, pad=4)
        ax.set_xlabel("Time after inoculation (h)", fontsize=8)
        style(ax)
        ax.text(-.04, 1.06, "abcd"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
        if k:
            plt.setp(ax.get_yticklabels(), visible=False)
    top[0].set_ylabel("Occupied area (cell equivalents)", fontsize=8)
    for q in (20, 30, 50):
        top[3].text(29.5, q * 1.04, f"{q}", fontsize=6.5, color="#6b6a64", ha="right", va="bottom")
    top[0].set_yticks([FLOOR, 1, 10, 100], ["\u22640.3", "1", "10", "100"])
    cond_handles = [plt.Line2D([], [], color=COND[c][1], lw=2, label=COND[c][0])
                    for c in ("FC1", "FC2", "FC3", "Tr14", "Tr37")]
    fig.legend(handles=cond_handles, loc="upper center", ncol=5, frameon=False, fontsize=7.5,
               bbox_to_anchor=(.5, .985), handlelength=1.6, columnspacing=1.4)

    bottom = [fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6]), fig.add_subplot(gs[1, 6:12])]
    for k, (ax, key, title) in enumerate(zip(bottom[:2], ["relative5_h", "absolute30_h"],
                                             ["Fivefold increase\n(per recording)",
                                              "30 cell equivalents\n(common scale)"])):
        status_key = key.replace("_h", "_status")
        for s in (1, 2, 3):
            y0, y1 = float(ends[f"FC1-{s}"][key]), float(ends[f"FC2-{s}"][key])
            ax.plot([0, 1], [y0, y1], color="#9a9890", lw=1.1, zorder=1)
            for x, y, rec, c in ((0, y0, f"FC1-{s}", "FC1"), (1, y1, f"FC2-{s}", "FC2")):
                censored = ends[rec].get(status_key, "observed") in ("at_start", "left_censored")
                if censored:
                    ax.scatter([x], [y], s=34, marker=SERIES_MARK[s], facecolor="white",
                               edgecolor=COND[c][1], linewidth=1.2, zorder=3)
                    ax.annotate("\u2264" + f"{y + 1e-9:.2f} h", (x, y), xytext=(-7, 1), textcoords="offset points",
                                fontsize=6.5, color="#2b2a27", va="center", ha="right")
                else:
                    ax.scatter([x], [y], s=34, marker=SERIES_MARK[s], color=COND[c][1], zorder=3,
                               edgecolor="white", linewidth=.8)
        ax.set_xticks([0, 1], ["First\nexposure", "After\ntransfer"], fontsize=8)
        ax.set_xlim(-.3, 1.3)
        ax.set_ylim(-1, 29)
        ax.set_title(title, fontsize=9, pad=4)
        style(ax)
        ax.text(-.04, 1.06, "ef"[k], transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    bottom[0].set_ylabel("Time to sustained threshold (h)", fontsize=8)
    bottom[1].sharey(bottom[0])
    plt.setp(bottom[1].get_yticklabels(), visible=False)
    handles = [plt.Line2D([], [], marker=SERIES_MARK[s], ls="", color="#52514e", ms=5.5, label=f"Series {s}")
               for s in (1, 2, 3)]
    cens = plt.Line2D([], [], marker="o", ls="", mfc="white", mec="#52514e", ms=5.5, label="Left-censored")
    bottom[0].legend(handles=handles + [cens], frameon=False, fontsize=6.8, loc="lower left",
                     handletextpad=.2, borderaxespad=.1)

    ax = bottom[2]
    for s in (1, 2, 3):
        v0 = float(dyn[f"FC1-{s}"]["occupancy_ratio_4_8h"])
        v1 = float(dyn[f"FC2-{s}"]["occupancy_ratio_4_8h"])
        dx = {1: -.1, 2: 0., 3: .1}[s]
        ax.plot([dx, 1 + dx], [v0, v1], color="#9a9890", lw=1.1, zorder=1)
        ax.scatter([dx], [v0], s=34, marker=SERIES_MARK[s], color=COND["FC1"][1], edgecolor="white",
                   linewidth=.8, zorder=3)
        ax.scatter([1 + dx], [v1], s=34, marker=SERIES_MARK[s], color=COND["FC2"][1], edgecolor="white",
                   linewidth=.8, zorder=3)
    singles = [("FC3-3", "FC3", "^"), ("Tr14", "Tr14", "D"), ("Tr37", "Tr37", "D")]
    for i, (rec, c, mk) in enumerate(singles):
        ax.scatter([2.4 + i], [float(dyn[rec]["occupancy_ratio_4_8h"])], s=34, marker=mk, color=COND[c][1],
                   edgecolor="white", linewidth=.8, zorder=3)
    ax.axvline(1.7, color="#c9c8c2", lw=.8, zorder=0)
    ax.text(.5, 6.6, "Paired series", ha="center", fontsize=7, color="#6b6a64")
    ax.text(3.4, 6.6, "Single recordings", ha="center", fontsize=7, color="#6b6a64")
    ax.axhline(1, color="#9a9890", lw=.7, ls=(0, (2, 2)), zorder=0)
    ax.set_xticks([0, 1, 2.4, 3.4, 4.4], ["First\nexposure", "After\ntransfer", "Second\ntransfer",
                                          "14-h\nhold", "37-h\nhold"], fontsize=7.2)
    ax.set_xlim(-.45, 4.9)
    ax.set_ylim(-.3, 7.2)
    ax.set_ylabel("Occupancy at 4–8 h / 0–2 h", fontsize=8)
    ax.set_title("Early occupancy ratio", fontsize=9, pad=4)
    style(ax)
    ax.text(-.04, 1.06, "g", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
