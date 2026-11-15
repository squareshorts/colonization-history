"""Main figures for the design results: missing cells (Fig. 4) and repeated contacts (Fig. 5).

  python paper/analysis/make_design_figures.py --results paper/results --figures paper/figures

Figure 4 reads fine_loss_detection.csv: joint mother-daughter detection with
the assumption-free lower bound for a moderate induced effect (predictor shift
0.85), under uniform loss that is independent of adhesion, applied to the
exposed arm only, to both arms equally, or to the sham arm at half the exposed
rate.

Figure 5 reads repeated_contact_information.csv (two contacts, oracle test) and
repeated_contact_design.csv (Tarone test, 2, 5 and 10 contacts) for
(p_H, p_L) = (0.10, 0.03).
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCEN = [("pulse_only_uniform", "Exposed arm only", "-", "o", "#2b2a27"),
        ("both_uniform_equal", "Both arms equally", (0, (5, 2)), "s", "#2a78d6"),
        ("sham_half_uniform", "Sham arm at half the exposed rate", (0, (1.2, 1.5)), "^", "#c05a28")]
KIN = [("logistic_exponential", "Logistic / exponential"),
       ("cloglog_weibull", "Complementary log–log / Weibull")]
HOLD = {0.: ("Immediately after transfer", "#8a8984", "o"),
        14.: ("After 14-h hold", "#2a78d6", "s"),
        37.: ("After 37-h hold", "#c05a28", "^")}


def read(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def fmt(n):
    """Two significant figures; scientific notation from 100,000 upward."""
    n = float(n)
    if n >= 1e5:
        e = int(np.floor(np.log10(n)))
        return f"{n / 10 ** e:.1f}×10$^{{{e}}}$"
    d = int(np.floor(np.log10(n))) - 1
    return f"{int(round(n, -d)):,}"


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=.18, lw=.6)


def fig4(results, out):
    rows = read(results / "fine_loss_detection.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
    for ax, (k, title) in zip(axs, KIN):
        for sc, lab, ls, mk, col in SCEN:
            rr = sorted((r for r in rows if r["kinetic"] == k and r["scenario"] == sc
                         and abs(float(r["effect"]) - .85) < 1e-9), key=lambda r: float(r["nominal_loss"]))
            ax.plot([100 * float(r["nominal_loss"]) for r in rr], [100 * float(r["bounded_joint"]) for r in rr],
                    ls=ls, marker=mk, ms=4.5, lw=1.6, color=col, mfc="white" if mk == "^" else col, label=lab)
        ax.set_title(title, fontsize=9, pad=4)
        ax.set_xlabel("Nominal exposed-arm loss (%)", fontsize=8)
        ax.set_xticks([0, 2, 4, 6, 8, 10])
        ax.set_ylim(-3, 103)
        style(ax)
    axs[0].set_ylabel("Joint mother–daughter detection (%)", fontsize=8)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7.5, loc="upper center", ncol=3, bbox_to_anchor=(.5, 1.06),
               title="Uniform loss, independent of adhesion, applied to", title_fontsize=7.5)
    fig.tight_layout(w_pad=1.5, rect=(0, 0, 1, .9))
    fig.savefig(out / "fine_loss_cliff.pdf", bbox_inches="tight")
    plt.close(fig)


def fig5(results, out):
    info = [r for r in read(results / "repeated_contact_information.csv") if abs(float(r["p_high"]) - .10) < 1e-9]
    design = [r for r in read(results / "repeated_contact_design.csv") if abs(float(r["p_high"]) - .10) < 1e-9]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1, 1.2], "wspace": .3})
    ax = axs[0]
    xs = [float(r["hold_h"]) for r in info]
    ys = [float(r["normal_oracle_n80"]) for r in info]
    ax.plot(xs, ys, color="#c9c8c2", lw=1, zorder=1)
    for x, y in zip(xs, ys):
        lab, col, mk = HOLD[x]
        ax.scatter([x], [y], s=40, marker=mk, color=col, zorder=3, edgecolor="white", linewidth=.6)
        ax.annotate(fmt(y), (x, y), xytext=(7, 5), textcoords="offset points", fontsize=7.5)
    ax.set_yscale("log")
    ax.set_xticks([0, 14, 37])
    ax.set_xlim(-4, 48)
    ax.set_ylim(3e3, 1e8)
    ax.set_xlabel("Time held in suspension (h)", fontsize=8)
    ax.set_ylabel("Cells required for 80% power", fontsize=8)
    ax.set_title("Two contacts per cell, known parameters", fontsize=9, pad=4)
    ax.text(-.12, 1.04, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    style(ax)

    ax = axs[1]
    for hold, (lab, col, mk) in HOLD.items():
        rr = sorted((r for r in design if float(r["hold_h"]) == hold and int(r["contacts_per_cell"]) in (2, 5, 10)),
                    key=lambda r: int(r["contacts_per_cell"]))
        k = [int(r["contacts_per_cell"]) for r in rr]
        n = [float(r["tarone_n80"]) for r in rr]
        ax.plot(k, n, color=col, lw=1.5, marker=mk, ms=5, mec="white", mew=.6, label=lab)
        for kk, nn in zip(k, n):
            ax.annotate(fmt(nn), (kk, nn), xytext=(5, 4), textcoords="offset points", fontsize=7)
    ax.set_yscale("log")
    ax.set_xticks([2, 5, 10])
    ax.set_xlim(1.3, 12.5)
    ax.set_ylim(1e2, 1e8)
    ax.set_xlabel("Contacts per cell", fontsize=8)
    ax.set_ylabel("Cells required for 80% power", fontsize=8)
    ax.set_title("Estimated parameters (Tarone test)", fontsize=9, pad=4)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.text(-.12, 1.04, "b", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
    style(ax)
    fig.savefig(out / "repeated_contact_requirements.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--figures", type=Path, required=True)
    a = ap.parse_args()
    a.figures.mkdir(parents=True, exist_ok=True)
    fig4(a.results, a.figures)
    fig5(a.results, a.figures)


if __name__ == "__main__":
    main()
