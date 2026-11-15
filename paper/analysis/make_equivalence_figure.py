"""Main figure: identical colonization-lag curves from selection and a fading cellular response.

  python paper/analysis/make_equivalence_figure.py --out paper/figures/lag_carrier_equivalence.pdf

The selection curve uses equation (2) with f0 and c chosen so that the curve
passes through the published lags (Lee et al. 2018), which therefore calibrate
the parameters and do not validate the model. The fading-response curve is
computed separately from a Hill response to an exponentially decaying state,
with nu = 2 and tau, K/m0 set by the identities of equation (3); the two curves
coincide to numerical precision.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import expit, logit

F0, C = .990359, .1456747
LN, L0 = 17.56, 4.68
PUBLISHED = ((0., 4.68), (14., 5.46), (37., 13.41))


def selection_lag(t):
    f = expit(logit(F0) - C * t)
    return LN - (LN - L0) * f / F0


def fading_response_lag(t, nu=2.):
    tau = nu / C
    m0 = 1.
    K = m0 * ((1 - F0) / F0) ** (1 / nu)
    m = m0 * np.exp(-t / tau)
    h = 1. / (1. + (K / m) ** nu)
    return LN - (LN - L0) * h / F0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    t = np.linspace(0, 60, 601)
    sel, hill = selection_lag(t), fading_response_lag(t)
    assert np.max(np.abs(sel - hill)) < 1e-9
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.axhline(LN, ls=(0, (4, 3)), lw=1, color="#8a8984", zorder=0)
    ax.text(59, LN + .35, "Lag before surface exposure (17.56 h)", ha="right", va="bottom", fontsize=7.5,
            color="#6b6a64")
    ax.plot(t, sel, lw=3.2, color="#2a78d6", alpha=.45, label="Selection among persistent phenotypes")
    ax.plot(t, hill, lw=1.3, color="#1d3f6e", ls=(0, (5, 2.5)), label="Fading response within each cell")
    ax.scatter([p[0] for p in PUBLISHED], [p[1] for p in PUBLISHED], s=46, facecolor="white",
               edgecolor="#c05a28", linewidth=1.5, zorder=3, label="Published fitted lags (used for calibration)")
    ax.set_xlabel("Time held in suspension before transfer (h)", fontsize=9)
    ax.set_ylabel("Colonization lag (h)", fontsize=9)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 19.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=.18, lw=.6)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
