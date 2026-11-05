"""Manuscript figures from completed result tables (no simulation)."""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from model import MECHANISMS, LABELS, SHORT
from run_experiments import COLORS, STYLE

R = Path("results")
FIG = R / "figures"
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "font.family": "DejaVu Sans", "axes.titlesize": 9,
                     "axes.titleweight": "bold", "axes.titlelocation": "left"})
GROUP_COLORS = {"FC1": "#5b6470", "FC2": "#1a5c86", "FC3": "#5c5597", "Tr14": "#5fa7c9", "Tr37": "#a9cfe2"}
GROUP_LABELS = {"FC1": "FC1 fresh", "FC2": "FC2 exposed, 0 h", "FC3": "FC3 twice exposed",
                "Tr14": "exposed, held 14 h", "Tr37": "exposed, held 37 h"}


def read(p):
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def fig1():
    fig = plt.figure(figsize=(7.2, 3.0), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[0]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    ax.set_title("a  Preparation, gap and challenge")
    box = lambda x, y, w, h, t, c: (ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                                                 fc=c, ec="none", alpha=.9)),
                                    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=7.5))
    box(0.0, 3.6, 2.7, 1.5, "transient contacts\nend in complete\ndetachment", "#e3ecf2")
    box(3.3, 3.6, 2.5, 1.5, "gap $g$:\nhistories decay,\nresidents lost", "#f2ece3")
    box(6.4, 3.6, 3.5, 1.5, "challenge with the\nsame dose and\nschedule in all arms", "#e8e6f2")
    for x in (2.78, 5.88):
        ax.annotate("", xy=(x + .45, 4.35), xytext=(x, 4.35), arrowprops=dict(arrowstyle="->", color="#444"))
    arms = [("fresh cells", "fresh surface", "FF"), ("exposed cells", "fresh surface", "EF"),
            ("fresh cells", "visited surface", "FV"), ("exposed cells", "visited surface", "EV")]
    for i, (c, s, k) in enumerate(arms):
        x = 0.3 + 2.45 * i
        ax.add_patch(FancyBboxPatch((x, 0.5), 2.1, 2.3, boxstyle="round,pad=0.05",
                                    fc="#ffffff", ec="#9aa7b1", lw=.8))
        ax.text(x + 1.05, 2.25, k, ha="center", fontweight="bold")
        ax.text(x + 1.05, 1.35, f"{c}\n{s}", ha="center", va="center", fontsize=7)
    ax = fig.add_subplot(gs[1]); ax.set_title("b  Contrasts moved by each hypothesis")
    pred = {"passive": (0, 0), "carryover": (0, 1), "selection": (1, 0), "cell": (1, 0),
            "surface": (0, 1), "both": (1, 1)}
    extra = {"carryover": "residents at $t=0$", "selection": "no acquisition", "cell": "acquired on surface",
             "surface": "deposited by contacts", "both": "", "passive": ""}
    for i, m in enumerate(MECHANISMS):
        y = len(MECHANISMS) - 1 - i
        for j in range(2):
            on = pred[m][j]
            ax.scatter(j, y, s=170, color=COLORS[m] if on else "#ffffff", edgecolor=COLORS[m], lw=1.2)
        ax.text(1.45, y, extra[m], va="center", fontsize=7, color="#444")
    ax.set_yticks(range(6), [LABELS[m] for m in MECHANISMS][::-1], fontsize=7.5)
    ax.set_xticks([0, 1], ["EF $-$ FF", "FV $-$ FF"])
    ax.set_xlim(-.5, 3.2); ax.set_ylim(-.6, 5.6)
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    fig.savefig(FIG / "fig1_design.png"); plt.close(fig)


def fig4(sigs):
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.7), layout="constrained")
    def sig(m, arm, k):
        return next(s[k] for s in sigs if s["mechanism"] == m and s["arm"] == arm and s["gap"] == 0)
    x = np.arange(6)
    ratio = [sig(m, "exposed_fresh", "attach_late") / max(sig(m, "exposed_fresh", "attach_early"), 1e-9) /
             (sig(m, "fresh_fresh", "attach_late") / max(sig(m, "fresh_fresh", "attach_early"), 1e-9))
             for m in MECHANISMS]
    axs[0].axhline(1, color="#999", lw=.6)
    axs[0].vlines(x, 1, ratio, color=[COLORS[m] for m in MECHANISMS], lw=2)
    axs[0].scatter(x, ratio, color=[COLORS[m] for m in MECHANISMS], zorder=3, s=22)
    axs[0].set(ylabel="Late/early attachment,\nexposed relative to fresh", title="a  Advantage across influx")
    axs[1].bar(x, [100 * sig(m, "fresh_conditioned", "residual_fraction_final") for m in MECHANISMS],
               color=[COLORS[m] for m in MECHANISMS])
    axs[1].set(ylabel="Final cells from residents (%)", title="b  Resident lineages")
    fq = [100 * (sig(m, "exposed_fresh", "founder_quick_fraction") - sig(m, "fresh_fresh", "founder_quick_fraction"))
          for m in MECHANISMS]
    axs[2].bar(x, fq, color=[COLORS[m] for m in MECHANISMS])
    axs[2].axhline(0, color="#999", lw=.6)
    axs[2].set(ylabel="Change in founders leaving\nwithin 1 time unit (pp)", title="c  Early founder loss")
    for ax in axs:
        ax.set_xticks(x, [SHORT[m] for m in MECHANISMS], rotation=40, ha="right")
    fig.savefig(FIG / "fig4_event_signatures.png"); plt.close(fig)


def fig5():
    acc = read(R / "tables" / "recovery_accuracy.csv")
    det = json.loads((R / "logs" / "recovery_details.json").read_text())
    fig = plt.figure(figsize=(7.2, 3.5), layout="constrained")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1])
    ax = fig.add_subplot(gs[0])
    styles = {"full": ("-", "o"), "combined": ("--", "s"), "exposed": (":", "^")}
    names = {"full": "full 2$\\times$2", "combined": "FF + EV", "exposed": "FF + EF"}
    for tier, col in (("tracking", "#1a5c86"), ("counts", "#b5642a")):
        for arms in ("full", "combined", "exposed"):
            rr = sorted([r for r in acc if r["tier"] == tier and r["arms"] == arms and r["gaps"] == "0-4-16"],
                        key=lambda r: r["replicates"])
            ax.errorbar([r["replicates"] for r in rr], [100 * r["accuracy"] for r in rr],
                        yerr=[196 * r["se"] for r in rr], color=col, ls=styles[arms][0],
                        marker=styles[arms][1], ms=3.5, capsize=1.5, lw=1.1,
                        label=f"{tier}, {names[arms]}")
        rr = sorted([r for r in acc if r["tier"] == tier and r["arms"] == "baseline"], key=lambda r: r["replicates"])
        ax.plot([r["replicates"] for r in rr], [100 * r["accuracy"] for r in rr], color=col, lw=.8, alpha=.5,
                label=f"{tier}, FF only (control)")
    ax.axhline(100 / 6, color="#999", lw=.6)
    ax.text(8.1, 100 / 6 + 1, "chance", fontsize=6.5, color="#777", ha="right")
    ax.set(xlabel="Replicate surfaces per arm", ylabel="Held-out accuracy (%)", xticks=[2, 4, 8],
           title="a  Recovery by design", ylim=(10, 80))
    ax.legend(fontsize=5.8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    for k, (key, title) in enumerate((("counts|n8|gaps0-4-16|full", "b  Counts, full design"),
                                      ("tracking|n8|gaps0-4-16|full", "c  Tracking, full design"))):
        ax = fig.add_subplot(gs[k + 1])
        cm = np.array(det[key]["confusion"], float)
        cm = cm / cm.sum(axis=1, keepdims=True)
        ax.imshow(cm, cmap="Blues", vmin=0, vmax=1)
        for i in range(6):
            for j in range(6):
                ax.text(j, i, f"{100*cm[i,j]:.0f}", ha="center", va="center", fontsize=6,
                        color="white" if cm[i, j] > .55 else "#203040")
        ax.set_xticks(range(6), [SHORT[m] for m in MECHANISMS], rotation=50, ha="right", fontsize=6.5)
        ax.set_yticks(range(6), [SHORT[m] for m in MECHANISMS] if k == 0 else [], fontsize=6.5)
        ax.set_xlabel("Assigned", fontsize=7)
        if k == 0:
            ax.set_ylabel("Generating", fontsize=7)
        ax.set_title(title)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.savefig(FIG / "fig5_recovery.png"); plt.close(fig)


def fig6():
    E = json.loads((R / "empirical" / "empirical_fits.json").read_text())
    val = read(R / "empirical" / "occupancy_validation.csv")
    fig, axs = plt.subplots(2, 2, figsize=(6.4, 5.0), layout="constrained")
    axs = axs.ravel()
    ax = axs[0]
    seen = set()
    for n, o in E["occupancy"].items():
        g = o["group"]
        if g not in GROUP_COLORS:
            continue
        t, c = np.array(o["time_h"]), np.array(o["cells"])
        k = 10
        cs = np.convolve(c, np.ones(k) / k, mode="valid")
        ax.plot(t[:len(cs)], np.maximum(cs, .3), color=GROUP_COLORS[g], lw=1,
                label=GROUP_LABELS[g] if g not in seen else None)
        seen.add(g)
    ax.set_yscale("log")
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    ax.set(xlabel="Time after inoculation (h)", ylabel="Attached cell equivalents", title="a  Surface occupancy")
    ax = axs[1]
    for r in val:
        if r["group"] in GROUP_COLORS and r["lag_h"] != "":
            ax.scatter(r["arrivals_per_h_first3h"], r["lag_h"], color=GROUP_COLORS[r["group"]], s=22,
                       edgecolor="white", lw=.5)
    ax.set_xscale("log")
    ax.set(xlabel="Arrivals per hour, first 3 h", ylabel="Lag to exponential growth (h)", title="b  Dose and lag")
    ax = axs[2]
    grid = np.array(E["groups"]["FC1"]["km_grid_s"], float)
    sel = E["families"]["selection"]
    import sys
    sys.path.insert(0, "empirical")
    import analyze as A
    a = np.logspace(np.log10(6), np.log10(3000), 80)
    for g, v in E["groups"].items():
        ax.scatter(grid, v["km"], color=GROUP_COLORS[g], s=10, zorder=3)
        S = A.S_frailty(a, sel["lambda_per_s"], sel["sigma"], sel["gamma"], sel["acquisition_rate_per_s"],
                        sel["adhesiveness_shift"][g], 0.0)
        S0 = A.S_frailty([3.0], sel["lambda_per_s"], sel["sigma"], sel["gamma"], sel["acquisition_rate_per_s"],
                         sel["adhesiveness_shift"][g], 0.0)[0]
        ax.plot(a, S / S0, color=GROUP_COLORS[g], lw=1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.008, 1.05)
    ax.set(xlabel="Residence after attachment (s)", ylabel="Fraction still attached",
           title="c  Residence of attached arrivals")
    ax = axs[3]
    ci = E["adhesiveness_shift_ci"]
    gaps = {"FC2": 0, "Tr14": 14, "Tr37": 37}
    for g, x in gaps.items():
        y = sel["adhesiveness_shift"][g]
        ax.errorbar(x, y, yerr=[[y - ci[g][0]], [ci[g][1] - y]], color=GROUP_COLORS[g], marker="o", ms=5, capsize=2)
    y = sel["adhesiveness_shift"]["FC3"]
    ax.errorbar(-2.5, y, yerr=[[y - ci["FC3"][0]], [ci["FC3"][1] - y]], color=GROUP_COLORS["FC3"], marker="D",
                ms=4, capsize=2)
    ax.axhline(0, color="#999", lw=.6)
    ax.set(xlabel="Hours held off-surface before challenge", ylabel="Adhesiveness shift relative to FC1",
           title="d  Exposure effect after holding")
    ax.annotate("FC3", (-2.5, y), xytext=(-7.5, y + .12), fontsize=6.5, color=GROUP_COLORS["FC3"])
    ax.set_xlim(-9, 40)
    fig.savefig(FIG / "fig6_empirical.png"); plt.close(fig)


def km_ci(f, c, t=300.0):
    """Kaplan-Meier S(t) with Greenwood 95% interval (log-log scale)."""
    te = f * 3.0
    S, var = 1.0, 0.0
    for u in np.unique(te[(te <= t) & ~c]):
        n = np.sum(te >= u); d = np.sum((te == u) & ~c)
        if n > d > 0:
            S *= 1 - d / n; var += d / (n * (n - d))
    if S <= 0 or S >= 1:
        return S, S, S
    se = np.sqrt(var) / abs(np.log(S))
    lo, hi = S ** np.exp(1.96 * se), S ** np.exp(-1.96 * se)
    return S, lo, hi


def fig7():
    import sys
    sys.path.insert(0, "empirical")
    import analyze as A
    import analyze_v2 as V
    import design_power as D
    E1 = json.loads((R / "empirical" / "empirical_fits.json").read_text())
    E2 = json.loads((R / "empirical_v2" / "fits_v2.json").read_text())
    P = json.loads((R / "empirical_v2" / "design_power.json").read_text())
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 5.2), layout="constrained")
    axs = axs.ravel()
    ax = axs[0]
    bins = {"FC1": [0, 3, 8, 20, 30], "FC2": [0, 1, 2, 3, 5, 8], "FC3": [0, 1, 2, 3],
            "Tr14": [0, 3, 6, 11], "Tr37": [0, 3, 10, 20]}
    for g, bb in bins.items():
        recs = V.group(V.GROUPS[g], "extended")
        t = np.array([r["start"] * 3 / 3600 for r in recs])
        f = np.array([r["frames"] for r in recs]); c = np.array([r["cens"] for r in recs])
        xs, ys, lo, hi = [], [], [], []
        for a, b in zip(bb[:-1], bb[1:]):
            s = (t >= a) & (t < b)
            if s.sum() >= 30:
                S, l, h = km_ci(f[s], c[s])
                xs.append(0.5 * (a + b)); ys.append(100 * S); lo.append(100 * (S - l)); hi.append(100 * (h - S))
        ax.errorbar(xs, ys, yerr=[lo, hi], color=GROUP_COLORS[g], marker="o", ms=3.5, lw=1, capsize=1.5,
                    label=GROUP_LABELS[g])
    ax.set(xlabel="Arrival time after inoculation (h)", ylabel="Still attached at 5 min (%)",
           title="a  Retention drifts with time")
    ax.legend(fontsize=6, frameon=False)
    ax = axs[1]
    specs = [("first 3 h, windowed", E1["delta_loglik_history_minus_selection"] * -1, None),
             ("first 3 h, continuous", E2["P1_first3h"]["dll_offset_minus_head"], E2["P1_first3h"].get("bootstrap_dll_q")),
             ("full records", E2["P1_extended"]["dll_offset_minus_head"], E2["P1_extended"].get("bootstrap_dll_q")),
             ("full, time-matched", E2["time_matched"]["dll_offset_minus_head"], E2["time_matched"].get("bootstrap_dll_q")),
             ("full, occupancy-adjusted", E2["occupancy_adjusted"]["dll_offset_minus_head"],
              E2["occupancy_adjusted"].get("bootstrap_dll_q"))]
    ys = np.arange(len(specs))[::-1]
    for y, (name, v, q) in zip(ys, specs):
        if q:
            ax.plot([q[0], q[2]], [y, y], color="#1a5c86", lw=1.5)
        ax.scatter(v, y, color="#1a5c86", zorder=3, s=20)
    ax.axvline(0, color="#999", lw=.6)
    ax.set_xlim(-5, 90)
    ax.set_yticks(ys, [s[0] for s in specs], fontsize=7)
    ax.set(xlabel="log L(offset) $-$ log L(head start)", title="b  Offset vs head start")
    ax = axs[2]
    o = E1["families"]["selection"]
    rng = np.random.default_rng(5)
    grid = np.array(E1["groups"]["FC1"]["km_grid_s"], float)
    for g in ("FC1", "FC2"):
        mu = o["adhesiveness_shift"][g]
        sims = []
        for r in range(40):
            f, c = D.simulate(E1["groups"][g]["attached"], o["lambda_per_s"], o["sigma"], o["gamma"],
                              o["acquisition_rate_per_s"], mu, 0.0, rng)
            sims.append(A.km(f, c, grid))
        sims = np.array(sims)
        ax.fill_between(grid, np.quantile(sims, .025, 0), np.quantile(sims, .975, 0), color=GROUP_COLORS[g], alpha=.25, lw=0)
        ax.scatter(grid, E1["groups"][g]["km"], color=GROUP_COLORS[g], s=12, zorder=3, label=GROUP_LABELS[g])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.01, 1.05)
    ax.set(xlabel="Residence after attachment (s)", ylabel="Fraction still attached", title="c  Posterior predictive check")
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    ax = axs[3]
    for truth, col in (("offset", "#1a5c86"), ("head", "#b5642a")):
        rr = sorted([p for p in P if p["truth"] == truth], key=lambda p: p["attached_per_group"])
        n = [p["attached_per_group"] for p in rr]
        ax.plot(n, [100 * p["power_correct"] for p in rr], color=col, marker="o", ms=3.5, lw=1.1,
                label=f"{truth} generates, correct family preferred")
        ax.plot(n, [100 * p["power_LR_gt_3"] for p in rr], color=col, marker="s", ms=3.5, lw=1.1, ls="--",
                label=f"{truth} generates, log LR > 3")
    ax.axvline(E1["groups"]["FC1"]["attached"], color="#999", lw=.6, ls=":")
    ax.text(E1["groups"]["FC1"]["attached"] * 1.07, 88, "FC1 in\ndeposit", fontsize=6, color="#777")
    ax.set_xscale("log")
    ax.set_xticks([250, 1000, 3000], ["250", "1000", "3000"])
    ax.minorticks_off()
    ax.set(xlabel="Attached arrivals per group", ylabel="Simulated experiments (%)", ylim=(0, 104),
           title="d  Arrivals needed")
    ax.legend(fontsize=5.6, frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.42))
    fig.savefig(FIG / "fig7_robustness_design.png"); plt.close(fig)


if __name__ == "__main__":
    from run_experiments import plot
    summary = read(R / "tables" / "factorial_summary.csv")
    contrasts = read(R / "tables" / "paired_history_contrasts.csv")
    traces = read(R / "logs" / "ensemble_trajectories.csv")
    sigs = read(R / "tables" / "event_signatures.csv")
    plot(summary, contrasts, traces, sigs, R)
    fig1(); fig4(sigs); fig5()
    if (R / "empirical_v2" / "design_power.json").exists() and (R / "empirical_v2" / "fits_v2.json").exists():
        fig7()
    if (R / "empirical" / "empirical_fits.json").exists():
        try:
            fig6()
        except KeyError as e:
            print("fig6 skipped; empirical fit incomplete:", e)
    print("figures written")
