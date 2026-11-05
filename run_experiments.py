"""Calibrated factorial comparison of six colonization hypotheses.

Every history mechanism uses baseline rates calibrated so its fresh-cell /
fresh-surface occupancy trajectory matches the heterogeneous passive null;
confounder strengths are calibrated to reproduce the zero-gap history effects.
Any remaining difference between hypotheses therefore lies in the transfer
arms, the gap dependence, or event-level observables.
"""
import argparse, csv, json, platform, time
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import (Config, MECHANISMS, ARMS, ACTIVE, LABELS, SHORT, ARM_LABELS, VERSION,
                   simulate, OBS_TIMES)
from calibrate import calibrated_overrides

METRICS = ("established", "final_attached", "integrated_occupancy", "retained_founders",
           "retained_descendants", "retained_residual", "divisions", "departures", "rejected")
PAIRS = (("cell_transfer", "exposed_fresh", "fresh_fresh"),
         ("surface_transfer", "fresh_conditioned", "fresh_fresh"),
         ("both_transfer", "exposed_conditioned", "fresh_fresh"))
GAPS = (0., 2., 8., 24.)
COLORS = {"passive": "#5b6470", "carryover": "#c2a03a", "selection": "#5fa7c9",
          "cell": "#1a5c86", "surface": "#b5642a", "both": "#5c5597"}
STYLE = {"passive": "-", "carryover": "--", "selection": "--", "cell": "-", "surface": "-", "both": "-"}
PRIMARY_SEED, STEP_SEED, EFFECT_SEED, BOOT_SEED, BOOT = 11003, 22007, 44017, 405, 2000


def write_csv(p, rows):
    rows = list(rows)
    if rows:
        with Path(p).open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)


def summarize(r):
    out = []
    for k in METRICS:
        x = np.asarray(r[k], float)
        se = x.std(ddof=1) / np.sqrt(len(x))
        out.append({"metric": k, "mean": float(x.mean()), "mc_se": float(se),
                    "ci_low": float(x.mean() - 1.96 * se), "ci_high": float(x.mean() + 1.96 * se)})
    return out


def signatures(r, cfg):
    """Event-level observables available to lineage tracking."""
    half = cfg.arrivals // 2
    fd = r["founder_departures"]
    return {"attach_early": float(r["attached_early"].mean() / half),
            "attach_late": float(r["attached_late"].mean() / half),
            "founder_quick_fraction": float(r["founder_quick"].sum() / max(fd.sum(), 1)),
            "residual_fraction_final": float(r["retained_residual"].sum() / max(r["final_attached"].sum(), 1)),
            "initial_residual": float(r["initial_residual"].mean())}


def differences(a, b, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    out = []
    for key in ("established", "final_attached", "integrated_occupancy"):
        x = np.asarray(a[key], float) - np.asarray(b[key], float)
        boots = x[rng.integers(len(x), size=(BOOT, len(x)))].mean(axis=1)
        lo, hi = np.quantile(boots, [.025, .975])
        out.append({"metric": key, "difference": float(x.mean()), "ci_low": float(lo), "ci_high": float(hi)})
    return out


def plot(summary, contrasts, traces, sigs, out):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 300, "font.family": "DejaVu Sans"})
    # Figure 2: matched baselines and zero-gap transfer effects
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), layout="constrained")
    for m in MECHANISMS:
        rr = [t for t in traces if t["mechanism"] == m and t["arm"] == "fresh_fresh"]
        axs[0].plot([t["time"] for t in rr], [t["mean_attached"] for t in rr], color=COLORS[m],
                    ls=STYLE[m], lw=1.4, label=LABELS[m])
    axs[0].set(xlabel="Challenge time", ylabel="Mean attached cells", title="a  Fresh cells, fresh surface")
    axs[0].legend(fontsize=7, frameon=False)
    for ax, c, title in ((axs[1], "cell_transfer", "b  Exposed cells, fresh surface"),
                         (axs[2], "surface_transfer", "c  Fresh cells, visited surface")):
        y, lo, hi = [], [], []
        for m in MECHANISMS:
            r = next(x for x in contrasts if x["mechanism"] == m and x["contrast"] == c
                     and x["gap"] == 0 and x["metric"] == "integrated_occupancy")
            y.append(r["difference"]); lo.append(r["difference"] - r["ci_low"]); hi.append(r["ci_high"] - r["difference"])
        ax.bar(range(6), y, yerr=[lo, hi], color=[COLORS[m] for m in MECHANISMS], capsize=2,
               hatch=None)
        ax.set_xticks(range(6), [SHORT[m] for m in MECHANISMS], rotation=35, ha="right")
        ax.axhline(0, color="#999", lw=.6)
        ax.set(ylabel="Change in integrated occupancy", title=title)
    fig.savefig(out / "figures" / "fig2_matched_baseline_transfer.png"); plt.close(fig)

    # Figure 3: gap dependence
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), layout="constrained")
    for ax, c, ms, title in ((axs[0], "cell_transfer", ("selection", "cell", "both"), "a  Exposed-cell effect"),
                             (axs[1], "surface_transfer", ("carryover", "surface", "both"), "b  Visited-surface effect")):
        base = {m: next(x["difference"] for x in contrasts if x["mechanism"] == m and x["contrast"] == c
                        and x["gap"] == 0 and x["metric"] == "integrated_occupancy") for m in ms}
        for m in ms:
            rr = sorted([x for x in contrasts if x["mechanism"] == m and x["contrast"] == c
                         and x["metric"] == "integrated_occupancy"], key=lambda x: x["gap"])
            g = [x["gap"] for x in rr]
            ax.errorbar(g, [x["difference"] / base[m] for x in rr],
                        yerr=[[(x["difference"] - x["ci_low"]) / base[m] for x in rr],
                              [(x["ci_high"] - x["difference"]) / base[m] for x in rr]],
                        color=COLORS[m], ls=STYLE[m], marker="o", ms=4, capsize=2, label=LABELS[m])
        ax.axhline(0, color="#999", lw=.6)
        ax.set(xlabel="Gap between preparation and challenge", ylabel="Effect relative to zero gap", title=title)
        ax.legend(fontsize=7, frameon=False)
    fig.savefig(out / "figures" / "fig3_gap_dependence.png"); plt.close(fig)

    # Figure 4: event-level signatures at zero gap
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), layout="constrained")
    def sig(m, arm, k):
        return next(s[k] for s in sigs if s["mechanism"] == m and s["arm"] == arm and s["gap"] == 0)
    x = np.arange(6)
    ratio = [sig(m, "exposed_fresh", "attach_late") / max(sig(m, "exposed_fresh", "attach_early"), 1e-9) /
             (sig(m, "fresh_fresh", "attach_late") / max(sig(m, "fresh_fresh", "attach_early"), 1e-9)) for m in MECHANISMS]
    axs[0].bar(x, ratio, color=[COLORS[m] for m in MECHANISMS])
    axs[0].axhline(1, color="#999", lw=.6)
    axs[0].set_ylim(min(ratio) * 0.9, max(ratio) * 1.05)
    axs[0].set(ylabel="Late/early attachment, exposed vs fresh", title="a  Decay of exposed advantage")
    axs[1].bar(x, [100 * sig(m, "fresh_conditioned", "residual_fraction_final") for m in MECHANISMS],
               color=[COLORS[m] for m in MECHANISMS])
    axs[1].set(ylabel="Final cells from residual founders (%)", title="b  Lineages without a landing event")
    fq = [100 * (sig(m, "exposed_fresh", "founder_quick_fraction") - sig(m, "fresh_fresh", "founder_quick_fraction"))
          for m in MECHANISMS]
    axs[2].bar(x, fq, color=[COLORS[m] for m in MECHANISMS])
    axs[2].axhline(0, color="#999", lw=.6)
    axs[2].set(ylabel="Change in departures < 1 time unit (pp)", title="c  Early-residence detachment")
    for ax in axs:
        ax.set_xticks(x, [SHORT[m] for m in MECHANISMS], rotation=35, ha="right")
    fig.savefig(out / "figures" / "fig4_event_signatures.png"); plt.close(fig)


def run(output=Path("results"), n=400, check_n=250):
    start = time.perf_counter()
    out = Path(output)
    for d in ("tables", "figures", "logs"):
        (out / d).mkdir(parents=True, exist_ok=True)
    cfg = Config()
    cal = json.loads((out / "logs" / "calibration.json").read_text())
    ov = calibrated_overrides(cal)
    raw, summary, contrasts, traces, sigs = [], [], [], [], []
    runs = 0
    for m in MECHANISMS:
        for gap in GAPS:
            block = {}
            for a in ARMS:
                r = simulate(cfg, m, a, gap, n=n, seed=PRIMARY_SEED, record=(gap == 0), overrides=ov[m])
                block[a] = r; runs += n
                summary.extend({"mechanism": m, "gap": gap, "arm": a, "n": n, **s} for s in summarize(r))
                sigs.append({"mechanism": m, "gap": gap, "arm": a, **signatures(r, cfg)})
                for i in range(n):
                    raw.append({"mechanism": m, "gap": gap, "arm": a, "replicate": i,
                                "establishment_time": float(r["establishment_time"][i]) if r["established"][i] else "",
                                **{k: float(r[k][i]) for k in METRICS},
                                "retained_lineages": int(r["retained_lineages"][i])})
                if gap == 0:
                    traces.extend({"mechanism": m, "arm": a, **t} for t in r["trace"])
            for name, a, b in PAIRS:
                contrasts.extend({"mechanism": m, "gap": gap, "contrast": name, "n_pairs": n, **d}
                                 for d in differences(block[a], block[b]))
        print("factorial", m, flush=True)
    write_csv(out / "tables" / "factorial_trials.csv", raw)
    write_csv(out / "tables" / "factorial_summary.csv", summary)
    write_csv(out / "tables" / "paired_history_contrasts.csv", contrasts)
    write_csv(out / "tables" / "event_signatures.csv", sigs)
    write_csv(out / "logs" / "ensemble_trajectories.csv", traces)

    steps = []
    for dt in (.05, .025, .0125):
        c = replace(cfg, dt=dt)
        for m in MECHANISMS:
            for a in ("fresh_fresh", "exposed_conditioned"):
                r = simulate(c, m, a, n=check_n, seed=STEP_SEED, overrides=ov[m]); runs += check_n
                steps.extend({"dt": dt, "mechanism": m, "arm": a, "n": check_n, **s} for s in summarize(r))
    write_csv(out / "tables" / "time_step_sensitivity.csv", steps)
    print("step sensitivity done", flush=True)

    effects = []
    for factor in (.5, 1., 1.5):
        for m in MECHANISMS:
            # Scale only the couplings that this hypothesis switches on.
            o = dict(ov[m])
            on = ACTIVE[m]
            if "cell" in on:
                for k in ("cell_attachment_effect", "cell_retention_effect"):
                    o[k] = getattr(cfg, k) * factor
            if "surface" in on:
                for k in ("surface_attachment_effect", "surface_retention_effect"):
                    o[k] = getattr(cfg, k) * factor
            if "selection" in on:
                o["selection_shift"] = cal["selection_shift"] * factor
            if "carryover" in on:
                o["carryover_cells"] = cal["carryover_cells"] * factor
            block = {a: simulate(cfg, m, a, n=check_n, seed=EFFECT_SEED, overrides=o) for a in ARMS}
            runs += 4 * check_n
            for name, a, b in PAIRS:
                effects.extend({"effect_multiplier": factor, "mechanism": m, "contrast": name,
                                "n_pairs": check_n, **d} for d in differences(block[a], block[b]))
    write_csv(out / "tables" / "effect_strength_sensitivity.csv", effects)
    plot(summary, contrasts, traces, sigs, out)
    audit = {"version": VERSION, "status": "synthetic; calibrated to a common passive baseline",
             "python": platform.python_version(), "numpy": np.__version__,
             "matplotlib": matplotlib.__version__, "primary_n_per_arm": n,
             "primary_trajectories": len(raw), "sensitivity_n_per_arm": check_n,
             "total_trajectories": runs, "primary_seed": PRIMARY_SEED, "step_seed": STEP_SEED,
             "effect_seed": EFFECT_SEED, "bootstrap_resamples": BOOT, "primary_dt": cfg.dt,
             "max_cell_accounting_error": 0, "elapsed_seconds": time.perf_counter() - start}
    (out / "logs" / "run_audit.json").write_text(json.dumps(audit, indent=2))
    (out / "logs" / "config.json").write_text(json.dumps(asdict(cfg), indent=2))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=400)
    p.add_argument("--check-n", type=int, default=250)
    p.add_argument("--output", type=Path, default=Path("results"))
    a = p.parse_args()
    run(a.output, a.n, a.check_n)
