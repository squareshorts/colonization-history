"""LaTeX tables for the manuscript and supplement, generated from result files."""
import csv, json
from dataclasses import asdict
from pathlib import Path
import numpy as np
from model import Config, MECHANISMS, LABELS
from recovery import PRIORS

R = Path("results")
OUT = Path("manuscript/generated")
OUT.mkdir(parents=True, exist_ok=True)
SYM = {"attachment_probability": "$p_0$", "detachment_rate": "$d_0$", "division_rate": "$b$",
       "patch_heterogeneity": "$\\sigma_u$", "cell_heterogeneity": "$\\sigma_z$",
       "memory_acquisition_rate": "$k_m$", "cell_memory_decay_time": "$\\tau_m$",
       "cell_attachment_effect": "$\\beta_m$", "cell_retention_effect": "$\\gamma_m$",
       "surface_decay_time": "$\\tau_s$", "surface_attachment_effect": "$\\beta_s$",
       "surface_retention_effect": "$\\gamma_s$", "contact_trace": "$q$",
       "resident_deposition_rate": "$k_s$", "carryover_cells": "$\\lambda_c$",
       "storage_loss_rate": "$\\ell$", "selection_shift": "$\\delta$",
       "selection_relaxation_time": "$\\tau_z$"}
KIND = {"u": "uniform", "log": "log-uniform", "logit": "logit-uniform"}


def f(x, d=3):
    return f"{x:.{d}g}"


def priors_table():
    cfg = asdict(Config())
    rows = [f"{SYM[k]} & {k.replace('_', ' ')} & {f(cfg[k])} & {KIND[kind]} [{f(lo)}, {f(hi)}] \\\\"
            for k, (kind, lo, hi) in PRIORS.items()]
    t = ("\\begin{tabular}{@{}llll@{}}\n\\toprule\nSymbol & Parameter & Default & Recovery prior \\\\\n\\midrule\n"
         + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_priors.tex").write_text(t)


def calibration_table():
    c = json.loads((R / "logs" / "calibration.json").read_text())
    rows = []
    for m in ("cell", "surface", "both"):
        b = c["baselines"][m]
        rows.append(f"{LABELS[m]} & {b['attachment_probability']:.3f} & {b['detachment_rate']:.3f} & "
                    f"{b['division_rate']:.3f} & {b['rms_log_error']:.3f} \\\\")
    cfg = Config()
    rows.insert(0, f"Heterogeneous passive (target) & {cfg.attachment_probability:.3f} & {cfg.detachment_rate:.3f} & "
                   f"{cfg.division_rate:.3f} & -- \\\\")
    t = ("\\begin{tabular}{@{}lrrrr@{}}\n\\toprule\nHypothesis & $p_0$ & $d_0$ & $b$ & RMS error \\\\\n\\midrule\n"
         + "\n".join(rows) + "\n\\midrule\n"
         + f"\\multicolumn{{5}}{{@{{}}l}}{{Carryover $\\lambda_c={c['carryover_cells']:.2f}$; "
           f"selection $\\delta={c['selection_shift']:.3f}$}} \\\\\n\\bottomrule\n\\end{{tabular}}\n")
    (OUT / "table_calibration.tex").write_text(t)


def dt_table():
    rows = list(csv.DictReader(open(R / "tables" / "time_step_sensitivity.csv")))
    out = []
    for m in MECHANISMS:
        cells = []
        for dt in ("0.05", "0.025", "0.0125"):
            r = next(x for x in rows if x["mechanism"] == m and x["arm"] == "exposed_conditioned"
                     and x["dt"] == dt and x["metric"] == "integrated_occupancy")
            cells.append(f"{float(r['mean']):.1f} $\\pm$ {1.96*float(r['mc_se']):.1f}")
        out.append(f"{LABELS[m]} & " + " & ".join(cells) + " \\\\")
    t = ("\\begin{tabular}{@{}lrrr@{}}\n\\toprule\nHypothesis & $\\Delta t=0.05$ & $0.025$ & $0.0125$ \\\\\n\\midrule\n"
         + "\n".join(out) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_dt.tex").write_text(t)


def effect_table():
    rows = list(csv.DictReader(open(R / "tables" / "effect_strength_sensitivity.csv")))
    out = []
    for m in MECHANISMS[1:]:
        c = "surface_transfer" if m in ("carryover", "surface") else "cell_transfer"
        cells = []
        for fac in ("0.5", "1.0", "1.5"):
            r = next(x for x in rows if x["mechanism"] == m and x["contrast"] == c and x["effect_multiplier"] == fac
                     and x["metric"] == "integrated_occupancy")
            cells.append(f"{float(r['difference']):.1f} [{float(r['ci_low']):.1f}, {float(r['ci_high']):.1f}]")
        out.append(f"{LABELS[m]} & {'FV' if c == 'surface_transfer' else 'EF'} $-$ FF & " + " & ".join(cells) + " \\\\")
    t = ("\\begin{tabular}{@{}llrrr@{}}\n\\toprule\nHypothesis & Contrast & $\\times0.5$ & $\\times1$ & $\\times1.5$ \\\\\n\\midrule\n"
         + "\n".join(out) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_effect.tex").write_text(t)


def recovery_table():
    rows = list(csv.DictReader(open(R / "tables" / "recovery_accuracy.csv")))
    out = []
    names = {"full": "full $2\\times2$", "combined": "FF + EV", "exposed": "FF + EF", "baseline": "FF only"}
    for tier in ("counts", "tracking"):
        for arms in ("full", "combined", "exposed", "baseline"):
            for gaps in (("0", "0-16", "0-4-16") if arms != "baseline" else ("none",)):
                cells = []
                for n in ("2", "4", "8"):
                    r = next(x for x in rows if x["tier"] == tier and x["arms"] == arms and x["gaps"] == gaps
                             and x["replicates"] == n)
                    cells.append(f"{100*float(r['accuracy']):.1f}")
                g = gaps.replace("-", ", ") if gaps != "none" else "--"
                out.append(f"{tier} & {names[arms]} & {g} & " + " & ".join(cells) + " \\\\")
    t = ("\\begin{tabular}{@{}lllrrr@{}}\n\\toprule\nTier & Arms & Gaps & $n=2$ & $n=4$ & $n=8$ \\\\\n\\midrule\n"
         + "\n".join(out) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_recovery.tex").write_text(t)


def empirical_tables():
    p = R / "empirical" / "empirical_fits.json"
    if not p.exists():
        return
    E = json.loads(p.read_text())
    if "families" not in E:
        return
    val = {r["dataset"]: r for r in csv.DictReader(open(R / "empirical" / "occupancy_validation.csv"))}
    out = []
    for g, v in E["groups"].items():
        for n in v["datasets"]:
            r = val[n]
            out.append(f"{g} & {n.replace('WT YFP-CFP ', '').replace('_', ' ')} & {float(r['record_h']):.1f} & "
                       f"{float(r['arrivals_per_h_first3h']):.0f} & {float(r['initial_cells']):.1f} & "
                       f"{r['lag_h'] if r['lag_h'] else '--'} \\\\")
    t = ("\\begin{tabular}{@{}llrrrr@{}}\n\\toprule\nGroup & Dataset & Record (h) & Arrivals h$^{-1}$ & "
         "Initial cells & Lag (h) \\\\\n\\midrule\n" + "\n".join(out) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_datasets.tex").write_text(t)
    out = []
    for g, v in E["groups"].items():
        km = dict(zip(v["km_grid_s"], v["km"]))
        out.append(f"{g} & {v['arrivals']} & {v['attached']} & {100*v['attachment_fraction']:.1f} & "
                   f"{100*km[60]:.1f} & {100*km[300]:.1f} & {100*km[1800]:.1f} \\\\")
    t = ("\\begin{tabular}{@{}lrrrrrr@{}}\n\\toprule\nGroup & Arrivals & Attached & Attached (\\%) & "
         "$S(60\\,\\mathrm{s})$ & $S(300\\,\\mathrm{s})$ & $S(1800\\,\\mathrm{s})$ \\\\\n\\midrule\n"
         + "\n".join(out) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_residence.tex").write_text(t)
    fam_names = {"null": "no exposure effect", "selection_no_acquisition": "offset, no acquisition",
                 "selection": "adhesiveness offset", "history": "head start on acquisition",
                 "selection_and_history": "offset and head start"}
    best = min(v["aic"] for v in E["families"].values())
    out = [f"{fam_names[k]} & {v['k']} & {v['loglik']:.1f} & {v['aic'] - best:.1f} \\\\" for k, v in E["families"].items()]
    t = ("\\begin{tabular}{@{}lrrr@{}}\n\\toprule\nFamily & Parameters & $\\log L$ & $\\Delta$AIC \\\\\n\\midrule\n"
         + "\n".join(out) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_families.tex").write_text(t)
    ms, mh, sat = E["mixture_selection_fit"], E["mixture_history_fit"], E["mixture_saturated"]
    best = E["families"]["selection"]
    lines = [("two-type reweighting", ms), ("single-type head start", mh),
             ("two-type mixture per group", sat), ("lognormal adhesiveness offset", best)]
    body = "\n".join(f"{n} & {v['k']} & {v['loglik']:.1f} & {v['aic']:.1f} \\\\" for n, v in lines)
    t = ("\\begin{tabular}{@{}lrrr@{}}\n\\toprule\nFamily & Parameters & $\\log L$ & AIC \\\\\n\\midrule\n"
         + body + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_mixture.tex").write_text(t)
    out = []
    for k, v in E["mutants_fc1"].items():
        km = dict(zip(E["groups"]["FC1"]["km_grid_s"], v["km"]))
        out.append(f"{k} & {v['arrivals']} & {v['attached']} & {100*v['attachment_fraction']:.1f} & "
                   f"{100*km[60]:.1f} & {100*km[300]:.1f} \\\\")
    t = ("\\begin{tabular}{@{}lrrrrr@{}}\n\\toprule\nStrain & Arrivals & Attached & Attached (\\%) & "
         "$S(60\\,\\mathrm{s})$ & $S(300\\,\\mathrm{s})$ \\\\\n\\midrule\n" + "\n".join(out)
         + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_mutants.tex").write_text(t)




def v2_tables():
    p = R / "empirical_v2" / "fits_v2.json"
    if not p.exists():
        return
    E2 = json.loads(p.read_text())
    if "time_matched" not in E2:
        return
    E1 = json.loads((R / "empirical" / "empirical_fits.json").read_text())
    def q(x):
        return f"[{x[0]:.1f}, {x[2]:.1f}]" if x else "--"
    def fc2ci(r):
        c = r.get("offset_ci", {}).get("FC2")
        return f" [{c[0]:.2f}, {c[1]:.2f}]" if c else ""
    rows = [("First 3 h, 1-h windows (v1)", sum(v["attached"] for v in E1["groups"].values()),
             -E1["delta_loglik_history_minus_selection"], None,
             f"{E1['families']['selection']['adhesiveness_shift']['FC2']:.2f} "
             f"[{E1['adhesiveness_shift_ci']['FC2'][0]:.2f}, {E1['adhesiveness_shift_ci']['FC2'][1]:.2f}]")]
    for name, key in (("First 3 h, continuous tracking", "P1_first3h"), ("All arrivals before crowding", "P1_extended"),
                      ("All arrivals, arrival-time effect", "time_matched"), ("All arrivals, occupancy covariate", "occupancy_adjusted")):
        r = E2[key]
        n = sum((r.get("n_attached") or r.get("n")).values())
        rows.append((name, n, r["dll_offset_minus_head"], r.get("bootstrap_dll_q"),
                     f"{r['fits']['offset']['offset']['FC2']:.2f}" + fc2ci(r)))
    body = "\n".join(f"{a} & {n} & {d:.1f} & {q(b)} & {c} \\\\" for a, n, d, b, c in rows)
    t = ("\\begin{tabular}{@{}lrrrr@{}}\n\\toprule\nAnalysis & Attached arrivals & $\\Delta\\log L$ & 95\\% interval & FC2 offset \\\\\n"
         "\\midrule\n" + body + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_robust.tex").write_text(t)
    tm = E2["time_matched"]["fits"]["offset"]
    oc = E2["occupancy_adjusted"]["fits"]["offset"]
    labels = ["1--3 h", "3--8 h", "8--20 h", "$>$20 h"]
    body = "\n".join(f"arrival {l} & {v:.2f} \\\\" for l, v in zip(labels, tm["time_effect"]))
    body += f"\nlog occupancy at arrival (separate fit) & {oc['time_effect'][0]:.3f} \\\\"
    body += "\n\\midrule\n" + "\n".join(f"offset {g} (arrival-time fit) & {v:.2f} \\\\" for g, v in tm["offset"].items() if g != "FC1")
    t = ("\\begin{tabular}{@{}lr@{}}\n\\toprule\nTerm (log-hazard reduction relative to reference) & Estimate \\\\\n\\midrule\n"
         + body + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_time.tex").write_text(t)
    pp = E2["P2_P3"]
    t = ("\\begin{tabular}{@{}lr@{}}\n\\toprule\nQuantity & Value \\\\\n\\midrule\n"
         + "\n".join(f"cAMP readings, {g} & {n} \\\\" for g, n in pp["n_read"].items())
         + f"\n$\\beta$ (per s.d.\\ of log YFP/CFP) & {pp['offset_with_cAMP']['beta']:.3f} \\\\"
         + f"\nLRT for $\\beta$ & $\\chi^2_1={pp['lrt_beta']['stat']:.2f}$, $p={pp['lrt_beta']['p']:.2f}$ \\\\"
         + f"\nFC2 offset without / with cAMP & {pp['offset_without_cAMP']['offset']['FC2']:.2f} / {pp['offset_with_cAMP']['offset']['FC2']:.2f} \\\\"
         + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_camp.tex").write_text(t)
    P = json.loads((R / "empirical_v2" / "design_power.json").read_text())
    body = "\n".join(f"{x['truth']} & {x['attached_per_group']} & {100*x['power_correct']:.0f} & {100*x['power_LR_gt_3']:.0f} \\\\" for x in P)
    t = ("\\begin{tabular}{@{}lrrr@{}}\n\\toprule\nGenerating family & Attached per group & Correct (\\%) & Correct, $\\Delta\\log L>3$ (\\%) \\\\\n"
         "\\midrule\n" + body + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_power.tex").write_text(t)


def channel_table():
    import glob, os
    rows = []
    for d in sorted(glob.glob("empirical/observables/*/tracks_full.csv")):
        r = [x for x in csv.DictReader(open(d)) if x["yfp"] != ""]
        if not r:
            continue
        y = np.array([float(x["yfp"]) for x in r]); c = np.array([float(x["cfp"]) for x in r])
        name = os.path.basename(os.path.dirname(d)).replace("WT YFP-CFP ", "")
        rows.append(f"{name} & {len(r)} & {np.median(y):.1f} & {np.median(c):.1f} & {100*np.mean(c <= 1):.0f} \\\\")
    t = ("\\begin{tabular}{@{}lrrrr@{}}\n\\toprule\nDataset & Readings & Median YFP & Median CFP & CFP $\\le 1$ (\\%) \\\\\n"
         "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
    (OUT / "table_channels.tex").write_text(t)


if __name__ == "__main__":
    priors_table(); calibration_table(); dt_table(); effect_table(); recovery_table(); empirical_tables()
    v2_tables(); channel_table()
    print(sorted(p.name for p in OUT.iterdir()))
