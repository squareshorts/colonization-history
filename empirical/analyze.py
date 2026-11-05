"""Residence-time and occupancy analysis of the Lee et al. (2018) transfer experiments.

Groups follow the published design (P. aeruginosa PA14, tandem flow cells):
  FC1   surface-naive cells on a virgin surface              (fresh / fresh)
  FC2   cells harvested from FC1 into a virgin flow cell     (exposed / fresh, gap 0)
  FC3   cells harvested from FC2 into a third flow cell      (twice exposed / fresh)
  Tr14, Tr37  FC1 cells held 14 h or 37 h in culture tubes, regrown ~3 h,
        then introduced to a virgin flow cell              (exposed / fresh, gap > 0)

A track counts as an attached arrival if it is visible in at least two
consecutive 3-s frames. Its residence is interval-censored between the last
frame it was seen and the next frame, or right-censored at the end of the
tracking window. Two hypothesis families with equal parameter counts are fitted
jointly across groups:

  selection (frailty): S_g(a) = w_g e^{-l1 a} + (1 - w_g) e^{-l2 a}
      shared component hazards, group-specific mixing weight.
  history (acquisition): h(a) = l_inf + (l_0 - l_inf) e^{-k a}; group g enters
      at an equivalent residence age a0_g (a0 = 0 for FC1), so exposed cells
      follow the fresh-cell hazard curve shifted in age.

The same likelihood is used for both; models are compared by log likelihood
(equal parameter counts) and AIC.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
OBS = ROOT / "observables"
FRAME_S = 3.0
GROUPS = {
    "FC1": ["WT YFP-CFP (T=29.5h;FC1)", "WT YFP-CFP (T=26.75h;FC1)",
            "WT YFP-CFP (T1=32.5h,T2=14.75h;FC1)", "WT YFP-CFP (t1=0-1h;FC1)"],
    "FC2": ["WT YFP-CFP (T=29.5h;FC2)", "WT YFP-CFP (T=26.75h;FC2)",
            "WT YFP-CFP (T1=32.5h,T2=14.75h;FC2)"],
    "FC3": ["WT YFP-CFP (T1=32.5h,T2=14.75h;FC3)"],
    "Tr14": ["WT YFP-CFP (Tr=14h)"],
    "Tr37": ["WT YFP-CFP (Tr=37h)"],
}
MUTANTS = {"WT (no reporter)": ["WT"], "pilA": ["pilA"], "pilJ": ["pilJ"],
           "pilT": ["pilT"], "pilU": ["pilU"], "pilU YFP-CFP": ["pilU YFP-CFP"]}
GAP_H = {"FC2": 0.0, "Tr14": 14.0, "Tr37": 37.0}
WINDOWS = (0, 1, 2)          # first three tracking hours after inoculation
MIN_FRAMES = 2


def load_tracks(dataset, windows=WINDOWS):
    rows = list(csv.DictReader((OBS / dataset / "tracks.csv").open()))
    out = []
    for r in rows:
        if int(r["window"]) not in windows or r["origin"] != "arrival":
            continue
        d = int(r["end"]) - int(r["start"]) + 1
        out.append((d, r["right_censored"] == "True"))
    return np.array([o[0] for o in out], int), np.array([o[1] for o in out], bool)


def group_data(names, windows=WINDOWS):
    D, C, arrivals = [], [], 0
    for n in names:
        if not (OBS / n / "tracks.csv").exists():
            continue
        d, c = load_tracks(n, windows)
        arrivals += len(d)
        D.append(d); C.append(c)
    if not D:
        return None
    D, C = np.concatenate(D), np.concatenate(C)
    keep = D >= MIN_FRAMES
    return {"frames": D[keep], "censored": C[keep], "arrivals": int(arrivals),
            "attached": int(keep.sum()), "datasets": [n for n in names if (OBS / n / "tracks.csv").exists()]}


# ---- survival families -------------------------------------------------------
def S_mix(a, w, l1, l2):
    return w * np.exp(-l1 * a) + (1 - w) * np.exp(-l2 * a)


def H_acq(a, l0, linf, k):
    return linf * a + (l0 - linf) * (1 - np.exp(-k * a)) / k


def loglik_group(S, frames, cens):
    """Conditional on survival to MIN_FRAMES frames (attachment)."""
    a_last = frames * FRAME_S            # last seen at age (frames-1)*3; gone by frames*3
    a_seen = (frames - 1) * FRAME_S
    a0 = (MIN_FRAMES - 1) * FRAME_S
    s0 = S(a0)
    ev = ~cens
    ll = np.sum(np.log(np.clip(S(a_seen[cens]) / s0, 1e-300, None)))
    ll += np.sum(np.log(np.clip((S(a_seen[ev]) - S(a_last[ev])) / s0, 1e-300, None)))
    return ll


def fit_selection(data, order):
    def unpack(x):
        l1, l2 = np.exp(x[0]), np.exp(x[0]) * np.exp(-np.exp(x[1]))      # l2 < l1
        w = 1 / (1 + np.exp(-x[2:]))
        return l1, l2, dict(zip(order, w))
    def nll(x):
        l1, l2, w = unpack(x)
        return -sum(loglik_group(lambda a, g=g: S_mix(a, w[g], l1, l2), *data[g]) for g in order)
    best = None
    for s in ([np.log(0.1), 2.0], [np.log(0.03), 1.0], [np.log(0.3), 3.0]):
        x0 = np.array(s + [0.0] * len(order))
        r = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-6, "fatol": 1e-8})
        r = minimize(nll, r.x, method="BFGS")
        if best is None or r.fun < best.fun:
            best = r
    l1, l2, w = unpack(best.x)
    return {"loglik": -float(best.fun), "k": len(best.x), "l1_per_s": float(l1), "l2_per_s": float(l2),
            "weights_fast": {g: float(v) for g, v in w.items()}}, best.x


def fit_history(data, order):
    exposed = [g for g in order if g != "FC1"]
    def unpack(x):
        l0, linf, k = np.exp(x[0]), np.exp(x[1]), np.exp(x[2])
        a0 = {"FC1": 0.0, **{g: float(np.exp(v)) for g, v in zip(exposed, x[3:])}}
        return l0, linf, k, a0
    def nll(x):
        l0, linf, k, a0 = unpack(x)
        tot = 0.0
        for g in order:
            o = a0[g]
            S = lambda a, o=o: np.exp(-(H_acq(a + o, l0, linf, k) - H_acq(o, l0, linf, k)))
            tot += loglik_group(S, *data[g])
        return -tot
    best = None
    for s in ([np.log(0.2), np.log(0.0005), np.log(0.02)], [np.log(0.05), np.log(0.001), np.log(0.005)],
              [np.log(0.5), np.log(0.0002), np.log(0.1)]):
        x0 = np.array(s + [np.log(60.0)] * len(exposed))
        r = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-6, "fatol": 1e-8})
        r = minimize(nll, r.x, method="BFGS")
        if best is None or r.fun < best.fun:
            best = r
    l0, linf, k, a0 = unpack(best.x)
    m0 = {g: float(1 - np.exp(-k * v)) for g, v in a0.items()}
    return {"loglik": -float(best.fun), "k": len(best.x), "l0_per_s": float(l0), "linf_per_s": float(linf),
            "acquisition_rate_per_s": float(k), "age_offset_s": a0, "equivalent_memory": m0}, best.x


def fit_saturated(data, order):
    """Separate two-component mixtures per group: goodness-of-fit reference."""
    tot, k = 0.0, 0
    for g in order:
        r, _ = fit_selection({g: data[g]}, [g])
        tot += r["loglik"]; k += r["k"]
    return {"loglik": tot, "k": k}


# ---- model-faithful families: lognormal adhesiveness heterogeneity + acquisition
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(24)
GH_W = GH_W / GH_W.sum()
from scipy.special import expi


def cum_hazard(a, lam, gam, k, m0):
    """lam * integral_0^a exp(-gam m(s)) ds with m(s) = 1 - (1 - m0) exp(-k s)."""
    c = gam * (1 - m0)
    if c < 1e-8:
        return lam * np.exp(-gam) * a
    return lam * np.exp(-gam) / k * (expi(c) - expi(c * np.exp(-k * a)))


def S_frailty(a, lam, sig, gam, k, mu, m0):
    L = cum_hazard(np.asarray(a, float), lam, gam, k, m0)
    rel = np.exp(-(mu + sig * GH_X))                     # exp(-z)
    return np.exp(-np.outer(L, rel)) @ GH_W


FAMILIES = {
    # name: (acquisition on, group shift in adhesiveness, group head start in history)
    "null": (True, False, False),
    "selection_no_acquisition": (False, True, False),
    "selection": (True, True, False),
    "history": (True, False, True),
    "selection_and_history": (True, True, True),
}


def fit_family(data, order, family, starts=4, seed=1):
    acq, shift, head = FAMILIES[family]
    ex = [g for g in order if g != "FC1"]
    nb = 4 if acq else 2
    n = nb + (len(ex) if shift else 0) + (len(ex) if head else 0)
    def unpack(x):
        lam, sig = np.exp(x[0]), np.exp(x[1])
        gam, k = (np.exp(x[2]), np.exp(x[3])) if acq else (0.0, 1.0)
        i = nb
        mu = {"FC1": 0.0}; m0 = {"FC1": 0.0}
        for g in ex:
            mu[g] = 0.0; m0[g] = 0.0
        if shift:
            for g in ex:
                mu[g] = float(x[i]); i += 1
        if head:
            for g in ex:
                m0[g] = float(1 / (1 + np.exp(-x[i]))); i += 1
        return lam, sig, gam, k, mu, m0
    def nll(x):
        lam, sig, gam, k, mu, m0 = unpack(x)
        if not (1e-6 < sig < 20 and gam < 50 and 1e-6 < k < 10):
            return 1e12
        tot = 0.0
        for g in order:
            S = lambda a, g=g: S_frailty(a, lam, sig, gam, k, mu[g], m0[g])
            tot += loglik_group(S, *data[g])
        return -tot if np.isfinite(tot) else 1e12
    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        x0 = np.zeros(n)
        x0[0] = np.log(0.02) + rng.normal(0, .5); x0[1] = np.log(1.5) + rng.normal(0, .3)
        if acq:
            x0[2] = np.log(2.0) + rng.normal(0, .5); x0[3] = np.log(0.02) + rng.normal(0, .5)
        x0[nb:] = rng.normal(0, .3, n - nb)
        r = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 6000, "xatol": 1e-5, "fatol": 1e-7})
        r = minimize(nll, r.x, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    lam, sig, gam, k, mu, m0 = unpack(best.x)
    return {"family": family, "loglik": -float(best.fun), "k": n, "aic": 2 * n + 2 * float(best.fun),
            "lambda_per_s": float(lam), "sigma": float(sig), "gamma": float(gam),
            "acquisition_rate_per_s": float(k) if acq else None,
            "adhesiveness_shift": mu if shift else None, "history_head_start": m0 if head else None}, best.x


def km(frames, cens, grid_s):
    """Kaplan-Meier survival of attached arrivals at the given ages."""
    t_event = frames * FRAME_S
    order = np.argsort(t_event)
    t, c = t_event[order], cens[order]
    at_risk = len(t)
    S, out, j = 1.0, [], 0
    times = np.unique(t)
    surv = {}
    for u in times:
        d = np.sum((t == u) & ~c)
        n = np.sum(t >= u)
        if n > 0:
            S *= 1 - d / n
        surv[u] = S
    keys = np.array(sorted(surv))
    vals = np.array([surv[k] for k in keys])
    return [float(vals[keys <= g][-1]) if np.any(keys <= g) else 1.0 for g in grid_s]


def occupancy(dataset):
    rows = list(csv.DictReader((OBS / dataset / "counts.csv").open()))
    meta = json.loads((OBS / dataset / "meta.json").read_text())
    t = np.array([float(r["time_h"]) for r in rows])
    area = np.array([float(r["area_px"]) for r in rows])
    return t, area / meta["median_component_area_first2h"], meta


def arrival_rate(dataset, windows=WINDOWS):
    d, _ = load_tracks(dataset, windows)
    return len(d) / len(windows)


def bootstrap(data, order, B=200, seed=8):
    rng = np.random.default_rng(seed)
    rows = []
    for b in range(B):
        bd = {}
        for g in order:
            f, c = data[g]
            i = rng.integers(len(f), size=len(f))
            bd[g] = (f[i], c[i])
        s, _ = fit_family(bd, order, "selection", starts=2, seed=b)
        h, _ = fit_family(bd, order, "history", starts=2, seed=b)
        rows.append({"dll": h["loglik"] - s["loglik"], "m0": h["history_head_start"],
                     "mu": s["adhesiveness_shift"]})
    return rows


def decay_time(values):
    """Exponential relaxation time of the exposed-cell effect from gaps 0, 14, 37 h."""
    g = np.array([GAP_H[k] for k in ("FC2", "Tr14", "Tr37")])
    y = np.clip(np.array([values[k] for k in ("FC2", "Tr14", "Tr37")]), 1e-9, None)
    slope, _ = np.polyfit(g, np.log(y), 1)
    return float(-1 / slope) if slope < 0 else float("inf")


def run(B=60, out=None):
    order = [g for g in GROUPS if group_data(GROUPS[g]) is not None]
    gd = {g: group_data(GROUPS[g]) for g in order}
    data = {g: (gd[g]["frames"], gd[g]["censored"]) for g in order}
    fits = {f: fit_family(data, order, f, starts=6)[0] for f in FAMILIES}
    mix_sel, _ = fit_selection(data, order)
    mix_his, _ = fit_history(data, order)
    sat = fit_saturated(data, order)
    for r in (mix_sel, mix_his, sat):
        r["aic"] = 2 * r["k"] - 2 * r["loglik"]
    loo = {}
    for g in ("FC1", "FC2"):
        for n in gd[g]["datasets"]:
            keep = [x for x in gd[g]["datasets"] if x != n]
            if not keep:
                continue
            dd = dict(data)
            sub = group_data(keep)
            dd[g] = (sub["frames"], sub["censored"])
            a = fit_family(dd, order, "selection", starts=3)[0]
            b = fit_family(dd, order, "history", starts=3)[0]
            loo[n] = b["loglik"] - a["loglik"]
    grid = [6, 15, 30, 60, 120, 300, 600, 1200, 1800, 3000]
    summary = {}
    for g in order:
        f, c = data[g]
        summary[g] = {"datasets": gd[g]["datasets"], "arrivals": gd[g]["arrivals"],
                      "attached": gd[g]["attached"],
                      "attachment_fraction": gd[g]["attached"] / max(gd[g]["arrivals"], 1),
                      "km_grid_s": grid, "km": km(f, c, grid),
                      "retained_5min_given_arrival": gd[g]["attached"] / gd[g]["arrivals"] * km(f, c, [300])[0],
                      "arrivals_per_hour": [arrival_rate(n) for n in gd[g]["datasets"]]}
    boots = bootstrap(data, order, B=B)
    dll = np.array([b["dll"] for b in boots])
    sel, his, both = fits["selection"], fits["history"], fits["selection_and_history"]
    from scipy.stats import chi2
    dfx = len(order) - 1
    res = {"groups": summary, "families": fits,
           "delta_loglik_history_minus_selection": his["loglik"] - sel["loglik"],
           "lrt_history_beyond_selection": {"stat": 2 * (both["loglik"] - sel["loglik"]), "df": dfx,
                                            "p": float(chi2.sf(2 * (both["loglik"] - sel["loglik"]), dfx))},
           "lrt_selection_beyond_history": {"stat": 2 * (both["loglik"] - his["loglik"]), "df": dfx,
                                            "p": float(chi2.sf(2 * (both["loglik"] - his["loglik"]), dfx))},
           "leave_one_dataset_out_delta_loglik": loo,
           "mixture_selection_fit": mix_sel, "mixture_history_fit": mix_his, "mixture_saturated": sat,
           "bootstrap_B": B,
           "bootstrap_delta_loglik_quantiles": np.quantile(dll, [0.025, 0.5, 0.975]).tolist(),
           "bootstrap_fraction_history_better": float(np.mean(dll > 0))}
    m0s = {g: np.array([b["m0"][g] for b in boots]) for g in his["history_head_start"]}
    res["history_head_start_ci"] = {g: np.quantile(v, [0.025, 0.975]).tolist() for g, v in m0s.items()}
    res["decay_time_h_history"] = decay_time(his["history_head_start"])
    taus = [decay_time(b["m0"]) for b in boots]
    res["decay_time_h_history_ci"] = np.quantile(np.clip(taus, 0, 1e4), [0.025, 0.975]).tolist()
    mus = {g: np.array([b["mu"][g] for b in boots]) for g in sel["adhesiveness_shift"]}
    res["adhesiveness_shift_ci"] = {g: np.quantile(v, [0.025, 0.975]).tolist() for g, v in mus.items()}
    res["decay_time_h_selection"] = decay_time(sel["adhesiveness_shift"])
    taus = [decay_time(b["mu"]) for b in boots]
    res["decay_time_h_selection_ci"] = np.quantile(np.clip(taus, 0, 1e4), [0.025, 0.975]).tolist()
    null = fits["null"]
    res["lrt_selection_vs_null"] = {"stat": 2 * (sel["loglik"] - null["loglik"]), "df": dfx,
                                    "p": float(chi2.sf(2 * (sel["loglik"] - null["loglik"]), dfx))}
    res["lrt_history_vs_null"] = {"stat": 2 * (his["loglik"] - null["loglik"]), "df": dfx,
                                  "p": float(chi2.sf(2 * (his["loglik"] - null["loglik"]), dfx))}
    mut = {}
    for name, ds in MUTANTS.items():
        g = group_data(ds)
        if g:
            mut[name] = {"arrivals": g["arrivals"], "attached": g["attached"],
                         "attachment_fraction": g["attached"] / max(g["arrivals"], 1),
                         "km": km(g["frames"], g["censored"], grid)}
    res["mutants_fc1"] = mut
    occ = {}
    for g in order:
        for n in gd[g]["datasets"]:
            t, cells, meta = occupancy(n)
            occ[n] = {"group": g, "hours": meta["hours"], "time_h": t[::20].tolist(), "cells": cells[::20].tolist()}
    for name, ds in MUTANTS.items():
        for n in ds:
            if (OBS / n / "counts.csv").exists():
                t, cells, meta = occupancy(n)
                occ[n] = {"group": name, "hours": meta["hours"], "time_h": t[::20].tolist(), "cells": cells[::20].tolist()}
    res["occupancy"] = occ
    out = Path(out) if out else ROOT.parent / "results" / "empirical"
    out.mkdir(parents=True, exist_ok=True)
    (out / "empirical_fits.json").write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    r = run(B, sys.argv[2] if len(sys.argv) > 2 else None)
    print(json.dumps({k: v for k, v in r.items() if k not in ("occupancy",)}, indent=1)[:6000])
