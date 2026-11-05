"""Analysis plan v2: full-record residence fits and cAMP-linked tests (see ANALYSIS_PLAN_v2.md)."""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, kruskal, mannwhitneyu, spearmanr
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as A

OBS = A.OBS
FRAME_S = 3.0
GROUPS = A.GROUPS
CROWD_CELLS = 100.0


def crowd_frame(ds):
    t, cells, meta = A.occupancy(ds)
    k = 20
    cs = np.convolve(cells, np.ones(k) / k, mode="same")
    i = np.flatnonzero(cs > CROWD_CELLS)
    return int(t[i[0]] * 3600 / FRAME_S) if len(i) else 10 ** 9


def load(ds, window="extended"):
    p = OBS / ds / "tracks_full.csv"
    meta = json.loads((OBS / ds / "meta_full.json").read_text())
    rows = list(csv.DictReader(p.open()))
    limit = 3 * 1200 if window == "first3h" else crowd_frame(ds)
    out = []
    for r in rows:
        if r["origin"] != "arrival":
            continue
        s, e = int(r["start"]), int(r["end"])
        if s >= limit:
            continue
        d = e - s + 1
        if d < A.MIN_FRAMES:
            continue
        rec = {"frames": d, "cens": r["right_censored"] == "True", "start": s, "track": int(r["track"])}
        if r["yfp"] != "" and meta["fluorescence_used"]:
            y, c = float(r["yfp"]), float(r["cfp"])
            rec.update({"read_age_s": (int(r["read_frame"]) - s) * FRAME_S,
                        "log_ratio": float(np.log(max(y, 1.0) / max(c, 1.0)))})
        out.append(rec)
    return out, rows, meta


def group(ds_list, window):
    recs = []
    for ds in ds_list:
        if (OBS / ds / "tracks_full.csv").exists():
            r, _, _ = load(ds, window)
            for x in r:
                x["dataset"] = ds
            recs += r
    return recs


# ---- likelihood on compressed, weighted rows ---------------------------------------
def arrays(recs, truncate=False, covariate=False):
    """Collapse records to unique rows (exact) with integer weights."""
    f = np.array([r["frames"] for r in recs], float)
    c = np.array([r["cens"] for r in recs], float)
    a0 = (np.array([r["read_age_s"] for r in recs]) if truncate
          else np.full(len(recs), (A.MIN_FRAMES - 1) * FRAME_S))
    x = np.array([r.get("log_ratio", 0.0) for r in recs]) if covariate else np.zeros(len(recs))
    key = np.column_stack([f, c, a0, x])
    u, inv, w = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return {"f": u[:, 0], "c": u[:, 1].astype(bool), "a0": u[:, 2], "x": u[:, 3], "w": w.astype(float),
            "inv": inv.ravel(), "n": len(recs)}


def loglik(D, lam, sig, gam, k, mu, m0, beta=0.0, xmean=0.0, xsd=1.0, w=None):
    w = D["w"] if w is None else w
    a_seen, a_last = (D["f"] - 1) * FRAME_S, D["f"] * FRAME_S
    a0 = np.minimum(D["a0"], a_seen)
    xs = (D["x"] - xmean) / xsd if beta else 0.0
    rel = np.exp(-(mu + beta * xs))
    rel = (np.atleast_1d(rel)[:, None] if np.ndim(rel) else np.full((1, 1), rel)) * np.exp(-sig * A.GH_X)[None, :]
    def S(a):
        L = A.cum_hazard(a, lam, gam, k, m0)
        return np.exp(-L[:, None] * rel) @ A.GH_W
    s0, ss = S(a0), S(a_seen)
    ll = np.where(D["c"], np.log(np.clip(ss / s0, 1e-300, None)),
                  np.log(np.clip((ss - S(a_last)) / s0, 1e-300, None)))
    return float(np.dot(w, ll))


V1_START = None


def fit(data, order, exposure="offset", covariate=False, starts=3, seed=0, xstats=(0, 1), x0=None, weights=None):
    ex = [g for g in order if g != "FC1"]
    nb = 4
    n = nb + (len(ex) if exposure in ("offset", "head") else 0) + (1 if covariate else 0)
    def unpack(x):
        lam, sig, gam, k = np.exp(x[:4])
        mu = {g: 0.0 for g in order}; m0 = {g: 0.0 for g in order}
        i = nb
        if exposure == "offset":
            for g in ex:
                mu[g] = float(x[i]); i += 1
        elif exposure == "head":
            for g in ex:
                m0[g] = float(1 / (1 + np.exp(-np.clip(x[i], -30, 30)))); i += 1
        beta = float(x[i]) if covariate else 0.0
        return lam, sig, gam, k, mu, m0, beta
    def nll(x):
        lam, sig, gam, k, mu, m0, beta = unpack(x)
        if not (1e-6 < sig < 20 and gam < 50 and 1e-6 < k < 10):
            return 1e12
        tot = sum(loglik(data[g], lam, sig, gam, k, mu[g], m0[g], beta, *xstats,
                         w=None if weights is None else weights[g]) for g in order)
        return -tot if np.isfinite(tot) else 1e12
    cands = []
    if x0 is not None:
        cands = [np.asarray(x0, float)]
    else:
        rng = np.random.default_rng(seed)
        base = np.array([np.log(0.055), np.log(2.0), np.log(0.9), np.log(0.0125)])
        for s in range(starts):
            v = np.zeros(n); v[:4] = base + (rng.normal(0, .3, 4) if s else 0); v[4:] = rng.normal(0, .2, n - 4) if s else 0.1
            cands.append(v)
    best = None
    for v in cands:
        r = minimize(nll, v, method="L-BFGS-B")
        if x0 is None:
            r = minimize(nll, r.x, method="Nelder-Mead", options={"maxiter": 1500, "xatol": 1e-5, "fatol": 1e-6})
            r = minimize(nll, r.x, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    lam, sig, gam, k, mu, m0, beta = unpack(best.x)
    return {"loglik": -float(best.fun), "k": n, "aic": 2 * n + 2 * float(best.fun), "lambda": float(lam),
            "sigma": float(sig), "gamma": float(gam), "acq_rate": float(k), "offset": mu, "head_start": m0,
            "beta": beta, "x": best.x.tolist()}


def lrt(big, small, df):
    st = 2 * (big["loglik"] - small["loglik"])
    return {"stat": st, "df": df, "p": float(chi2.sf(max(st, 0), df))}


def boot_weights(data, rng):
    """Resampling individuals with replacement = multinomial counts over original records."""
    out = {}
    for g, D in data.items():
        cnt = np.bincount(rng.integers(D["n"], size=D["n"]), minlength=D["n"])
        out[g] = np.bincount(D["inv"], weights=cnt, minlength=len(D["w"])).astype(float)
    return out


def p1(window, B):
    order = [g for g in GROUPS if group(GROUPS[g], window)]
    recs = {g: group(GROUPS[g], window) for g in order}
    data = {g: arrays(v) for g, v in recs.items()}
    fits = {"null": fit(data, order, "none"), "offset": fit(data, order, "offset"),
            "head": fit(data, order, "head")}
    dfx = len(order) - 1
    res = {"n_attached": {g: len(v) for g, v in recs.items()},
           "datasets": {g: sorted({r["dataset"] for r in v}) for g, v in recs.items()},
           "fits": fits, "lrt_offset_vs_null": lrt(fits["offset"], fits["null"], dfx),
           "lrt_head_vs_null": lrt(fits["head"], fits["null"], dfx),
           "dll_offset_minus_head": fits["offset"]["loglik"] - fits["head"]["loglik"]}
    rng = np.random.default_rng(11)
    d, offs = [], []
    for b in range(B):
        w = boot_weights(data, rng)
        o = fit(data, order, "offset", x0=fits["offset"]["x"], weights=w)
        h = fit(data, order, "head", x0=fits["head"]["x"], weights=w)
        d.append(o["loglik"] - h["loglik"]); offs.append(o["offset"])
    res["bootstrap_B"] = B
    if B:
        res["bootstrap_dll_q"] = np.quantile(d, [0.025, 0.5, 0.975]).tolist()
        res["bootstrap_frac_offset_better"] = float(np.mean(np.array(d) > 0))
        res["offset_ci"] = {g: np.quantile([o[g] for o in offs], [0.025, 0.975]).tolist() for g in order}
    return res


def p2p3(B):
    recs = {}
    for g in GROUPS:
        r = [x for x in group(GROUPS[g], "extended") if "log_ratio" in x]
        if r:
            recs[g] = r
    order = [g for g in GROUPS if g in recs]
    allx = np.concatenate([[r["log_ratio"] for r in v] for v in recs.values()])
    xs = (float(allx.mean()), float(allx.std()))
    data = {g: arrays(v, truncate=True, covariate=True) for g, v in recs.items()}
    base = fit(data, order, "offset", covariate=False, xstats=xs)
    cov = fit(data, order, "offset", covariate=True, xstats=xs)
    res = {"n_read": {g: len(v) for g, v in recs.items()}, "x_mean_sd": xs,
           "offset_without_cAMP": base, "offset_with_cAMP": cov, "lrt_beta": lrt(cov, base, 1),
           "mediated_fraction_FC2": 1 - cov["offset"]["FC2"] / base["offset"]["FC2"]}
    rng = np.random.default_rng(12)
    med, betas = [], []
    for b in range(B):
        w = boot_weights(data, rng)
        f0 = fit(data, order, "offset", xstats=xs, x0=base["x"], weights=w)
        f1 = fit(data, order, "offset", covariate=True, xstats=xs, x0=cov["x"], weights=w)
        betas.append(f1["beta"])
        if abs(f0["offset"]["FC2"]) > 1e-6:
            med.append(1 - f1["offset"]["FC2"] / f0["offset"]["FC2"])
    res["bootstrap_B"] = B
    if B:
        res["beta_ci"] = np.quantile(betas, [0.025, 0.975]).tolist()
        res["mediated_fraction_ci"] = np.quantile(med, [0.025, 0.975]).tolist() if med else None
    vals = {g: np.array([r["log_ratio"] for r in v]) for g, v in recs.items()}
    res["S1_median_log_ratio"] = {g: float(np.median(v)) for g, v in vals.items()}
    res["S1_kruskal"] = [float(x) for x in kruskal(*vals.values())]
    res["S1_fc2_vs_fc1"] = [float(x) for x in mannwhitneyu(vals["FC2"], vals["FC1"])]
    return res


def s2():
    out = {}
    for g, ds_list in GROUPS.items():
        mothers, daughters = [], []
        for ds in ds_list:
            p = OBS / ds / "tracks_full.csv"
            if not p.exists():
                continue
            rows = {int(r["track"]): r for r in csv.DictReader(p.open())}
            for r in rows.values():
                if r["origin"] == "split" and int(r["parent"]) in rows:
                    m = rows[int(r["parent"])]
                    if r["right_censored"] == "True":
                        continue
                    mothers.append(int(r["start"]) - int(m["start"]))
                    daughters.append(int(r["end"]) - int(r["start"]) + 1)
        if len(mothers) > 10:
            rho = spearmanr(mothers, daughters)
            out[g] = {"n": len(mothers), "rho": float(rho[0]), "p": float(rho[1])}
    return out


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    out = Path(__file__).resolve().parent.parent / "results" / "empirical_v2"
    out.mkdir(parents=True, exist_ok=True)
    res = {"plan": "ANALYSIS_PLAN_v2.md", "P1_extended": p1("extended", B), "P1_first3h": p1("first3h", B),
           "P2_P3": p2p3(B), "S2": s2()}
    (out / "fits_v2.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1)[:8000])


# ---- time-matched comparison (exploratory, after the plan; see deviations log) ----------
TIME_BINS_H = (0.0, 1.0, 3.0, 8.0, 20.0, 1e9)


def occupancy_lookup(ds):
    t, cells, _ = A.occupancy(ds)
    k = 20
    cs = np.convolve(cells, np.ones(k) / k, mode="same")
    return t, cs


def time_records(g):
    out = []
    for ds in GROUPS[g]:
        if not (OBS / ds / "tracks_full.csv").exists():
            continue
        t, cs = occupancy_lookup(ds)
        recs, _, _ = load(ds, "extended")
        for r in recs:
            th = r["start"] * FRAME_S / 3600
            r["bin"] = int(np.searchsorted(TIME_BINS_H, th, side="right") - 1)
            r["locc"] = round(float(np.log1p(cs[min(np.searchsorted(t, th), len(cs) - 1)])), 1)
            r["dataset"] = ds
            out.append(r)
    return out


def tarrays(recs, mode):
    f = np.array([r["frames"] for r in recs], float)
    c = np.array([r["cens"] for r in recs], float)
    cov = np.array([r["bin"] if mode == "time" else r["locc"] for r in recs], float)
    key = np.column_stack([f, c, cov])
    u, inv, w = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    return {"f": u[:, 0], "c": u[:, 1].astype(bool), "a0": np.full(len(u), (A.MIN_FRAMES - 1) * FRAME_S),
            "x": u[:, 2], "w": w.astype(float), "inv": inv.ravel(), "n": len(recs)}


def tfit(data, order, exposure, mode, x0=None, weights=None, starts=3):
    ex = [g for g in order if g != "FC1"]
    nt = len(TIME_BINS_H) - 2 if mode == "time" else 1
    n = 4 + nt + (len(ex) if exposure in ("offset", "head") else 0)
    def unpack(x):
        lam, sig, gam, k = np.exp(x[:4])
        tau = x[4:4 + nt]
        mu = {g: 0.0 for g in order}; m0 = {g: 0.0 for g in order}
        i = 4 + nt
        for g in ex:
            if exposure == "offset":
                mu[g] = float(x[i]); i += 1
            elif exposure == "head":
                m0[g] = float(1 / (1 + np.exp(-np.clip(x[i], -30, 30)))); i += 1
        return lam, sig, gam, k, tau, mu, m0
    def ll_group(D, lam, sig, gam, k, tau, mu, m0, w):
        if mode == "time":
            shift = np.concatenate([[0.0], tau])[D["x"].astype(int)]
        else:
            shift = tau[0] * D["x"]
        a_seen, a_last = (D["f"] - 1) * FRAME_S, D["f"] * FRAME_S
        rel = np.exp(-(mu + shift))[:, None] * np.exp(-sig * A.GH_X)[None, :]
        S = lambda a: np.exp(-A.cum_hazard(a, lam, gam, k, m0)[:, None] * rel) @ A.GH_W
        s0, ss = S(D["a0"]), S(a_seen)
        ll = np.where(D["c"], np.log(np.clip(ss / s0, 1e-300, None)),
                      np.log(np.clip((ss - S(a_last)) / s0, 1e-300, None)))
        return float(np.dot(D["w"] if w is None else w, ll))
    def nll(x):
        lam, sig, gam, k, tau, mu, m0 = unpack(x)
        if not (1e-6 < sig < 20 and gam < 50 and 1e-6 < k < 10):
            return 1e12
        v = sum(ll_group(data[g], lam, sig, gam, k, tau, mu[g], m0[g], None if weights is None else weights[g])
                for g in order)
        return -v if np.isfinite(v) else 1e12
    cands = [np.asarray(x0, float)] if x0 is not None else []
    rng = np.random.default_rng(1)
    for s in range(0 if x0 is not None else starts):
        v = np.zeros(n); v[:4] = np.array([np.log(0.03), np.log(1.8), np.log(1.1), np.log(0.011)]) + (rng.normal(0, .3, 4) if s else 0)
        v[4:] = rng.normal(0, .2, n - 4) if s else 0.05
        cands.append(v)
    best = None
    for v in cands:
        r = minimize(nll, v, method="L-BFGS-B")
        if x0 is None:
            r = minimize(nll, r.x, method="Nelder-Mead", options={"maxiter": 2000, "xatol": 1e-5, "fatol": 1e-6})
            r = minimize(nll, r.x, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    lam, sig, gam, k, tau, mu, m0 = unpack(best.x)
    return {"loglik": -float(best.fun), "k": n, "aic": 2 * n + 2 * float(best.fun), "sigma": float(sig),
            "gamma": float(gam), "acq_rate": float(k), "lambda": float(lam), "time_effect": [float(t) for t in tau],
            "offset": mu, "head_start": m0, "x": best.x.tolist()}


def time_matched(B=200, mode="time"):
    order = list(GROUPS)
    recs = {g: time_records(g) for g in order}
    data = {g: tarrays(v, mode) for g, v in recs.items()}
    fits = {e: tfit(data, order, e, mode) for e in ("none", "offset", "head")}
    dfx = len(order) - 1
    res = {"mode": mode, "bins_h": TIME_BINS_H[:-1], "n": {g: len(v) for g, v in recs.items()}, "fits": fits,
           "lrt_offset_vs_null": lrt(fits["offset"], fits["none"], dfx),
           "lrt_head_vs_null": lrt(fits["head"], fits["none"], dfx),
           "dll_offset_minus_head": fits["offset"]["loglik"] - fits["head"]["loglik"]}
    rng = np.random.default_rng(21)
    d, offs = [], []
    for b in range(B):
        w = boot_weights(data, rng)
        o = tfit(data, order, "offset", mode, x0=fits["offset"]["x"], weights=w)
        h = tfit(data, order, "head", mode, x0=fits["head"]["x"], weights=w)
        d.append(o["loglik"] - h["loglik"]); offs.append(o["offset"])
    if B:
        res["bootstrap_B"] = B
        res["bootstrap_dll_q"] = np.quantile(d, [0.025, 0.5, 0.975]).tolist()
        res["bootstrap_frac_offset_better"] = float(np.mean(np.array(d) > 0))
        res["offset_ci"] = {g: np.quantile([o[g] for o in offs], [0.025, 0.975]).tolist() for g in order}
    return res
