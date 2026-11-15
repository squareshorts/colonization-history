"""Cell-level residence after attachment in the first three hours after inoculation.

Inputs are the windowed tracks (tracks.csv) previously extracted from the
Lee et al. (2018) binary image stacks. An arrival is a track that appears after
the first frame; it is attached when visible in at least two consecutive 3-s
frames. Windows 0-2 cover the first three hours after inoculation in every flow
cell, so first-exposure (FC1) and transferred (FC2) arrivals are time-matched.

Three analyses:

1. Five-minute retention per arrival, P(visible 300 s after arrival), by
   Kaplan-Meier over discrete frames; arrival-level bootstrap intervals and
   per-field estimates.
2. Composition change versus graded shift. Residence after attachment follows a
   mixture of exponential detachment types. Composition change keeps the
   type-specific hazards fixed and lets the type weights differ between FC1 and
   FC2 (pre-existing types selected, or cells switched between fixed states).
   A graded shift keeps the weights fixed and multiplies every FC2 hazard by
   exp(-delta). With two types both families have four parameters, so their
   log likelihoods are directly comparable; the model containing both is the
   common alternative for likelihood-ratio tests. A parametric bootstrap
   measures how often the generating family is recovered at the archived
   sample sizes.
3. Design thresholds. The archived five-minute retention probabilities are
   used as sham (FC1) and exposed (FC2) success probabilities in the expected
   lower-bound inequalities and in a founder-indexed Monte Carlo design.

Usage from a repository checkout:
  python paper/analysis/residence_composition.py --repo . --out paper/results
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2, t as student_t

FRAME_S = 3.0
WINDOWS = (0, 1, 2)
GROUPS = {
    "FC1": ["WT YFP-CFP (T=29.5h;FC1)", "WT YFP-CFP (T=26.75h;FC1)",
            "WT YFP-CFP (T1=32.5h,T2=14.75h;FC1)", "WT YFP-CFP (t1=0-1h;FC1)"],
    "FC2": ["WT YFP-CFP (T=29.5h;FC2)", "WT YFP-CFP (T=26.75h;FC2)",
            "WT YFP-CFP (T1=32.5h,T2=14.75h;FC2)"],
    "FC3": ["WT YFP-CFP (T1=32.5h,T2=14.75h;FC3)"],
    "Tr14": ["WT YFP-CFP (Tr=14h)"],
    "Tr37": ["WT YFP-CFP (Tr=37h)"],
}
FIVE_MIN_FRAMES = 101  # visible at age 300 s <=> at least 101 frames


# ---------------------------------------------------------------- data
def load_dataset(repo, name, windows=WINDOWS):
    loc = Path(repo) / "empirical" / "observables" / name
    meta = json.loads((loc / "meta.json").read_text())
    wf = meta["windows_frames"]
    rows = []
    with (loc / "tracks.csv").open(newline="") as f:
        for r in csv.DictReader(f):
            w = int(r["window"])
            if w not in windows or r["origin"] != "arrival":
                continue
            s, e = int(r["start"]), int(r["end"])
            limit = wf[w][1] - s  # frames observable from arrival to window end
            rows.append({"frames": e - s + 1, "cens": r["right_censored"] == "True",
                         "limit": limit, "dataset": name})
    return rows


def load_group(repo, names):
    out = []
    for n in names:
        if (Path(repo) / "empirical" / "observables" / n / "tracks.csv").exists():
            out += load_dataset(repo, n)
    return out


# ---------------------------------------------------------------- Kaplan-Meier
def km_visible_at(frames, cens, k=FIVE_MIN_FRAMES):
    """P(D >= k) for discrete visible-frame counts D with right censoring."""
    frames, cens = np.asarray(frames), np.asarray(cens, bool)
    s = 1.0
    for j in range(1, k):
        at_risk = np.sum(frames >= j)
        events = np.sum((frames == j) & ~cens)
        if at_risk == 0:
            break
        s *= 1 - events / at_risk
    return float(s)


def retention_summary(rows, rng, boot=2000):
    f = np.array([r["frames"] for r in rows])
    c = np.array([r["cens"] for r in rows])
    attached = f >= 2
    est = {"arrivals": int(len(f)), "attached": int(attached.sum()),
           "attached_fraction": float(attached.mean()),
           "visible_60s_given_attached": km_visible_at(f[attached], c[attached], 21),
           "visible_300s_given_attached": km_visible_at(f[attached], c[attached]),
           "retained_300s_per_arrival": km_visible_at(f, c)}
    bs = []
    ds = np.array([r["dataset"] for r in rows])
    idx_by = [np.flatnonzero(ds == d) for d in np.unique(ds)]
    for _ in range(boot):
        ii = np.concatenate([rng.choice(ix, len(ix)) for ix in idx_by])
        bs.append(km_visible_at(f[ii], c[ii]))
    est["retained_300s_per_arrival_ci95"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
    est["per_field_retained_300s"] = {
        d: km_visible_at(f[ds == d], c[ds == d]) for d in np.unique(ds)}
    est["per_field_arrivals"] = {d: int(np.sum(ds == d)) for d in np.unique(ds)}
    return est


# ---------------------------------------------------------------- mixture likelihood
def compress(rows):
    """Attached arrivals collapsed to unique (frames, censored) with counts."""
    f = np.array([r["frames"] for r in rows if r["frames"] >= 2], float)
    c = np.array([r["cens"] for r in rows if r["frames"] >= 2], bool)
    key, w = np.unique(np.column_stack([f, c]), axis=0, return_counts=True)
    return {"f": key[:, 0], "c": key[:, 1].astype(bool), "w": w.astype(float)}


def mix_loglik(D, weights, rates):
    """Interval-censored log likelihood conditional on attachment (age >= 3 s)."""
    weights, rates = np.asarray(weights), np.asarray(rates)
    def S(a):
        return np.exp(-np.outer(a, rates)) @ weights
    a_seen = (D["f"] - 1) * FRAME_S
    a_gone = D["f"] * FRAME_S
    s0 = S(np.array([FRAME_S]))[0]
    ss = S(a_seen)
    ll = np.where(D["c"], np.log(np.clip(ss / s0, 1e-300, None)),
                  np.log(np.clip((ss - S(a_gone)) / s0, 1e-300, None)))
    return float(D["w"] @ ll)


def softmax(z):
    z = np.concatenate([[0.], z])
    e = np.exp(z - z.max())
    return e / e.sum()


def rates_from(x):
    """Strictly decreasing positive rates: first rate free, then log-gaps."""
    r = [np.exp(x[0])]
    for g in x[1:]:
        r.append(r[-1] * np.exp(-np.exp(g)))
    return np.array(r)


def fit_family(data, family, m=2, starts=8, seed=0):
    """Fit composition ('reweight'), graded shift ('shift') or both ('both')."""
    g1, g2 = data["FC1"], data["FC2"]
    nr = m
    nw = m - 1

    def unpack(x):
        rates = rates_from(x[:nr])
        i = nr
        if family == "reweight":
            w1 = softmax(x[i:i + nw]); w2 = softmax(x[i + nw:i + 2 * nw]); delta = 0.
        elif family == "shift":
            w1 = w2 = softmax(x[i:i + nw]); delta = x[i + nw]
        elif family == "both":
            w1 = softmax(x[i:i + nw]); w2 = softmax(x[i + nw:i + 2 * nw]); delta = x[i + 2 * nw]
        elif family == "null":
            w1 = w2 = softmax(x[i:i + nw]); delta = 0.
        else:
            raise ValueError(family)
        return rates, w1, w2, delta

    n = nr + {"reweight": 2 * nw, "shift": nw + 1, "both": 2 * nw + 1, "null": nw}[family]

    def nll(x):
        rates, w1, w2, delta = unpack(x)
        v = mix_loglik(g1, w1, rates) + mix_loglik(g2, w2, rates * np.exp(-delta))
        return -v if np.isfinite(v) else 1e15

    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        x0 = np.zeros(n)
        x0[0] = np.log(0.08) + (rng.normal(0, .5) if s else 0)
        x0[1:nr] = np.log(1.5) + (rng.normal(0, .4, nr - 1) if s else 0)
        x0[nr:] = rng.normal(0, .6, n - nr) if s else 0.
        r = minimize(nll, x0, method="Nelder-Mead",
                     options={"maxiter": 6000, "xatol": 1e-7, "fatol": 1e-8})
        r = minimize(nll, r.x, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    rates, w1, w2, delta = unpack(best.x)
    return {"family": family, "components": m, "k": n, "loglik": -float(best.fun),
            "aic": 2 * n + 2 * float(best.fun), "rates_per_s": rates.tolist(),
            "weights_FC1": w1.tolist(), "weights_FC2": w2.tolist(), "delta": float(delta)}


# ---------------------------------------------------------------- simulation
def simulate_group(rng, weights, rates, limits):
    """Attached residence frames for arrivals with given follow-up limits."""
    weights, rates = np.asarray(weights), np.asarray(rates)
    cond = weights * np.exp(-rates * FRAME_S)
    cond /= cond.sum()
    comp = rng.choice(len(rates), size=len(limits), p=cond)
    T = FRAME_S + rng.exponential(1 / rates[comp])
    d = np.floor(T / FRAME_S).astype(int) + 1
    cens = d >= limits
    d = np.where(cens, limits, d)
    return [{"frames": int(a), "cens": bool(b)} for a, b in zip(d, cens)]


def bootstrap_discrimination(rng, data_rows, fits, sims=100):
    lim = {g: np.array([r["limit"] for r in data_rows[g] if r["frames"] >= 2]) for g in ("FC1", "FC2")}
    out = {}
    for gen in ("reweight", "shift"):
        f = fits[gen]
        wins, margins = 0, []
        for s in range(sims):
            sim = {"FC1": compress(simulate_group(rng, f["weights_FC1"], f["rates_per_s"], lim["FC1"])),
                   "FC2": compress(simulate_group(rng, f["weights_FC2"],
                                                  np.array(f["rates_per_s"]) * np.exp(-f["delta"]),
                                                  lim["FC2"]))}
            a = fit_family(sim, "reweight", starts=3, seed=s)["loglik"]
            b = fit_family(sim, "shift", starts=3, seed=s)["loglik"]
            m = (a - b) if gen == "reweight" else (b - a)
            margins.append(m)
            wins += m > 0
        out[f"generated_by_{gen}"] = {
            "simulations": sims, "generating_family_preferred": wins / sims,
            "preferred_by_more_than_2_loglik": float(np.mean(np.array(margins) > 2)),
            "margin_quantiles_5_50_95": [float(q) for q in np.quantile(margins, [.05, .5, .95])]}
    return out


# ---------------------------------------------------------------- design calibration
def threshold_rows(p0, p1, pulse_losses=(0., .05, .10)):
    rows = []
    for e1 in pulse_losses:
        rows.append({"p_sham": p0, "p_exposed": p1, "pulse_loss": e1,
                     "max_sham_loss_for_positive_expected_lower_bound":
                         ((1 - e1) * p1 - p0) / (1 - p0)})
    rows.append({"p_sham": p0, "p_exposed": p1, "pulse_loss": "equal",
                 "max_sham_loss_for_positive_expected_lower_bound": (p1 - p0) / (1 + p1 - p0)})
    return rows


def design_power(rng, p0, p1, loss, blocks, founders=48, studies=4000, alpha=.025):
    """Single-endpoint founder-indexed design with equal uniform loss in both arms."""
    n1 = rng.binomial(founders, 1 - loss, size=(studies, blocks))
    n0 = rng.binomial(founders, 1 - loss, size=(studies, blocks))
    s1 = rng.binomial(n1, p1)
    s0 = rng.binomial(n0, p0)
    low = s1 / founders - (s0 + founders - n0) / founders
    raw = s1 / np.maximum(n1, 1) - s0 / np.maximum(n0, 1)
    def test(d):
        se = d.std(axis=1, ddof=1) / np.sqrt(blocks)
        return d.mean(axis=1) - student_t.ppf(1 - alpha, blocks - 1) * se > 0
    return float(test(low).mean()), float(test(raw).mean())


def blocks_for_power(rng, p0, p1, loss, target=.8, founders=48):
    for b in (8, 16, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024):
        if design_power(rng, p0, p1, loss, b, founders, studies=2000)[0] >= target:
            return b
    return None


def write_csv(path, rows):
    path = Path(path)
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sims", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20261115)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rows = {g: load_group(a.repo, names) for g, names in GROUPS.items()}

    retention = {g: retention_summary(r, rng) for g, r in rows.items()}

    data = {g: compress(rows[g]) for g in ("FC1", "FC2")}
    fits = {}
    for m in (2, 3):
        for fam in ("null", "reweight", "shift", "both"):
            fits[f"{fam}_{m}"] = fit_family(data, fam, m=m)
    tests = {}
    for m in (2, 3):
        both, rw, sh, nu = (fits[f"{x}_{m}"] for x in ("both", "reweight", "shift", "null"))
        tests[f"{m}_types"] = {
            "reweight_minus_shift_loglik": rw["loglik"] - sh["loglik"],
            "shift_beyond_reweight": {"stat": 2 * (both["loglik"] - rw["loglik"]), "df": both["k"] - rw["k"],
                                      "p": float(chi2.sf(max(2 * (both["loglik"] - rw["loglik"]), 0),
                                                         both["k"] - rw["k"]))},
            "reweight_beyond_shift": {"stat": 2 * (both["loglik"] - sh["loglik"]), "df": both["k"] - sh["k"],
                                      "p": float(chi2.sf(max(2 * (both["loglik"] - sh["loglik"]), 0),
                                                         both["k"] - sh["k"]))},
            "any_difference_vs_null": {"stat": 2 * (both["loglik"] - nu["loglik"]), "df": both["k"] - nu["k"],
                                       "p": float(chi2.sf(max(2 * (both["loglik"] - nu["loglik"]), 0),
                                                          both["k"] - nu["k"]))},
        }
    disc = bootstrap_discrimination(rng, rows, {"reweight": fits["reweight_2"], "shift": fits["shift_2"]},
                                    sims=a.sims)

    p0 = retention["FC1"]["retained_300s_per_arrival"]
    p1 = retention["FC2"]["retained_300s_per_arrival"]
    thr = threshold_rows(p0, p1)
    design = []
    for loss in (0., .01, .02, .03, .05):
        pb, pr = design_power(rng, p0, p1, loss, blocks=32)
        design.append({"p_sham": p0, "p_exposed": p1, "equal_loss": loss, "blocks": 32, "founders_per_arm": 48,
                       "bounded_power": pb, "recovered_only_power": pr,
                       "blocks_for_80pct_bounded": blocks_for_power(rng, p0, p1, loss)})

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "residence_composition.json").write_text(json.dumps(
        {"window": "first three hours after inoculation (tracking windows 0-2)",
         "retention": retention, "fits": fits, "tests": tests,
         "parametric_bootstrap": disc}, indent=2) + "\n")
    write_csv(out / "archive_calibrated_thresholds.csv", thr)
    write_csv(out / "archive_calibrated_design.csv", design)
    print(json.dumps({"retention": {g: (v["arrivals"], v["attached"],
                                        round(v["visible_300s_given_attached"], 4),
                                        round(v["retained_300s_per_arrival"], 4),
                                        [round(x, 4) for x in v["retained_300s_per_arrival_ci95"]])
                                    for g, v in retention.items()},
                      "tests": tests, "bootstrap": disc, "thresholds": thr, "design": design}, indent=1))


if __name__ == "__main__":
    main()
