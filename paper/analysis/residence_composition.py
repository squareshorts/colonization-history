"""Persistence of tracked image components in the first three hours after inoculation.

Inputs are the windowed tracks (tracks.csv) extracted from the Lee et al. (2018)
binary image stacks by empirical/extract.py. That tracker links connected
components by nearest centroid within 20 px, bridges gaps of up to two frames,
labels a second component claiming the same predecessor as a split, and ends a
track when no successor appears within the gap allowance. Merges and exits
through the field boundary therefore also end tracks. The measured outcome is
the persistence of a linked image component, which approximates the residence
of an attached cell.

An arrival is a track whose origin is 'arrival' and which starts after the
first frame of its window; it counts as attached when its track spans at least
`min_frames` frames (default two). Windows 0-2 cover the first three hours after
inoculation in every flow cell.

Primary comparison: the three matched experimental sequences (FC1 and FC2 of
series 1-3). Sensitivity analyses: adding the separate first-hour FC1 recording,
each sequence alone, stricter attachment criteria, and counting split-labelled
tracks as arrivals.

Families for persistence after attachment, a mixture of m exponential types:
  composition  type rates shared, type proportions differ between FC1 and FC2
  rate_shift   proportions shared, every FC2 rate multiplied by exp(-delta)
  age_advance  every FC2 component follows the FC1 hazard advanced by A seconds,
               S2(a) = S1(a + A) / S1(A); with two types this is the composition
               family reparameterized, with more types it is a one-parameter
               restriction of it
  both         composition and rate_shift together
  null         no difference

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
MATCHED = [("WT YFP-CFP (T=29.5h;FC1)", "WT YFP-CFP (T=29.5h;FC2)"),
           ("WT YFP-CFP (T=26.75h;FC1)", "WT YFP-CFP (T=26.75h;FC2)"),
           ("WT YFP-CFP (T1=32.5h,T2=14.75h;FC1)", "WT YFP-CFP (T1=32.5h,T2=14.75h;FC2)")]
EXTRA_FC1 = "WT YFP-CFP (t1=0-1h;FC1)"
GROUPS = {
    "FC1": [p[0] for p in MATCHED],
    "FC2": [p[1] for p in MATCHED],
    "FC1_with_first_hour_recording": [p[0] for p in MATCHED] + [EXTRA_FC1],
    "FC3": ["WT YFP-CFP (T1=32.5h,T2=14.75h;FC3)"],
    "Tr14": ["WT YFP-CFP (Tr=14h)"],
    "Tr37": ["WT YFP-CFP (Tr=37h)"],
}
FIVE_MIN_FRAMES = 101  # visible at age 300 s <=> track spans at least 101 frames


# ---------------------------------------------------------------- data
def load_dataset(repo, name, windows=WINDOWS, origins=("arrival",)):
    loc = Path(repo) / "empirical" / "observables" / name
    meta = json.loads((loc / "meta.json").read_text())
    wf = meta["windows_frames"]
    rows = []
    with (loc / "tracks.csv").open(newline="") as f:
        for r in csv.DictReader(f):
            w = int(r["window"])
            if w not in windows or r["origin"] not in origins:
                continue
            s, e = int(r["start"]), int(r["end"])
            rows.append({"frames": e - s + 1, "cens": r["right_censored"] == "True",
                         "limit": wf[w][1] - s, "dataset": name})
    return rows


def load_group(repo, names, origins=("arrival",)):
    out = []
    for n in names:
        if (Path(repo) / "empirical" / "observables" / n / "tracks.csv").exists():
            out += load_dataset(repo, n, origins=origins)
    return out


# ---------------------------------------------------------------- Kaplan-Meier
def km_visible_at(frames, cens, k=FIVE_MIN_FRAMES):
    """P(D >= k) for discrete track spans D with right censoring."""
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
    ds = np.array([r["dataset"] for r in rows])
    est = {"arrivals": int(len(f)), "attached": int(attached.sum()),
           "attached_fraction": float(attached.mean()),
           "visible_60s_given_attached": km_visible_at(f[attached], c[attached], 21),
           "visible_300s_given_attached": km_visible_at(f[attached], c[attached]),
           "retained_300s_per_arrival": km_visible_at(f, c)}
    idx_by = [np.flatnonzero(ds == d) for d in np.unique(ds)]
    bs = [km_visible_at(f[ii], c[ii]) for ii in
          (np.concatenate([rng.choice(ix, len(ix)) for ix in idx_by]) for _ in range(boot))]
    est["retained_300s_per_arrival_ci95_within_fields"] = [float(np.quantile(bs, .025)),
                                                           float(np.quantile(bs, .975))]
    est["per_field_retained_300s"] = {d: km_visible_at(f[ds == d], c[ds == d]) for d in np.unique(ds)}
    est["per_field_arrivals"] = {d: int(np.sum(ds == d)) for d in np.unique(ds)}
    return est


def paired_fields(repo):
    """Five-minute retention per arrival in each matched sequence."""
    rows, diffs = [], []
    for s, (a, b) in enumerate(MATCHED, 1):
        r0, r1 = load_dataset(repo, a), load_dataset(repo, b)
        p0 = km_visible_at([x["frames"] for x in r0], [x["cens"] for x in r0])
        p1 = km_visible_at([x["frames"] for x in r1], [x["cens"] for x in r1])
        rows.append({"series": s, "first_exposure": p0, "after_transfer": p1, "difference": p1 - p0,
                     "arrivals_first_exposure": len(r0), "arrivals_after_transfer": len(r1)})
        diffs.append(p1 - p0)
    d = np.array(diffs)
    half = student_t.ppf(.975, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
    return {"pairs": rows, "mean_difference": float(d.mean()),
            "t95_interval_between_sequences": [float(d.mean() - half), float(d.mean() + half)]}


# ---------------------------------------------------------------- mixture likelihood
def compress(rows, min_frames=2):
    f = np.array([r["frames"] for r in rows if r["frames"] >= min_frames], float)
    c = np.array([r["cens"] for r in rows if r["frames"] >= min_frames], bool)
    key, w = np.unique(np.column_stack([f, c]), axis=0, return_counts=True)
    return {"f": key[:, 0], "c": key[:, 1].astype(bool), "w": w.astype(float),
            "a0": (min_frames - 1) * FRAME_S}


def mix_loglik(D, weights, rates):
    """Interval-censored log likelihood conditional on persistence to the attachment age."""
    weights, rates = np.asarray(weights), np.asarray(rates)
    def S(a):
        return np.exp(-np.outer(a, rates)) @ weights
    a_seen = (D["f"] - 1) * FRAME_S
    a_gone = D["f"] * FRAME_S
    s0 = S(np.array([D.get("a0", FRAME_S)]))[0]
    ss = S(a_seen)
    ll = np.where(D["c"], np.log(np.clip(ss / s0, 1e-300, None)),
                  np.log(np.clip((ss - S(a_gone)) / s0, 1e-300, None)))
    return float(D["w"] @ ll)


def softmax(z):
    z = np.concatenate([[0.], z])
    e = np.exp(z - z.max())
    return e / e.sum()


def rates_from(x):
    r = [np.exp(x[0])]
    for g in x[1:]:
        r.append(r[-1] * np.exp(-np.exp(g)))
    return np.array(r)


N_EXTRA = {"composition": lambda nw: 2 * nw, "rate_shift": lambda nw: nw + 1,
           "age_advance": lambda nw: nw + 1, "both": lambda nw: 2 * nw + 1, "null": lambda nw: nw}


def fit_family(data, family, m=2, starts=8, seed=0):
    g1, g2 = data["FC1"], data["FC2"]
    nr, nw = m, m - 1

    def unpack(x):
        rates = rates_from(x[:nr])
        i = nr
        delta, age = 0., 0.
        if family == "composition":
            w1 = softmax(x[i:i + nw]); w2 = softmax(x[i + nw:i + 2 * nw])
        elif family == "rate_shift":
            w1 = w2 = softmax(x[i:i + nw]); delta = x[i + nw]
        elif family == "age_advance":
            w1 = softmax(x[i:i + nw]); age = float(np.exp(np.clip(x[i + nw], -20, 12)))
            lw = np.log(np.maximum(w1, 1e-300)) - rates * age
            w2 = np.exp(lw - lw.max()); w2 = w2 / w2.sum()
        elif family == "both":
            w1 = softmax(x[i:i + nw]); w2 = softmax(x[i + nw:i + 2 * nw]); delta = x[i + 2 * nw]
        elif family == "null":
            w1 = w2 = softmax(x[i:i + nw])
        else:
            raise ValueError(family)
        return rates, w1, w2, delta, age

    n = nr + N_EXTRA[family](nw)

    def nll(x):
        rates, w1, w2, delta, _ = unpack(x)
        v = mix_loglik(g1, w1, rates) + mix_loglik(g2, w2, rates * np.exp(-delta))
        return -v if np.isfinite(v) else 1e15

    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        x0 = np.zeros(n)
        x0[0] = np.log(0.08) + (rng.normal(0, .5) if s else 0)
        x0[1:nr] = np.log(1.5) + (rng.normal(0, .4, nr - 1) if s else 0)
        x0[nr:] = rng.normal(0, .6, n - nr) if s else 0.
        if family == "age_advance":
            x0[-1] = np.log(15.) + (rng.normal(0, 1.) if s else 0)
        r = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 8000, "xatol": 1e-7, "fatol": 1e-8})
        r = minimize(nll, r.x, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    rates, w1, w2, delta, age = unpack(best.x)
    return {"family": family, "components": m, "k": n, "loglik": -float(best.fun),
            "aic": 2 * n + 2 * float(best.fun), "rates_per_s": rates.tolist(),
            "weights_FC1": list(map(float, w1)), "weights_FC2": list(map(float, w2)),
            "delta": float(delta), "age_advance_s": age}


def compare(data, ms=(2, 3), families=("null", "composition", "rate_shift", "age_advance", "both")):
    out = {}
    for m in ms:
        fits = {fam: fit_family(data, fam, m=m) for fam in families}
        c, r = fits["composition"], fits["rate_shift"]
        out[f"{m}_types"] = {
            "fits": fits,
            "composition_minus_rate_shift_loglik": c["loglik"] - r["loglik"],
            "aic_rate_shift_minus_composition": r["aic"] - c["aic"],
            "composition_minus_age_advance_loglik": c["loglik"] - fits["age_advance"]["loglik"],
            "lrt_rate_shift_beyond_composition": lrt(fits["both"], c),
            "lrt_composition_beyond_rate_shift": lrt(fits["both"], r),
        }
    return out


def lrt(big, small):
    st = max(2 * (big["loglik"] - small["loglik"]), 0.)
    df = big["k"] - small["k"]
    return {"stat": st, "df": df, "p": float(chi2.sf(st, df))}


def age_advance_equivalence(fit2):
    """Age advance that reproduces a two-type composition fit exactly."""
    lf, ls = fit2["rates_per_s"]
    w0, w1 = fit2["weights_FC1"][1], fit2["weights_FC2"][1]
    lg = lambda p: np.log(p / (1 - p))
    A = (lg(w1) - lg(w0)) / (lf - ls)
    a = np.linspace(FRAME_S, 3600, 2000)
    S0 = lambda x: (1 - w0) * np.exp(-lf * x) + w0 * np.exp(-ls * x)
    S1 = lambda x: (1 - w1) * np.exp(-lf * x) + w1 * np.exp(-ls * x)
    err = float(np.max(np.abs(S1(a) / S1(FRAME_S) - S0(a + A) / S0(A + FRAME_S))))
    return {"age_advance_s": float(A), "max_abs_survival_difference": err}


def simulate_group(rng, weights, rates, limits):
    """Track spans for attached arrivals from a mixture, with follow-up limits."""
    weights, rates = np.asarray(weights), np.asarray(rates)
    cond = weights * np.exp(-rates * FRAME_S)
    cond /= cond.sum()
    comp = rng.choice(len(rates), size=len(limits), p=cond)
    T = FRAME_S + rng.exponential(1 / rates[comp])
    d = np.floor(T / FRAME_S).astype(int) + 1
    cens = d >= limits
    d = np.where(cens, limits, d)
    return [{"frames": int(a), "cens": bool(b)} for a, b in zip(d, cens)]


# ---------------------------------------------------------------- design calibration
def threshold_rows(p0, p1, label, pulse_losses=(0., .05, .10)):
    rows = [{"calibration": label, "p_sham": p0, "p_exposed": p1, "pulse_loss": e1,
             "max_sham_loss_for_positive_expected_lower_bound": ((1 - e1) * p1 - p0) / (1 - p0)}
            for e1 in pulse_losses]
    rows.append({"calibration": label, "p_sham": p0, "p_exposed": p1, "pulse_loss": "equal",
                 "max_sham_loss_for_positive_expected_lower_bound": (p1 - p0) / (1 + p1 - p0)})
    return rows


def design_power(rng, p0, p1, loss, blocks, founders=48, studies=4000, alpha=.025):
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
    ap.add_argument("--seed", type=int, default=20261115)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rows = {g: load_group(a.repo, names) for g, names in GROUPS.items()}
    retention = {g: retention_summary(r, rng) for g, r in rows.items()}
    pairs = paired_fields(a.repo)

    matched = {"FC1": compress(rows["FC1"]), "FC2": compress(rows["FC2"])}
    pooled = {"FC1": compress(rows["FC1_with_first_hour_recording"]), "FC2": compress(rows["FC2"])}
    models = {"matched": compare(matched), "with_first_hour_recording": compare(pooled)}
    models["matched"]["equivalence_2_types"] = age_advance_equivalence(
        models["matched"]["2_types"]["fits"]["composition"])
    per_sequence = {}
    for s, (n1, n2) in enumerate(MATCHED, 1):
        d = {"FC1": compress(load_dataset(a.repo, n1)), "FC2": compress(load_dataset(a.repo, n2))}
        c = fit_family(d, "composition", m=2)
        r = fit_family(d, "rate_shift", m=2)
        per_sequence[f"series_{s}"] = {"attached_FC1": int(d["FC1"]["w"].sum()),
                                       "attached_FC2": int(d["FC2"]["w"].sum()),
                                       "composition_minus_rate_shift_loglik": c["loglik"] - r["loglik"]}
    tracking = []
    for label, kw, origins in (("attachment span >= 2 frames (primary)", 2, ("arrival",)),
                               ("attachment span >= 3 frames", 3, ("arrival",)),
                               ("attachment span >= 5 frames", 5, ("arrival",)),
                               ("split-labelled tracks counted as arrivals", 2, ("arrival", "split"))):
        r1 = load_group(a.repo, GROUPS["FC1"], origins)
        r2 = load_group(a.repo, GROUPS["FC2"], origins)
        d = {"FC1": compress(r1, kw), "FC2": compress(r2, kw)}
        row = {"analysis": label,
               "retained_300s_FC1": km_visible_at([x["frames"] for x in r1], [x["cens"] for x in r1]),
               "retained_300s_FC2": km_visible_at([x["frames"] for x in r2], [x["cens"] for x in r2]),
               "attached_FC1": int(d["FC1"]["w"].sum()), "attached_FC2": int(d["FC2"]["w"].sum())}
        for m in (2, 3):
            row[f"composition_minus_rate_shift_loglik_{m}_types"] = (
                fit_family(d, "composition", m=m)["loglik"] - fit_family(d, "rate_shift", m=m)["loglik"])
        tracking.append(row)

    p1 = retention["FC2"]["retained_300s_per_arrival"]
    cal = {"matched first-exposure fields": retention["FC1"]["retained_300s_per_arrival"],
           "first-exposure fields with first-hour recording":
               retention["FC1_with_first_hour_recording"]["retained_300s_per_arrival"]}
    thr, design = [], []
    for label, p0 in cal.items():
        thr += threshold_rows(p0, p1, label)
        for loss in (0., .01, .02, .03, .04, .05):
            pb, pr = design_power(rng, p0, p1, loss, blocks=32)
            design.append({"calibration": label, "p_sham": p0, "p_exposed": p1, "equal_loss": loss,
                           "blocks": 32, "founders_per_arm": 48, "bounded_power": pb,
                           "recovered_only_power": pr,
                           "blocks_for_80pct_bounded": blocks_for_power(rng, p0, p1, loss)})
    for pr_ in pairs["pairs"]:
        d = pr_["difference"]
        thr.append({"calibration": f"series {pr_['series']} fields", "p_sham": pr_["first_exposure"],
                    "p_exposed": pr_["after_transfer"], "pulse_loss": "equal",
                    "max_sham_loss_for_positive_expected_lower_bound": d / (1 + d) if d > 0 else 0.})

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "residence_composition.json").write_text(json.dumps(
        {"window": "first three hours after inoculation (tracking windows 0-2)",
         "outcome": "persistence of a linked image component",
         "retention": retention, "matched_pairs": pairs, "models": models,
         "per_sequence": per_sequence, "tracking_sensitivity": tracking}, indent=2) + "\n")
    write_csv(out / "archive_calibrated_thresholds.csv", thr)
    write_csv(out / "archive_calibrated_design.csv", design)
    write_csv(out / "residence_tracking_sensitivity.csv", tracking)
    summary = {
        "retention": {g: (v["arrivals"], v["attached"], round(v["retained_300s_per_arrival"], 4))
                      for g, v in retention.items()},
        "pairs": pairs,
        "matched": {k: {kk: vv for kk, vv in v.items() if kk != "fits"} if isinstance(v, dict) and "fits" in v
                    else v for k, v in models["matched"].items()},
        "pooled": {k: {kk: vv for kk, vv in v.items() if kk != "fits"} for k, v in
                   models["with_first_hour_recording"].items()},
        "per_sequence": per_sequence, "tracking": tracking, "thresholds": thr}
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
