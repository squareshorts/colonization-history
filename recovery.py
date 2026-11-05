"""Mechanism recovery from synthetic transfer experiments of realistic size.

Each synthetic experiment draws one parameter vector from broad priors under one
of six generative hypotheses, then simulates replicate surfaces in every arm of
the largest design (fresh/fresh plus three history arms at each gap). Smaller
designs are obtained by subsetting replicates, gaps and arms. A random-forest
classifier (likelihood-free model choice) is trained on summary statistics of
simulated experiments and scored on held-out experiments.

Experiments are retained so that the fresh-cell / fresh-surface arm has the same
distribution under every hypothesis (stratified matching on two baseline
summaries). A classifier given only the fresh/fresh arm is the negative control.
"""
import argparse, json, time
from dataclasses import asdict
from itertools import product
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from model import Config, MECHANISMS, ACTIVE, simulate_batch, OBS_TIMES

GAPS = (0.0, 4.0, 16.0)
N_MAX = 8
CELLS = [("fresh_fresh", 0.0)] + [(a, g) for g in GAPS for a in
                                   ("exposed_fresh", "fresh_conditioned", "exposed_conditioned")]
SAMPLE_TIMES = tuple(range(1, 18, 2))      # indices into OBS_TIMES -> t = 2, 4, ..., 18
RAW_KEYS = ("integrated_occupancy", "final_attached", "attached_early", "attached_late",
            "founder_departures", "founder_quick", "descendant_departures", "descendant_quick",
            "divisions", "initial_residual", "retained_residual")

# Priors: (kind, low, high); kind 'u' uniform, 'log' log-uniform, 'logit' logit-uniform.
PRIORS = {
    "attachment_probability": ("logit", 0.08, 0.50),
    "detachment_rate": ("log", 0.5, 3.0),
    "division_rate": ("log", 0.2, 1.0),
    "patch_heterogeneity": ("u", 0.0, 0.7),
    "cell_heterogeneity": ("u", 0.0, 0.8),
    "memory_acquisition_rate": ("log", 0.3, 3.0),
    "cell_memory_decay_time": ("log", 1.0, 16.0),
    "cell_attachment_effect": ("u", 0.5, 3.0),
    "cell_retention_effect": ("u", 0.5, 2.5),
    "surface_decay_time": ("log", 3.0, 48.0),
    "surface_attachment_effect": ("u", 0.5, 3.0),
    "surface_retention_effect": ("u", 0.5, 2.5),
    "contact_trace": ("u", 0.03, 0.15),
    "resident_deposition_rate": ("u", 0.05, 0.30),
    "carryover_cells": ("log", 1.0, 20.0),
    "storage_loss_rate": ("log", 0.01, 0.3),
    "selection_shift": ("u", 0.1, 1.0),
    "selection_relaxation_time": ("log", 4.0, 400.0),
}
OFF = {"cell": ("cell_attachment_effect", "cell_retention_effect"),
       "surface": ("surface_attachment_effect", "surface_retention_effect"),
       "carryover": ("carryover_cells",), "selection": ("selection_shift",)}


def draw_params(mechanism, n, rng):
    out = {}
    for k, (kind, lo, hi) in PRIORS.items():
        u = rng.random(n)
        if kind == "u":
            out[k] = lo + (hi - lo) * u
        elif kind == "log":
            out[k] = np.exp(np.log(lo) + (np.log(hi) - np.log(lo)) * u)
        else:
            a, b = np.log(lo / (1 - lo)), np.log(hi / (1 - hi))
            out[k] = 1 / (1 + np.exp(-(a + (b - a) * u)))
    for family, keys in OFF.items():
        if family not in ACTIVE[mechanism]:
            for k in keys:
                out[k] = np.zeros(n)
    return out


def simulate_chunk(job):
    """Simulate E experiments of one mechanism; return raw replicate outputs."""
    mechanism, E, seed = job
    cfg = Config()
    rng = np.random.default_rng(seed)
    theta = draw_params(mechanism, E, rng)
    C, N = len(CELLS), N_MAX
    rep = lambda v: np.repeat(v, C * N)
    params = {k: rep(v) for k, v in theta.items()}
    cell_idx = np.tile(np.repeat(np.arange(C), N), E)
    exposed = np.array([CELLS[c][0] in ("exposed_fresh", "exposed_conditioned") for c in cell_idx])
    conditioned = np.array([CELLS[c][0] in ("fresh_conditioned", "exposed_conditioned") for c in cell_idx])
    gap = np.array([CELLS[c][1] for c in cell_idx])
    r = simulate_batch(cfg, params, exposed, conditioned, gap, seed=seed + 1)
    raw = {k: np.asarray(r[k], np.float32).reshape(E, C, N) for k in RAW_KEYS}
    raw["counts"] = r["counts"].astype(np.float32).reshape(E, C, N, -1)
    raw["theta"] = np.column_stack([theta[k] for k in PRIORS]).astype(np.float32)
    return mechanism, raw


def simulate_pool(n_per_mechanism, chunk=200, workers=2, seed=60013, out=None):
    jobs, s = [], seed
    for m in MECHANISMS:
        left = n_per_mechanism
        while left > 0:
            e = min(chunk, left)
            jobs.append((m, e, s))
            s += 7919
            left -= e
    store = {m: [] for m in MECHANISMS}
    t0 = time.perf_counter()
    with Pool(workers) as pool:
        for i, (m, raw) in enumerate(pool.imap_unordered(simulate_chunk, jobs)):
            store[m].append(raw)
            if (i + 1) % 10 == 0:
                print(f"chunks {i+1}/{len(jobs)}  {time.perf_counter()-t0:.0f}s", flush=True)
    data = {}
    for k in RAW_KEYS + ("counts", "theta"):
        data[k] = np.concatenate([np.concatenate([c[k] for c in store[m]]) for m in MECHANISMS])
    data["label"] = np.concatenate([np.full(sum(len(c["theta"]) for c in store[m]), i)
                                    for i, m in enumerate(MECHANISMS)])
    if out:
        np.savez_compressed(out, **data)
    return data


def _safe_ratio(a, b):
    return np.where(b > 0, a / np.maximum(b, 1), np.nan)


def features(data, n=N_MAX, gaps=GAPS, arms="full", tier="counts"):
    """Summary statistics for a design. Replicates 0..n-1 are used in each arm."""
    keep = []
    for c, (arm, g) in enumerate(CELLS):
        if arm == "fresh_fresh":
            keep.append(c)
        elif g in gaps and (arms == "full" or (arms == "combined" and arm == "exposed_conditioned")
                            or (arms == "exposed" and arm == "exposed_fresh")):
            keep.append(c)
    if arms == "baseline":
        keep = [0]
    cnt = np.log1p(data["counts"][:, keep, :n][:, :, :, list(SAMPLE_TIMES)])
    occ = np.log1p(data["integrated_occupancy"][:, keep, :n])
    blocks = [cnt.mean(axis=2).reshape(len(cnt), -1),
              cnt[:, :, :, -1].std(axis=2),
              occ.mean(axis=2)]
    if len(keep) > 1:
        blocks.append(occ.mean(axis=2)[:, 1:] - occ.mean(axis=2)[:, :1])
    if tier == "tracking":
        half = Config().arrivals // 2
        g = lambda k: data[k][:, keep, :n]
        blocks += [g("attached_early").mean(axis=2) / half,
                   g("attached_late").mean(axis=2) / half,
                   np.nan_to_num(np.nanmean(_safe_ratio(g("founder_quick"), g("founder_departures")), axis=2), nan=-1),
                   np.nan_to_num(np.nanmean(_safe_ratio(g("descendant_quick"), g("descendant_departures")), axis=2), nan=-1),
                   np.log1p(g("divisions")).mean(axis=2),
                   g("initial_residual").mean(axis=2),
                   g("retained_residual").mean(axis=2)]
    return np.column_stack(blocks).astype(np.float32)


def matched_indices(data, rng, bins=8):
    """Equalize the fresh/fresh baseline distribution across hypotheses."""
    ff = np.log1p(data["counts"][:, 0, :, :])
    a = ff[:, :, 3].mean(axis=1)            # t = 4, end of influx
    b = ff[:, :, -1].mean(axis=1)           # t = 18
    qa = np.quantile(a, np.linspace(0, 1, bins + 1)[1:-1])
    qb = np.quantile(b, np.linspace(0, 1, bins + 1)[1:-1])
    cell = np.digitize(a, qa) * bins + np.digitize(b, qb)
    lab = data["label"]
    keep = []
    for c in np.unique(cell):
        members = [np.flatnonzero((cell == c) & (lab == k)) for k in range(len(MECHANISMS))]
        m = min(len(x) for x in members)
        for x in members:
            keep.extend(rng.choice(x, m, replace=False))
    return np.sort(np.array(keep, int))


def evaluate(X, y, train, test, seed=0, trees=400):
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier(n_estimators=trees, min_samples_leaf=2, max_features="sqrt",
                                 n_jobs=2, random_state=seed)
    clf.fit(X[train], y[train])
    prob = clf.predict_proba(X[test])
    pred = prob.argmax(axis=1)
    K = len(MECHANISMS)
    cm = np.zeros((K, K), int)
    np.add.at(cm, (y[test], pred), 1)
    acc = float((pred == y[test]).mean())
    # Binomial standard error of held-out accuracy
    se = float(np.sqrt(acc * (1 - acc) / len(test)))
    return {"accuracy": acc, "se": se, "confusion": cm.tolist(),
            "importance": clf.feature_importances_.tolist()}, prob


DESIGNS = ([{"tier": t, "n": n, "gaps": g, "arms": a}
            for t, n, g, a in product(("counts", "tracking"), (2, 4, 8),
                                      ((0.0,), (0.0, 16.0), (0.0, 4.0, 16.0)),
                                      ("full", "combined", "exposed"))]
           + [{"tier": t, "n": n, "gaps": (0.0,), "arms": "baseline"}
              for t, n in product(("counts", "tracking"), (2, 4, 8))])


def run(n_per_mechanism=4000, out=Path("results"), workers=2, reuse=False):
    out = Path(out)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    pool_path = out / "recovery_pool.npz"
    t0 = time.perf_counter()
    if reuse and pool_path.exists():
        data = dict(np.load(pool_path))
    else:
        data = simulate_pool(n_per_mechanism, workers=workers, out=pool_path)
    sim_seconds = time.perf_counter() - t0
    rng = np.random.default_rng(70001)
    idx = matched_indices(data, rng)
    perm = rng.permutation(idx)
    n_test = len(perm) // 4
    test, train = np.sort(perm[:n_test]), np.sort(perm[n_test:])
    y = data["label"]
    rows, details = [], {}
    for d in DESIGNS:
        X = features(data, d["n"], d["gaps"], d["arms"], d["tier"])
        res, prob = evaluate(X, y, train, test)
        key = f'{d["tier"]}|n{d["n"]}|gaps{"-".join(str(int(g)) for g in d["gaps"])}|{d["arms"]}'
        details[key] = {**d, **res, "n_features": X.shape[1]}
        rows.append({"tier": d["tier"], "replicates": d["n"],
                     "gaps": "-".join(str(int(g)) for g in d["gaps"]) if d["arms"] != "baseline" else "none",
                     "arms": d["arms"], "accuracy": res["accuracy"], "se": res["se"],
                     "n_features": X.shape[1]})
        print(key, round(res["accuracy"], 3), flush=True)
        if d["arms"] == "full" and d["n"] == 8 and len(d["gaps"]) == 3:
            np.save(out / "tables" / f'recovery_prob_{d["tier"]}.npy', prob)
    import csv
    with (out / "tables" / "recovery_accuracy.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    meta = {"n_per_mechanism_simulated": n_per_mechanism, "matched_experiments": int(len(idx)),
            "matched_per_mechanism": int(len(idx) // len(MECHANISMS)),
            "train": int(len(train)), "test": int(len(test)),
            "trajectories": int(len(y) * len(CELLS) * N_MAX), "cells": CELLS, "gaps": GAPS,
            "replicates_max": N_MAX, "priors": PRIORS, "simulation_seconds": sim_seconds,
            "test_labels": y[test].tolist(), "test_index": test.tolist()}
    (out / "logs" / "recovery_meta.json").write_text(json.dumps(meta, indent=1))
    (out / "logs" / "recovery_details.json").write_text(json.dumps(details, indent=1))
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=4000)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--output", type=Path, default=Path("results"))
    p.add_argument("--reuse", action="store_true")
    a = p.parse_args()
    run(a.n, a.output, a.workers, a.reuse)
