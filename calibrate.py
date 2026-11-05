"""Baseline and effect-size calibration for the illustrative factorial comparison.

1. Each history mechanism receives its own baseline attachment probability,
   detachment rate and division rate, chosen so its fresh-cell / fresh-surface
   mean occupancy trajectory matches the heterogeneous passive null. Without this step,
   history acquired during the challenge itself separates the mechanisms by
   baseline level alone.
2. Confounder strengths are then chosen so that, at zero gap, residual-cell
   carryover reproduces the surface-history transfer effect and adhesive
   selection reproduces the cell-history transfer effect.

Calibration uses common random numbers so objectives are smooth in parameters.
"""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, brentq
from model import Config, mechanism_params, simulate_batch, OBS_TIMES

CAL_N = 1000
STARTS = ((-0.9, 0.5, -0.92), (-0.5, 1.0, -1.39), (0.0, 1.5, -0.69))
CAL_SEED = 50021


def mean_curve(cfg, mechanism, exposed, conditioned, gap=0.0, overrides=None, n=CAL_N, seed=CAL_SEED):
    p = mechanism_params(cfg, mechanism)
    if overrides:
        p.update(overrides)
    r = simulate_batch(cfg, p, np.full(n, exposed), np.full(n, conditioned),
                       np.full(n, gap), seed=seed)
    return r["counts"].mean(axis=0), r["integrated_occupancy"].mean()


def _baseline(x):
    return {"attachment_probability": float(1 / (1 + np.exp(-x[0]))),
            "detachment_rate": float(np.exp(x[1])), "division_rate": float(np.exp(x[2]))}


def calibrate_baseline(cfg, mechanism, target):
    def loss(x):
        c, _ = mean_curve(cfg, mechanism, False, False, overrides=_baseline(x))
        return float(np.mean((np.log1p(c) - np.log1p(target)) ** 2))
    best, evals = None, 0
    for x0 in STARTS:
        res = minimize(loss, np.array(x0), method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-6, "maxiter": 100})
        evals += res.nfev
        if best is None or res.fun < best.fun:
            best = res
    return {**_baseline(best.x), "rms_log_error": float(np.sqrt(best.fun)), "evaluations": int(evals)}


def transfer_effect(cfg, mechanism, exposed, conditioned, overrides):
    _, a = mean_curve(cfg, mechanism, exposed, conditioned, overrides=overrides)
    _, b = mean_curve(cfg, mechanism, False, False, overrides=overrides)
    return a - b


def run(out=Path("results/logs/calibration.json"), cfg=None, reuse_baselines=False):
    cfg = cfg or Config()
    target, target_occ = mean_curve(cfg, "passive", False, False)
    if reuse_baselines:
        baselines = json.loads(out.read_text())["baselines"]
    else:
        baselines = {}
        for m in ("cell", "surface", "both"):
            baselines[m] = calibrate_baseline(cfg, m, target)
            print("baseline", m, baselines[m], flush=True)
    keys = ("attachment_probability", "detachment_rate", "division_rate")
    bo = lambda m: {k: baselines[m][k] for k in keys}
    cell_eff = transfer_effect(cfg, "cell", True, False, bo("cell"))
    surf_eff = transfer_effect(cfg, "surface", False, True, bo("surface"))
    sel = brentq(lambda s: transfer_effect(cfg, "selection", True, False,
                                           {"selection_shift": s}) - cell_eff, 0.01, 4.0, xtol=1e-3)
    car = brentq(lambda c: transfer_effect(cfg, "carryover", False, True,
                                           {"carryover_cells": c}) - surf_eff, 0.01, 40.0, xtol=1e-3)
    result = {"target_passive_fresh_curve": target.tolist(), "target_integrated_occupancy": float(target_occ),
              "obs_times": list(OBS_TIMES), "baselines": baselines,
              "cell_transfer_effect_gap0": float(cell_eff),
              "surface_transfer_effect_gap0": float(surf_eff),
              "selection_shift": float(sel), "carryover_cells": float(car),
              "n": CAL_N, "seed": CAL_SEED}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != "target_passive_fresh_curve"}, indent=2))
    return result


def calibrated_overrides(cal):
    """Per-mechanism parameter overrides implied by a calibration record."""
    o = {m: {} for m in ("passive", "carryover", "selection", "cell", "surface", "both")}
    for m, v in cal["baselines"].items():
        o[m] = {k: v[k] for k in ("attachment_probability", "detachment_rate", "division_rate")}
    o["selection"] = {"selection_shift": cal["selection_shift"]}
    o["carryover"] = {"carryover_cells": cal["carryover_cells"]}
    return o


if __name__ == "__main__":
    import sys
    run(reuse_baselines="--reuse-baselines" in sys.argv)
