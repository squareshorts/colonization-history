"""Prospective power: arrivals needed to separate an adhesiveness offset from a head start.

Residence times are simulated from the simulator's residence law (individual hazard
lambda exp(-z) exp(-gamma m(a)), lognormal z, acquisition m(a)) at the parameters fitted
to the deposited experiments, for two groups (fresh, exposed) with N attached arrivals
each, observed at 3-s frames within 1-h windows. Both families are refitted to every
simulated dataset; power is the fraction of datasets in which the generating family has
the higher likelihood, and the fraction in which the LR exceeds 3 log units.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as A

FRAME_S, WINDOW_S = 3.0, 3600.0


def simulate(n, lam, sig, gam, k, mu, m0, rng):
    out_f, out_c = [], []
    while sum(len(x) for x in out_f) < n:
        m = 4 * n
        z = mu + sig * rng.normal(size=m)
        E = rng.exponential(size=m)
        target = E * np.exp(z)                       # need cum_hazard(T) = target
        lo, hi = np.zeros(m), np.full(m, 1e7)
        for _ in range(60):                         # vectorised bisection
            mid = 0.5 * (lo + hi)
            c = gam * (1 - m0)
            if c < 1e-8:
                L = lam * np.exp(-gam) * mid
            else:
                L = lam * np.exp(-gam) / k * (A.expi(c) - A.expi(c * np.exp(-k * mid)))
            hi = np.where(L > target, mid, hi)
            lo = np.where(L > target, lo, mid)
        T = 0.5 * (lo + hi)
        win = WINDOW_S - rng.uniform(0, WINDOW_S, m)
        frames = np.floor(np.minimum(T, win) / FRAME_S).astype(int) + 1
        cens = T > win
        keep = frames >= A.MIN_FRAMES
        out_f.append(frames[keep]); out_c.append(cens[keep])
    f, c = np.concatenate(out_f)[:n], np.concatenate(out_c)[:n]
    return f, c


def run(Ns=(250, 1000, 3000), reps=20, seed=3):
    E = json.loads((A.ROOT.parent / "results" / "empirical" / "empirical_fits.json").read_text())
    off, head = E["families"]["selection"], E["families"]["history"]
    truths = {
        "offset": dict(lam=off["lambda_per_s"], sig=off["sigma"], gam=off["gamma"], k=off["acquisition_rate_per_s"],
                       mu=off["adhesiveness_shift"]["FC2"], m0=0.0),
        "head": dict(lam=head["lambda_per_s"], sig=head["sigma"], gam=head["gamma"], k=head["acquisition_rate_per_s"],
                     mu=0.0, m0=head["history_head_start"]["FC2"]),
    }
    rng = np.random.default_rng(seed)
    res = []
    for truth, p in truths.items():
        for n in Ns:
            wins, strong = 0, 0
            for r in range(reps):
                fc1 = simulate(n, p["lam"], p["sig"], p["gam"], p["k"], 0.0, 0.0, rng)
                fc2 = simulate(n, p["lam"], p["sig"], p["gam"], p["k"], p["mu"], p["m0"], rng)
                data = {"FC1": fc1, "FC2": fc2}
                o, _ = A.fit_family(data, ["FC1", "FC2"], "selection", starts=2, seed=r)
                h, _ = A.fit_family(data, ["FC1", "FC2"], "history", starts=2, seed=r)
                d = (o["loglik"] - h["loglik"]) * (1 if truth == "offset" else -1)
                wins += d > 0; strong += d > 3
            res.append({"truth": truth, "attached_per_group": n, "reps": reps,
                        "power_correct": wins / reps, "power_LR_gt_3": strong / reps})
            print(res[-1], flush=True)
    out = A.ROOT.parent / "results" / "empirical_v2"
    out.mkdir(parents=True, exist_ok=True)
    (out / "design_power.json").write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    run()
