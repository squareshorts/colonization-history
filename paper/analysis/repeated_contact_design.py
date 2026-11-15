"""Practical repeated-contact designs: more contacts per cell and estimated parameters.

A persistent two-type population (fraction f of high-adhesion cells with
attachment probability p_H, otherwise p_L) gives each cell a binomial mixture
for its number of successes in k independent contacts. A homogeneous
population with the same mean gives a single binomial. The number of successes
per cell is sufficient for both, so k = 2 reproduces the ordered-outcome oracle
calculation of model_core.py.

Two tests are evaluated:
  oracle  likelihood ratio with known parameters (normal approximation, as in
          the article);
  Tarone  score test for binomial overdispersion (Tarone 1979) with the
          attachment probability estimated from the data; its 5% critical value
          is calibrated by Monte Carlo under the homogeneous null.

Imaging effort uses the arrival rate in first-exposure fields during the first
three hours after inoculation (797 arrivals in 10 field-hours).

Usage:
  python paper/analysis/repeated_contact_design.py --out paper/results
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.special import expit, logit
from scipy.stats import binom, norm

F0, C = .990359, .145675
ARRIVALS_PER_FIELD_HOUR = 797 / 10.


def high_fraction(hold_h):
    return float(expit(logit(F0) - C * hold_h))


def count_distributions(k, f, p_hi, p_lo):
    j = np.arange(k + 1)
    mix = f * binom.pmf(j, k, p_hi) + (1 - f) * binom.pmf(j, k, p_lo)
    p = f * p_hi + (1 - f) * p_lo
    hom = binom.pmf(j, k, p)
    return mix, hom


def oracle_n80(mix, hom, alpha=.05, power=.8):
    keep = (mix > 0) & (hom > 0)
    ll = np.log(mix[keep] / hom[keep])
    a0, a1 = float(hom[keep] @ ll), float(mix[keep] @ ll)
    v0 = float(hom[keep] @ (ll - a0) ** 2)
    v1 = float(mix[keep] @ (ll - a1) ** 2)
    return float(((norm.ppf(1 - alpha) * np.sqrt(v0) + norm.ppf(power) * np.sqrt(v1)) / (a1 - a0)) ** 2)


def tarone_z(counts, k):
    """Tarone's overdispersion statistic from per-cell success counts tabulated by j."""
    counts = np.atleast_2d(counts)
    j = np.arange(k + 1)
    n = counts.sum(axis=1)
    phat = (counts @ j) / (n * k)
    phat = np.clip(phat, 1e-12, 1 - 1e-12)
    S = np.sum(counts * (j[None, :] - k * phat[:, None]) ** 2, axis=1) / (phat * (1 - phat))
    return (S - n * k) / np.sqrt(2 * n * k * (k - 1))


def tarone_power(rng, k, n, mix, hom, studies=4000):
    z0 = tarone_z(rng.multinomial(n, hom, size=studies), k)
    crit = float(np.quantile(z0, .95))
    z1 = tarone_z(rng.multinomial(n, mix, size=studies), k)
    return float(np.mean(z1 > crit)), crit


def tarone_n80(rng, k, mix, hom, guess, studies=4000):
    """Sample size for 80% power from a probit fit of Monte Carlo power on sqrt(n).

    Power is estimated on a geometric grid from guess/3 to 3*guess, where guess
    is the oracle sample size, and the 0.8 crossing of the fitted probit curve
    is returned with the grid used.
    """
    grid = np.unique(np.geomspace(max(guess / 2.5, 20), guess * 2.5, 15).astype(int))
    pw = np.array([tarone_power(rng, k, int(n), mix, hom, studies=studies)[0] for n in grid])
    ok = (pw > .03) & (pw < .99)
    zz = norm.ppf(pw[ok])
    b, a = np.polyfit(np.sqrt(grid[ok]), zz, 1)
    return int(round(((norm.ppf(.8) - a) / b) ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=20261116)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rows = []
    for p_hi, p_lo in ((.10, .03), (.15, .04)):
        for hold in (0., 14., 37.):
            f = high_fraction(hold)
            for k in (2, 3, 5, 10):
                mix, hom = count_distributions(k, f, p_hi, p_lo)
                n_or = oracle_n80(mix, hom)
                n_ta = tarone_n80(rng, k, mix, hom, n_or)
                check = tarone_power(rng, k, n_ta, mix, hom, studies=6000)[0] if n_ta else None
                contacts = None if n_ta is None else n_ta * k
                rows.append({"p_high": p_hi, "p_low": p_lo, "hold_h": hold, "high_fraction": round(f, 6),
                             "contacts_per_cell": k, "oracle_n80": round(n_or),
                             "tarone_n80": n_ta, "tarone_power_check": check,
                             "tarone_total_contacts": contacts,
                             "field_hours_at_archive_arrival_rate":
                                 None if contacts is None else round(contacts / ARRIVALS_PER_FIELD_HOUR, 1)})
                print(rows[-1])
    path = a.out / "repeated_contact_design.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
