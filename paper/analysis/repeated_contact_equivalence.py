"""Repeated-contact identifiability: stable mixture vs posterior-predictive history.

The history-conditioned probabilities are obtained by Bayes factorization of the
same stable two-phenotype mixture. Numerical agreement is an algebraic check,
not independent evidence of contact-induced physiology.

Usage (repository root):
    python paper/analysis/repeated_contact_equivalence.py
    python paper/analysis/repeated_contact_equivalence.py --out paper/results/repeated_contact_equivalence.csv
"""
import argparse
import csv
from itertools import product
from pathlib import Path

P_HIGH = 0.10
P_LOW = 0.03
# Rounded high-type fractions used in the main manuscript's illustrative model.
FRACTIONS = {'immediate': 0.99036, '14h': 0.93040, '37h': 0.31910}


def sequence_probability_static(bits, high_fraction, p_high=P_HIGH, p_low=P_LOW):
    successes = sum(bits)
    failures = len(bits) - successes
    return (high_fraction * p_high**successes * (1 - p_high)**failures
            + (1 - high_fraction) * p_low**successes * (1 - p_low)**failures)


def conditional_probability_history(bits, high_fraction, p_high=P_HIGH, p_low=P_LOW):
    """Posterior predictive P(next=1 | history), not an independently fitted mechanism."""
    successes = sum(bits)
    failures = len(bits) - successes
    numerator = (high_fraction * p_high**(successes + 1) * (1 - p_high)**failures
                 + (1 - high_fraction) * p_low**(successes + 1) * (1 - p_low)**failures)
    denominator = sequence_probability_static(bits, high_fraction, p_high, p_low)
    return numerator / denominator


def sequence_probability_history(bits, high_fraction, p_high=P_HIGH, p_low=P_LOW):
    joint = 1.0
    for position, outcome in enumerate(bits):
        p = conditional_probability_history(bits[:position], high_fraction, p_high, p_low)
        joint *= p if outcome else (1 - p)
    return joint


def result_rows():
    for hold, f in FRACTIONS.items():
        mean = f * P_HIGH + (1 - f) * P_LOW
        covariance = f * (1 - f) * (P_HIGH - P_LOW)**2
        p_after_success = conditional_probability_history((1,), f)
        p_after_failure = conditional_probability_history((0,), f)
        p_after_00 = conditional_probability_history((0, 0), f)
        p_after_10 = conditional_probability_history((1, 0), f)
        for k in (2, 5, 10):
            worst = max(abs(sequence_probability_static(bits, f)
                            - sequence_probability_history(bits, f))
                        for bits in product((0, 1), repeat=k))
            yield {'hold': hold, 'high_fraction': f, 'p_high': P_HIGH,
                   'p_low': P_LOW, 'mean_p': mean, 'cov_y1_y2': covariance,
                   'p_y2_given_y1_1': p_after_success,
                   'p_y2_given_y1_0': p_after_failure,
                   'p_y3_given_00': p_after_00,
                   'p_y3_given_10': p_after_10,
                   'contacts_per_cell': k, 'max_abs_sequence_difference': worst}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=None,
                        help='Optional CSV of exact-equivalence diagnostics')
    args = parser.parse_args()
    rows = list(result_rows())
    for row in rows:
        print(f"{row['hold']:>9} k={row['contacts_per_cell']:>2} "
              f"mean={row['mean_p']:.6f} cov={row['cov_y1_y2']:.8g} "
              f"max|delta|={row['max_abs_sequence_difference']:.3g}")
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == '__main__':
    main()
