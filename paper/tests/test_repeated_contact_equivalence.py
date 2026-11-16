"""Exact-sequence equivalence and restricted-history counterexample."""
import sys
import unittest
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
import repeated_contact_equivalence as eq


class RepeatedContactEquivalence(unittest.TestCase):
    def test_exact_sequences_and_normalization(self):
        for f in eq.FRACTIONS.values():
            for k in (2, 5, 10):
                bits = list(product((0, 1), repeat=k))
                p = [eq.sequence_probability_static(row, f) for row in bits]
                q = [eq.sequence_probability_history(row, f) for row in bits]
                self.assertAlmostEqual(sum(p), 1.0, places=12)
                self.assertLess(max(abs(a - b) for a, b in zip(p, q)), 1e-12)

    def test_within_cell_covariance(self):
        f = eq.FRACTIONS['37h']
        p11 = eq.sequence_probability_static((1, 1), f)
        mean = f * eq.P_HIGH + (1 - f) * eq.P_LOW
        self.assertAlmostEqual(p11 - mean**2,
                               f * (1 - f) * (eq.P_HIGH - eq.P_LOW)**2,
                               places=13)

    def test_three_contact_first_order_restriction(self):
        f = eq.FRACTIONS['37h']
        p00 = eq.conditional_probability_history((0, 0), f)
        p10 = eq.conditional_probability_history((1, 0), f)
        self.assertAlmostEqual(p00, .05012, delta=0.00005)
        self.assertAlmostEqual(p10, .07142, delta=0.00005)
        self.assertGreater(p10 - p00, .02)
        # Both histories have the same last outcome; a first-order process
        # cannot assign two distinct conditional probabilities to them.


if __name__ == '__main__':
    unittest.main()
