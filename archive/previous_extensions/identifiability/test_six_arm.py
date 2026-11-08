import unittest
import random
from six_arm import ordered_contrasts

class SixArmTests(unittest.TestCase):
    def test_telescoping_identity_for_nonlinear_values(self):
        rng = random.Random(20261107)
        for _ in range(1000):
            f = {i: rng.gauss(0, 3) for i in range(16)}
            d = ordered_contrasts(f)
            self.assertAlmostEqual(d['total'], f[15]-f[0], places=10)
    def test_requires_interventions(self):
        with self.assertRaises(ValueError):
            ordered_contrasts({0:0,3:1,12:2,15:3})

if __name__ == '__main__': unittest.main()
