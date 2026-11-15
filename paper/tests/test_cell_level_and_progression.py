"""Tests for matched progression, cell-level residence and repeated-contact design."""
import csv
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "paper" / "analysis"))

import matched_progression as mp  # noqa: E402
import repeated_contact_design as rc  # noqa: E402
import residence_composition as res  # noqa: E402
from model_core import normal_power_sample_size, paired_probabilities  # noqa: E402

RESULTS = ROOT / "paper" / "results"


def read(name):
    with (RESULTS / name).open(newline="") as f:
        return list(csv.DictReader(f))


class Progression(unittest.TestCase):
    def test_censoring_rules(self):
        row = {"record": "x", "record_end_h": 10., "t20_h": 2., "t20_status": "observed",
               "t30_h": None, "t30_status": "after_end"}
        self.assertEqual(mp.progression(row, 20, 30), (8., "lower_bound"))
        row["t20_status"] = "at_start"
        self.assertEqual(mp.progression(row, 20, 30)[1], "lower_threshold_at_start")

    def test_early_occupancy_ratio(self):
        t = np.arange(0, 10, .1)
        y = np.where(t < 2, 4., 1.)
        self.assertAlmostEqual(mp.early_occupancy_ratio(t, y), .25)

    def test_expansion_rate_recovers_exponential(self):
        t = np.arange(0, 20, 1 / 120)
        y = 10 * np.exp(.4 * t)
        rate, _ = mp.max_expansion_rate(t, y)
        self.assertAlmostEqual(rate, .4, places=6)

    def test_saved_decay_bound_verdicts(self):
        rows = read("decay_bound_by_readout.csv")
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(r["verdict"].startswith("violated") for r in rows))
        pub = rows[0]
        self.assertAlmostEqual(float(pub["mixture_lower_bound_r37"]), .9394 ** (37 / 14), places=3)

    def test_saved_dynamics_reproduce(self):
        saved = {r["record"]: r for r in read("occupancy_dynamics.csv")}
        new = {r["record"]: r for r in mp.dynamics_table(ROOT)}
        for rec in saved:
            self.assertAlmostEqual(float(saved[rec]["occupancy_ratio_4_8h"]), new[rec]["occupancy_ratio_4_8h"],
                                   places=4)
        fc1 = [float(saved[f"FC1-{s}"]["occupancy_ratio_4_8h"]) for s in (1, 2, 3)]
        fc2 = [float(saved[f"FC2-{s}"]["occupancy_ratio_4_8h"]) for s in (1, 2, 3)]
        self.assertTrue(max(fc1) < .05 < min(fc2))


class Residence(unittest.TestCase):
    def test_km_without_censoring_is_empirical(self):
        f = np.array([1, 2, 5, 101, 150, 30])
        self.assertAlmostEqual(res.km_visible_at(f, np.zeros(6, bool)), 2 / 6)

    def test_likelihood_normalizes(self):
        # Probabilities of all frame counts up to a large censoring point sum to one.
        frames = np.arange(2, 3001, dtype=float)
        D = {"f": np.append(frames, 3001.), "c": np.append(np.zeros(len(frames), bool), True),
             "w": np.ones(len(frames) + 1)}
        w, r = np.array([.7, .3]), np.array([1 / 30., 1 / 900.])
        def S(a):
            return np.exp(-np.outer(a, r)) @ w
        a_seen = (D["f"] - 1) * 3.
        p = np.where(D["c"], S(a_seen), S(a_seen) - S(D["f"] * 3.)) / S(np.array([3.]))[0]
        self.assertAlmostEqual(float(p.sum()), 1., places=10)

    def test_simulation_and_fit_recover_composition(self):
        rng = np.random.default_rng(3)
        lim = np.full(6000, 1200)
        r = np.array([1 / 30., 1 / 900.])
        data = {"FC1": res.compress(res.simulate_group(rng, [.75, .25], r, lim)),
                "FC2": res.compress(res.simulate_group(rng, [.6, .4], r, lim))}
        rw = res.fit_family(data, "composition", starts=3)
        sh = res.fit_family(data, "rate_shift", starts=3)
        self.assertGreater(rw["loglik"], sh["loglik"])
        self.assertAlmostEqual(rw["weights_FC2"][1], .4, delta=.05)

    def test_saved_counts_and_pooling_sensitivity(self):
        d = json.loads((RESULTS / "residence_composition.json").read_text())
        self.assertEqual(d["retention"]["FC1"]["arrivals"], 481)
        self.assertEqual(d["retention"]["FC1"]["attached"], 100)
        self.assertEqual(d["retention"]["FC1_with_first_hour_recording"]["arrivals"], 797)
        self.assertEqual(d["retention"]["FC2"]["arrivals"], 3048)
        m = d["models"]
        # Two-type preference reverses when the first-hour recording is removed.
        self.assertLess(m["matched"]["2_types"]["composition_minus_rate_shift_loglik"], 0)
        self.assertGreater(m["with_first_hour_recording"]["2_types"]["composition_minus_rate_shift_loglik"], 0)
        self.assertGreater(m["matched"]["3_types"]["composition_minus_rate_shift_loglik"], 0)

    def test_age_advance_reproduces_two_type_composition(self):
        d = json.loads((RESULTS / "residence_composition.json").read_text())
        eq = d["models"]["matched"]["equivalence_2_types"]
        self.assertLess(eq["max_abs_survival_difference"], 1e-12)
        m2 = d["models"]["matched"]["2_types"]
        self.assertAlmostEqual(m2["composition_minus_age_advance_loglik"], 0., places=4)

    def test_paired_sequences(self):
        d = json.loads((RESULTS / "residence_composition.json").read_text())
        diffs = [p["difference"] for p in d["matched_pairs"]["pairs"]]
        self.assertEqual(len(diffs), 3)
        lo, hi = d["matched_pairs"]["t95_interval_between_sequences"]
        self.assertLess(lo, 0.)
        self.assertGreater(hi, 0.)

    def test_threshold_formula(self):
        rows = res.threshold_rows(.05, .08, "test")
        self.assertAlmostEqual(rows[-1]["max_sham_loss_for_positive_expected_lower_bound"], .03 / 1.03)


class RepeatedContact(unittest.TestCase):
    def test_two_contacts_reproduce_article(self):
        for hold in (0., 14., 37.):
            f, _, sel, hom = paired_probabilities(hold, .10, .03)
            old = normal_power_sample_size(sel, hom)[0]
            mix, h = rc.count_distributions(2, rc.high_fraction(hold), .10, .03)
            self.assertAlmostEqual(rc.oracle_n80(mix, h) / old, 1., places=9)

    def test_more_contacts_need_fewer_cells(self):
        f = rc.high_fraction(37.)
        n = [rc.oracle_n80(*rc.count_distributions(k, f, .10, .03)) for k in (2, 3, 5, 10)]
        self.assertTrue(all(a > b for a, b in zip(n, n[1:])))

    def test_tarone_null_calibration(self):
        rng = np.random.default_rng(5)
        _, hom = rc.count_distributions(5, .5, .10, .03)
        z = rc.tarone_z(rng.multinomial(4000, hom, size=3000), 5)
        self.assertLess(abs(float(np.mean(z))), .1)


if __name__ == "__main__":
    unittest.main()
