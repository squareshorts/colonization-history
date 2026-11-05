import unittest
from dataclasses import replace
import numpy as np
from model import (Config, MECHANISMS, ARMS, conditioning_history, mechanism_params,
                   simulate, simulate_batch)


class Tests(unittest.TestCase):
    def setUp(self):
        self.c = replace(Config(), horizon=6., arrivals=12, arrival_duration=2.)

    def test_accounting_and_capacity(self):
        for m in MECHANISMS:
            r = simulate(self.c, m, "exposed_conditioned", n=16, seed=99)
            np.testing.assert_array_equal(
                r["retained_founders"] + r["retained_descendants"] + r["retained_residual"],
                r["final_attached"])
            np.testing.assert_array_equal(
                r["initial_residual"] + r["arrivals"] + r["divisions"],
                r["final_attached"] + r["rejected"] + r["departures"])
            self.assertTrue(np.all(r["final_attached"] <= self.c.patches * self.c.capacity_per_patch))

    def test_passive_invariant_to_history(self):
        a = simulate(self.c, "passive", ARMS[0], n=15, seed=23)
        for arm in ARMS[1:]:
            b = simulate(self.c, "passive", arm, gap=8, n=15, seed=23)
            for k in ("final_attached", "established", "divisions", "rejected", "departures"):
                np.testing.assert_array_equal(a[k], b[k])

    def test_cell_only_invariant_to_surface_preparation(self):
        for m in ("cell", "selection"):
            a = simulate(self.c, m, "exposed_fresh", n=15, seed=25)
            b = simulate(self.c, m, "exposed_conditioned", n=15, seed=25)
            np.testing.assert_array_equal(a["final_attached"], b["final_attached"])

    def test_surface_only_invariant_to_cell_preparation(self):
        for m in ("surface", "carryover"):
            a = simulate(self.c, m, "fresh_conditioned", n=15, seed=26)
            b = simulate(self.c, m, "exposed_conditioned", n=15, seed=26)
            np.testing.assert_array_equal(a["final_attached"], b["final_attached"])

    def test_confounders_equal_null_on_untreated_factor(self):
        p = simulate(self.c, "passive", "fresh_fresh", n=15, seed=29)
        for m, arm in (("carryover", "exposed_fresh"), ("selection", "fresh_conditioned")):
            r = simulate(self.c, m, arm, n=15, seed=29)
            np.testing.assert_array_equal(p["final_attached"], r["final_attached"])

    def test_no_arrivals_no_population_without_carryover(self):
        c = replace(self.c, arrivals=0)
        for m in ("passive", "selection", "cell", "surface", "both"):
            r = simulate(c, m, "exposed_conditioned", n=10)
            self.assertTrue(np.all(r["final_attached"] == 0))
            self.assertTrue(np.all(r["divisions"] == 0))
        r = simulate(replace(c, carryover_cells=30.), "carryover", "fresh_conditioned", n=10)
        self.assertGreater(r["initial_residual"].sum(), 0)
        self.assertTrue(np.all(r["retained_founders"] == 0))

    def test_no_division_no_descendants(self):
        c = replace(self.c, division_rate=0)
        r = simulate(c, "both", "exposed_conditioned", n=10)
        self.assertTrue(np.all(r["retained_descendants"] == 0))

    def test_long_gap_erases_preparations(self):
        for m in ("cell", "surface", "both", "selection", "carryover"):
            a = simulate(self.c, m, "exposed_conditioned", gap=1e6, n=15, seed=27)
            b = simulate(self.c, m, "fresh_fresh", n=15, seed=27)
            np.testing.assert_array_equal(a["final_attached"], b["final_attached"])

    def test_zero_couplings_collapse_to_null(self):
        c = replace(self.c, cell_attachment_effect=0, cell_retention_effect=0,
                    surface_attachment_effect=0, surface_retention_effect=0)
        a = simulate(c, "passive", "fresh_fresh", n=15, seed=28)
        b = simulate(c, "both", "exposed_conditioned", n=15, seed=28)
        np.testing.assert_array_equal(a["final_attached"], b["final_attached"])

    def test_prescribed_history_leaves_zero_cells(self):
        h = conditioning_history(mechanism_params(self.c, "both"))
        self.assertEqual(h["attached_cells"], 0)
        self.assertGreater(float(h["cell_memory"]), 0)
        self.assertGreater(float(h["surface_trace"]), 0)

    def test_vectorized_conditioning_matches_scalar(self):
        p = mechanism_params(self.c, "both")
        taus = np.array([2., 4., 9.])
        vec = conditioning_history({**p, "cell_memory_decay_time": taus})["cell_memory"]
        for k, tau in enumerate(taus):
            s = conditioning_history({**p, "cell_memory_decay_time": tau})["cell_memory"]
            self.assertAlmostEqual(float(vec[k]), float(s), places=12)

    def test_residual_count_is_poisson(self):
        c = replace(self.c, carryover_cells=3., arrivals=0, horizon=0.05, dt=0.05,
                    arrival_duration=0.05, division_rate=0., detachment_rate=0.)
        r = simulate(c, "carryover", "fresh_conditioned", n=4000, seed=3)
        x = r["initial_residual"]
        self.assertLess(abs(x.mean() - 3.), 0.12)
        self.assertLess(abs(x.var() - 3.), 0.4)

    def test_per_replicate_parameters_are_independent(self):
        n = 20
        p = mechanism_params(self.c, "both")
        p["division_rate"] = np.where(np.arange(n) < 10, 0., 0.9)
        r = simulate_batch(self.c, p, np.ones(n, bool), np.ones(n, bool), np.zeros(n), seed=5)
        self.assertTrue(np.all(r["divisions"][:10] == 0))
        self.assertGreater(r["divisions"][10:].sum(), 0)

    def test_overrides_cannot_activate_inactive_couplings(self):
        o = {"cell_attachment_effect": 3., "surface_retention_effect": 3., "selection_shift": 2.,
             "carryover_cells": 20.}
        a = simulate(self.c, "passive", "fresh_fresh", n=15, seed=33)
        b = simulate(self.c, "passive", "exposed_conditioned", n=15, seed=33, overrides=o)
        np.testing.assert_array_equal(a["final_attached"], b["final_attached"])

    def test_zero_persistence_still_requires_population(self):
        c = replace(self.c, arrivals=0, establishment_persistence=0)
        r = simulate(c, "passive", "fresh_fresh", n=5)
        self.assertTrue(np.all(~r["established"]))

    def test_equal_incoming_dose(self):
        for m in MECHANISMS:
            for arm in ARMS:
                r = simulate(self.c, m, arm, n=3, seed=10)
                self.assertTrue(np.all(r["arrivals"] == self.c.arrivals))

    def test_reproducibility(self):
        a = simulate(self.c, "both", "exposed_conditioned", n=12, seed=31)
        b = simulate(self.c, "both", "exposed_conditioned", n=12, seed=31)
        np.testing.assert_array_equal(a["establishment_time"], b["establishment_time"])
        np.testing.assert_array_equal(a["counts"], b["counts"])

    def test_threshold_is_observation_not_absorbing_state(self):
        c = replace(self.c, establishment_threshold=1, establishment_persistence=0,
                    detachment_rate=4., division_rate=0, arrivals=50, arrival_duration=1.)
        r = simulate(c, "passive", "fresh_fresh", n=15, seed=32)
        self.assertTrue(np.any(r["established"] & (r["final_attached"] == 0)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
