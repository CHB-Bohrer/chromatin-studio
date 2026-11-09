import unittest
import numpy as np
from chromatin_studio.chromatin_contacts import compute_chromatin_contacts as compute


class ChromatinContactTests(unittest.TestCase):
    def setUp(self):
        self.x = np.array([[0., 0, 0], [9, 9, 9], [8, 8, 8], [.6, .2, .1]])
        self.marks = {"A": np.array([True, False, False, True]),
                      "B": np.array([True, False, False, True])}
        self.rules = [dict(a="A", b="A", strength=2., radius=1.5, min_separation=3)]

    def test_force_is_negative_energy_gradient_and_balanced(self):
        r = compute(self.x, self.marks, self.rules)
        for i in (0, 3):
            for axis in range(3):
                plus, minus = self.x.copy(), self.x.copy()
                plus[i, axis] += 1e-6
                minus[i, axis] -= 1e-6
                gradient = (compute(plus, self.marks, self.rules)["energy"] - compute(minus, self.marks, self.rules)["energy"]) / 2e-6
                self.assertAlmostEqual(r["forces"][i, axis], -gradient, places=7)
        np.testing.assert_allclose(r["forces"].sum(axis=0), 0)
        self.assertGreater(r["forces"][0, 0], 0)
        self.assertEqual(r["report"]["total_contacts"], 1)
        self.assertTrue(np.all(r["stiffness"] >= 0))

    def test_cutoff_coincidence_and_sequence_exclusion(self):
        for distance, expected in [(0, -2), (1.5, 0), (2, 0)]:
            self.x[3] = [distance, 0, 0]
            result = compute(self.x, self.marks, self.rules)
            self.assertAlmostEqual(result["energy"], expected)
            np.testing.assert_allclose(result["forces"], 0)
        self.x[3] = [.5, 0, 0]
        self.rules[0]["min_separation"] = 4
        self.assertEqual(compute(self.x, self.marks, self.rules)["energy"], 0)

    def test_deepest_current_well_wins_without_mark_double_counting(self):
        short = dict(a="B", b="B", strength=20., radius=.7, min_separation=3)
        combined = compute(self.x, self.marks, self.rules + [short])
        separate = [compute(self.x, self.marks, [r]) for r in self.rules + [short]]
        best = min(separate, key=lambda r: r["energy"])
        self.assertEqual(combined["energy"], best["energy"])
        np.testing.assert_allclose(combined["forces"], best["forces"])
        self.assertEqual(combined["report"]["total_contacts"], 1)
        self.assertEqual(combined["report"]["rules"][0]["contacts"], 1)

    def test_current_marks_and_symmetric_state_matching(self):
        self.marks["A"][:] = False
        self.marks["B"][:] = False
        rule = dict(a="unmodified", b="A", strength=1, radius=1, min_separation=3)
        self.assertEqual(compute(self.x, self.marks, [rule])["energy"], 0)
        self.marks["A"][0] = True
        self.assertLess(compute(self.x, self.marks, [rule])["energy"], 0)
        self.marks["B"][3] = True
        self.assertEqual(compute(self.x, self.marks, [rule])["energy"], 0)

    def test_empty_and_disabled_rules_are_zero(self):
        for rules in ([], [dict(self.rules[0], strength=0)]):
            result = compute(self.x, self.marks, rules)
            np.testing.assert_allclose(result["forces"], 0)
            np.testing.assert_allclose(result["stiffness"], 0)
            self.assertEqual(result["energy"], 0)


if __name__ == "__main__":
    unittest.main()
