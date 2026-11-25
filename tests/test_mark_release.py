"""Optional mark-dependent release (proteins[].mark_off_rate_multiplier)."""
import copy
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.schema import ConfigError, checked_config, default_config


def model(multiplier=None):
    config = default_config()
    config.update(marks=[{"id": "K9", "name": "H3K9me3", "color": "#aa5566", "decay_rate": 0,
                          "spread_rate": 0, "contact_radius": 1.0},
                         {"id": "K27", "name": "H3K27me3", "color": "#5566aa", "decay_rate": 0,
                          "spread_rate": 0, "contact_radius": 1.0}],
                  interactions=[], bodies=[], genes=[])
    config["cohesin"] = {"enabled": False, "count": 0, "step_rate": 0, "off_rate": 0, "spring_k": 0}
    config["chromatin"].update(beads=12, barriers=[], initial_marks=[
        {"mark": "K9", "start": 0, "end": 3}, {"mark": "K27", "start": 2, "end": 5}])
    reader = {"id": "reader", "name": "Reader", "color": "#d6a34a", "count": 2, "diffusion": 0.5, "radius": 0.28,
              "valence": 1, "representation": "bead", "idr_beads": 4, "sticker_fraction": 1.0,
              "chromatin_affinity": {"unmodified": 1.0}, "on_rate": 0.0, "off_rate": 0.4,
              "writer_mark": None, "eraser_mark": None, "catalytic_rate": 0.0}
    if multiplier is not None:
        reader["mark_off_rate_multiplier"] = multiplier
    config["proteins"] = [reader]
    return checked_config(config)


class MarkReleaseTests(unittest.TestCase):
    def test_absent_field_keeps_every_release_rate(self):
        sim = Simulation(model())
        np.testing.assert_allclose(sim._chromatin_off_rates("reader"), 0.4)

    def test_listed_marks_set_factors_and_overlap_takes_the_smallest(self):
        sim = Simulation(model({"K9": 0.25, "K27": 0.5, "unmodified": 2.0}))
        rates = sim._chromatin_off_rates("reader")
        np.testing.assert_allclose(rates[:2], 0.1)        # K9 only
        np.testing.assert_allclose(rates[2:4], 0.1)       # K9 and K27: most stabilizing factor
        np.testing.assert_allclose(rates[4:6], 0.2)       # K27 only
        np.testing.assert_allclose(rates[6:], 0.8)        # no listed mark
        for bead in range(12):
            self.assertAlmostEqual(0.4 * sim._mark_release_factor("reader", bead), rates[bead])

    def test_bond_release_follows_the_current_mark(self):
        sim = Simulation(model({"K9": 0.25}))
        particle = sim.protein_particles[0]
        bond = {"a": 0, "b": particle["index"], "kind": "chromatin", "off_rate": 0.4, "k": 12.0, "rest": 0.5}
        self.assertAlmostEqual(sim._bond_off_rate(bond), 0.1)
        sim.marks["K9"][0] = False                         # mark removed while bound
        self.assertAlmostEqual(sim._bond_off_rate(bond), 0.4)

    def test_field_adds_no_random_draws(self):
        first, second = Simulation(model()), Simulation(model({"unmodified": 1.0}))
        first.step(40)
        second.step(40)
        np.testing.assert_array_equal(first.positions, second.positions)

    def test_schema_rejects_unknown_marks_and_out_of_range_factors(self):
        for bad in ({"H4K16ac": 0.5}, {"K9": -1}, {"K9": 1001}):
            config = copy.deepcopy(model())
            config["proteins"][0]["mark_off_rate_multiplier"] = bad
            with self.assertRaises(ConfigError):
                checked_config(config)


if __name__ == "__main__":
    unittest.main()
