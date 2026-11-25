"""Optional per-bead limit on protein bonds (chromatin.binding_capacity)."""
import copy
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.schema import ConfigError, checked_config, default_config


def model(capacity=None, count=30):
    config = default_config()
    config.update(marks=[], interactions=[], bodies=[], genes=[])
    config["cohesin"] = {"enabled": False, "count": 0, "step_rate": 0, "off_rate": 0, "spring_k": 0}
    config["simulation"].update(confinement_radius=3.0, dt=0.005)
    config["chromatin"].update(beads=8, barriers=[], initial_marks=[])
    if capacity is not None:
        config["chromatin"]["binding_capacity"] = capacity
    config["proteins"] = [{"id": "binder", "name": "Binder", "color": "#d6a34a", "count": count, "diffusion": 0.5,
                           "radius": 0.28, "valence": 1, "representation": "bead", "idr_beads": 4,
                           "sticker_fraction": 1.0, "chromatin_affinity": {"unmodified": 1.0}, "on_rate": 20.0,
                           "off_rate": 0.01, "writer_mark": None, "eraser_mark": None, "catalytic_rate": 0.0}]
    return checked_config(config)


def bonds_per_bead(sim):
    counts = np.zeros(sim.n_chromatin, dtype=int)
    for bond in sim.bindings:
        if bond["kind"] == "chromatin":
            counts[bond["a"]] += 1
    return counts


class BindingCapacityTests(unittest.TestCase):
    def test_absent_capacity_allows_several_bonds_per_bead(self):
        sim = Simulation(model())
        sim.step(200)
        self.assertGreater(bonds_per_bead(sim).max(), 1)

    def test_capacity_is_never_exceeded(self):
        for capacity in (1, 2):
            sim = Simulation(model(capacity))
            for _ in range(20):
                sim.step(10)
                self.assertLessEqual(bonds_per_bead(sim).max(), capacity)
            self.assertEqual(bonds_per_bead(sim).max(), capacity)      # saturating conditions

    def test_large_capacity_reproduces_the_unlimited_trajectory(self):
        first, second = Simulation(model()), Simulation(model(64))
        first.step(60)
        second.step(60)
        np.testing.assert_array_equal(first.positions, second.positions)
        self.assertEqual(len(first.bindings), len(second.bindings))

    def test_schema_rejects_invalid_capacity(self):
        for bad in (0, 1.5, 65, "one"):
            config = copy.deepcopy(model())
            config["chromatin"]["binding_capacity"] = bad
            with self.assertRaises(ConfigError):
                checked_config(config)


if __name__ == "__main__":
    unittest.main()
