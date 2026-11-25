"""Whole-nucleus envelope, tethering elements and their occupancy dependence."""

import copy
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.nucleus import NucleusEnvelope
from chromatin_studio.openmm_native import unsupported_features
from chromatin_studio.schema import default_config, validate_config
from chromatin_studio.studies import _ablate
from tests.test_engine import small_config
from tests.test_nucleosomes import nucleosomes


def tether(**kwargs):
    value = {"id": "lad", "name": "LAD", "color": "#7a8fa6", "start": 0, "end": 7, "enabled": True,
             "anchor": "envelope", "attach_rate": 5.0, "detach_rate": 0.0, "capture_distance": 50.0,
             "spring_k": 2.0, "rest_length": 0.5}
    value.update(kwargs)
    return value


def nucleus_config(**kwargs):
    config = small_config()
    config["nucleus"] = {"enabled": True, "center": [-20.0, 0.0, 0.0], "radii": [24.0, 30.0, 30.0],
                         "envelope_stiffness": 30.0, "tethers": [tether()]}
    config["nucleus"].update(kwargs)
    return config


class EnvelopeTests(unittest.TestCase):
    def test_sphere_reduces_to_reference_wall(self):
        wall = NucleusEnvelope({"radii": [5, 5, 5], "center": [1, 0, 0], "envelope_stiffness": 7})
        points = np.array([[8.0, 0, 0], [1, 6, 0], [1, 0, 1]])
        force, energy = wall.force_energy(points)
        self.assertAlmostEqual(energy, 0.5 * 7 * (2 ** 2 + 1 ** 2))
        np.testing.assert_allclose(force[0], [-14, 0, 0])
        np.testing.assert_allclose(force[2], [0, 0, 0])

    def test_ellipsoid_force_is_exact_gradient(self):
        wall = NucleusEnvelope({"radii": [3, 5, 2], "center": [0.5, -1, 0], "envelope_stiffness": 11})
        rng = np.random.default_rng(3)
        points = rng.normal(scale=5, size=(40, 3))
        force, _ = wall.force_energy(points)
        h = 1e-6
        for i in range(len(points)):
            for axis in range(3):
                plus, minus = points[i].copy(), points[i].copy()
                plus[axis] += h
                minus[axis] -= h
                numeric = -(wall.force_energy(plus)[1] - wall.force_energy(minus)[1]) / (2 * h)
                self.assertAlmostEqual(force[i, axis], numeric, places=4)
        self.assertTrue(np.allclose(wall.rho(wall.surface_point(points)), 1))

    def test_simulation_force_matches_energy_with_tethers(self):
        sim = Simulation(nucleus_config())
        sim.step(30)
        self.assertTrue(sim.tethers)
        x0 = sim.positions.copy()
        force = sim.forces()
        h = 1e-6
        for bead in range(sim.n_chromatin):
            for axis in range(3):
                sim.positions[:] = x0
                sim.positions[bead, axis] += h
                up = sim.potential_energy()
                sim.positions[:] = x0
                sim.positions[bead, axis] -= h
                down = sim.potential_energy()
                self.assertAlmostEqual(force[bead, axis], -(up - down) / (2 * h), places=4)
        sim.positions[:] = x0


class TetherTests(unittest.TestCase):
    def test_attach_detach_events_and_snapshot(self):
        config = nucleus_config()
        config["nucleus"]["tethers"] = [tether(detach_rate=2.0)]
        sim = Simulation(config)
        sim.step(400)
        types = {e["type"] for e in sim.genomic_events}
        self.assertTrue({"tether_attach", "tether_detach"} <= types)
        snap = sim.snapshot()
        self.assertTrue(snap["nucleus"]["enabled"])
        self.assertEqual(len(snap["nucleus"]["attachments"]), len(sim.tethers))
        self.assertEqual(sum(b["kind"] == "tether" for b in snap["bonds"]), len(sim.tethers))
        attached = sum((b - a) * v for a, b, v in snap["genome_tracks"]["tether_runs"])
        self.assertEqual(attached, len(sim.tethers))

    def test_snapshot_is_read_only(self):
        sim = Simulation(nucleus_config())
        sim.step(20)
        state = sim.rng.bit_generator.state
        before = copy.deepcopy(sim.tethers)
        sim.snapshot()
        sim.nucleus_state()
        self.assertEqual(state, sim.rng.bit_generator.state)
        self.assertEqual(before, sim.tethers)

    def test_out_of_capture_range_never_attaches(self):
        config = nucleus_config()
        config["nucleus"]["tethers"] = [tether(capture_distance=0.1)]
        config["nucleus"]["radii"] = [200, 200, 200]
        config["nucleus"]["center"] = [0, 0, 0]
        sim = Simulation(config)
        sim.step(100)
        self.assertEqual(sim.tethers, [])

    def test_attachment_hazard_depends_on_occupancy_and_mark(self):
        config = nucleus_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(initial_regions=[{"start": 0, "end": 3, "occupied": False}])
        config["nucleus"]["tethers"] = [tether(occupancy_attach_multiplier={"occupied": 1, "vacant": 0},
                                               occupancy_detach_multiplier={"occupied": 1, "vacant": 7})]
        sim = Simulation(config)
        spec = sim.tether_specs[0]
        self.assertEqual(sim._tether_hazard(spec, 0, "attach"), 0)
        self.assertEqual(sim._tether_hazard(spec, 5, "attach"), 5)
        self.assertEqual(sim._tether_hazard(spec, 0, "detach"), 0)
        spec["detach_rate"] = 1.0
        self.assertEqual(sim._tether_hazard(spec, 1, "detach"), 7)
        sim.step(200)
        self.assertTrue(all(t["bead"] >= 4 for t in sim.tethers))
        spec["required_mark"] = None

    def test_disabled_nucleus_preserves_legacy_trajectory(self):
        legacy = Simulation(small_config())
        config = small_config()
        config["nucleus"] = {"enabled": False, "radii": [10, 10, 10], "tethers": [tether()]}
        disabled = Simulation(config)
        legacy.step(50)
        disabled.step(50)
        np.testing.assert_array_equal(legacy.positions, disabled.positions)

    def test_ablations_and_native_refusal(self):
        config = nucleus_config()
        self.assertFalse(_ablate(config, "no_tethering")["nucleus"]["tethers"][0]["enabled"])
        self.assertFalse(_ablate(config, "no_nucleus")["nucleus"]["enabled"])
        self.assertIn("nuclear envelope and tethers", unsupported_features(config))
        sim = Simulation(_ablate(config, "no_tethering"))
        sim.step(50)
        self.assertEqual(sim.tethers, [])


class SchemaTests(unittest.TestCase):
    def preset(self):
        config = default_config()
        r = config["simulation"]["confinement_radius"]
        config["nucleus"] = {"enabled": True, "center": [-3 * r, 0, 0], "radii": [3.5 * r] * 3,
                             "envelope_stiffness": 30, "tethers": [tether(end=10, capture_distance=1.5)]}
        return config

    def test_valid_and_invalid_nucleus(self):
        self.assertTrue(validate_config(self.preset())["valid"], validate_config(self.preset())["errors"])
        bad = self.preset()
        bad["nucleus"]["center"] = [-4 * bad["simulation"]["confinement_radius"], 0, 0]
        self.assertFalse(validate_config(bad)["valid"])
        bad = self.preset()
        bad["nucleus"]["tethers"][0]["anchor"] = "body"
        bad["nucleus"]["tethers"][0]["body_id"] = "missing"
        self.assertFalse(validate_config(bad)["valid"])
        bad = self.preset()
        bad["nucleus"]["tethers"][0]["occupancy_attach_multiplier"] = {"open": 1}
        self.assertFalse(validate_config(bad)["valid"])

    def test_release_multiplier_range_allows_measured_nucleosome_effects(self):
        config = default_config()
        config["proteins"][0]["nucleosome_off_rate_multiplier"] = {"occupied": 500, "vacant": 1}
        self.assertTrue(validate_config(config)["valid"], validate_config(config)["errors"])
        config["proteins"][0]["nucleosome_off_rate_multiplier"] = {"occupied": 5000, "vacant": 1}
        self.assertFalse(validate_config(config)["valid"])


if __name__ == "__main__":
    unittest.main()
