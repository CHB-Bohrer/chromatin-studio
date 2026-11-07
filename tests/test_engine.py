"""Physics and kinetic checks for the actual reference simulation."""

import copy
import json
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.schema import default_config, preset_config


def small_config():
    return {
        "seed": 23, "simulation": {"dt": 0.002, "confinement_radius": 30},
        "chromatin": {"beads": 8, "bond_length": 1.0, "bond_k": 30.0,
                      "diffusion": 0.1, "persistence_length": 2.0,
                      "repulsion": 0.0, "initial_marks": [], "barriers": []},
        "marks": [], "proteins": [], "interactions": [], "bodies": [],
        "cohesin": {"enabled": False}, "genes": [],
    }


def protein(species="P", count=1, **kwargs):
    value = {"id": species, "name": species, "count": count, "diffusion": 0.5,
             "radius": 0.25, "valence": 2, "representation": "bead",
             "idr_beads": 4, "sticker_fraction": 1.0,
             "chromatin_affinity": {"unmodified": 0}, "on_rate": 0,
             "off_rate": 0, "writer_mark": None, "eraser_mark": None,
             "catalytic_rate": 0}
    value.update(kwargs)
    return value


def gene(**kwargs):
    value = {"id": "g", "name": "Gene", "start": 2, "end": 5,
             "on_rate": 0.3, "off_rate": 0.2, "initiation_rate": 1,
             "elongation_rate": 2, "processing_rate": 0.2, "rna_decay_rate": 0,
             "activator_species": None, "required_mark": None,
             "repressive_mark": None, "speckle_gain": 3,
             "persistence_multiplier": 0.4}
    value.update(kwargs)
    return value


def straight_chromatin(sim):
    sim.positions[:sim.n_chromatin] = np.column_stack(
        (np.arange(sim.n_chromatin) * sim.bond_length, np.zeros(sim.n_chromatin), np.zeros(sim.n_chromatin)))


class MechanicalTests(unittest.TestCase):
    def test_strong_narrow_body_field_requires_safe_substeps_or_clear_error(self):
        c = small_config()
        c["simulation"]["dt"] = 0.01
        c["proteins"] = [protein(diffusion=5, valence=0)]
        c["bodies"] = [{"kind": "speckle", "center": [0, 0, 0], "radius": 0.2,
                         "strength": 20, "chromatin_affinity": 0,
                         "protein_affinity": {"P": 20}, "enabled": True}]
        s = Simulation(c)
        s.positions[8] = [0.2, 0, 0]
        with self.assertRaisesRegex(ValueError, "1000 integration substeps"):
            s.integration_substeps()
        s.dt = 0.001
        subdivisions = s.integration_substeps()
        self.assertGreater(subdivisions, 100)
        deterministic_move = s.diffusion[8] * np.linalg.norm(s.get_forces()[8]) * s.dt / subdivisions
        self.assertLess(deterministic_move, 0.2)

    def test_force_is_negative_energy_gradient_for_all_force_types(self):
        c = small_config()
        c["chromatin"]["repulsion"] = 3
        c["proteins"] = [protein(count=2, representation="idr", idr_beads=3)]
        c["genes"] = [gene()]
        c["bodies"] = [{"kind": "speckle", "center": [0, 0, 0], "radius": 3,
                         "strength": 2, "chromatin_affinity": -0.7,
                         "protein_affinity": {"P": 2}, "enabled": True}]
        s = Simulation(c)
        straight_chromatin(s)
        s.positions[2] += [0.1, 0.25, -0.3]
        s.positions[4] += [0.1, -0.3, 0.2]
        s.positions[8:] = [[1.2, 0.5, 0], [1.7, 0.5, 0.1], [2.1, 0.8, 0],
                           [3, 0.8, 0.4], [3.4, 1.1, 0.5], [3.8, 1.4, 0.5]]
        s.genes[0]["polymerases"] = [2.1]
        s.bindings = [{"a": 1, "b": 8, "kind": "chromatin", "off_rate": 0, "k": 12, "rest": 0.75},
                      {"a": 10, "b": 11, "kind": "protein", "off_rate": 0, "k": 4, "rest": 0.5}]
        s.cohesin = [[1, 6]]
        s.cohesin_config["spring_k"] = 1.3
        force = s.get_forces()
        numerical = np.zeros_like(force)
        h = 1e-6
        for i in range(len(s.positions)):
            for d in range(3):
                s.positions[i, d] += h
                plus = s.potential_energy()
                s.positions[i, d] -= 2 * h
                minus = s.potential_energy()
                s.positions[i, d] += h
                numerical[i, d] = -(plus - minus) / (2 * h)
        np.testing.assert_allclose(force, numerical, atol=1e-6, rtol=2e-6)

    def test_confinement_gradient_and_no_coordinate_clipping(self):
        s = Simulation(small_config())
        straight_chromatin(s)
        s.positions += [35, 0, 0]
        force = s.get_forces()
        self.assertLess(float(force[:, 0].sum()), 0)
        self.assertGreater(float(s.positions[0, 0]), s.radius)
        with self.assertRaises(FloatingPointError):
            s.validate_positions(np.full_like(s.positions, 1e9))

    def test_binding_applies_equal_opposite_force_to_chromatin(self):
        c = small_config()
        c["proteins"] = [protein()]
        s = Simulation(c)
        straight_chromatin(s)
        s.positions[8] = [2, 1.5, 0]
        before = s.get_forces()
        s.bindings = [{"a": 2, "b": 8, "kind": "chromatin", "off_rate": 0, "k": 12, "rest": 0.75}]
        delta = s.get_forces() - before
        np.testing.assert_allclose(delta[2], -delta[8], atol=1e-12)
        np.testing.assert_allclose(delta[2], [0, 9, 0], atol=1e-12)
        np.testing.assert_allclose(delta.sum(axis=0), 0, atol=1e-12)

    def test_internal_forces_conserve_total_force(self):
        c = small_config()
        c["chromatin"]["repulsion"] = 5
        c["proteins"] = [protein(representation="idr")]
        s = Simulation(c)
        np.testing.assert_allclose(s.get_forces().sum(axis=0), 0, atol=1e-10)

    def test_free_particle_diffusion_and_idr_com_diffusion(self):
        for representation, beads in [("bead", 1), ("idr", 4)]:
            c = small_config()
            c["simulation"]["confinement_radius"] = 200
            c["chromatin"]["diffusion"] = 0
            c["proteins"] = [protein(count=80, representation=representation, idr_beads=beads,
                                      diffusion=0.4, valence=0)]
            s = Simulation(c)
            self.assertTrue(np.all(s.diffusion[8:] == 0.4 * beads))
            variances = []
            for _ in range(60):
                before = s.positions[8:].reshape(80, beads, 3).mean(axis=1)
                s.step(1)
                after = s.positions[8:].reshape(80, beads, 3).mean(axis=1)
                variances.append(np.sum((after - before) ** 2, axis=1))
            measured = float(np.mean(variances)) / (6 * s.dt)
            self.assertAlmostEqual(measured, 0.4, delta=0.035)

    def test_gene_occupancy_changes_local_persistence(self):
        c = small_config()
        c["genes"] = [gene()]
        s = Simulation(c)
        before = s.local_persistence()
        s.genes[0]["polymerases"] = [3.0]
        during = s.local_persistence()
        np.testing.assert_allclose(during[2:6], before[2:6] * 0.4)
        np.testing.assert_allclose(during[:2], before[:2])
        s.genes[0]["polymerases"] = []
        np.testing.assert_array_equal(s.local_persistence(), before)


class ReactionTests(unittest.TestCase):
    def test_initial_cohesin_honors_maximum_valid_requested_count(self):
        c = small_config()
        c["seed"] = 0
        c["cohesin"] = {"enabled": True, "count": 4, "step_rate": 1, "off_rate": 0, "spring_k": 1}
        s = Simulation(c)
        self.assertEqual(len(s.cohesin), 4)
        self.assertEqual(sorted(sum(s.cohesin, [])), list(range(8)))
        self.assertTrue(all(right == left + 1 for left, right in s.cohesin))

    def test_multivalent_binding_cannot_exceed_molecule_capacity(self):
        c = small_config()
        c["proteins"] = [protein(count=6, valence=2)]
        c["interactions"] = [{"a": "P", "b": "P", "on_rate": 1e6, "off_rate": 0, "strength": 4, "radius": 2}]
        s = Simulation(c)
        s.positions[8:] = np.column_stack((np.linspace(-0.5, 0.5, 6), np.zeros(6), np.zeros(6)))
        s._advance_binding(1)
        degree = np.zeros(6, dtype=int)
        for b in s.bindings:
            degree[s._particle_info[b["a"]]["molecule"]] += 1
            degree[s._particle_info[b["b"]]["molecule"]] += 1
        self.assertTrue(np.all(degree <= 2))
        self.assertGreaterEqual(s._largest_cluster(), 3)
        self.assertEqual(len({tuple(sorted((b["a"], b["b"]))) for b in s.bindings}), len(s.bindings))

    def test_idr_stickers_and_molecule_valence_both_limit_links(self):
        c = small_config()
        c["proteins"] = [protein(count=2, representation="idr", idr_beads=4, valence=5,
                                  sticker_fraction=0.5, on_rate=1e6,
                                  chromatin_affinity={"unmodified": 1})]
        s = Simulation(c)
        straight_chromatin(s)
        s.positions[8:] = [2, 0.7, 0]
        s._advance_binding(1)
        self.assertEqual(len(s.bindings), 4)
        self.assertEqual(len({b["b"] for b in s.bindings}), 4)
        self.assertTrue(all(s._particle_info[b["b"]]["sticker"] for b in s.bindings))

    def test_writer_requires_chromatin_binding_and_eraser_removes_mark(self):
        c = small_config()
        c["marks"] = [{"id": "M", "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}]
        c["proteins"] = [protein("writer", writer_mark="M", catalytic_rate=1e5),
                          protein("eraser", eraser_mark="M", catalytic_rate=1e5)]
        s = Simulation(c)
        s._advance_marks(1)
        self.assertFalse(np.any(s.marks["M"]))
        s.bindings = [{"a": 3, "b": 8, "kind": "chromatin", "off_rate": 0, "k": 12, "rest": 0.75}]
        s._advance_marks(1)
        self.assertEqual(np.flatnonzero(s.marks["M"]).tolist(), [3])
        s.bindings[0]["b"] = 9
        s._advance_marks(1)
        self.assertFalse(np.any(s.marks["M"]))

    def test_contact_spreading_is_spatial_and_affinity_changes(self):
        c = small_config()
        c["marks"] = [{"id": "M", "decay_rate": 0, "spread_rate": 1e5, "contact_radius": 0.4}]
        c["chromatin"]["initial_marks"] = [{"start": 0, "end": 0, "mark": "M"}]
        c["proteins"] = [protein(chromatin_affinity={"unmodified": 1, "M": 0})]
        s = Simulation(c)
        straight_chromatin(s)
        s.positions[7] = [0.2, 0, 0]
        s._advance_marks(1)
        self.assertEqual(np.flatnonzero(s.marks["M"]).tolist(), [0, 7])
        self.assertEqual(s._chromatin_affinity("P")[7], 0)
        self.assertEqual(s._chromatin_affinity("P")[6], 1)

    def test_directional_barriers_and_cohesin_collisions(self):
        c = small_config()
        c["cohesin"] = {"enabled": True, "count": 0, "step_rate": 1e3, "off_rate": 0, "spring_k": 1}
        c["chromatin"]["barriers"] = [{"bead": 1, "direction": "left"}, {"bead": 6, "direction": "right"}]
        s = Simulation(c)
        s.cohesin = [[2, 3]]
        s._advance_cohesin(1)
        self.assertEqual(s.cohesin, [[1, 6]])
        s.barriers = {}
        s.cohesin = [[1, 2], [5, 6]]
        s._advance_cohesin(1)
        self.assertEqual(s.cohesin[0][0], 0)
        self.assertEqual(s.cohesin[1][1], 7)
        self.assertLess(s.cohesin[0][1], s.cohesin[1][0])
        self.assertEqual(len(set(sum(s.cohesin, []))), 4)

    def test_mark_gates_and_locus_bound_activation(self):
        c = small_config()
        c["marks"] = [{"id": m, "decay_rate": 0, "spread_rate": 0, "contact_radius": 1} for m in ("active", "silent")]
        c["proteins"] = [protein("activator")]
        c["genes"] = [gene(required_mark="active", repressive_mark="silent", activator_species="activator")]
        s = Simulation(c)
        spec = s.genes[0]["config"]
        self.assertEqual(s._gene_on_rate(spec), 0)
        s.marks["active"][2] = True
        self.assertEqual(s._gene_on_rate(spec), 0.3)
        s.bindings = [{"a": 3, "b": 8, "kind": "chromatin", "off_rate": 0, "k": 12, "rest": 0.75}]
        self.assertEqual(s._gene_on_rate(spec), 0.6)
        s.marks["silent"][2] = True
        self.assertEqual(s._gene_on_rate(spec), 0)

    def test_polymerases_move_and_rna_is_conserved_without_decay(self):
        c = small_config()
        c["genes"] = [gene(on_rate=0, off_rate=0, initiation_rate=0, processing_rate=1e6)]
        s = Simulation(c)
        g = s.genes[0]
        g["polymerases"] = [2.0, 3.0]
        g["initiated"] = 2
        s._advance_genes(0.5)
        self.assertEqual(g["polymerases"], [3.0, 4.0])
        s._advance_genes(2)
        self.assertEqual(g["completed"], 2)
        self.assertEqual(g["mature"], 2)
        self.assertEqual(len(g["polymerases"]) + g["unprocessed"] + g["mature"], g["initiated"])

    def test_speckle_contact_increases_processing_only(self):
        c = small_config()
        c["genes"] = [gene(on_rate=0, off_rate=0, initiation_rate=0)]
        contact = copy.deepcopy(c)
        contact["bodies"] = [{"kind": "speckle", "center": [0, 0, 0], "radius": 50,
                               "strength": 0, "chromatin_affinity": 0, "protein_affinity": {}, "enabled": True}]
        outside, inside = Simulation(c), Simulation(contact)
        for s in (outside, inside):
            s.genes[0]["unprocessed"] = 20000
            s._advance_genes(1)
        self.assertGreater(inside.genes[0]["mature"], 2.5 * outside.genes[0]["mature"])
        self.assertEqual(inside.genes[0]["initiated"], outside.genes[0]["initiated"])


class IntegrationTests(unittest.TestCase):
    def test_seed_reproducibility_and_json_snapshot(self):
        a, b = Simulation(default_config()), Simulation(default_config())
        a.step(30)
        b.step(30)
        left, right = a.snapshot(), b.snapshot()
        # Session identity is observation metadata, not seeded physical state.
        left["genomic_events"].pop("run_id")
        right["genomic_events"].pop("run_id")
        self.assertEqual(left, right)
        json.dumps(a.snapshot(), allow_nan=False)
        self.assertAlmostEqual(a.time, 0.06)

    def test_default_and_presets_are_finite_for_1000_steps(self):
        for name in ("nuclear", "assembly", "epigenetic", "extrusion"):
            with self.subTest(preset=name):
                s = Simulation(preset_config(name))
                s.step(1000)
                self.assertTrue(np.all(np.isfinite(s.positions)))
                self.assertLess(np.max(np.linalg.norm(s.positions, axis=1)), s.radius + 2)
                self.assertEqual(s.step_count, 1000)
                self.assertTrue(np.isfinite(s.potential_energy()))

    def test_backend_hook_advances_shared_reactions_and_time(self):
        c = small_config()
        calls = []

        def dynamics(sim):
            calls.append(sim.time)

        s = Simulation(c, backend="test", dynamics=dynamics)
        before = s.positions.copy()
        s.step(3)
        np.testing.assert_array_equal(s.positions, before)
        self.assertEqual(len(calls), 3)
        self.assertAlmostEqual(s.time, 3 * s.dt)
        self.assertEqual(s.engine, "test")
        with self.assertRaises(ValueError):
            Simulation(c, backend="polychrom")


if __name__ == "__main__":
    unittest.main()
