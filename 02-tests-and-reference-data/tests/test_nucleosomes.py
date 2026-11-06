"""Occupancy-dependent binding, stiffness, remodeling and histone mark memory."""

import copy
import json
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from tests.test_engine import gene, protein, small_config
from tests.test_interaction_state import link
from tests.test_regulation import FixedRandom, region, regulatory_system


def nucleosomes(**updates):
    config = {"enabled": True, "initial_occupancy": 1, "initial_regions": [],
              "eviction_rate": 0, "reassembly_rate": 0,
              "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": 0.5}
    config.update(updates)
    return config


def marked_config():
    config = small_config()
    config["marks"] = [{"id": "A", "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}]
    config["chromatin"]["initial_marks"] = [{"start": 0, "end": 7, "mark": "A"}]
    config["chromatin"]["nucleosomes"] = nucleosomes()
    return config


class NucleosomeTests(unittest.TestCase):
    def test_initial_intervals_override_in_order_and_vacancy_clears_initial_marks(self):
        config = marked_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 1, "end": 3, "occupied": False},
                             {"start": 2, "end": 2, "occupied": True}],
            occupied_persistence_multiplier=1.2)
        sim = Simulation(config)
        expected = [1, 0, 1, 0, 1, 1, 1, 1]
        np.testing.assert_array_equal(sim.nucleosome_occupied, expected)
        np.testing.assert_array_equal(sim.marks["A"], expected)
        self.assertEqual(sim.nucleosome_state()["occupied_fraction"], 0.75)
        self.assertEqual(sim.genome_tracks()["nucleosome_occupancy_runs"],
                         [[0, 1, 1], [1, 2, 0], [2, 3, 1], [3, 4, 0], [4, 8, 1]])
        np.testing.assert_allclose(sim.local_persistence(), [2.4, 1, 2.4, 1, 2.4, 2.4, 2.4, 2.4])

    def test_occupancy_affinity_and_stiffness_multiply_existing_region_and_transcription_effects(self):
        config = marked_config()
        config["proteins"] = [protein("TF", chromatin_affinity={"unmodified": 2, "A": 5},
                                       nucleosome_affinity={"occupied": 0.2, "vacant": 3})]
        config["genes"] = [gene(persistence_multiplier=0.4)]
        config["regulatory_regions"] = [region(start=2, end=2, coactivator_species=None,
                                                tf_affinity_multiplier=4)]
        config["chromatin"]["nucleosomes"]["initial_regions"] = [{"start": 2, "end": 2, "occupied": False}]
        sim = Simulation(config)
        self.assertEqual(sim._chromatin_affinity("TF")[2], 24)  # Baseline2 × region4 × vacant3.
        self.assertEqual(sim._chromatin_affinity("TF")[3], 1)   # Mark5 × occupied0.2.
        sim.genes[0]["polymerases"] = [2.2]
        np.testing.assert_allclose(sim.local_persistence(), [2, 2, 0.4, 0.8, 0.8, 0.8, 2, 2])
        self.assertEqual(sim.genome_tracks()["persistence_runs"],
                         [[0, 2, 2.0], [2, 3, 0.4], [3, 6, 0.8], [6, 8, 2.0]])
        # An accessibility multiplier does not create affinity from zero.
        sim.species["TF"]["chromatin_affinity"]["unmodified"] = 0
        self.assertEqual(sim._chromatin_affinity("TF")[2], 0)

    def test_actual_binding_candidates_follow_occupied_and_vacant_preferences(self):
        config = small_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 3, "end": 3, "occupied": False}])
        config["proteins"] = [protein("TF", on_rate=1e6, chromatin_affinity={"unmodified": 1},
                                       nucleosome_affinity={"occupied": 0, "vacant": 1})]
        sim = Simulation(config)
        sim.positions[:8] = np.column_stack((np.arange(8) * 0.1, np.zeros(8), np.zeros(8)))
        sim.positions[8] = [0.3, 0.1, 0]
        sim._advance_binding(1)
        self.assertEqual([(bond["a"], bond["kind"]) for bond in sim.bindings], [(3, "chromatin")])

    def test_vacancy_changes_bending_force_consistently_with_model_energy(self):
        config = small_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(vacant_persistence_multiplier=0.2)
        sim = Simulation(config)
        sim.positions[:] = np.column_stack((np.arange(8), np.zeros(8), np.zeros(8)))
        sim.positions[3, 1] = 0.4
        occupied_force = sim.get_forces()
        sim.nucleosome_occupied[3] = False
        vacant_force = sim.get_forces()
        self.assertFalse(np.array_equal(vacant_force, occupied_force))
        numerical = np.zeros_like(vacant_force)
        h = 1e-6
        for bead in range(8):
            for axis in range(3):
                sim.positions[bead, axis] += h
                plus = sim.potential_energy()
                sim.positions[bead, axis] -= 2 * h
                minus = sim.potential_energy()
                sim.positions[bead, axis] += h
                numerical[bead, axis] = -(plus - minus) / (2 * h)
        np.testing.assert_allclose(vacant_force, numerical, atol=1e-7, rtol=1e-7)

    def test_remodeling_rates_deduplicate_actual_direct_and_recruited_anchors(self):
        template = regulatory_system()
        config = copy.deepcopy(template.config)
        config["chromatin"]["nucleosomes"] = nucleosomes(eviction_rate=0.5, reassembly_rate=0.25)
        config["proteins"][0]["nucleosome_eviction_rate"] = 1
        config["proteins"][1].update(nucleosome_eviction_rate=2, nucleosome_reassembly_rate=3)
        config["regulatory_regions"].append(region("overlap"))
        config["genes"][0]["nucleosome_eviction_rate"] = 0.7
        sim = Simulation(config)
        sim.bindings = copy.deepcopy(template.bindings)
        sim.bindings.append(link(0, 16, "chromatin"))  # CO molecule also directly attached at TF anchor.
        sim.genes[0]["polymerases"] = [0.3, 0.9, 7.99]
        eviction, reassembly = sim._nucleosome_rates()
        np.testing.assert_allclose(eviction, [4.9, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.2])
        np.testing.assert_allclose(reassembly, [3.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])

    def test_eviction_and_reassembly_within_interval_clear_marks_even_with_same_endpoint(self):
        config = marked_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(eviction_rate=1, reassembly_rate=1)
        sim = Simulation(config)
        # P(final vacant)=.4323 < .55 < P(any eviction)=.6321.
        sim.rng = FixedRandom(0.55)
        sim._advance_nucleosomes(1)
        self.assertTrue(sim.nucleosome_occupied.all())
        self.assertFalse(sim.marks["A"].any())
        self.assertEqual(sim.nucleosome_state()["endpoint_evictions"], 0)
        self.assertEqual(sim.nucleosome_state()["endpoint_reassemblies"], 0)
        self.assertEqual(sim.nucleosome_state()["mark_reset_intervals"], 8)

    def test_reassembly_has_no_historical_marks_and_existing_protein_bond_is_retained(self):
        config = marked_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(eviction_rate=1)
        config["proteins"] = [protein("reader")]
        sim = Simulation(config)
        sim.bindings = [link(0, 8, "chromatin")]
        sim.rng = FixedRandom(0.5)
        sim._advance_nucleosomes(1)
        self.assertFalse(sim.nucleosome_occupied.any())
        self.assertFalse(sim.marks["A"].any())
        self.assertEqual(len(sim.bindings), 1)
        sim.nucleosome_config.update(eviction_rate=0, reassembly_rate=1)
        sim._advance_nucleosomes(1)
        self.assertTrue(sim.nucleosome_occupied.all())
        self.assertFalse(sim.marks["A"].any())
        self.assertEqual(sim.nucleosome_state()["endpoint_evictions"], 8)
        self.assertEqual(sim.nucleosome_state()["endpoint_reassemblies"], 8)

    def test_vacant_beads_cannot_receive_or_donate_histone_marks(self):
        config = small_config()
        config["marks"] = [{"id": "A", "decay_rate": 0, "spread_rate": 1, "contact_radius": 1}]
        config["chromatin"]["initial_marks"] = [{"start": 0, "end": 0, "mark": "A"}]
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 1, "end": 1, "occupied": False}])
        config["proteins"] = [protein("writer", writer_mark="A", catalytic_rate=1e6)]
        sim = Simulation(config)
        sim.positions[:8] = np.column_stack((np.arange(8) * 10, np.zeros(8), np.zeros(8)))
        sim.positions[1] = [0.1, 0, 0]
        sim.positions[2] = [0.2, 0, 0]
        sim.bindings = [link(1, 8, "chromatin")]
        sim.rng = FixedRandom(0)
        sim._advance_marks(1)
        self.assertFalse(sim.marks["A"][1])
        self.assertTrue(sim.marks["A"][2])
        sim.marks["A"][:] = False
        sim.marks["A"][1] = True  # Force invalid external state to test donor masking.
        sim.bindings = []
        sim._advance_marks(1)
        self.assertFalse(sim.marks["A"].any())

    def test_disabled_and_neutral_enabled_models_preserve_seeded_legacy_dynamics(self):
        base = small_config()
        base["genes"] = [gene()]
        base["proteins"] = [protein("P", count=2, on_rate=2, chromatin_affinity={"unmodified": 1})]
        absent = Simulation(base)
        disabled_config = copy.deepcopy(base)
        disabled_config["chromatin"]["nucleosomes"] = nucleosomes(enabled=False, eviction_rate=20,
                                                                 initial_occupancy=0)
        disabled = Simulation(disabled_config)
        neutral_config = copy.deepcopy(base)
        neutral_config["chromatin"]["nucleosomes"] = nucleosomes()
        neutral = Simulation(neutral_config)
        for sim in (absent, disabled, neutral):
            sim.step(80)
        for other in (disabled, neutral):
            np.testing.assert_array_equal(absent.positions, other.positions)
            self.assertEqual(absent.genes, other.genes)
            self.assertEqual(absent.bindings, other.bindings)
            self.assertEqual(absent.rng.bit_generator.state, other.rng.bit_generator.state)
        self.assertEqual(disabled.nucleosome_state()["occupied"], [])
        self.assertIsNone(disabled.nucleosome_state()["occupied_fraction"])
        self.assertEqual(disabled.genome_tracks()["nucleosome_occupancy_runs"], [])

    def test_stochastic_initial_occupancy_and_reports_are_reproducible_and_read_only(self):
        config = marked_config()
        config["chromatin"]["nucleosomes"]["initial_occupancy"] = 0.5
        left, right = Simulation(config), Simulation(config)
        np.testing.assert_array_equal(left.nucleosome_occupied, right.nucleosome_occupied)
        before = copy.deepcopy(left.rng.bit_generator.state)
        first = left.snapshot()
        self.assertEqual(first, left.snapshot())
        self.assertEqual(before, left.rng.bit_generator.state)
        json.dumps(first, allow_nan=False)
        first["nucleosomes"]["occupied"][0] = 5
        first["genome_tracks"]["persistence_runs"][0][2] = 999
        self.assertNotEqual(left.nucleosome_state()["occupied"][0], 5)
        self.assertNotEqual(left.genome_tracks()["persistence_runs"][0][2], 999)


if __name__ == "__main__":
    unittest.main()
