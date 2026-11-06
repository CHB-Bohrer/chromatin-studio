"""Confinement, optional body hypotheses, and continuous promoter occupancy."""

import copy
import json
import math
import unittest
from unittest.mock import patch

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.transcription import sample_promoter_interval
from tests.test_engine import gene, protein, small_config


def body(identifier="speckle", **changes):
    result = {"id": identifier, "kind": "speckle", "center": [0, 0, 0], "radius": 2,
              "strength": 0, "chromatin_affinity": 0, "protein_affinity": {}, "enabled": True}
    result.update(changes)
    return result


class BirthRandom:
    """Two deterministic birth times; no other reactions are enabled."""
    def poisson(self, mean):
        return 2

    def random(self, size=None):
        return np.asarray([0.25, 0.75]) if size is not None else 0.5


class ConfinementTests(unittest.TestCase):
    def test_wall_force_and_energy_use_configured_stiffness(self):
        config = small_config()
        config["proteins"] = [protein()]
        config["simulation"]["confinement_stiffness"] = 7.3
        confined = Simulation(config)
        config["simulation"]["confinement_stiffness"] = 0
        free = Simulation(config)
        confined.positions[-1] = [33, 4, 0]
        free.positions[:] = confined.positions
        r = np.linalg.norm(confined.positions[-1])
        expected = -7.3 * (r - confined.radius) * confined.positions[-1] / r
        np.testing.assert_allclose((confined.get_forces() - free.get_forces())[-1], expected)
        self.assertAlmostEqual(confined.potential_energy() - free.potential_energy(),
                               0.5 * 7.3 * (r - confined.radius) ** 2)
        for axis in range(3):
            original = confined.positions[-1, axis]
            confined.positions[-1, axis] = original + 1e-5
            plus = confined.potential_energy()
            confined.positions[-1, axis] = original - 1e-5
            minus = confined.potential_energy()
            confined.positions[-1, axis] = original
            self.assertAlmostEqual(-(plus - minus) / 2e-5, confined.get_forces()[-1, axis], places=5)
        free.positions += 100000
        free.validate_positions()  # Unconfined diffusion has no hidden radius rejection.

    def test_legacy_wall_default_and_zero_wall_rna_drift(self):
        config = small_config()
        config["genes"] = [gene()]
        config["transcription"] = {"enabled": True, "rna_diffusion": 0.4}
        sim = Simulation(config)
        self.assertEqual(sim.confinement_stiffness, 30)
        rna = sim.explicit_transcription._new_nascent(sim.genes[0], 0)
        rna.update(state="unprocessed", position=np.array([35., 0, 0]))
        class NoNoise:
            def normal(self, size):
                return np.zeros(size)
        sim.rng = NoNoise()
        sim.confinement_stiffness = 0
        sim.explicit_transcription._brownian(rna, 0.01)
        np.testing.assert_array_equal(rna["position"], [35, 0, 0])
        sim.confinement_stiffness = 5
        sim.explicit_transcription._brownian(rna, 0.01)
        np.testing.assert_allclose(rna["position"], [34.9, 0, 0])

    def test_diagnostics_count_released_rna_separately_and_are_read_only(self):
        config = small_config()
        config["simulation"]["confinement_stiffness"] = 2
        config["genes"] = [gene()]
        config["transcription"] = {"enabled": True}
        sim = Simulation(config)
        sim.positions[0] = [31, 0, 0]
        sim.explicit_transcription._new_nascent(sim.genes[0], 0)
        released = sim.explicit_transcription._new_nascent(sim.genes[0], 0)
        released.update(state="unprocessed", position=np.array([32., 0, 0]), bead_position=None)
        state = copy.deepcopy(sim.rng.bit_generator.state)
        before = sim.positions.copy()
        report = sim.snapshot()["confinement"]
        self.assertEqual(report["particles"]["outside_count"], 1)
        self.assertEqual(report["rna"]["total"], 1)
        self.assertEqual(report["rna"]["max_penetration"], 2)
        self.assertEqual(report["energy"], 5)
        self.assertEqual(report, sim.confinement_state())
        self.assertEqual(state, sim.rng.bit_generator.state)
        np.testing.assert_array_equal(sim.positions, before)
        sim.confinement_stiffness = 0
        report = sim.confinement_state()
        self.assertFalse(report["enabled"])
        self.assertEqual(report["energy"], 0)
        self.assertEqual(report["particles"]["outside_count"], 1)


class BodyResponseTests(unittest.TestCase):
    def make_sim(self, responses):
        config = small_config()
        config["genes"] = [gene(on_rate=0.4, initiation_rate=3, body_responses=responses)]
        config["bodies"] = [body(), body("nucleolus", kind="nucleolus")]
        sim = Simulation(config)
        sim.positions[2] = [0, 0, 0]
        return sim

    def test_center_distance_gaussian_and_repressive_response(self):
        sim = self.make_sim({"speckle": {"activation_multiplier": 4, "initiation_multiplier": 9},
                             "nucleolus": {"activation_multiplier": 0.5}})
        report = sim.regulation_state()["genes"][0]
        self.assertAlmostEqual(report["effective_on_rate"], 0.8)
        self.assertAlmostEqual(report["effective_initiation_rate"], 27)
        self.assertEqual(report["body_responses"][0]["weight"], 1)
        sim.positions[2] = [2, 0, 0]
        report = sim.regulation_state()["genes"][0]
        w = math.exp(-0.5)
        self.assertAlmostEqual(report["effective_on_rate"], 0.4 * 2 ** w)
        self.assertAlmostEqual(report["effective_initiation_rate"], 3 * 9 ** w)
        self.assertEqual(report["causal_terms"]["on_rate"]["effective"], report["effective_on_rate"])
        self.assertEqual(report["body_responses"][0]["distance"], 2)
        sim.positions[2] = [200, 0, 0]
        report = sim.regulation_state()["genes"][0]
        self.assertEqual(report["effective_on_rate"], 0.4)
        self.assertEqual(report["effective_initiation_rate"], 3)

    def test_mark_gate_and_disabled_body_do_not_gain_activation(self):
        sim = self.make_sim({"speckle": {"activation_multiplier": 4}})
        sim.marks["R"] = np.zeros(8, dtype=bool)
        sim.marks["R"][2] = True
        sim.genes[0]["config"]["repressive_mark"] = "R"
        report = sim.regulation_state()["genes"][0]
        self.assertEqual(report["effective_on_rate"], 0)
        self.assertFalse(report["causal_terms"]["on_rate"]["mark_gate"])
        config = sim.config
        config["bodies"][0]["enabled"] = False
        disabled = Simulation(config)
        self.assertEqual(disabled.regulation_state()["genes"][0]["environment_activation_multiplier"], 1)

    def test_neutral_response_has_identical_dynamics_and_read_only_reports(self):
        absent = self.make_sim({})
        neutral = self.make_sim({"speckle": {"activation_multiplier": 1, "initiation_multiplier": 1}})
        for _ in range(5):
            absent._advance_genes(0.4)
            neutral._advance_genes(0.4)
        self.assertEqual(absent.rng.bit_generator.state, neutral.rng.bit_generator.state)
        for key in ("on", "polymerases", "unprocessed", "mature", "initiated", "completed"):
            self.assertEqual(absent.genes[0][key], neutral.genes[0][key])
        before = copy.deepcopy(neutral.rng.bit_generator.state)
        report = json.dumps(neutral.snapshot(), sort_keys=True)
        self.assertEqual(report, json.dumps(neutral.snapshot(), sort_keys=True))
        self.assertEqual(before, neutral.rng.bit_generator.state)

    def test_overlapping_body_rates_fail_before_large_event_allocation(self):
        sim = self.make_sim({"speckle": {"initiation_multiplier": 20},
                             "nucleolus": {"initiation_multiplier": 20}})
        sim.genes[0]["on"] = True
        sim.genes[0]["config"].update(off_rate=0, initiation_rate=100)
        with self.assertRaisesRegex(ValueError, "reduce dt.*body response"):
            sim._advance_genes(1)
        self.assertEqual(sim.genes[0]["initiated"], 0)
        self.assertEqual(sim.genes[0]["polymerases"], [])


class PromoterTimingTests(unittest.TestCase):
    def test_on_time_integral_and_birth_mean_match_telegraph_solution(self):
        rng = np.random.default_rng(2349)
        kon, koff, production, dt, size = 0.7, 1.1, 3, 1.2, 10000
        stationary = kon / (kon + koff)
        expected = stationary * (dt - (-math.expm1(-(kon + koff) * dt)) / (kon + koff))
        paths = [sample_promoter_interval(rng, False, kon, koff, production, dt) for _ in range(size)]
        self.assertAlmostEqual(np.mean([p["on_time"] for p in paths]), expected, delta=0.015)
        self.assertAlmostEqual(np.mean([len(p["initiation_times"]) for p in paths]), production * expected, delta=0.035)
        self.assertTrue(all(0 <= t < dt for p in paths for t in p["initiation_times"]))

    def test_new_births_only_elongate_for_remaining_time_in_both_modes(self):
        for explicit in (False, True):
            with self.subTest(explicit=explicit):
                config = small_config()
                config["genes"] = [gene(start=1, end=1, on_rate=0, off_rate=0, initiation_rate=2,
                                         elongation_rate=2, processing_rate=0)]
                config["transcription"] = {"enabled": explicit, "max_transcripts": 10,
                                            "rna_diffusion": 0, "mature_rna_diffusion": 0}
                sim = Simulation(config)
                sim.genes[0]["on"] = True
                sim.rng = BirthRandom()
                sim._advance_genes(1)
                state = sim.genes[0]
                self.assertEqual(state["polymerases"], [1.5])
                self.assertEqual(state["completed"], 1)
                self.assertEqual(state["unprocessed"], 1)
                self.assertEqual(state["initiated"], 2)
                self.assertEqual(state["last_transcription_interval"]["on_time"], 1)
                if explicit:
                    records = list(sim.explicit_transcription.transcripts.values())
                    self.assertEqual([r["initiated_time"] for r in records], [0.25, 0.75])
                    self.assertEqual(records[0]["released_time"], 0.75)

    def test_count_processing_and_decay_match_exact_probabilities(self):
        config = small_config()
        config["genes"] = [gene(on_rate=0, off_rate=0, initiation_rate=0, processing_rate=1,
                                 rna_decay_rate=0.7)]
        sim = Simulation(config)
        size = 50000
        sim.genes[0]["unprocessed"] = size
        sim._advance_genes(1)
        self.assertAlmostEqual(sim.genes[0]["unprocessed"] / size, math.exp(-1.7), delta=0.008)
        self.assertAlmostEqual(sim.genes[0]["mature"] / size, math.exp(-0.7) * (1 - math.exp(-1)), delta=0.008)

    def test_explicit_capacity_is_checked_at_actual_birth_time(self):
        config = small_config()
        config["genes"] = [gene(on_rate=0, off_rate=0, initiation_rate=2,
                                 elongation_rate=0, processing_rate=0)]
        config["transcription"] = {"enabled": True, "max_transcripts": 1}
        sim = Simulation(config)
        explicit = sim.explicit_transcription
        old = explicit._new_nascent(sim.genes[0], 0)
        old["state"] = "unprocessed"
        sim.genes[0]["on"] = True
        sim.rng = BirthRandom()
        def degrade_at_half(record, spec, dt, start):
            if start + dt >= 0.5:
                del explicit.transcripts[record["id"]]
                explicit.degraded[spec["id"]] += 1
        with patch.object(explicit, "_advance_free", side_effect=degrade_at_half):
            sim._advance_genes(1)
        self.assertEqual(explicit.blocked["g"], 1)
        self.assertEqual([r["initiated_time"] for r in explicit.transcripts.values()], [0.75])
        self.assertTrue(any("capacity" in warning for warning in sim.warnings))


if __name__ == "__main__":
    unittest.main()
