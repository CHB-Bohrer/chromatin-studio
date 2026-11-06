"""Individual RNA identity, kinetics, mass balance, and spatial diffusion."""

import json
import math
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.transcription import ExplicitTranscription, disabled_transcription_state
from tests.test_engine import gene, protein, small_config, straight_chromatin
from tests.test_regulation import FixedRandom


def model(genes=None, **options):
    config = small_config()
    config["chromatin"]["bp_per_bead"] = 200
    config["genes"] = genes or [gene(start=1, end=3, on_rate=0, off_rate=0,
                                    initiation_rate=0, elongation_rate=1,
                                    processing_rate=0, rna_decay_rate=0)]
    config["transcription"] = {"enabled": True, "max_transcripts": 5000,
                               "rna_diffusion": 0, "mature_rna_diffusion": 0, **options}
    sim = Simulation(config)
    straight_chromatin(sim)
    explicit = ExplicitTranscription(sim)
    return sim, explicit


def rates(sim, initiation=0):
    return {g["config"]["id"]: {"effective_on_rate": 0,
                                 "effective_initiation_rate": initiation} for g in sim.genes}


def release(explicit, gene_state, count=1, position=None):
    records = []
    for _ in range(count):
        r = explicit._new_nascent(gene_state, 0)
        r.update(state="unprocessed", bead_position=None, released_time=0,
                 length_nt=600, position=np.asarray(position or [0, 0, 0], dtype=float))
        gene_state["completed"] += 1
        records.append(r)
    explicit._synchronize()
    return records


class ExplicitTranscriptionTests(unittest.TestCase):
    def test_identity_growth_attachment_release_and_mass_balance(self):
        sim, explicit = model()
        record = explicit._new_nascent(sim.genes[0], 0)
        explicit.advance(1.25, rates(sim))
        state = explicit.snapshot()
        r = state["transcripts"][0]
        self.assertEqual(r["id"], record["id"])
        self.assertEqual(r["polymerase_id"], "pol-1")
        self.assertEqual(r["length_nt"], 250)
        self.assertEqual(r["bead_position"], 2.25)
        np.testing.assert_allclose(r["position"], [2.25, 0, 0])
        sim.positions[2:4, 1] = 2
        self.assertEqual(explicit.snapshot()["transcripts"][0]["position"], [2.25, 2, 0])
        sim.time = 1.25
        explicit.advance(2, rates(sim))
        r = explicit.snapshot()["transcripts"][0]
        self.assertEqual(r["state"], "unprocessed")
        self.assertEqual(r["released_time"], 3)
        self.assertEqual(r["length_nt"], 600)
        self.assertIsNone(r["bead_position"])
        self.assertEqual(r["position"], [3, 2, 0])
        self.assertEqual(sim.genes[0]["polymerases"], [])
        self.assertEqual(sim.genes[0]["completed"], 1)
        self.assertEqual(sim.genes[0]["initiated"], 1)
        json.dumps(explicit.snapshot())

    def test_occupancy_changes_speed_at_each_crossed_boundary_and_can_stall(self):
        sim, explicit = model()
        sim.nucleosomes_enabled = True
        sim.nucleosome_occupied = np.ones(sim.n_chromatin, dtype=bool)
        sim.nucleosome_occupied[1] = False
        sim.genes[0]["config"]["occupied_elongation_multiplier"] = 0.5
        record = explicit._new_nascent(sim.genes[0], 0)
        explicit.advance(2, rates(sim))
        self.assertEqual(record["bead_position"], 2.5)
        sim.genes[0]["config"]["occupied_elongation_multiplier"] = 0
        explicit.advance(10, rates(sim))
        self.assertEqual(record["bead_position"], 2.5)
        sim.nucleosome_occupied[2] = False
        explicit.advance(0.5, rates(sim))
        self.assertEqual(record["bead_position"], 3)

    def test_pause_halts_growth_and_resume_restores_it(self):
        sim, explicit = model()
        record = explicit._new_nascent(sim.genes[0], 0)
        record["paused"] = True
        sim.genes[0]["config"]["resume_rate"] = 0
        explicit.advance(1, rates(sim))
        self.assertEqual(record["length_nt"], 0)
        self.assertEqual(explicit.snapshot()["genes"][0]["paused_polymerases"], 1)
        self.assertEqual(explicit.snapshot()["transcripts"][0]["elongation_rate"], 0)
        sim.genes[0]["config"]["resume_rate"] = 1000
        explicit.advance(1, rates(sim))
        self.assertGreater(record["length_nt"], 195)
        self.assertFalse(record["paused"])

    def test_capacity_rejects_attempts_without_erasing_existing_transcripts(self):
        sim, explicit = model(max_transcripts=2)
        sim.genes[0]["on"] = True
        explicit.advance(1, rates(sim, initiation=100))
        state = explicit.snapshot()
        self.assertEqual(state["totals"]["live"], 2)
        self.assertGreater(state["totals"]["blocked_initiations"], 50)
        first_ids = [r["id"] for r in state["transcripts"]]
        explicit.advance(1, rates(sim, initiation=100))
        self.assertEqual([r["id"] for r in explicit.snapshot()["transcripts"]], first_ids)
        self.assertEqual(sim.genes[0]["initiated"], 2)

    def test_chronological_births_do_not_always_favor_first_gene(self):
        winners = []
        for seed in range(40):
            a = gene(id="a", start=1, end=2, on_rate=0, off_rate=0)
            b = gene(id="b", start=4, end=5, on_rate=0, off_rate=0)
            sim, explicit = model([a, b], max_transcripts=1)
            sim.rng = np.random.default_rng(seed)
            for g in sim.genes:
                g["on"] = True
            explicit.advance(1, rates(sim, initiation=100))
            winners.append(explicit.snapshot()["transcripts"][0]["gene_id"])
        self.assertGreater(winners.count("a"), 8)
        self.assertGreater(winners.count("b"), 8)

    def test_processing_depends_on_rna_location_not_gene_contact(self):
        sim, explicit = model()
        sim.positions += [3, 0, 0]
        sim.bodies = [{"kind": "speckle", "center": [0, 0, 0], "radius": 0.5}]
        spec = sim.genes[0]["config"]
        spec.update(processing_rate=1, speckle_gain=9)
        record = release(explicit, sim.genes[0], position=[0, 0, 0])[0]
        self.assertFalse(sim._gene_speckle_contact(spec))
        self.assertEqual(explicit._processing_rate(record, spec), 10)
        record["position"] = np.asarray([10, 0, 0], dtype=float)
        sim.positions[1:4] = 0
        self.assertTrue(sim._gene_speckle_contact(spec))
        self.assertEqual(explicit._processing_rate(record, spec), 1)

    def test_competing_processing_decay_matches_two_state_solution(self):
        sim, explicit = model()
        sim.genes[0]["config"].update(processing_rate=1.2, rna_decay_rate=0.7)
        release(explicit, sim.genes[0], count=4000)
        explicit.advance(0.8, rates(sim))
        totals = explicit.snapshot()["totals"]
        survival = math.exp(-0.7 * 0.8)
        unprocessed = math.exp(-(1.2 + 0.7) * 0.8)
        self.assertAlmostEqual(totals["unprocessed"] / 4000, unprocessed, delta=0.025)
        self.assertAlmostEqual(totals["mature"] / 4000, survival - unprocessed, delta=0.025)
        self.assertAlmostEqual(totals["degraded"] / 4000, 1 - survival, delta=0.025)
        self.assertEqual(totals["initiated"], totals["live"] + totals["degraded"])
        self.assertEqual(totals["completed"], totals["unprocessed"] + totals["mature"] + totals["degraded"])

    def test_released_rna_diffusion_and_release_residual_time(self):
        sim, explicit = model(rna_diffusion=0.4)
        released = release(explicit, sim.genes[0], count=2000)
        explicit.advance(0.01, rates(sim))
        msd = np.mean([np.dot(r["position"], r["position"]) for r in released])
        self.assertAlmostEqual(msd / (6 * 0.01), 0.4, delta=0.025)
        sim, explicit = model(rna_diffusion=0.4)
        sim.genes[0]["config"]["elongation_rate"] = 100
        for _ in range(2000):
            explicit._new_nascent(sim.genes[0], 0)
        explicit.advance(0.04, rates(sim))
        records = explicit.snapshot()["transcripts"]
        self.assertTrue(all(abs(r["released_time"] - 0.03) < 1e-12 for r in records))
        displacement = np.asarray([r["position"] for r in records]) - [3, 0, 0]
        measured = np.mean(np.sum(displacement ** 2, axis=1)) / (6 * 0.01)
        self.assertAlmostEqual(measured, 0.4, delta=0.025)

    def test_mature_diffusion_has_its_own_coefficient(self):
        sim, explicit = model(rna_diffusion=0.4, mature_rna_diffusion=0.1)
        released = release(explicit, sim.genes[0], count=2000)
        for r in released:
            r["state"] = "mature"
        explicit.advance(0.01, rates(sim))
        msd = np.mean([np.dot(r["position"], r["position"]) for r in released])
        self.assertAlmostEqual(msd / (6 * 0.01), 0.1, delta=0.008)

    def test_disabled_helper_contains_no_invented_molecules(self):
        self.assertEqual(disabled_transcription_state()["transcripts"], [])
        self.assertFalse(disabled_transcription_state()["enabled"])

    def test_engine_hook_creates_identity_and_tracks_gene_polymerases(self):
        sim, _ = model()
        sim.genes[0]["on"] = True
        sim.genes[0]["config"]["initiation_rate"] = 100
        sim._advance_genes(0.1)
        snapshot = sim.snapshot()
        records = snapshot["transcription"]["transcripts"]
        self.assertTrue(records)
        expected = [1.1 - record["initiated_time"] for record in records]
        np.testing.assert_allclose(snapshot["genes"][0]["polymerase_positions"], expected)
        sim.genes[0]["config"]["initiation_rate"] = 0
        sim.time = 0.1
        sim._advance_genes(0.5)
        snapshot = sim.snapshot()
        np.testing.assert_allclose([r["length_nt"] for r in snapshot["transcription"]["transcripts"]],
                                   [200 * (0.6 - r["initiated_time"]) for r in records])
        np.testing.assert_allclose(snapshot["genes"][0]["polymerase_positions"], np.asarray(expected) + 0.5)

    def test_newly_released_rna_processes_only_for_remaining_interval(self):
        sim, explicit = model()
        spec = sim.genes[0]["config"]
        spec.update(elongation_rate=4, processing_rate=1)
        for _ in range(4000):
            explicit._new_nascent(sim.genes[0], 0)
        explicit.advance(1, rates(sim))
        totals = explicit.snapshot()["totals"]
        self.assertAlmostEqual(totals["mature"] / 4000, 1 - math.exp(-0.25), delta=0.025)
        self.assertTrue(all(r["processed_time"] is None or 0.75 <= r["processed_time"] <= 1
                            for r in explicit.snapshot()["transcripts"]))

    def test_soft_wall_drift_is_finite_and_does_not_clip_to_nuclear_surface(self):
        sim, explicit = model(rna_diffusion=0.1)
        records = release(explicit, sim.genes[0], count=2000, position=[35, 0, 0])
        explicit.advance(0.001, rates(sim))
        positions = np.asarray([r["position"] for r in records])
        self.assertTrue(np.all(np.isfinite(positions)))
        self.assertGreater(float(positions[:, 0].min()), sim.radius)
        self.assertAlmostEqual(float(positions[:, 0].mean()), 35 - 0.1 * 30 * 5 * 0.001, delta=0.001)

    def test_explicit_polymerase_eviction_changes_only_local_accessibility_and_stiffness(self):
        config = small_config()
        config["marks"] = [{"id": "A", "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}]
        config["chromatin"]["initial_marks"] = [{"start": 0, "end": 7, "mark": "A"}]
        config["chromatin"]["nucleosomes"] = {"enabled": True, "initial_occupancy": 1,
            "eviction_rate": 0, "reassembly_rate": 0,
            "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": 0.5}
        config["proteins"] = [protein("TF", chromatin_affinity={"unmodified": 1, "A": 1},
                                      nucleosome_affinity={"occupied": 0.1, "vacant": 1})]
        config["genes"] = [gene(start=2, end=5, on_rate=0, off_rate=0, initiation_rate=0,
                                 elongation_rate=0, resume_rate=0, persistence_multiplier=1,
                                 nucleosome_eviction_rate=2)]
        config["transcription"] = {"enabled": True, "max_transcripts": 10,
                                   "rna_diffusion": 0, "mature_rna_diffusion": 0}
        sim = Simulation(config)
        record = sim.explicit_transcription._new_nascent(sim.genes[0], 0)
        record["paused"] = True
        sim._advance_genes(0.1)
        self.assertEqual(sim.genes[0]["polymerases"], [2])
        before_affinity, before_stiffness = sim._chromatin_affinity("TF"), sim.local_persistence()
        sim.rng = FixedRandom(0.5)
        sim._advance_nucleosomes(1)
        np.testing.assert_array_equal(sim.nucleosome_occupied, [1, 1, 0, 1, 1, 1, 1, 1])
        np.testing.assert_array_equal(sim.marks["A"], sim.nucleosome_occupied)
        after_affinity, after_stiffness = sim._chromatin_affinity("TF"), sim.local_persistence()
        self.assertEqual(after_affinity[2], 10 * before_affinity[2])
        self.assertEqual(after_stiffness[2], 0.5 * before_stiffness[2])
        np.testing.assert_array_equal(np.delete(after_affinity, 2), np.delete(before_affinity, 2))
        np.testing.assert_array_equal(np.delete(after_stiffness, 2), np.delete(before_stiffness, 2))


if __name__ == "__main__":
    unittest.main()
