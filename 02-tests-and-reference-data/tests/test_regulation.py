"""Explicit TF occupancy, one-hop recruitment, catalytic targeting, and rates."""

import copy
import json
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from tests.test_engine import gene, protein, small_config
from tests.test_interaction_state import link


def region(identifier="enhancer", **overrides):
    value = {"id": identifier, "name": identifier, "color": "#227788",
             "kind": "enhancer", "start": 0, "end": 1, "target_gene": "g",
             "tf_species": "TF", "coactivator_species": "CO",
             "tf_affinity_multiplier": 2, "activation_gain": 2,
             "initiation_gain": 4, "coactivator_gain": 2,
             "contact_radius": 0.5, "enabled": True}
    value.update(overrides)
    return value


def regulatory_system(explicit_rna=False):
    config = small_config()
    config["genes"] = [gene(start=5, end=7, on_rate=0.2, initiation_rate=3,
                            off_rate=0, activator_species="TF")]
    config["marks"] = [{"id": key, "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}
                       for key in ("A", "R")]
    config["proteins"] = [protein("TF", count=2, representation="idr", idr_beads=3, valence=6),
                          protein("CO", count=2, representation="idr", idr_beads=3,
                                  valence=6, writer_mark="A", catalytic_rate=1),
                          protein("relay")]
    config["regulatory_regions"] = [region(), region("promoter", kind="promoter", start=5, end=6,
                                                    tf_species=None, coactivator_species=None)]
    config["transcription"] = {"enabled": explicit_rna, "max_transcripts": 5000}
    sim = Simulation(config)
    # The enhancer is close to bead 6 of the promoter interval, not gene.start=5.
    sim.positions[:8] = [[0, 0, 0], [20, 0, 0], [25, 0, 0], [30, 0, 0],
                         [35, 0, 0], [5, 0, 0], [0.5, 0, 0], [40, 0, 0]]
    sim.bindings = [link(0, 8, "chromatin"), link(0, 9, "chromatin"),
                    link(6, 11, "chromatin"),
                    link(9, 14), link(10, 15),
                    link(10, 20), link(20, 17), link(12, 18)]
    return sim


class FixedRandom:
    def __init__(self, value):
        self.value = value
        self.poisson_means = []

    def random(self, size=None):
        return self.value if size is None else np.full(size, self.value)

    def poisson(self, mean):
        self.poisson_means.append(mean)
        return 0

    def exponential(self, scale):
        return -np.log1p(-self.value) * scale


class RegulationTests(unittest.TestCase):
    def test_no_bound_tf_preserves_basal_stochastic_activation_in_both_rna_modes(self):
        for explicit in (False, True):
            with self.subTest(explicit_rna=explicit):
                sim = regulatory_system(explicit_rna=explicit)
                sim.bindings = []
                rates = sim.regulation_state()["genes"][0]
                self.assertEqual(rates["effective_on_rate"], 0.2)
                self.assertEqual(rates["effective_initiation_rate"], 3)
                self.assertFalse(sim.genes[0]["on"])
                sim.rng = FixedRandom(0)
                sim._advance_genes(0.2)
                self.assertTrue(sim.genes[0]["on"])
                self.assertAlmostEqual(sim.rng.poisson_means[-1], 0.6)

    def test_frozen_promoter_endpoint_statistics_match_analytic_two_state_solution(self):
        # An ensemble of independent promoter copies tests the observed
        # transition fraction against the CTMC solution, in both RNA paths.
        size, on_rate, off_rate, dt = 4096, 0.7, 1.3, 0.8
        for explicit in (False, True):
            for initially_on in (False, True):
                with self.subTest(explicit_rna=explicit, initially_on=initially_on):
                    config = small_config()
                    config["genes"] = [gene(id=f"g{i}", on_rate=on_rate, off_rate=off_rate,
                                            initiation_rate=0, processing_rate=0)
                                       for i in range(size)]
                    config["transcription"] = {"enabled": explicit, "max_transcripts": 5000}
                    sim = Simulation(config)
                    for state in sim.genes:
                        state["on"] = initially_on
                    sim.rng = np.random.default_rng(9501 + int(initially_on))
                    sim._advance_genes(dt)
                    switch_rate = off_rate if initially_on else on_rate
                    expected = switch_rate / (on_rate + off_rate) * (
                        1 - np.exp(-(on_rate + off_rate) * dt))
                    observed = np.mean([state["on"] != initially_on for state in sim.genes])
                    standard_error = np.sqrt(expected * (1 - expected) / size)
                    self.assertLess(abs(observed - expected), 6 * standard_error)

    def test_initiation_uses_actual_on_time_including_final_off_intervals(self):
        class OneSwitchRandom(FixedRandom):
            def __init__(self):
                super().__init__(0.5)
                self.waits = iter([0.05, np.inf])

            def exponential(self, scale):
                return next(self.waits)

        for explicit in (False, True):
            for initially_on in (False, True):
                with self.subTest(explicit_rna=explicit, initially_on=initially_on):
                    sim = regulatory_system(explicit_rna=explicit)
                    sim.bindings = []
                    sim.genes[0]["config"]["off_rate"] = 1
                    sim.genes[0]["on"] = initially_on
                    sim.rng = OneSwitchRandom()
                    sim._advance_genes(0.2)
                    expected_on_time = 0.05 if initially_on else 0.15
                    self.assertEqual(sim.genes[0]["on"], not initially_on)
                    self.assertEqual(len(sim.rng.poisson_means), 1)
                    self.assertAlmostEqual(sim.rng.poisson_means[0], 3 * expected_on_time)
                    self.assertAlmostEqual(sim.genes[0]["last_transcription_interval"]["on_time"], expected_on_time)
                    self.assertAlmostEqual(sim.events[-1]["time"], 0.05)

    def test_tf_motif_multiplier_is_separate_from_marks_with_max_overlap_and_zero(self):
        config = small_config()
        config["genes"] = [gene()]
        config["marks"] = [{"id": "A", "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}]
        config["proteins"] = [protein("TF", chromatin_affinity={"unmodified": 2, "A": 5}),
                              protein("other", chromatin_affinity={"unmodified": 3})]
        config["regulatory_regions"] = [region("zero", start=0, end=2, tf_affinity_multiplier=0),
                                        region("positive", start=2, end=4, tf_affinity_multiplier=3),
                                        region("disabled", start=0, end=7, tf_affinity_multiplier=100,
                                               enabled=False)]
        sim = Simulation(config)
        sim.marks["A"][3] = True
        np.testing.assert_array_equal(sim._chromatin_affinity("TF"), [0, 0, 6, 15, 6, 2, 2, 2])
        np.testing.assert_array_equal(sim._chromatin_affinity("other"), np.full(8, 3))
        self.assertEqual(sim.genome_tracks()["affinity_runs"]["TF"],
                         [[0, 2, 0.0], [2, 3, 6.0], [3, 4, 15.0], [4, 5, 6.0], [5, 8, 2.0]])

    def test_recruitment_requires_local_direct_tf_anchor_and_one_actual_protein_hop(self):
        sim = regulatory_system()
        report = sim.regulation_state()["regions"][0]
        self.assertEqual(report["tf_molecules"], [0])
        self.assertEqual(report["tf_bound_beads"], [0])
        self.assertEqual(report["recruited_molecules"], [2])
        self.assertEqual(report["recruited_anchors"], [{"molecule": 2, "beads": [0]}])
        self.assertEqual(report["tf_count"], 1)
        self.assertEqual(report["recruited_count"], 1)
        # Molecule 3 is connected through a relay, and to TF 1 outside region.
        # Moving it close cannot create a recruitment bond.
        sim.positions[17:20] = sim.positions[0]
        sim.bindings = [bond for bond in sim.bindings if not (
            bond["kind"] == "protein" and bond["b"] in (14, 15))]
        self.assertEqual(sim.regulation_state()["regions"][0]["recruited_count"], 0)

    def test_recruited_writer_deduplicates_domains_regions_and_own_direct_attachment(self):
        for directly_attached in (False, True):
            with self.subTest(directly_attached=directly_attached):
                sim = regulatory_system()
                # An overlapping enhancer and two IDR domain links must still
                # contribute one catalytic hazard for molecule 2 at bead 0.
                duplicate = region("overlapping")
                sim.regulatory_regions.append(duplicate)
                sim._enabled_regions.append(duplicate)
                if directly_attached:
                    sim.bindings.append(link(0, 16, "chromatin"))
                sim.rng = FixedRandom(0.75)
                sim._advance_marks(1)
                self.assertFalse(sim.marks["A"].any())  # 1-exp(-1)<.75<1-exp(-2).
                sim.rng = FixedRandom(0.5)
                sim._advance_marks(1)
                np.testing.assert_array_equal(sim.marks["A"], [True, False, False, False,
                                                             False, False, False, False])

    def test_recruited_eraser_targets_tf_anchor_without_chromatin_contact(self):
        sim = regulatory_system()
        sim.species["CO"]["writer_mark"] = None
        sim.species["CO"]["eraser_mark"] = "A"
        sim.marks["A"][:] = True
        sim.positions[5:7] = [10, 10, 10]  # No enhancer-promoter contact.
        self.assertFalse(sim.regulation_state()["regions"][0]["contact"])
        sim.rng = FixedRandom(0.5)
        sim._advance_marks(1)
        self.assertFalse(sim.marks["A"][0])
        self.assertTrue(sim.marks["A"][1])

    def test_region_rates_use_minimum_promoter_distance_mark_gate_and_replace_legacy_gain(self):
        sim = regulatory_system()
        report = sim.regulation_state()
        enhancer, promoter = report["regions"]
        self.assertTrue(enhancer["contact"])
        self.assertEqual(enhancer["distance"], 0.5)
        self.assertEqual(enhancer["target_promoter"], {"start": 5, "end": 6, "region_id": "promoter"})
        self.assertTrue(promoter["contact"])
        self.assertEqual(promoter["distance"], 0)
        self.assertEqual(enhancer["signal"], 3)
        self.assertEqual(enhancer["activation_contribution"], 6)
        self.assertEqual(enhancer["initiation_contribution"], 12)
        rates = report["genes"][0]
        self.assertEqual(rates["mode"], "regions")
        self.assertEqual(rates["legacy_activator_count"], 0)
        self.assertAlmostEqual(rates["effective_on_rate"], 1.4)
        self.assertEqual(rates["effective_initiation_rate"], 39)
        spec = sim.genes[0]["config"]
        spec["required_mark"] = "A"
        self.assertEqual(sim._gene_on_rate(spec), 0)
        sim.marks["A"][5] = True
        self.assertAlmostEqual(sim._gene_on_rate(spec), 1.4)
        spec["repressive_mark"] = "R"
        sim.marks["R"][5] = True
        self.assertEqual(sim._gene_on_rate(spec), 0)
        self.assertEqual(sim.regulation_state()["genes"][0]["effective_initiation_rate"], 39)

    def test_actual_initiation_uses_regional_rate_and_lost_contact_does_not_force_off(self):
        sim = regulatory_system()
        sim.rng = FixedRandom(0.99)
        sim.genes[0]["on"] = True
        sim._advance_genes(0.2)
        self.assertAlmostEqual(sim.rng.poisson_means[-1], 7.8)
        sim.positions[5:7] = [10, 10, 10]
        sim._advance_genes(0.2)
        self.assertTrue(sim.genes[0]["on"])
        self.assertAlmostEqual(sim.rng.poisson_means[-1], 0.6)
        self.assertEqual(sim.regulation_state()["genes"][0]["activation_multiplier"], 1)

    def test_disabled_and_absent_regions_have_identical_seeded_legacy_dynamics(self):
        config = small_config()
        config["genes"] = [gene(activator_species="TF")]
        config["proteins"] = [protein("TF", count=2, on_rate=2,
                                      chromatin_affinity={"unmodified": 1})]
        disabled = copy.deepcopy(config)
        disabled["regulatory_regions"] = [region(enabled=False, tf_affinity_multiplier=0)]
        absent, inactive = Simulation(config), Simulation(disabled)
        absent.step(80)
        inactive.step(80)
        np.testing.assert_array_equal(absent.positions, inactive.positions)
        self.assertEqual(absent.bindings, inactive.bindings)
        self.assertEqual(absent.genes, inactive.genes)
        self.assertEqual(absent.rng.bit_generator.state, inactive.rng.bit_generator.state)
        self.assertEqual(list(absent.events), list(inactive.events))
        self.assertEqual(absent.regulation_state()["genes"], inactive.regulation_state()["genes"])
        self.assertEqual(inactive.regulation_state()["genes"][0]["mode"], "legacy")

    def test_reports_refresh_from_current_state_without_mutation_or_rng_draws(self):
        sim = regulatory_system()
        rng, positions, bindings = copy.deepcopy(sim.rng.bit_generator.state), sim.positions.copy(), copy.deepcopy(sim.bindings)
        report = sim.snapshot()
        self.assertEqual(report["regulation"]["genes"][0], report["genes"][0]["regulation"])
        json.dumps(report, allow_nan=False)
        self.assertEqual(sim.rng.bit_generator.state, rng)
        np.testing.assert_array_equal(sim.positions, positions)
        self.assertEqual(sim.bindings, bindings)
        report["regulation"]["regions"][0]["recruited_anchors"][0]["beads"].append(999)
        self.assertEqual(sim.regulation_state()["regions"][0]["recruited_anchors"][0]["beads"], [0])
        sim.marks["A"][0] = True
        self.assertEqual(sim.regulation_state()["regions"][0]["mark_fractions"]["A"], 0.5)
        sim.positions[5:7] = [10, 10, 10]
        self.assertFalse(sim.regulation_state()["regions"][0]["contact"])

    def test_chromatin_mark_state_changes_force_energy_report_and_stability(self):
        config = small_config()
        config["simulation"]["dt"] = 0.1
        config["marks"] = [{"id": "A", "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}]
        config["chromatin"]["state_interactions"] = [
            {"a": "A", "b": "A", "strength": 20, "radius": 1.5, "min_separation": 3}]
        sim = Simulation(config)
        sim.positions[:] = np.column_stack((np.arange(8), np.zeros(8), np.zeros(8)))
        sim.positions[4] = [0.6, 0.5, 0]
        baseline_force, baseline_energy = sim.get_forces(), sim.potential_energy()
        baseline_substeps = sim.integration_substeps()
        self.assertEqual(sim.snapshot()["chromatin_contacts"]["total_contacts"], 0)
        sim.marks["A"][[0, 4]] = True
        report = sim.snapshot()["chromatin_contacts"]
        self.assertEqual(report["total_contacts"], 1)
        self.assertAlmostEqual(sim.potential_energy(), baseline_energy + report["energy"])
        force = sim.get_forces()
        self.assertGreater(force[0, 0] - baseline_force[0, 0], 0)
        np.testing.assert_allclose((force - baseline_force).sum(axis=0), 0, atol=1e-12)
        self.assertGreater(sim.integration_substeps(), baseline_substeps)
        numerical = np.zeros_like(force)
        h = 1e-6
        for bead in range(8):
            for axis in range(3):
                sim.positions[bead, axis] += h
                plus = sim.potential_energy()
                sim.positions[bead, axis] -= 2 * h
                minus = sim.potential_energy()
                sim.positions[bead, axis] += h
                numerical[bead, axis] = -(plus - minus) / (2 * h)
        np.testing.assert_allclose(force, numerical, atol=2e-7, rtol=2e-7)
        sim.marks["A"][:] = False
        np.testing.assert_allclose(sim.get_forces(), baseline_force, atol=1e-12)
        self.assertEqual(sim.snapshot()["chromatin_contacts"]["total_contacts"], 0)

    def test_openmm_uploads_new_mark_dependent_force_on_next_step(self):
        from chromatin_studio.openmm_backend import PolychromDynamics, backend_status
        if not backend_status()["available"]:
            self.skipTest("optional polychrom/OpenMM unavailable")
        config = small_config()
        config["marks"] = [{"id": "A", "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}]
        config["chromatin"]["state_interactions"] = [
            {"a": "A", "b": "A", "strength": 0.4, "radius": 1.5, "min_separation": 3}]
        sim = Simulation(config)
        sim.positions[:] = np.column_stack((np.arange(8), np.zeros(8), np.zeros(8)))
        sim.positions[4] = [0.6, 0.5, 0]
        sim.diffusion[:] = 0  # Isolate force upload from stochastic displacement.
        dynamics = PolychromDynamics(sim, seed=391, platform="Reference")
        expected = sim.get_forces()
        dynamics(sim)
        uploaded = np.asarray(dynamics.integrator.getPerDofVariableByName("model_force"), dtype=float)
        np.testing.assert_allclose(uploaded, expected)
        sim.marks["A"][[0, 4]] = True
        marked = sim.get_forces()
        self.assertFalse(np.array_equal(marked, expected))
        dynamics(sim)
        uploaded = np.asarray(dynamics.integrator.getPerDofVariableByName("model_force"), dtype=float)
        np.testing.assert_allclose(uploaded, marked)


if __name__ == "__main__":
    unittest.main()
