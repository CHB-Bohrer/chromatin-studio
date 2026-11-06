"""Actual event observations, occupancy-dependent residence, and remodeling."""

import copy
import math
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.schema import preset_config, validate_config
from tests.test_engine import protein, small_config
from tests.test_interaction_state import link
from tests.test_nucleosomes import marked_config, nucleosomes
from tests.test_regulation import FixedRandom, region, regulatory_system


class OneJumpRandom:
    def __init__(self):
        self.waits = iter([.01, math.inf])
        self.scales = []

    def exponential(self, scale):
        self.scales.append(scale)
        return next(self.waits)

    def random(self):
        return .1


class LiveNucleosomeTests(unittest.TestCase):
    def test_current_occupancy_changes_actual_dissociation_survival(self):
        config = small_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 1, "end": 1, "occupied": False}])
        config["proteins"] = [protein("TF", count=2000, off_rate=.8,
            nucleosome_off_rate_multiplier={"occupied": 2, "vacant": .25})]
        sim = Simulation(config)
        sim.bindings = [dict(link(i % 2, 8 + i, "chromatin"), off_rate=.8) for i in range(2000)]
        sim.rng = np.random.default_rng(814)
        sim._advance_binding(.6)
        for bead, rate in ((0, 1.6), (1, .2)):
            observed = sum(bond["a"] == bead for bond in sim.bindings) / 1000
            expected = math.exp(-rate * .6)
            self.assertLess(abs(observed - expected), 5 * math.sqrt(expected * (1 - expected) / 1000))
        self.assertTrue(all(event["type"] == "chromatin_unbind" for event in sim.genomic_events))
        self.assertTrue(all(0 <= event["time"] < .6 for event in sim.genomic_events))
        # Changing occupancy changes the current hazard of an existing bond.
        bond = dict(link(0, 8, "chromatin"), off_rate=.8)
        self.assertEqual(sim._bond_off_rate(bond), 1.6)
        sim.nucleosome_occupied[0] = False
        self.assertEqual(sim._bond_off_rate(bond), .2)

    def test_protein_pair_dissociation_is_not_scaled_and_zero_is_allowed(self):
        config = small_config()
        config["chromatin"]["nucleosomes"] = nucleosomes()
        config["proteins"] = [protein(count=2, off_rate=1,
            nucleosome_off_rate_multiplier={"occupied": 0, "vacant": 20})]
        sim = Simulation(config)
        self.assertEqual(sim._bond_off_rate(dict(link(0, 8, "chromatin"), off_rate=1)), 0)
        self.assertEqual(sim._bond_off_rate(dict(link(8, 9), off_rate=1)), 1)

    def test_binding_and_turnover_events_report_only_observed_information(self):
        config = marked_config()
        config["proteins"] = [protein("TF", on_rate=1e6, off_rate=1,
                                       chromatin_affinity={"unmodified": 1})]
        sim = Simulation(config)
        sim.positions[:8] = np.column_stack((np.arange(8) * 10, np.zeros(8), np.zeros(8)))
        sim.positions[8] = [.1, 0, 0]
        sim._advance_binding(.1)
        event = sim.genomic_event_state()["events"][0]
        self.assertEqual((event["type"], event["bead"], event["species"], event["particle"]),
                         ("chromatin_bind", 0, "TF", 8))
        self.assertEqual(event["on_rate"], 1e6)
        self.assertEqual(event["timing"], "sampled_encounter")
        sim.nucleosome_config.update(eviction_rate=1, reassembly_rate=1)
        sim.rng = FixedRandom(.55)
        sim._advance_nucleosomes(1)
        resets = [e for e in sim.genomic_events if e["type"] == "nucleosome_reset"]
        self.assertEqual(len(resets), 8)
        self.assertTrue(all(e["timing"] == "interval_endpoint" and e["time"] == 1 for e in resets))
        self.assertTrue(sim.nucleosome_occupied.all())
        self.assertFalse(sim.marks["A"].any())

    def test_event_retention_and_snapshot_are_read_only_and_run_identity_is_unique(self):
        sim = Simulation(small_config())
        before = copy.deepcopy(sim.rng.bit_generator.state)
        for i in range(2100):
            sim._genomic_event("test", i % 8, .1)
        report = sim.genomic_event_state()
        self.assertEqual((report["capacity"], report["total_emitted"], report["dropped"]), (2048, 2100, 52))
        self.assertEqual((report["first_retained_id"], report["last_id"]), (53, 2100))
        self.assertEqual([e["id"] for e in report["events"]], list(range(53, 2101)))
        report["events"][0]["bead"] = -1
        self.assertGreaterEqual(sim.genomic_events[0]["bead"], 0)
        self.assertEqual(before, sim.rng.bit_generator.state)
        self.assertNotEqual(sim.run_id, Simulation(small_config()).run_id)

    def test_reporting_and_neutral_new_controls_preserve_seeded_dynamics(self):
        config = marked_config()
        config["proteins"] = [protein("TF", count=2, on_rate=2, off_rate=.5,
                                       chromatin_affinity={"unmodified": 1})]
        left = Simulation(config)
        neutral = copy.deepcopy(config)
        neutral["proteins"][0].update(nucleosome_off_rate_multiplier={"occupied": 1, "vacant": 1},
                                      nucleosome_sliding_rate=0)
        right = Simulation(neutral)
        right._genomic_event = lambda *args, **kwargs: None
        left.step(50)
        right.step(50)
        np.testing.assert_array_equal(left.positions, right.positions)
        self.assertEqual(left.bindings, right.bindings)
        self.assertEqual(left.rng.bit_generator.state, right.rng.bit_generator.state)
        before = copy.deepcopy(left.rng.bit_generator.state)
        tracks = left.genome_tracks()
        self.assertEqual(tracks["on_rate_runs"]["TF"], [[0, 8, 2.0]])
        self.assertEqual(tracks["off_rate_runs"]["TF"], [[0, 8, .5]])
        self.assertEqual(before, left.rng.bit_generator.state)

    def test_sliding_moves_marks_and_preserves_nucleosome_and_mark_counts(self):
        config = marked_config()
        config["chromatin"]["initial_marks"] = [{"start": 2, "end": 2, "mark": "A"}]
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 3, "end": 3, "occupied": False}])
        config["proteins"] = [protein("remodeler", nucleosome_sliding_rate=2)]
        sim = Simulation(config)
        sim.bindings = [link(2, 8, "chromatin")]
        sim.rng = OneJumpRandom()
        sim._advance_nucleosome_sliding(.1)
        self.assertEqual(sim.nucleosome_occupied.sum(), 7)
        self.assertFalse(sim.nucleosome_occupied[2])
        self.assertTrue(sim.nucleosome_occupied[3])
        self.assertEqual(np.flatnonzero(sim.marks["A"]).tolist(), [3])
        self.assertEqual(sim.nucleosome_state()["slides"], 1)
        event = sim.genomic_events[-1]
        self.assertEqual((event["type"], event["bead"], event["to_bead"], event["anchor_bead"]),
                         ("nucleosome_slide", 2, 3, 2))
        self.assertEqual(event["time"], .01)
        self.assertEqual(sim.rng.scales[0], .5)

    def test_sliding_requires_actual_anchor_and_handles_edges_without_wraparound(self):
        config = marked_config()
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 7, "end": 7, "occupied": False}])
        config["proteins"] = [protein("remodeler", nucleosome_sliding_rate=2)]
        sim = Simulation(config)
        before = copy.deepcopy(sim.rng.bit_generator.state)
        sim._advance_nucleosome_sliding(10)
        self.assertEqual(before, sim.rng.bit_generator.state)
        sim.bindings = [link(0, 8, "chromatin")]
        sim._advance_nucleosome_sliding(10)
        self.assertEqual(before, sim.rng.bit_generator.state)
        self.assertEqual(sim.nucleosome_state()["slides"], 0)

    def test_recruited_sliding_deduplicates_repeated_regions_and_domain_links(self):
        template = regulatory_system()
        config = copy.deepcopy(template.config)
        config["chromatin"]["nucleosomes"] = nucleosomes(
            initial_regions=[{"start": 1, "end": 1, "occupied": False},
                             {"start": 7, "end": 7, "occupied": False}])
        config["proteins"][1]["nucleosome_sliding_rate"] = 2
        config["regulatory_regions"].append(region("duplicate"))
        sim = Simulation(config)
        sim.bindings = copy.deepcopy(template.bindings)
        # CO molecule2 is recruited twice by TF domain links, and also bound directly.
        sim.bindings.append(link(0, 14, "chromatin"))
        anchors = [(m, b) for m, b in sim._remodeling_targets() if sim.molecules[m]["species"] == "CO"]
        self.assertEqual(anchors, [(2, 0)])
        sim.rng = OneJumpRandom()
        sim._advance_nucleosome_sliding(.1)
        self.assertEqual(sim.rng.scales[0], .5)
        self.assertEqual(sim.nucleosome_occupied.sum(), 6)

    def test_new_controls_validate_and_demo_contains_real_open_regions(self):
        config = preset_config("nucleosome_access")
        report = validate_config(config)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["genomic_span_bp"], 25600)
        sim = Simulation(config)
        self.assertEqual(np.count_nonzero(~sim.nucleosome_occupied), 11)
        self.assertTrue(all(not sim.nucleosome_occupied[b] for b in (10, 36, 72)))
        for field, invalid in (("nucleosome_off_rate_multiplier", {"other": 1}),
                               ("nucleosome_off_rate_multiplier", {"vacant": -1}),
                               ("nucleosome_off_rate_multiplier", {"occupied": True}),
                               ("nucleosome_sliding_rate", -1), ("nucleosome_sliding_rate", 21)):
            changed = copy.deepcopy(config)
            changed["proteins"][0][field] = invalid
            self.assertFalse(validate_config(changed)["valid"], field)


if __name__ == "__main__":
    unittest.main()
