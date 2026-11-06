"""Molecule-level topology and direct-versus-indirect chromatin binding reports."""

import copy
import json
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from tests.test_engine import small_config, protein


def link(left, right, kind="protein"):
    return {"a": left, "b": right, "kind": kind, "off_rate": 0, "k": 1, "rest": 1}


def interacting_system():
    config = small_config()
    config["proteins"] = [protein("P", count=4, representation="idr", idr_beads=3, valence=6),
                          protein("Q", count=2, valence=4), protein("absent", count=0)]
    config["interactions"] = [{"a": a, "b": b, "on_rate": 1, "off_rate": 1,
                               "radius": 2, "strength": 1}
                              for a, b in (("P", "P"), ("Q", "P"), ("Q", "absent"))]
    sim = Simulation(config)
    # Molecules 0..3 are P with particles 8..19; molecules 4..5 are Q with 20..21.
    sim.bindings = [link(0, 8, "chromatin"), link(0, 9, "chromatin"),
                    link(7, 10, "chromatin"), link(1, 11, "chromatin"),
                    link(8, 11), link(9, 12), link(10, 14),
                    link(15, 20), link(16, 20), link(17, 21)]
    return sim


class InteractionStateTests(unittest.TestCase):
    def test_domain_links_deduplicate_molecule_pairs_and_direct_binding(self):
        sim = interacting_system()
        state = sim.snapshot()["interaction_state"]
        self.assertEqual(state["totals"], {"molecules": 6, "direct_bound": 2,
                         "indirect_bound": 2, "free": 2, "protein_links": 6,
                         "molecule_pairs": 4, "components": 2, "assemblies": 2})
        species = {item["id"]: item for item in state["species"]}
        self.assertEqual(species["P"], {"id": "P", "total": 4, "direct_bound": 2,
                                       "indirect_bound": 1, "free": 1})
        self.assertEqual(species["Q"], {"id": "Q", "total": 2, "direct_bound": 0,
                                       "indirect_bound": 1, "free": 1})
        self.assertEqual(species["absent"]["total"], 0)
        self.assertEqual([item["binding_status"] for item in state["molecules"]],
                         ["direct", "direct", "indirect", "free", "indirect", "free"])
        self.assertEqual(state["molecules"][0]["chromatin_beads"], [0, 7])
        self.assertEqual(state["molecules"][0]["used_valence"], 6)
        self.assertEqual(state["molecules"][0]["particle_ids"], [8, 9, 10])
        pair = next(item for item in state["pairs"] if item["a"] == item["b"] == "P")
        self.assertEqual(pair["links"], {"total": 3, "free_free": 0,
                                         "bound_free": 1, "bound_bound": 2})
        self.assertEqual(pair["molecule_pairs"], {"total": 2, "free_free": 0,
                                                  "bound_free": 1, "bound_bound": 1})

    def test_transitive_tether_is_distinct_from_direct_pair_context(self):
        state = interacting_system().interaction_state()
        pair = next(item for item in state["pairs"] if (item["a"], item["b"]) == ("P", "Q"))
        # The 2--4 edge has BOTH partners indirectly tethered. The 3--5 edge
        # belongs to a detached complex. Neither edge has a directly bound end.
        self.assertEqual(pair["links"], {"total": 3, "free_free": 3,
                                         "bound_free": 0, "bound_bound": 0})
        self.assertEqual(pair["molecule_pairs"]["free_free"], 2)
        self.assertEqual(state["components"][0], {
            "id": 0, "molecule_ids": [0, 1, 2, 4], "size": 4,
            "species_counts": {"P": 3, "Q": 1}, "direct_bound": 2,
            "indirect_bound": 2, "chromatin_bound": True, "protein_links": 5,
            "unique_molecule_pairs": 3, "chromatin_beads": [0, 1, 7]})
        self.assertFalse(state["components"][1]["chromatin_bound"])
        self.assertEqual(state["components"][1]["id"], 3)
        self.assertEqual(state["molecule_links"], [
            {"a": 0, "b": 1, "links": 2, "context": "bound_bound"},
            {"a": 0, "b": 2, "links": 1, "context": "bound_free"},
            {"a": 2, "b": 4, "links": 2, "context": "free_free"},
            {"a": 3, "b": 5, "links": 1, "context": "free_free"}])

    def test_zero_configured_pairs_and_existing_unconfigured_pairs_are_visible(self):
        sim = interacting_system()
        sim.bindings.append(link(20, 21))  # Existing Q--Q link absent from rules.
        state = sim.interaction_state()
        pairs = {(item["a"], item["b"]): item for item in state["pairs"]}
        self.assertTrue(pairs[("Q", "absent")]["configured"])
        self.assertTrue(all(value == 0 for value in pairs[("Q", "absent")]["links"].values()))
        self.assertFalse(pairs[("Q", "Q")]["configured"])
        self.assertEqual(pairs[("Q", "Q")]["molecule_pairs"]["free_free"], 1)
        self.assertEqual(state["totals"]["components"], 1)
        self.assertEqual(state["totals"]["indirect_bound"], 4)
        self.assertEqual(state["totals"]["free"], 0)

    def test_shared_chromatin_is_not_a_protein_edge_and_unbinding_updates_topology(self):
        sim = interacting_system()
        sim.bindings = [link(0, 8, "chromatin"), link(0, 11, "chromatin")]
        state = sim.interaction_state()
        self.assertEqual(state["totals"]["components"], 6)
        self.assertEqual(state["totals"]["assemblies"], 0)
        self.assertEqual(state["totals"]["direct_bound"], 2)
        self.assertEqual(state["totals"]["indirect_bound"], 0)
        sim.bindings.extend([link(8, 14), link(15, 20)])
        self.assertEqual(sim.interaction_state()["totals"]["indirect_bound"], 2)
        sim.bindings = [bond for bond in sim.bindings if bond["kind"] == "protein"]
        detached = sim.interaction_state()
        self.assertEqual(detached["totals"]["free"], 6)
        self.assertEqual(detached["totals"]["direct_bound"], 0)
        self.assertEqual(detached["totals"]["indirect_bound"], 0)
        self.assertEqual(detached["totals"]["assemblies"], 1)

    def test_empty_system_and_zero_population_species_have_valid_reports(self):
        sim = Simulation(small_config())
        state = sim.interaction_state()
        for key in ("species", "molecules", "pairs", "components", "molecule_links"):
            self.assertEqual(state[key], [])
        self.assertTrue(all(value == 0 for value in state["totals"].values()))
        config = small_config()
        config["proteins"] = [protein("P", count=0)]
        config["interactions"] = [{"a": "P", "b": "P", "on_rate": 0,
                                   "off_rate": 0, "radius": 2, "strength": 1}]
        state = Simulation(config).interaction_state()
        self.assertEqual(len(state["pairs"]), 1)
        self.assertEqual(state["species"], [{"id": "P", "total": 0,
                                          "direct_bound": 0, "indirect_bound": 0, "free": 0}])
        self.assertEqual(state["pairs"][0]["molecule_pairs"]["total"], 0)

    def test_order_independence_detached_output_and_read_only_snapshot(self):
        sim = interacting_system()
        # Cycle formation must not duplicate members or depend on union order.
        sim.bindings.append(link(12, 16))
        rng = copy.deepcopy(sim.rng.bit_generator.state)
        positions = sim.positions.copy()
        original_bindings = copy.deepcopy(sim.bindings)
        original_events = list(sim.events)
        counters = (sim.time, sim.step_count)
        first = sim.interaction_state()
        second = sim.snapshot()["interaction_state"]
        self.assertEqual(first, second)
        json.dumps(second, allow_nan=False)
        self.assertEqual(sim.rng.bit_generator.state, rng)
        np.testing.assert_array_equal(sim.positions, positions)
        self.assertEqual(sim.bindings, original_bindings)
        self.assertEqual(list(sim.events), original_events)
        self.assertEqual((sim.time, sim.step_count), counters)
        sim.bindings.reverse()
        for bond in sim.bindings:
            if bond["kind"] == "protein":
                bond["a"], bond["b"] = bond["b"], bond["a"]
        self.assertEqual(sim.interaction_state(), second)
        first["molecules"][0]["particle_ids"].append(999)
        first["components"][0]["molecule_ids"].append(999)
        first["species"][0]["total"] = 999
        self.assertEqual(sim.interaction_state(), second)


if __name__ == "__main__":
    unittest.main()
