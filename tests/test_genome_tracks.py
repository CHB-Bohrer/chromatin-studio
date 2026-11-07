"""Authoritative affinity/occupancy track format and read-only snapshot checks."""

import copy
import json
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.schema import preset_config
from tests.test_engine import small_config, protein


def decode_runs(runs, bead_count):
    values = np.zeros(bead_count)
    end_of_previous = 0
    for start, end, value in runs:
        if start != end_of_previous or not start < end <= bead_count:
            raise AssertionError("Track does not cover the bead axis contiguously.")
        values[start:end] = value
        end_of_previous = end
    if end_of_previous != bead_count:
        raise AssertionError("Track does not cover the final bead.")
    return values


class GenomeTrackTests(unittest.TestCase):
    def test_affinity_runs_reuse_mixed_mark_max_zero_and_fallback_semantics(self):
        config = small_config()
        config["chromatin"]["bp_per_bead"] = 1000
        config["marks"] = [{"id": m, "decay_rate": 0, "spread_rate": 0, "contact_radius": 1}
                           for m in ("A", "B", "unmapped")]
        config["proteins"] = [protein("P", chromatin_affinity={"unmodified": 2, "A": 5, "B": 0}),
                              protein("Q", chromatin_affinity={"B": 3}),
                              protein("absent", count=0, chromatin_affinity={})]
        sim = Simulation(config)
        sim.marks["A"][[1, 3, 5]] = True
        sim.marks["B"][[2, 3, 4, 5]] = True
        sim.marks["unmapped"][[0, 4, 7]] = True
        tracks = sim.snapshot()["genome_tracks"]
        self.assertEqual(tracks["bp_per_bead"], 1000)
        self.assertEqual(tracks["span_bp"], 8000)
        self.assertEqual(tracks["affinity_runs"]["P"],
                         [[0, 1, 2.0], [1, 2, 5.0], [2, 3, 0.0], [3, 4, 5.0],
                          [4, 5, 0.0], [5, 6, 5.0], [6, 8, 2.0]])
        self.assertEqual(tracks["affinity_runs"]["Q"], [[0, 2, 0.0], [2, 6, 3.0], [6, 8, 0.0]])
        self.assertEqual(tracks["affinity_runs"]["absent"], [[0, 8, 0.0]])
        for species, runs in tracks["affinity_runs"].items():
            np.testing.assert_array_equal(decode_runs(runs, 8), sim._chromatin_affinity(species))
        # Tracks are recomputed from current marks rather than cached initial labels.
        sim.marks["B"][7] = True
        self.assertEqual(sim.genome_tracks()["affinity_runs"]["P"][-1], [7, 8, 0.0])

    def test_occupancy_counts_distinct_directly_bound_molecules_per_bead(self):
        config = small_config()
        config["proteins"] = [protein("P", count=2, representation="idr", idr_beads=3, valence=4),
                              protein("Q", count=1), protein("empty", count=0)]
        sim = Simulation(config)
        sim.bindings = [
            {"a": 0, "b": 8, "kind": "chromatin"},
            {"a": 0, "b": 9, "kind": "chromatin"},  # Same IDR molecule as bead 8.
            {"a": 7, "b": 10, "kind": "chromatin"},  # Same molecule on another locus.
            {"a": 0, "b": 11, "kind": "chromatin"},
            {"a": 3, "b": 12, "kind": "chromatin"},
            {"a": 0, "b": 14, "kind": "chromatin"},  # Different species.
            {"a": 10, "b": 13, "kind": "protein"},  # No extra direct occupancy.
        ]
        tracks = sim.genome_tracks()
        self.assertEqual(tracks["occupancy_runs"]["P"],
                         [[0, 1, 2], [1, 3, 0], [3, 4, 1], [4, 7, 0], [7, 8, 1]])
        self.assertEqual(tracks["occupancy_runs"]["Q"], [[0, 1, 1], [1, 8, 0]])
        self.assertEqual(tracks["occupancy_runs"]["empty"], [[0, 8, 0]])
        # Affinity can be zero while existing bonds still provide occupancy.
        self.assertEqual(tracks["affinity_runs"]["P"], [[0, 8, 0.0]])
        self.assertTrue(all(isinstance(run[2], int) for run in tracks["occupancy_runs"]["P"]))

    def test_actual_valence_limited_binding_is_reported_without_altering_links(self):
        config = small_config()
        config["proteins"] = [protein("P", count=3, valence=2, on_rate=1e6,
                                      chromatin_affinity={"unmodified": 1})]
        sim = Simulation(config)
        sim.positions[:8] = np.column_stack((np.arange(8) * 0.2, np.zeros(8), np.zeros(8)))
        sim.positions[8:] = [0.7, 0.2, 0]
        sim._advance_binding(1)
        self.assertEqual(len(sim.bindings), 6)
        counts = decode_runs(sim.genome_tracks()["occupancy_runs"]["P"], 8)
        expected = np.zeros(8)
        by_molecule = {}
        for bond in sim.bindings:
            molecule = sim._particle_info[bond["b"]]["molecule"]
            by_molecule.setdefault(molecule, set()).add(bond["a"])
        for sites in by_molecule.values():
            self.assertEqual(len(sites), 2)
            for bead in sites:
                expected[bead] += 1
        np.testing.assert_array_equal(counts, expected)
        self.assertEqual(int(counts.sum()), 6)

    def test_track_reads_do_not_mutate_rng_state_or_simulation(self):
        sim = Simulation(preset_config("nuclear"))
        sim.step(10)
        rng_before = copy.deepcopy(sim.rng.bit_generator.state)
        positions_before = sim.positions.copy()
        marks_before = {key: values.copy() for key, values in sim.marks.items()}
        bindings_before = copy.deepcopy(sim.bindings)
        genes_before = copy.deepcopy(sim.genes)
        events_before = list(sim.events)
        counters_before = (sim.time, sim.step_count)
        first, second = sim.snapshot(), sim.snapshot()
        self.assertEqual(first, second)
        json.dumps(first, allow_nan=False)
        self.assertEqual(sim.rng.bit_generator.state, rng_before)
        np.testing.assert_array_equal(sim.positions, positions_before)
        for key, values in marks_before.items():
            np.testing.assert_array_equal(sim.marks[key], values)
        self.assertEqual(sim.bindings, bindings_before)
        self.assertEqual(sim.genes, genes_before)
        self.assertEqual(list(sim.events), events_before)
        self.assertEqual((sim.time, sim.step_count), counters_before)
        species = next(iter(sim.species))
        first["genome_tracks"]["affinity_runs"][species][0][2] = 9999
        self.assertNotEqual(sim.genome_tracks()["affinity_runs"][species][0][2], 9999)

    def test_twenty_mb_tracks_have_full_half_open_coverage(self):
        sim = Simulation(preset_config("twenty_mb"))
        tracks = sim.genome_tracks()
        self.assertEqual(tracks["bp_per_bead"], 10_000)
        self.assertEqual(tracks["span_bp"], 20_000_000)
        self.assertEqual(set(tracks["affinity_runs"]), set(sim.species))
        self.assertEqual(set(tracks["occupancy_runs"]), set(sim.species))
        for kind in ("affinity_runs", "occupancy_runs"):
            for species, runs in tracks[kind].items():
                self.assertEqual(runs[0][0], 0)
                self.assertEqual(runs[-1][1], 2000)
                self.assertEqual(runs[-1][1] * tracks["bp_per_bead"], tracks["span_bp"])
                decoded = decode_runs(runs, 2000)
                for left, right in zip(runs[:-1], runs[1:]):
                    self.assertNotEqual(left[2], right[2], "Adjacent equal runs should be merged.")
                if kind == "affinity_runs":
                    np.testing.assert_array_equal(decoded, sim._chromatin_affinity(species))
                else:
                    self.assertEqual(runs, [[0, 2000, 0]])


if __name__ == "__main__":
    unittest.main()
