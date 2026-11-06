"""Neighbor-search equivalence and a finite 20 Mb integration check."""

import copy
import unittest
from unittest.mock import patch

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.neighbors import GrowingNeighbors, neighbor_pairs
from chromatin_studio.schema import preset_config
from tests.test_engine import small_config, protein


def dense_pairs(positions, cutoff, tree=None):
    """Independent all-pairs reference, used only for small test systems."""
    left, right = np.triu_indices(len(positions), 1)
    delta = positions[left] - positions[right]
    inside = np.sum(delta * delta, axis=1) < cutoff * cutoff
    return np.column_stack((left[inside], right[inside])) if cutoff > 0 else np.empty((0, 2), dtype=int)


class DenseTree:
    """Small-system reference for radius queries in reaction-equivalence tests."""

    def __init__(self, positions):
        self.positions = np.array(positions, copy=True)

    def query_ball_point(self, points, radii, return_sorted=True):
        return [np.flatnonzero(np.sum((self.positions - point) ** 2, axis=1) <= radius * radius).tolist()
                for point, radius in zip(points, radii)]


class NeighborEquivalenceTests(unittest.TestCase):
    def test_neighbor_pairs_match_dense_with_exact_boundary_and_coincident_points(self):
        positions = np.random.default_rng(123).normal(size=(90, 3))
        positions[:4] = [[0, 0, 0], [1, 0, 0], [0, 0, 0], [0, 2, 0]]
        for cutoff in (0, 0.2, 1, 2.5):
            np.testing.assert_array_equal(neighbor_pairs(positions, cutoff), dense_pairs(positions, cutoff))

    def test_spatial_hash_overlap_checks_and_backtracking_match_dense(self):
        rng = np.random.default_rng(66)
        positions = rng.normal(size=(80, 3)) * 3
        radii = rng.uniform(0.1, 0.7, len(positions))
        growing = GrowingNeighbors(0.8, 0.7)
        for i, (point, radius) in enumerate(zip(positions, radii)):
            growing.add(i, point, radius)
        for _ in range(120):
            point, radius = rng.normal(size=3) * 3, float(rng.uniform(0.1, 0.7))
            expected = bool(np.all(np.sum((positions - point) ** 2, axis=1) > (0.65 * (radii + radius)) ** 2))
            self.assertEqual(growing.is_clear(point, radius), expected)
        for i in range(40, 80):
            growing.remove(i)
        for i in range(40, 80):
            expected = bool(np.all(np.sum((positions[:40] - positions[i]) ** 2, axis=1) > (0.65 * (radii[:40] + radii[i])) ** 2))
            self.assertEqual(growing.is_clear(positions[i], radii[i]), expected)

    def test_all_force_and_energy_results_match_dense_search(self):
        config = small_config()
        config["chromatin"]["repulsion"] = 7
        config["proteins"] = [protein("A", count=3, radius=0.2),
                              protein("B", count=2, radius=0.5, representation="idr", idr_beads=3)]
        sim = Simulation(config)
        sim.positions[:] = np.random.default_rng(4).normal(size=sim.positions.shape) * 0.9
        forces, energy = sim._mechanics()
        with patch("chromatin_studio.engine.neighbor_pairs", dense_pairs):
            reference_forces, reference_energy = sim._mechanics()
        np.testing.assert_allclose(forces, reference_forces, atol=1e-12, rtol=1e-12)
        self.assertAlmostEqual(energy, reference_energy, places=12)
        # Check the actual heterogeneous overlap list and backbone exclusions too.
        pairs, sigma2 = sim._repulsion_pairs()
        all_pairs = dense_pairs(sim.positions, 2 * float(sim.radii.max()))
        excluded = set(sim._excluded_pair_codes.tolist())
        expected = []
        for a, b in all_pairs:
            if int(a) * len(sim.positions) + int(b) in excluded:
                continue
            if np.sum((sim.positions[a] - sim.positions[b]) ** 2) < (sim.radii[a] + sim.radii[b]) ** 2:
                expected.append([a, b])
        np.testing.assert_array_equal(pairs, np.asarray(expected))
        self.assertEqual(len(sigma2), len(expected))

    def test_seeded_binding_and_mark_reactions_match_dense_search(self):
        config = small_config()
        config["chromatin"]["beads"] = 32
        config["marks"] = [{"id": "M", "decay_rate": 0.3, "spread_rate": 2, "contact_radius": 0.8},
                           {"id": "Q", "decay_rate": 0.1, "spread_rate": 0.7, "contact_radius": 1.4}]
        config["chromatin"]["initial_marks"] = [{"start": 0, "end": 10, "mark": "M"},
                                                {"start": 18, "end": 25, "mark": "Q"}]
        config["proteins"] = [protein("A", count=5, on_rate=6, off_rate=0.1, writer_mark="M", catalytic_rate=3,
                                      chromatin_affinity={"unmodified": 1, "M": 4}),
                              protein("B", count=4, on_rate=4, off_rate=0.2, representation="idr", idr_beads=3,
                                      chromatin_affinity={"unmodified": 2, "Q": 0.5})]
        config["interactions"] = [{"a": "A", "b": "A", "radius": 1.6, "on_rate": 5, "off_rate": 0.1, "strength": 3},
                                  {"a": "A", "b": "B", "radius": 1.1, "on_rate": 8, "off_rate": 0.3, "strength": 2}]
        actual, reference = Simulation(config), Simulation(copy.deepcopy(config))
        points = np.random.default_rng(910).normal(size=actual.positions.shape) * 0.7
        actual.positions[:] = reference.positions[:] = points
        for _ in range(5):
            actual.advance_reactions(0.1)
            with patch("chromatin_studio.engine.cKDTree", DenseTree), patch("chromatin_studio.engine.neighbor_pairs", dense_pairs):
                reference.advance_reactions(0.1)
            self.assertEqual(actual.bindings, reference.bindings)
            left, right = actual.snapshot(), reference.snapshot()
            left["genomic_events"].pop("run_id")
            right["genomic_events"].pop("run_id")
            self.assertEqual(left, right)

    def test_twenty_mb_initialization_and_100_steps_are_finite(self):
        sim = Simulation(preset_config("twenty_mb"))
        self.assertEqual(sim.n_chromatin, 2000)
        self.assertEqual(sim.n_chromatin * sim.chrom["bp_per_bead"], 20_000_000)
        self.assertFalse(hasattr(sim, "_pair_i"))
        self.assertFalse(hasattr(sim, "_pair_j"))
        self.assertLess(len(sim._excluded_pair_codes), 2 * len(sim.positions))
        pairs, sigma2 = sim._repulsion_pairs()
        delta = sim.positions[pairs[:, 0]] - sim.positions[pairs[:, 1]]
        self.assertTrue(np.all(np.sum(delta * delta, axis=1) > 0.65 ** 2 * sigma2 - 1e-12))
        lengths = np.linalg.norm(np.diff(sim.positions[:sim.n_chromatin], axis=0), axis=1)
        np.testing.assert_allclose(lengths, sim.bond_length, atol=1e-12)
        sim.step(100)
        self.assertTrue(np.isfinite(sim.positions).all())
        self.assertTrue(np.isfinite(sim.potential_energy()))
        self.assertLess(float(np.linalg.norm(sim.positions, axis=1).max()), sim.radius + 2)
        self.assertEqual(sim.step_count, 100)


if __name__ == "__main__":
    unittest.main()
