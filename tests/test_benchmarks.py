"""Checks for independently derived Gaussian-chain reference benchmarks."""
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from chromatin_studio import benchmarks as bench
from chromatin_studio.engine import Simulation


class GaussianReferenceTests(unittest.TestCase):
    def test_open_chain_eigenvalues_and_translational_invariance(self):
        n = 19
        matrix = bench.graph_laplacian(n)
        expected = 2-2*np.cos(np.arange(n)*np.pi/n)
        np.testing.assert_allclose(np.linalg.eigvalsh(matrix), expected, atol=1e-14)
        np.testing.assert_allclose(matrix@np.ones(n), 0, atol=1e-14)
        looped = bench.graph_laplacian(n, [(2, 15), (5, 12)])
        self.assertTrue(np.all(np.linalg.eigvalsh(looped-matrix) >= -1e-12))
        for loops in ([(1, 2)], [(1, n)]):
            with self.assertRaises(ValueError):
                bench.graph_laplacian(n, loops)

    def test_two_bead_exact_result_and_zero_lag(self):
        times, diffusion, stiffness = np.array([0, .2, 2.]), .7, 1.3
        matrix = bench.graph_laplacian(2)
        expected = 3/(2*stiffness)*(1-np.exp(-2*diffusion*stiffness*times))
        np.testing.assert_allclose(bench.eigenmode_msd(matrix, times, diffusion, stiffness), expected)
        full = bench.eigenmode_msd(matrix, times, diffusion, stiffness, remove_com=False)
        np.testing.assert_allclose(full-expected, 3*diffusion*times)

    def test_finite_step_prediction_converges_and_rejects_instability(self):
        matrix = bench.graph_laplacian(32, [(3, 19)])
        times = np.array([.2, .8, 3.2])
        exact = bench.eigenmode_msd(matrix, times)
        coarse = bench.eigenmode_msd(matrix, times, dt=.02)
        fine = bench.eigenmode_msd(matrix, times, dt=.01)
        self.assertLess(np.max(abs(fine-exact)), np.max(abs(coarse-exact)))
        with self.assertRaises(ValueError):
            bench.eigenmode_msd(matrix, times, dt=1.)
        with self.assertRaises(ValueError):
            bench.eigenmode_msd(matrix, np.array([.123]), dt=.01)

    def test_independent_edge_forces_match_matrix_and_conserve_com(self):
        sim, matrix = bench._gaussian_fixture(12, 4, [(1, 8), (3, 10)], .7, 2.3, .01, 78)
        points = sim.positions.reshape(4, 12, 3)
        force = sim.get_forces().reshape(4, 12, 3)
        np.testing.assert_allclose(force, -2.3*np.einsum('ij,rjd->rid', matrix, points), atol=1e-13)
        np.testing.assert_allclose(force.sum(axis=1), 0, atol=1e-13)
        self.assertGreater(np.std(points), 0)

    def test_app_integrator_is_exercised_and_wrong_noise_is_detected(self):
        original = Simulation._brownian_move
        with patch.object(Simulation, '_brownian_move', side_effect=original) as move:
            measured = bench._sample_free(4096, 1., .01, 301, [1, 2, 4])
            self.assertEqual(move.call_count, 4)
        target = 6*np.array([.01, .02, .04])
        good = bench._curve('free', 'Free', 'free', [.01, .02, .04], measured, target, target)
        self.assertTrue(bench._curve_agreement(good)['passed'])
        def inflated_noise(sim):
            before = sim.positions.copy()
            original(sim)
            sim.positions[:] = before+1.4*(sim.positions-before)
        with patch.object(Simulation, '_brownian_move', side_effect=inflated_noise):
            wrong = bench._sample_free(4096, 1., .01, 301, [1, 2, 4])
        bad = bench._curve('free', 'Free', 'free', [.01, .02, .04], wrong, target, target)
        self.assertFalse(bench._curve_agreement(bad)['passed'])

    def test_full_reference_run_measures_expected_behavior_and_finite_json(self):
        report = bench.run_benchmarks()
        self.assertEqual(report['status'], 'passed')
        self.assertEqual(report['summary']['passed'], report['summary']['total'])
        self.assertEqual(len(report['curves']), 5)
        self.assertFalse(report['provenance']['experimental_data_used'])
        self.assertFalse(report['provenance']['full_paper_replication'])
        self.assertEqual(report['provenance']['integrator'], 'Simulation._brownian_move')
        self.assertLess(report['summary']['loop_suppression_observed'], .9)
        self.assertLess(abs(report['summary']['rouse_slope_observed']-.5), .1)
        json.dumps(report, allow_nan=False)
        for curve in report['curves']:
            times = [p['lag'] for p in curve['points']]
            self.assertEqual(times, sorted(set(times)))
            self.assertEqual(times[-1], report['parameters']['duration'])
            self.assertTrue(all(p['observed'] > 0 and p['standard_error'] > 0 for p in curve['points']))

    def test_cli_preserves_existing_report_and_parameters_are_bounded(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'report.json'
            path.write_text('{"keep": true}')
            with self.assertRaises(SystemExit) as error:
                bench.main(['--output', str(path)])
            self.assertEqual(error.exception.code, 2)
            self.assertEqual(json.loads(path.read_text()), {'keep': True})
        for options in ({'replicas': 1}, {'beads': True}, {'duration': math.inf}, {'dt': 0}, {'seed': -1}):
            with self.assertRaises(ValueError):
                bench.run_benchmarks(**options)


if __name__ == '__main__':
    unittest.main()
