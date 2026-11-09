"""Scientific acceptance suite must reject wrong physics and device failures."""
import unittest
import json
from unittest.mock import patch

from chromatin_studio import validation


class ValidationTests(unittest.TestCase):
    def test_numpy_matches_independent_analytical_expectations(self):
        report = validation.run_validation()
        self.assertTrue(report["passed"], report)
        self.assertGreaterEqual(len(report["checks"]), 27)
        self.assertIn("source_sha256", report)

    def test_incorrect_diffusion_does_not_pass_acceptance(self):
        original = validation.Simulation._brownian_move

        def wrong_noise(sim):
            before = sim.positions.copy()
            original(sim)
            sim.positions[:] = before + 2*(sim.positions-before)

        with patch.object(validation.Simulation, "_brownian_move", wrong_noise):
            result = validation._mechanics("numpy", "mixed", None)
        diffusion = [check for check in result["checks"] if check["name"].startswith("free diffusion")]
        self.assertTrue(all(not check["passed"] for check in diffusion))

    def test_unavailable_requested_device_is_failure_not_skipped_success(self):
        original = validation._mechanics

        def unavailable(platform, precision, device):
            if platform == "CUDA":
                raise RuntimeError("test device unavailable")
            return original(platform, precision, device)

        with patch.object(validation, "_mechanics", side_effect=unavailable):
            report = validation.run_validation("CUDA")
        self.assertFalse(report["passed"])
        self.assertIn("test device unavailable", report["errors"][0])
        self.assertEqual(len(report["mechanics"]), 1)

    def test_nonfinite_observation_is_not_a_pass(self):
        check = validation._check("bad", float("nan"), 0, 1)
        self.assertFalse(check["passed"])
        json.dumps(check, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
