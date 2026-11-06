"""Headless execution streams samples and records reproducible success/failure."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from chromatin_studio.batch import inspect_state, run_batch
from chromatin_studio.engine import Simulation
from chromatin_studio.schema import default_config, preset_config


def small_model():
    config = default_config()
    config["proteins"] = []
    config["interactions"] = []
    config["genes"] = []
    config["bodies"] = []
    config["cohesin"]["enabled"] = False
    return config


class BatchTests(unittest.TestCase):
    def test_samples_snapshots_final_state_and_no_overwrite(self):
        config = small_model()
        job = {"platform": "numpy", "replicas": 2, "steps": 7, "sample_every": 3, "snapshot_every": 4}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "replica"
            manifest = run_batch(config, job, replica=1, output=output)
            metrics = [json.loads(v) for v in (output / "metrics.jsonl").read_text().splitlines()]
            snapshots = [json.loads(v) for v in (output / "snapshots.jsonl").read_text().splitlines()]
            self.assertEqual([v["step"] for v in metrics], [0, 3, 6, 7])
            self.assertEqual([v["step"] for v in snapshots], [0, 4, 7])
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["seed"], config["seed"] + 1)
            self.assertEqual(manifest["completed_steps"], 7)
            self.assertIn("engine.py", manifest["source_sha256"])
            reference = Simulation({**config, "seed": config["seed"] + 1})
            reference.step(7)
            saved, expected = json.loads((output / "final-state.json").read_text()), reference.snapshot()
            saved["genomic_events"].pop("run_id")
            expected["genomic_events"].pop("run_id")
            self.assertEqual(saved, expected)
            with self.assertRaises(FileExistsError):
                run_batch(config, job, replica=1, output=output)

    def test_preflight_initializes_without_advancing_or_claiming_restart(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "preflight"
            manifest = run_batch(small_model(), {"platform": "numpy"}, replica=0, output=output, preflight=True)
            self.assertEqual(manifest["completed_steps"], 0)
            self.assertTrue(manifest["preflight_only"])
            self.assertFalse(manifest["restart_checkpoint"])
            self.assertFalse((output / "metrics.jsonl").exists())
            self.assertEqual(json.loads((output / "preflight-state.json").read_text())["step"], 0)

    def test_failed_requested_backend_is_recorded_without_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "failed"
            with patch("chromatin_studio.batch._create_simulation", side_effect=RuntimeError("No CUDA device")) as create:
                with self.assertRaisesRegex(RuntimeError, "No CUDA"):
                    run_batch(small_model(), {"platform": "CUDA"}, replica=0, output=output)
                self.assertEqual(create.call_count, 1)
            manifest = json.loads((output / "run-manifest.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["requested_job"]["platform"], "CUDA")
            self.assertFalse((output / "final-state.json").exists())

    def test_invariant_check_rejects_corrupted_rna_and_nonfinite_state(self):
        sim = Simulation(preset_config("nucleosome_transcription"))
        state = sim.snapshot()
        state["transcription"]["totals"]["initiated"] += 1
        with self.assertRaisesRegex(ValueError, "RNA mass balance"):
            inspect_state(sim, state)
        sim.positions[0, 0] = float("nan")
        with self.assertRaises((ValueError, FloatingPointError)):
            inspect_state(sim, sim.snapshot())

    def test_capacity_blocks_are_recorded_as_quantitative_manifest_warnings(self):
        config = preset_config("nucleosome_transcription")
        sim = Simulation(config)
        sim.explicit_transcription.blocked[config["genes"][0]["id"]] = 9
        with tempfile.TemporaryDirectory() as directory:
            with patch("chromatin_studio.batch._create_simulation", return_value=sim):
                manifest = run_batch(config, {"platform": "numpy"}, replica=0, output=Path(directory) / "capacity", preflight=True)
            self.assertEqual(manifest["rna_totals"]["blocked_initiations"], 9)
            self.assertIn("blocked 9 initiation", manifest["warnings"][0])

    def test_metrics_only_run_and_cli_invocation(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder / "model.json").write_text(json.dumps(small_model()))
            (folder / "job.json").write_text(json.dumps({"platform": "numpy", "replicas": 1, "steps": 2, "sample_every": 50, "snapshot_every": 0}))
            completed = subprocess.run([sys.executable, "-m", "chromatin_studio.batch", "--model", str(folder / "model.json"), "--job", str(folder / "job.json"), "--output", str(folder / "run")], capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(json.loads(completed.stdout)["steps"], 2)
            self.assertFalse((folder / "run" / "snapshots.jsonl").exists())
            self.assertEqual(len((folder / "run" / "metrics.jsonl").read_text().splitlines()), 2)


if __name__ == "__main__":
    unittest.main()
