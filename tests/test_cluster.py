"""Portable job validation, seed independence, and bundle integrity contracts."""

import copy
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from chromatin_studio.batch import verify_bundle
from chromatin_studio.cluster import build_bundle, normalize_job, plan_job, sha256
from chromatin_studio.schema import default_config


class ClusterTests(unittest.TestCase):
    def test_plan_uses_distinct_seeded_replicas_without_mutating_model(self):
        config = default_config()
        original = copy.deepcopy(config)
        plan = plan_job(config, {"platform": "CUDA", "replicas": 3, "steps": 123, "device_index": ""})
        self.assertEqual(config, original)
        self.assertEqual([v["seed"] for v in plan["replicas"]], [config["seed"] + i for i in range(3)])
        self.assertNotIn("device_index", plan["job"])
        self.assertEqual(plan["model"]["reduced_duration"], 123 * config["simulation"]["dt"])
        self.assertIn("#SBATCH --gres=gpu:1", plan["scripts"]["slurm"])
        self.assertIn("#SBATCH --array=0-2", plan["scripts"]["slurm"])
        self.assertLess(plan["scripts"]["shell"].index("verify_bundle('.')"), plan["scripts"]["shell"].index("python -m chromatin_studio.validation"))
        self.assertFalse(plan["resume_supported"])

    def test_large_output_requests_have_explicit_sampling_counts_and_warning(self):
        plan = plan_job(default_config(), {"steps": 100000, "sample_every": 2, "snapshot_every": 1, "replicas": 2})
        self.assertEqual(plan["output"]["full_snapshots_per_replica"], 100001)
        self.assertEqual(plan["output"]["metric_records_per_replica"], 50001)
        self.assertTrue(any("Large output request" in warning for warning in plan["warnings"]))

    def test_rejects_injection_bad_counts_and_seed_wrap(self):
        config = default_config()
        for field, value in (("partition", "gpu\n#SBATCH --export=ALL"), ("account", "$(touch x)"),
                             ("job_name", "--option"), ("walltime", "00:99:00"), ("walltime", "00:00:00"),
                             ("walltime", "1-24:00:00"), ("steps", True), ("replicas", 1.5),
                             ("steps", float("nan")), ("steps", 10**1000), ("snapshot_every", -1),
                             ("device_index", "0,1"), ("device_index", False), ("device_index", -1)):
            with self.subTest(field=field, value=str(value)[:80]), self.assertRaises(ValueError):
                normalize_job({field: value}, config)
        config["seed"] = 2**32 - 1
        with self.assertRaisesRegex(ValueError, "seed"):
            normalize_job({"replicas": 2}, config)
        with self.assertRaisesRegex(ValueError, "Unknown"):
            normalize_job({"unexpected": 1}, default_config())

    def test_cpu_jobs_never_request_gpu_and_devices_require_gpu(self):
        for platform in ("numpy", "CPU", "Reference"):
            plan = plan_job(default_config(), {"platform": platform})
            self.assertNotIn("--gres", plan["scripts"]["slurm"])
            self.assertNotIn("--device-index", plan["scripts"]["shell"])
            with self.assertRaisesRegex(ValueError, "only applicable"):
                normalize_job({"platform": platform, "device_index": 0}, default_config())
        normalized = normalize_job({"device_index": "02", "partition": "gpu-a", "account": "lab_1", "walltime": "2-03:04:05"}, default_config())
        self.assertEqual(normalized["device_index"], 2)

    def test_existing_particle_limits_still_apply(self):
        config = default_config()
        config["chromatin"]["beads"] = 100000
        with self.assertRaises(ValueError):
            plan_job(config, {"platform": "CUDA"})

    def test_bundle_is_portable_and_hashes_detect_tampering(self):
        raw = build_bundle(default_config(), {"platform": "numpy", "replicas": 1, "steps": 3})
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            names = archive.namelist()
            self.assertIn("chromatin_studio/batch.py", names)
            self.assertIn("chromatin_studio/validation.py", names)
            self.assertFalse(any(name.startswith("/") or ".." in Path(name).parts for name in names))
            self.assertFalse(any("__pycache__" in name for name in names))
            manifest = json.loads(archive.read("bundle-manifest.json"))
            for name, value in manifest["files"].items():
                self.assertEqual(sha256(archive.read(name)), value["sha256"])
            with tempfile.TemporaryDirectory() as directory:
                archive.extractall(directory)
                self.assertIsInstance(verify_bundle(directory), str)
                source = Path(directory) / "chromatin_studio" / "engine.py"
                source.write_text(source.read_text() + "\n# changed\n")
                with self.assertRaisesRegex(ValueError, "integrity"):
                    verify_bundle(directory)


if __name__ == "__main__":
    unittest.main()
