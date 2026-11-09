"""Physical mappings survive the UI/service/batch boundary without changing dynamics."""
import copy
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
import zipfile
from http.server import ThreadingHTTPServer

from chromatin_studio import server
from chromatin_studio.batch import run_batch
from chromatin_studio.calibration import fit_time_mapping
from chromatin_studio.cluster import build_bundle
from chromatin_studio.engine import Simulation
from chromatin_studio.schema import default_config, validate_config


def demo_payload(config):
    lags = [.01, .03, .1, .3, 1., 3., 10., 30., 100.]
    meta = {"dimension": 3, "observable": "locus", "drift_correction": "none"}
    return {"config": config,
            "simulation": {"lag": lags, "msd": [x**.5 for x in lags], **meta,
                           "source": "Generated test curve; not a measured trajectory"},
            "experiment": {"lag": [.2, .6, 2., 6., 20.],
                           "msd": [.01*x**.5 for x in [.1, .3, 1., 3., 10.]], **meta,
                           "source": "Synthetic test comparison; not an experiment"},
            "length_nm": 100., "length_source": "Synthetic test assumption",
            "fit_min_s": .2, "fit_max_s": 20., "source_kind": "synthetic", "notes": "Test only"}


class CalibrationIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.http = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.thread = threading.Thread(target=cls.http.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = "http://127.0.0.1:" + str(cls.http.server_port)

    @classmethod
    def tearDownClass(cls):
        cls.http.shutdown()
        cls.http.server_close()
        cls.thread.join()

    def setUp(self):
        server.WORKSPACE = server.Workspace()

    def request(self, path, body=None):
        req = urllib.request.Request(self.base + path,
              data=None if body is None else json.dumps(body).encode(),
              headers={"Content-Type": "application/json"})
        try:
            response = urllib.request.urlopen(req, timeout=10)
        except urllib.error.HTTPError as exc:
            response = exc
        with response:
            return response.status, json.load(response)

    def test_fitting_and_summary_do_not_change_live_run(self):
        self.request("/api/reset", {"config": default_config()})
        self.request("/api/step", {"steps": 3})
        _, before = self.request("/api/export")
        self.assertEqual(before["metadata"]["physical_mapping"]["status"], "uncalibrated")
        status, result = self.request("/api/calibration/fit", demo_payload(default_config()))
        self.assertEqual(status, 200, result)
        self.assertEqual(result["profile"]["status"], "demo")
        self.assertAlmostEqual(result["profile"]["time_s"], 2.)
        draft = default_config()
        draft["calibration"] = result["profile"]
        status, summary = self.request("/api/calibration/summary", {"config": draft})
        self.assertEqual(status, 200, summary)
        self.assertTrue(summary["baseline_matches"])
        self.assertAlmostEqual(summary["simulation"]["dt_s"], .004)
        _, after = self.request("/api/export")
        for key in ("config", "state", "history"):
            self.assertEqual(before[key], after[key])

    def test_tampered_profile_is_rejected_before_replacing_model(self):
        config = default_config()
        config["calibration"] = fit_time_mapping(demo_payload(config))["profile"]
        self.assertTrue(validate_config(config)["valid"])
        config["calibration"]["status"] = "fitted"
        self.assertFalse(validate_config(config)["valid"])
        status, _ = self.request("/api/reset", {"config": config})
        self.assertEqual(status, 400)
        self.assertIsNone(server.WORKSPACE.sim)

    def test_saved_reference_results_detect_changed_integrator(self):
        status, report = self.request("/api/calibration/benchmarks")
        self.assertEqual(status, 200)
        self.assertTrue(report["current_source_matches"])
        with patch("inspect.getsource", return_value="changed integrator"):
            status, stale = self.request("/api/calibration/benchmarks")
        self.assertEqual(status, 200)
        self.assertFalse(stale["current_source_matches"])

    def test_mapping_does_not_change_seeded_trajectory_and_survives_bundle(self):
        original = default_config()
        mapped = copy.deepcopy(original)
        mapped["calibration"] = fit_time_mapping(demo_payload(original))["profile"]
        a, b = Simulation(original), Simulation(mapped)
        a.step(10)
        b.step(10)
        sa, sb = a.snapshot(), b.snapshot()
        sa.pop("warnings"); sb.pop("warnings")
        sa["genomic_events"].pop("run_id"); sb["genomic_events"].pop("run_id")
        self.assertEqual(sa, sb)
        with zipfile.ZipFile(io.BytesIO(build_bundle(mapped, {"platform": "numpy"}))) as archive:
            self.assertEqual(json.loads(archive.read("model.json"))["calibration"], mapped["calibration"])
            units = json.loads(archive.read("physical-units.json"))
            self.assertEqual(units["status"], "demo")
            self.assertTrue(units["baseline_matches"])
            self.assertIn("chromatin_studio/calibration.py", archive.namelist())

    def test_batch_to_msd_cli_preserves_units_and_replica_baseline(self):
        config = default_config()
        config.update(proteins=[], interactions=[], genes=[], bodies=[])
        config["cohesin"]["enabled"] = False
        config["calibration"] = fit_time_mapping(demo_payload(config))["profile"]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/"replica"
            run = run_batch(config, {"platform": "numpy", "replicas": 2, "steps": 20,
                            "sample_every": 2, "snapshot_every": 2}, replica=1, output=output)
            self.assertTrue(run["physical_mapping"]["baseline_matches"])
            self.assertEqual(run["physical_mapping"], json.loads((output/"physical-units.json").read_text()))
            curve = Path(directory)/"curve.json"
            proc = subprocess.run([sys.executable, "-m", "chromatin_studio.calibration", "msd",
                    "--snapshots", str(output/"snapshots.jsonl"), "--beads", "10,20", "--dimension", "2",
                    "--observable", "relative_pair", "--output", str(curve)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            data = json.loads(curve.read_text())
            self.assertEqual(data["observable"], "relative_pair")
            self.assertEqual(data["dimension"], 2)
            self.assertEqual(data["lag_unit"], "reduced_time")
            self.assertEqual(data["baseline_fingerprint"], config["calibration"]["baseline_fingerprint"])
            self.assertEqual(data["analysis"]["retained_frames"], 11)


if __name__ == "__main__":
    unittest.main()
