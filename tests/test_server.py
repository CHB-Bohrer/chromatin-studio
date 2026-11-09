"""HTTP contract tests against a temporary loopback server."""
import json
import io
import zipfile
import threading
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

from chromatin_studio import server
from chromatin_studio.schema import default_config, preset_config, MAX_CHROMATIN_BEADS, MAX_PARTICLES


class ServerTests(unittest.TestCase):
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

    def request(self, path, data=None, headers=None):
        hdr = {"Content-Type": "application/json"}
        hdr.update(headers or {})
        req = urllib.request.Request(self.base + path, data=None if data is None else json.dumps(data).encode(), headers=hdr)
        try:
            response = urllib.request.urlopen(req, timeout=10)
        except urllib.error.HTTPError as exc:
            response = exc
        with response:
            return response.status, json.load(response)

    def test_step_and_export_preserve_applied_config(self):
        status, state = self.request("/api/reset", {"config": default_config(), "backend": "numpy"})
        self.assertEqual(status, 200)
        self.assertEqual(state["step"], 0)
        status, advanced = self.request("/api/step", {"steps": 5})
        self.assertEqual(status, 200)
        self.assertEqual(advanced["step"], 5)
        status, export = self.request("/api/export")
        self.assertEqual(export["config"], default_config())
        self.assertEqual(export["state"]["step"], 5)
        self.assertEqual(len(export["history"]), 2)

    def test_study_plan_and_reference_leave_live_run_unchanged(self):
        self.request("/api/reset", {"config": default_config()})
        self.request("/api/step", {"steps": 3})
        _, before = self.request("/api/export")
        payload = {"config": preset_config("manuscript_locus"), "options": {
            "engine": "reference", "conditions": ["uncoupled"], "cells": 64,
            "duration": 2, "bootstrap_samples": 0}}
        for path in ("/api/studies/plan", "/api/studies/run"):
            status, result = self.request(path, payload)
            self.assertEqual(status, 200, result)
        _, after = self.request("/api/export")
        for key in ("config", "state", "history"):
            self.assertEqual(before[key], after[key], key)
        self.assertFalse(server.STUDY_LOCK.locked())
        self.assertEqual(self.request("/api/studies/catalog")[0], 200)

    def test_study_lock_released_on_error_and_busy_run_rejected(self):
        payload = {"config": default_config(), "options": {}}
        with patch("chromatin_studio.studies.run_study", side_effect=ValueError("test failure")):
            status, result = self.request("/api/studies/run", payload)
            self.assertEqual(status, 400)
            self.assertIn("test failure", result["error"])
        self.assertFalse(server.STUDY_LOCK.locked())
        server.STUDY_LOCK.acquire()
        try:
            self.assertEqual(self.request("/api/studies/run", payload)[0], 409)
        finally:
            server.STUDY_LOCK.release()
        self.assertIsNone(server.WORKSPACE.sim)

    def test_invalid_reset_keeps_previous_working_model(self):
        self.request("/api/reset", {"config": default_config()})
        self.request("/api/step", {"steps": 3})
        c = default_config()
        c["genes"][0]["required_mark"] = "missing"
        status, result = self.request("/api/reset", {"config": c})
        self.assertEqual(status, 400)
        self.assertIn("details", result)
        self.assertEqual(self.request("/api/state")[1]["step"], 3)

    def test_cluster_plan_and_zip_export_leave_live_run_unchanged(self):
        self.request("/api/reset", {"config": default_config()})
        self.request("/api/step", {"steps": 3})
        _, before = self.request("/api/export")
        draft = preset_config("nucleosome_transcription")
        payload = {"config": draft, "job": {"platform": "CUDA", "replicas": 2}}
        status, plan = self.request("/api/cluster/plan", payload)
        self.assertEqual(status, 200, plan)
        self.assertEqual(plan["model"]["genomic_span_bp"], 80000)
        self.assertEqual(plan["replicas"], [{"index": 0, "seed": 42}, {"index": 1, "seed": 43}])
        request = urllib.request.Request(self.base + "/api/cluster/export", data=json.dumps(payload).encode(),
                                         headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(request) as response:
            self.assertEqual(response.headers["Content-Type"], "application/zip")
            raw = response.read()
        with zipfile.ZipFile(io.BytesIO(raw)) as archive:
            self.assertEqual(json.loads(archive.read("model.json")), draft)
            self.assertIn("chromatin_studio/validation.py", archive.namelist())
            self.assertIn(b"--gres=gpu:1", archive.read("submit.slurm"))
        _, after = self.request("/api/export")
        for key in ("config", "state", "history"):
            self.assertEqual(before[key], after[key])
        payload["job"]["partition"] = "bad\n#SBATCH --nodes=100"
        self.assertEqual(self.request("/api/cluster/plan", payload)[0], 400)

    def test_regulatory_preset_survives_api_step_and_export(self):
        config = preset_config("nucleosome_locus")
        _, caps = self.request("/api/capabilities")
        self.assertEqual(caps["limits"]["regulatory_regions"], 128)
        status, state = self.request("/api/reset", {"config": config, "backend": "numpy"})
        self.assertEqual(status, 200, state)
        self.assertEqual(len(state["regulation"]["regions"]), 4)
        self.assertEqual(state["regulation"]["genes"][0]["mode"], "regions")
        self.assertEqual(len(state["chromatin_contacts"]["rules"]), 1)
        status, state = self.request("/api/step", {"steps": 10})
        self.assertEqual(status, 200, state)
        _, export = self.request("/api/export")
        self.assertEqual(export["config"], config)
        self.assertEqual(export["state"]["regulation"], state["regulation"])
        self.assertEqual(export["state"]["chromatin_contacts"], state["chromatin_contacts"])
        self.assertEqual(export["state"]["genome_tracks"]["span_bp"], 80_000)

    def test_nucleosomes_and_individual_rna_survive_api_and_export(self):
        config = preset_config("nucleosome_transcription")
        config["simulation"]["dt"] = .01
        for gene in config["genes"]:
            gene.update(on_rate=20, off_rate=0, initiation_rate=20,
                        elongation_rate=100, end=gene["start"])
        _, capabilities = self.request("/api/capabilities")
        self.assertEqual(capabilities["limits"]["transcripts"], 5000)
        status, state = self.request("/api/reset", {"config": config, "backend": "numpy"})
        self.assertEqual(status, 200, state)
        self.assertTrue(state["nucleosomes"]["enabled"])
        self.assertEqual(len(state["nucleosomes"]["occupied"]), 400)
        status, state = self.request("/api/step", {"steps": 100})
        self.assertEqual(status, 200, state)
        self.assertTrue(state["transcription"]["enabled"])
        self.assertGreater(state["transcription"]["totals"]["completed"], 0)
        totals = state["transcription"]["totals"]
        self.assertEqual(totals["initiated"], totals["live"] + totals["degraded"])
        _, saved = self.request("/api/export")
        self.assertEqual(saved["state"]["transcription"], state["transcription"])
        self.assertEqual(saved["state"]["nucleosomes"], state["nucleosomes"])
        self.assertEqual(saved["config"], config)
        self.assertTrue(saved["state"]["genome_tracks"]["persistence_runs"])

    def test_cross_origin_mutation_is_rejected(self):
        status, _ = self.request("/api/step", {"steps": 3}, {"Origin": "https://unrelated.example"})
        self.assertEqual(status, 403)

    def test_unknown_engine_and_invalid_step_rejected(self):
        self.assertEqual(self.request("/api/reset", {"config": default_config(), "backend": "made-up"})[0], 400)
        for value in (0, 101, True, "10"):
            self.assertEqual(self.request("/api/step", {"steps": value})[0], 400)

    def test_validation_returns_resource_errors_without_initialization(self):
        c = default_config()
        c["proteins"][0]["count"] = 10000
        status, result = self.request("/api/validate", c)
        self.assertEqual(status, 200)
        self.assertFalse(result["valid"])
        self.assertIsNone(server.WORKSPACE.sim)

    def test_physical_binding_estimate_uses_applied_config_without_mutation(self):
        config = default_config()
        self.request("/api/reset", {"config": config})
        self.request("/api/step", {"steps": 3})
        _, before = self.request("/api/export")
        payload = {"a": "srrm2", "b": "activator", "length_nm": 100,
                   "time_s": 0.01, "target_kd_um": 1, "efficiency": 1}
        status, estimate = self.request("/api/kinetics/estimate", payload)
        self.assertEqual(status, 200, estimate)
        self.assertGreater(estimate["kon_limit_M_inv_s"], 0)
        self.assertEqual(estimate["kon_M_inv_s"], estimate["kon_limit_M_inv_s"])
        self.assertAlmostEqual(estimate["nominal_koff_s_inv"],
                               estimate["kon_M_inv_s"] * 1e-6)
        _, after = self.request("/api/export")
        self.assertEqual(before["state"], after["state"])
        self.assertEqual(before["history"], after["history"])
        self.assertEqual(before["config"], after["config"])

    def test_invalid_physical_estimate_does_not_initialize_or_reset(self):
        payload = {"a": "srrm2", "b": "activator", "length_nm": 100,
                   "time_s": 0.01, "target_kd_um": 1, "efficiency": 1}
        for field, value in (("length_nm", 0), ("time_s", True),
                             ("target_kd_um", -1), ("efficiency", 2),
                             ("a", "missing")):
            with self.subTest(field=field):
                status, report = self.request("/api/kinetics/estimate", {**payload, field: value})
                self.assertEqual(status, 400, report)
                self.assertIn("error", report)
                self.assertIsNone(server.WORKSPACE.sim)

    def test_twenty_mb_model_runs_and_exports_with_reported_limits(self):
        config = preset_config("twenty_mb")
        status, report = self.request("/api/validate", config)
        self.assertEqual(status, 200)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["genomic_span_bp"], 20_000_000)
        _, capabilities = self.request("/api/capabilities")
        self.assertEqual(capabilities["limits"]["particles"], MAX_PARTICLES)
        self.assertEqual(capabilities["limits"]["chromatin_beads"], MAX_CHROMATIN_BEADS)
        status, state = self.request("/api/reset", {"config": config, "backend": "numpy"})
        self.assertEqual(status, 200, state)
        status, advanced = self.request("/api/step", {"steps": 2})
        self.assertEqual(status, 200, advanced)
        self.assertEqual(advanced["step"], 2)
        self.assertEqual(advanced["metrics"]["particle_count"], report["particle_count"])
        tracks = advanced["genome_tracks"]
        self.assertEqual(tracks["span_bp"], 20_000_000)
        self.assertEqual(tracks["bp_per_bead"], 10_000)
        self.assertEqual(set(tracks["affinity_runs"]), {p["id"] for p in config["proteins"]})
        for family in ("affinity_runs", "occupancy_runs"):
            for runs in tracks[family].values():
                self.assertEqual(runs[0][0], 0)
                self.assertEqual(runs[-1][1], 2000)
        _, exported = self.request("/api/export")
        self.assertEqual(exported["config"], config)
        self.assertEqual(exported["state"]["genome_tracks"], tracks)
        interactions = advanced["interaction_state"]
        self.assertEqual(interactions["totals"]["molecules"], 194)
        self.assertEqual(len(interactions["pairs"]), len(config["interactions"]))
        self.assertEqual(exported["state"]["interaction_state"], interactions)
        self.assertEqual(len(exported["history"]), 2)


if __name__ == "__main__":
    unittest.main()
