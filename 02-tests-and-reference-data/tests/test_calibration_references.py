"""Public experimental references and their read-only HTTP contract."""
import json
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import numpy as np

from chromatin_studio import server
from chromatin_studio.calibration import _curve


class CalibrationReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.http = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.thread = threading.Thread(target=cls.http.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.http.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.http.shutdown()
        cls.http.server_close()
        cls.thread.join()

    def get(self, path):
        try:
            response = urllib.request.urlopen(self.base + path, timeout=10)
        except urllib.error.HTTPError as exc:
            response = exc
        with response:
            return response.status, json.load(response)

    def test_catalog_has_all_conditions_without_mutating_workspace(self):
        workspace = server.WORKSPACE
        before = (workspace.sim, workspace.config, list(workspace.history))
        status, catalog = self.get("/api/calibration/references")
        self.assertEqual(status, 200)
        defaults = [r for r in catalog["references"] if r["is_default_for_condition"]]
        self.assertEqual(len(defaults), 12)
        self.assertEqual(len({r["condition"] for r in defaults}), 12)
        self.assertEqual(catalog["default_reference_id"], "fbn2-rad21-2h-fit-corrected")
        for entry in catalog["references"]:
            status, response = self.get("/api/calibration/reference?id=" + entry["id"])
            self.assertEqual(status, 200)
            self.assertEqual(response["reference"], entry)
            curve = _curve(response["experiment"], "experiment")
            self.assertEqual(curve["dimension"], 3)
            self.assertEqual(curve["observable"], "relative_pair")
            self.assertEqual(len(curve["lag"]), 180)
            self.assertNotIn("time_s", response)
            self.assertNotIn("baseline_fingerprint", curve)
        self.assertEqual(before, (workspace.sim, workspace.config, list(workspace.history)))

    def test_baseline_acquisitions_are_disjoint_and_noise_is_explicit(self):
        curves = {}
        for split in ("fit", "holdout", "all"):
            for noise in ("raw", "corrected"):
                _, response = self.get(f"/api/calibration/reference?id=fbn2-rad21-2h-{split}-{noise}")
                curves[split, noise] = response["experiment"]
        fit = curves["fit", "raw"]["analysis"]
        holdout = curves["holdout", "raw"]["analysis"]
        all_tracks = curves["all", "raw"]["analysis"]
        self.assertFalse(set(fit["cluster_ids"]) & set(holdout["cluster_ids"]))
        self.assertEqual(set(fit["cluster_ids"]) | set(holdout["cluster_ids"]), set(all_tracks["cluster_ids"]))
        self.assertEqual(fit["tracks"] + holdout["tracks"], 537)
        for split in ("fit", "holdout", "all"):
            raw, corrected = curves[split, "raw"], curves[split, "corrected"]
            np.testing.assert_allclose(np.asarray(raw["msd"]) - corrected["msd"], 2 * .203**2)
            self.assertEqual(raw["counts"], corrected["counts"])
        fit_curve, hold_curve, full_curve = [curves[s, "raw"] for s in ("fit", "holdout", "all")]
        fit_counts, hold_counts = np.asarray(fit_curve["counts"]), np.asarray(hold_curve["counts"])
        np.testing.assert_array_equal(fit_counts + hold_counts, full_curve["counts"])
        np.testing.assert_allclose((fit_counts * fit_curve["msd"] + hold_counts * hold_curve["msd"]) /
                                   (fit_counts + hold_counts), full_curve["msd"])

    def test_reference_ids_do_not_select_arbitrary_files(self):
        for path, expected in [("", 400), ("?id=", 400), ("?id=x&id=y", 400),
                               ("?id=unknown", 404), ("?id=../../README.md", 404)]:
            self.assertEqual(self.get("/api/calibration/reference" + path)[0], expected)


if __name__ == "__main__":
    unittest.main()
