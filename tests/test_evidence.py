"""Evidence view backend: read-only report index, scoped paths, hidden machine paths."""
import http.client
import json
import tempfile
import threading
import unittest
import urllib.parse
from http.server import ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch

from chromatin_studio import server
from chromatin_studio.components import evidence


def report(study="polymer_wlc", **fields):
    value = {"format": 2, "study": study, "title": "Test study", "scale": "smoke", "analyzed_utc": "2026-11-21T00:00:00Z",
             "check_summary": {"total": 1, "passed": 1, "by_evidence": {}}, "figures": [{"id": "f"}],
             "checks": [{"id": "a", "label": "a", "evidence": "numerical_check", "observed": 1.0,
                         "expected": 1.0, "passed": True}],
             "summary": {"warm_start": r"C:\Users\someone\Dropbox\runs\final_positions.npy",
                         "posix": "/home/someone/runs/x.npy", "url": "https://doi.org/10.1/x"},
             "missing_jobs": ["a-r1"]}
    value.update(fields)
    return value


class EvidenceIndexTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        study = self.root / "2026-11-21" / "polymer_wlc"
        (study / "report-history").mkdir(parents=True)
        (study / "report.json").write_text(json.dumps(report()), encoding="utf-8")
        (study / "report-history" / "report-20261120T000000Z.json").write_text(json.dumps(report(format=1)),
                                                                              encoding="utf-8")
        (study / "plan.json").write_text("{}", encoding="utf-8")
        broken = self.root / "2026-11-21" / "broken"
        broken.mkdir()
        (broken / "report.json").write_text("{not json", encoding="utf-8")
        (self.root / "secret.json").write_text('{"secret": true}', encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def test_unconfigured_folder_lists_the_catalog_and_explains_how_to_connect(self):
        result = evidence.index(None)
        self.assertFalse(result["configured"])
        self.assertTrue(result["catalog"])
        self.assertIn("--evidence", result["hint"])

    def test_index_uses_relative_paths_and_reports_bad_files_without_failing(self):
        result = evidence.index(self.root)
        paths = {entry["path"]: entry for entry in result["reports"]}
        self.assertEqual(set(paths), {"2026-11-21/polymer_wlc/report.json", "2026-11-21/broken/report.json"})
        good = paths["2026-11-21/polymer_wlc/report.json"]
        self.assertEqual(good["history"], ["2026-11-21/polymer_wlc/report-history/report-20261120T000000Z.json"])
        self.assertEqual((good["figures"], good["missing_jobs"], good["registered"]), (1, 1, True))
        self.assertIn("error", paths["2026-11-21/broken/report.json"])
        self.assertNotIn(str(self.root), json.dumps(result))

    def test_report_hides_absolute_paths_but_keeps_links(self):
        data = evidence.report(self.root, "2026-11-21/polymer_wlc/report.json")
        self.assertEqual(data["summary"]["warm_start"], "…/final_positions.npy")
        self.assertEqual(data["summary"]["posix"], "…/x.npy")
        self.assertEqual(data["summary"]["url"], "https://doi.org/10.1/x")
        archived = evidence.report(self.root, "2026-11-21/polymer_wlc/report-history/report-20261120T000000Z.json")
        self.assertEqual(archived["format"], 1)

    def test_only_reports_inside_the_folder_can_be_opened(self):
        for bad in ("../secret.json", "secret.json", "2026-11-21/polymer_wlc/plan.json", "/etc/passwd",
                    "C:/Windows/win.ini", "2026-11-21\\polymer_wlc\\report.json", "", "a/../../report.json"):
            with self.subTest(path=bad), self.assertRaises(ValueError):
                evidence.report(self.root, bad)
        with self.assertRaises(FileNotFoundError):
            evidence.report(self.root, "2026-11-21/missing/report.json")
        with self.assertRaises(ValueError):
            evidence.report(None, "2026-11-21/polymer_wlc/report.json")


class EvidenceServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name).resolve()
        (cls.root / "run" / "tf_binding").mkdir(parents=True)
        (cls.root / "run" / "tf_binding" / "report.json").write_text(json.dumps(report("tf_binding")),
                                                                    encoding="utf-8")
        cls.patch = patch.dict(server.EVIDENCE, {"root": cls.root})
        cls.patch.start()
        cls.http = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.thread = threading.Thread(target=cls.http.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.http.shutdown()
        cls.http.server_close()
        cls.thread.join()
        cls.patch.stop()
        cls.temporary.cleanup()

    def get(self, path):
        # Plain HTTP connection: urllib builds an SSL context, which can fail
        # intermittently on a heavily loaded Windows machine.
        connection = http.client.HTTPConnection("127.0.0.1", self.http.server_port, timeout=10)
        try:
            connection.request("GET", path)
            response = connection.getresponse()
            return response.status, json.loads(response.read())
        finally:
            connection.close()

    def test_index_and_report_routes(self):
        status, index = self.get("/api/evidence")
        self.assertEqual(status, 200)
        self.assertEqual([r["path"] for r in index["reports"]], ["run/tf_binding/report.json"])
        status, data = self.get("/api/evidence/report?path=" + urllib.parse.quote("run/tf_binding/report.json"))
        self.assertEqual((status, data["study"]), (200, "tf_binding"))

    def test_bad_requests_are_refused_with_client_errors(self):
        self.assertEqual(self.get("/api/evidence/report")[0], 400)
        self.assertEqual(self.get("/api/evidence/report?path=..%2Fx%2Freport.json")[0], 400)
        self.assertEqual(self.get("/api/evidence/report?path=run%2Fnone%2Freport.json")[0], 404)


if __name__ == "__main__":
    unittest.main()
