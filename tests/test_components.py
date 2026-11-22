"""Component studies: module contract, runner behavior and honest check rules."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from chromatin_studio import components
from chromatin_studio.components import core, runner
from chromatin_studio.engine import Simulation

# Shrunk parameters keep the whole file to a few seconds; values are not evidence.
TINY = {
    "polymer_wlc": {"steps": 600, "burn_in_steps": 200},
    "protein_diffusion": {"steps": 400},
    "nucleosome_turnover": {"beads": 64},
    "telegraph_transcription": {"samples": 4, "burn_in": 50.0},
    "cohesin_extrusion": {"residences": 3, "beads": 200},
    "tf_binding": {"steps": 1500, "burn_in_steps": 300},
    "nucleosome_tf_kinetics": {"time": 10.0, "burn_in": 2.0},
    "hp1_mark_density": {"time": 10.0},
    "fbn2_cohesin": {"steps": 750, "burn_in_steps": 250, "lattice_prerun": 200.0},
}


class RegistryTests(unittest.TestCase):
    def test_every_registered_study_meets_the_module_contract(self):
        self.assertEqual(set(TINY), set(components.REGISTRY))
        for study_id in components.REGISTRY:
            study = components.load(study_id).STUDY
            self.assertTrue(set(study["evidence"]) <= set(core.EVIDENCE_LEVELS), study_id)
            self.assertTrue(study["sources"], study_id)
            for record in study["sources"]:
                self.assertTrue(record["citation"] and record["measured"], study_id)
            for scale in study["scales"].values():
                self.assertIn("replicas", scale)
                self.assertIn("seed", scale)
        self.assertEqual(len(components.catalog()), len(components.REGISTRY))

    def test_unknown_study_and_parameter_are_rejected(self):
        with self.assertRaises(ValueError):
            components.load("not_a_study")
        with self.assertRaises(ValueError):
            runner.parameters("polymer_wlc", "smoke", {"no_such_parameter": 1})
        with self.assertRaises(ValueError):
            runner.parameters("polymer_wlc", "no_such_scale")


class RunnerTests(unittest.TestCase):
    def test_every_study_runs_end_to_end_at_tiny_scale(self):
        for study_id, overrides in TINY.items():
            with self.subTest(study=study_id), tempfile.TemporaryDirectory() as folder:
                report = runner.run(study_id, folder, scale="smoke", replicas=1, overrides=overrides)
                self.assertEqual(report["failed_jobs"], [], report.get("job_outcomes"))
                self.assertGreater(report["results_used"], 0)
                self.assertTrue(report["source_matches_current"])
                for check in report["checks"]:
                    self.assertIn(check["evidence"], core.EVIDENCE_LEVELS)
                    self.assertIn("passed", check)
                for record in report.get("figures", []):
                    self.assertLessEqual(len(core.colored_groups(record)), core.MAX_COLORED_GROUPS,
                                         (study_id, record["id"]))
                json.dumps(report, allow_nan=False)

    def test_resume_keeps_finished_results_and_refuses_a_changed_plan(self):
        with tempfile.TemporaryDirectory() as folder:
            first = runner.run("nucleosome_turnover", folder, scale="smoke", replicas=1, overrides={"beads": 64})
            self.assertEqual([o["status"] for o in first["job_outcomes"]], ["complete"])
            raw = Path(folder) / "raw" / "default-r0.json"
            before = raw.read_bytes()
            second = runner.run("nucleosome_turnover", folder, scale="smoke", replicas=1, overrides={"beads": 64})
            self.assertEqual([o["status"] for o in second["job_outcomes"]], ["kept"])
            self.assertEqual(raw.read_bytes(), before)
            with self.assertRaises(ValueError):
                runner.run("nucleosome_turnover", folder, scale="smoke", replicas=1, overrides={"beads": 80})

    def test_partial_results_are_analyzed_and_labelled(self):
        with tempfile.TemporaryDirectory() as folder:
            runner.run("fbn2_cohesin", folder, scale="smoke", replicas=1, overrides=TINY["fbn2_cohesin"])
            raw = Path(folder) / "raw"
            finished = raw / "cohesin-r0.json"
            record = json.loads(finished.read_text(encoding="utf-8"))
            record["complete"] = False
            (raw / "cohesin-r0.partial.json").write_text(json.dumps(record), encoding="utf-8")
            finished.unlink()
            report = runner.analyze(folder)
            self.assertEqual(report["partial_results"], 1)
            self.assertTrue(any("partial" in line for line in report["findings"]))


class EvidenceTests(unittest.TestCase):
    def test_check_rules(self):
        passed = core.check("a", "a", 1.02, 1.0, evidence=core.ANALYTICAL, tolerance=0.05)
        failed = core.check("b", "b", 1.2, 1.0, evidence=core.ANALYTICAL, tolerance=0.05)
        noisy = core.check("c", "c", 1.2, 1.0, evidence=core.ANALYTICAL, tolerance=0.05, sem=0.1)
        outside = core.check("d", "d", 0.2, 0.05, evidence=core.EXPERIMENTAL, low=0.03, high=0.065)
        missing = core.check("e", "e", None, 1.0, evidence=core.NUMERICAL, tolerance=0.05)
        self.assertEqual([c["passed"] for c in (passed, failed, noisy, outside, missing)],
                         [True, False, True, False, False])
        self.assertEqual(core.summarize([passed, failed, outside])["by_evidence"][core.EXPERIMENTAL],
                         {"total": 1, "passed": 0})
        with self.assertRaises(ValueError):
            core.check("f", "f", 1, 1, evidence="validated")

    def test_wrong_physics_does_not_pass(self):
        module = components.load("nucleosome_turnover")
        params = runner.parameters("nucleosome_turnover", "smoke")
        wrong = {"complete": True, "case": {"id": "default"},
                 "observables": {"occupied_fraction": {str(t): 0.5 for t in params["times"]}}}
        analysis = module.analyze(params, [wrong, wrong])
        self.assertTrue(analysis["checks"])
        self.assertFalse(any(check["passed"] for check in analysis["checks"]))

    def test_chemistry_only_stepping_freezes_coordinates_and_advances_time(self):
        config = core.bare_config(16, seed=3)
        config["cohesin"] = {"enabled": True, "count": 2, "step_rate": 1.0, "off_rate": 0.1, "spring_k": 2.5}
        sim = Simulation(core.validated(config))
        before = sim.positions.copy()
        core.chemistry_only_step(sim, 10, 0.5)
        self.assertTrue(np.array_equal(sim.positions, before))
        self.assertAlmostEqual(sim.time, 5.0)
        self.assertEqual(sim.step_count, 10)

    def test_bare_config_has_no_other_mechanism(self):
        config = core.validated(core.bare_config(12))
        sim = Simulation(config)
        self.assertEqual((len(sim.molecules), len(sim.genes), len(sim.cohesin), len(sim.marks), len(sim.bodies)),
                         (0, 0, 0, 0, 0))


class StudyHelperTests(unittest.TestCase):
    def test_hp1_episodes_merge_short_gaps_and_count_search_time_at_risk(self):
        from chromatin_studio.components.hp1_mark_density import merge_episodes
        # Gap 0.5 merges; gap 3 starts a new episode. Search time: 1 before the first
        # arrival, 3 - 1 (merge window) before the second, 10 - 7 - 1 after the last.
        self.assertEqual(merge_episodes([(1, 2), (2.5, 3), (6, 7)], None, 1.0, 0.0, 10.0), ([2.0, 1.0], 2, 5.0))
        # A bond still open at the end that would merge leaves its episode unfinished.
        self.assertEqual(merge_episodes([(1, 2)], 2.4, 1.0, 0.0, 5.0), ([], 1, 1.0))
        # A bond still open after a long gap is a new arrival; the finished episode counts.
        self.assertEqual(merge_episodes([(1, 2)], 4.0, 1.0, 0.0, 5.0), ([1.0], 2, 2.0))
        # No binding at all: the whole run is search time with no arrival.
        self.assertEqual(merge_episodes([], None, 1.0, 0.0, 5.0), ([], 0, 5.0))


if __name__ == "__main__":
    unittest.main()
