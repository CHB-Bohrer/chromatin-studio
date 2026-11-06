import copy
import json
import math
import unittest

import numpy as np

from chromatin_studio.schema import default_config, preset_config
from chromatin_studio.studies import (DEFAULTS, analyze_gene_pair, build_study, correlation,
                                     detection_phi, four_promoter_statistics,
                                     run_study, simulate_reference, study_catalog)


class ReferenceStudyTests(unittest.TestCase):
    def options(self, **changes):
        result = copy.deepcopy(DEFAULTS["reference"])
        result.update(cells=4096, duration=30, environment_gain=0, competition_strength=0,
                      extrinsic_cv=0, bootstrap_samples=0)
        result.update(changes)
        return result

    def test_matched_homolog_contrast_and_independent_cell_bootstrap(self):
        # A=B within each homolog, with independent homologs, is exactly cis=1,
        # trans=0 in the balanced four-cell population. Resampling must keep
        # all four allele observations from a cell together.
        values = np.array([[0, 0, 0, 0], [0, 0, 1, 1], [1, 1, 0, 0], [1, 1, 1, 1]] * 40)
        result = four_promoter_statistics(values, bootstrap_samples=100, seed=17)
        self.assertAlmostEqual(result["phi_cis"], 1)
        self.assertAlmostEqual(result["phi_trans"], 0)
        self.assertAlmostEqual(result["delta_phi"], 1)
        self.assertEqual(result["bootstrap_fraction"], 1)
        self.assertEqual(result["bootstrap_unit"], "independent cell")
        self.assertGreater(result["ci95"][0], .7)
        # Homolog-exchange symmetry must not change a paired contrast.
        swapped = four_promoter_statistics(values[:, [2, 3, 0, 1]], bootstrap_samples=0)
        self.assertEqual(result["delta_phi"], swapped["delta_phi"])

    def test_undefined_correlation_is_not_reported_as_zero(self):
        self.assertIsNone(correlation([1, 1, 1], [0, 1, 0]))
        result = four_promoter_statistics(np.zeros((64, 4)), bootstrap_samples=20)
        self.assertIsNone(result["delta_phi"])
        self.assertIsNone(result["bootstrap_fraction"])
        self.assertEqual(result["valid_bootstraps"], 0)
        json.dumps(result, allow_nan=False)

    def test_independent_detection_matches_exact_joint_probability(self):
        # Joint P11=.30, P10=.20, P01=.10, P00=.40.
        pa, pb, p11 = .5, .4, .3
        latent = (p11-pa*pb) / math.sqrt(pa*(1-pa)*pb*(1-pb))
        qa, qb = .3, .6
        observed = (qa*qb*p11-(qa*pa)*(qb*pb)) / math.sqrt(qa*pa*(1-qa*pa)*qb*pb*(1-qb*pb))
        self.assertAlmostEqual(detection_phi(latent, pa, pb, qa, qb), observed)
        self.assertAlmostEqual(detection_phi(latent, pa, pb, 1), latent)
        self.assertIsNone(detection_phi(latent, pa, pb, 0))
        self.assertNotEqual(detection_phi(latent, .1, .1, .3)/latent,
                            detection_phi(latent, .8, .8, .3)/latent)

    def test_exact_uncoupled_telegraph_moments_and_zero_correlation(self):
        options = self.options(on_rate=.5, off_rate=1, initiation_rate=2, clearance_rate=1)
        result = simulate_reference(options, seed=802)
        counts = result["nascent_counts"].ravel()
        on = result["promoter_on"].ravel()
        target_on = .5 / 1.5
        target_mean = 2 * target_on
        target_variance = target_mean + 2**2 * .5 * 1 / (1 * 1.5**2 * 2.5)
        self.assertLess(abs(float(on.mean())-target_on), 6*math.sqrt(target_on*(1-target_on)/len(on)))
        self.assertLess(abs(float(counts.mean())-target_mean), 6*math.sqrt(target_variance/len(counts)))
        self.assertLess(abs(float(counts.var())-target_variance), .12)
        stats = four_promoter_statistics(result["nascent_positive"], bootstrap_samples=0)
        self.assertLess(abs(stats["delta_phi"]), .06)
        self.assertLess(abs(stats["phi_cis"]), .06)
        self.assertLess(abs(stats["phi_trans"]), .06)

    def test_mutual_inhibition_matches_stationary_promoter_distribution(self):
        options = self.options(competition_strength=20, initiation_rate=6, clearance_rate=5)
        result = simulate_reference(options, seed=902)
        ratio = options["on_rate"] / options["off_rate"]
        double = ratio**2 / (1+options["competition_strength"])
        z = 1+2*ratio+double
        p, p11 = (ratio+double)/z, double/z
        target_phi = (p11-p*p)/(p*(1-p))
        promoter = four_promoter_statistics(result["promoter_on"], bootstrap_samples=0)
        self.assertLess(abs(promoter["phi_cis"]-target_phi), .04)
        self.assertLess(abs(promoter["phi_trans"]), .06)
        nascent = four_promoter_statistics(result["nascent_positive"], bootstrap_samples=0)
        self.assertLess(nascent["delta_phi"], -.08)

    def test_shared_environment_can_coordinate_without_competition(self):
        options = self.options(cells=2048, duration=35, environment_gain=12,
                               environment_on_rate=.04, environment_off_rate=.04,
                               on_rate=.15, initiation_rate=4, clearance_rate=2)
        result = simulate_reference(options, seed=417)
        stats = four_promoter_statistics(result["nascent_positive"], bootstrap_samples=0)
        self.assertGreater(stats["delta_phi"], .1)
        self.assertLess(abs(stats["phi_trans"]), .06)

    def test_observation_loss_is_separate_from_promoter_and_latent_states(self):
        result = run_study(None, {"conditions": ["uncoupled"], "cells": 256,
                                 "duration": 15, "detection_efficiency": 0, "bootstrap_samples": 20})
        metrics = result["conditions"][0]["metrics"]
        self.assertIsNotNone(metrics["latent"]["phi_cis"])
        self.assertIsNotNone(metrics["promoter"]["phi_cis"])
        self.assertIsNone(metrics["observed"]["phi_cis"])
        self.assertEqual(metrics["observed_nascent_probability"], [0, 0, 0, 0])


class PolymerStudyTests(unittest.TestCase):
    def test_ablation_plan_is_valid_and_does_not_mutate_its_source(self):
        config = preset_config("nucleosome_transcription")
        config["proteins"][0]["nucleosome_off_rate_multiplier"] = {"occupied": 2, "vacant": .5}
        config["proteins"][0]["nucleosome_sliding_rate"] = .8
        original = copy.deepcopy(config)
        conditions = [c["id"] for c in study_catalog()["conditions"] if "polymer" in c["engines"]]
        plan = build_study(config, {"engine": "polymer", "conditions": conditions})
        self.assertEqual(config, original)
        rows = {c["id"]: c["config"] for c in plan["conditions"]}
        self.assertEqual(rows["baseline"], config)
        self.assertEqual(rows["no_cohesin"]["cohesin"]["enabled"], False)
        self.assertTrue(all(r["activation_gain"] == r["initiation_gain"] == 0 for r in rows["no_regulatory_gain"]["regulatory_regions"]))
        self.assertTrue(all(p["on_rate"] == 0 for p in rows["no_protein_binding"]["interactions"]))
        self.assertEqual(rows["no_protein_binding"]["proteins"], config["proteins"])
        self.assertTrue(all(b["strength"] == 0 for b in rows["no_body_affinity"]["bodies"]))
        self.assertEqual(rows["no_body_affinity"]["genes"], config["genes"])
        self.assertEqual(rows["no_speckle_processing"]["bodies"], config["bodies"])
        self.assertEqual(rows["no_confinement"]["simulation"]["confinement_stiffness"], 0)
        self.assertEqual(rows["no_nucleosome_access"]["chromatin"], config["chromatin"])
        self.assertTrue(all(p["nucleosome_affinity"] == {"occupied": 1, "vacant": 1}
                            for p in rows["no_nucleosome_access"]["proteins"]))
        self.assertTrue(all(p["nucleosome_off_rate_multiplier"] == {"occupied": 1, "vacant": 1}
                            for p in rows["no_nucleosome_access"]["proteins"]))
        frozen = rows["no_nucleosome_turnover"]
        self.assertEqual(frozen["chromatin"]["nucleosomes"]["initial_regions"], config["chromatin"]["nucleosomes"]["initial_regions"])
        self.assertEqual(frozen["chromatin"]["nucleosomes"]["eviction_rate"], 0)
        self.assertTrue(all(p["nucleosome_eviction_rate"] == p["nucleosome_reassembly_rate"] == 0 for p in frozen["proteins"]))
        self.assertTrue(all(p["nucleosome_sliding_rate"] == 0 for p in frozen["proteins"]))
        self.assertTrue(all(g["nucleosome_eviction_rate"] == 0 for g in frozen["genes"]))
        rows["baseline"]["genes"][0]["on_rate"] = 0
        self.assertEqual(config, original)

    def test_contact_threshold_and_replica_uncertainty(self):
        samples = []
        for replica in range(3):
            for frame in range(4):
                active = frame % 2
                samples.append({"replica": replica, "time": frame,
                                "genes": {"a": {"nascent": active, "on": bool(1-active), "position": [0, 0, 0], "start_bp": 0},
                                          "b": {"nascent": 2*active, "on": bool(active), "position": [3, 4, 0], "start_bp": 10000}}})
        result = analyze_gene_pair(samples, "a", "b", contact_radius=5, bootstrap_samples=20)
        self.assertAlmostEqual(result["nascent_phi"], 1)
        self.assertAlmostEqual(result["promoter_phi"], -1)
        self.assertAlmostEqual(result["count_correlation"], 1)
        self.assertEqual(result["contact_fraction"], 1)
        self.assertEqual(result["mean_distance"], 5)
        self.assertEqual(result["genomic_separation_bp"], 10000)
        self.assertEqual(result["n_replicas"], 3)
        self.assertEqual(result["bootstrap_ci95"], [1, 1])
        self.assertNotIn("phi_trans", result)
        self.assertNotIn("delta_phi", result)
        self.assertEqual(analyze_gene_pair(samples, "a", "b", contact_radius=4.99, bootstrap_samples=0)["contact_fraction"], 0)

    def test_small_polymer_pilot_is_reproducible_and_keeps_input_unchanged(self):
        config = default_config()
        config["chromatin"].update(beads=16, initial_marks=[], barriers=[])
        config["proteins"], config["interactions"], config["bodies"] = [], [], []
        config["cohesin"]["enabled"] = False
        template = config["genes"][0]
        template.update(start=2, end=4, activator_species=None, required_mark=None, repressive_mark=None,
                        on_rate=10, initiation_rate=5)
        config["genes"].append(dict(copy.deepcopy(template), id="gene_b", name="Gene B", start=10, end=12))
        before = copy.deepcopy(config)
        options = {"engine": "polymer", "conditions": ["baseline"], "replicas": 2,
                   "warmup_steps": 2, "samples_per_replica": 5, "sample_every": 3, "bootstrap_samples": 0}
        first, second = run_study(config, options), run_study(config, options)
        self.assertEqual(first["conditions"], second["conditions"])
        self.assertEqual(config, before)
        self.assertEqual(len(first["conditions"][0]["samples"]), 10)
        json.dumps(first, allow_nan=False)

    def test_rejects_invalid_options_and_keeps_large_plans_exportable(self):
        for options in ({"cells": True}, {"engine": "other"}, {"conditions": ["no_cohesin"]},
                        {"conditions": ["baseline", "baseline"]}, {"detection_efficiency": float("nan")},
                        {"unknown": 1}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                build_study(None, options)
        with self.assertRaises(ValueError):
            build_study(default_config(), {"engine": "polymer"})
        plan = build_study(None, {"cells": 4096, "duration": 200})
        self.assertEqual(plan["options"]["cells"], 4096)
        with self.assertRaisesRegex(ValueError, "budget"):
            run_study(None, plan["options"])


if __name__ == "__main__":
    unittest.main()
