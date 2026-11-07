import copy
import unittest

from chromatin_studio.schema import (default_config, preset_config, PRESETS,
                                    validate_config, MAX_CHROMATIN_BEADS, MAX_PARTICLES,
                                    MAX_GENES, MAX_COHESIN)


class SchemaTests(unittest.TestCase):
    def test_all_presets_are_valid_and_finite_in_size(self):
        for p in PRESETS:
            with self.subTest(preset=p["id"]):
                result = validate_config(preset_config(p["id"]))
                self.assertTrue(result["valid"], result["errors"])
                self.assertLessEqual(result["particle_count"], MAX_PARTICLES)

    def test_rejects_dangling_species_and_mark_references(self):
        c = default_config()
        c["proteins"] = c["proteins"][1:]
        c["marks"] = c["marks"][1:]
        report = validate_config(c)
        self.assertFalse(report["valid"])
        self.assertTrue(any("unknown identifier" in e for e in report["errors"]))

    def test_rejects_nonfinite_and_invalid_region(self):
        c = default_config()
        c["simulation"]["dt"] = float("nan")
        c["genes"][0]["start"] = 1000
        c["genes"][0]["end"] = -1
        report = validate_config(c)
        self.assertFalse(report["valid"])
        self.assertGreaterEqual(len(report["errors"]), 3)

    def test_particle_limit_counts_explicit_idr_beads(self):
        c = default_config()
        p = c["proteins"][0]
        p.update(count=400, representation="idr", idr_beads=16)
        report = validate_config(c)
        self.assertFalse(report["valid"])
        self.assertGreater(report["particle_count"], MAX_PARTICLES)
        self.assertTrue(any("at most 6000" in error for error in report["errors"]))

    def test_twenty_mb_span_and_coarse_resolution_warning(self):
        c = preset_config("twenty_mb")
        report = validate_config(c)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(c["chromatin"]["beads"], 2000)
        self.assertEqual(c["chromatin"]["bp_per_bead"], 10000)
        self.assertEqual(report["genomic_span_bp"], 20_000_000)
        self.assertEqual(report["particle_count"], 2314)
        self.assertEqual(c["cohesin"]["count"], 24)
        self.assertEqual(len(c["proteins"]), 6)
        self.assertEqual(len(c["genes"]), 6)
        self.assertTrue(any("Coarse genomic resolution" in w for w in report["warnings"]))
        self.assertTrue(any("does not calibrate or rescale" in w for w in report["warnings"]))

    def test_bp_per_bead_requires_positive_integer(self):
        for value in (1000.5, 0, -1, True, "10000", float("nan")):
            with self.subTest(value=value):
                c = default_config()
                c["chromatin"]["bp_per_bead"] = value
                report = validate_config(c)
                self.assertFalse(report["valid"])
                self.assertIsNone(report["genomic_span_bp"])
                self.assertTrue(any("chromatin.bp_per_bead" in e for e in report["errors"]))

    def test_chromatin_limit_and_inclusive_region_boundaries(self):
        c = preset_config("twenty_mb")
        c["chromatin"]["beads"] = MAX_CHROMATIN_BEADS
        c["genes"][0].update(start=MAX_CHROMATIN_BEADS-1, end=MAX_CHROMATIN_BEADS-1)
        self.assertTrue(validate_config(c)["valid"])
        c["genes"][0]["end"] = MAX_CHROMATIN_BEADS
        self.assertFalse(validate_config(c)["valid"])
        c["genes"][0]["end"] = MAX_CHROMATIN_BEADS-1
        c["chromatin"]["beads"] = MAX_CHROMATIN_BEADS+1
        self.assertFalse(validate_config(c)["valid"])

    def test_exact_particle_limit_is_allowed(self):
        c = default_config()
        c["chromatin"]["beads"] = MAX_CHROMATIN_BEADS
        for species in c["proteins"]:
            species["count"] = 500 if species["representation"] == "bead" else 0
        report = validate_config(c)
        self.assertEqual(report["particle_count"], MAX_PARTICLES)
        self.assertTrue(report["valid"], report["errors"])
        c["proteins"][2]["count"] = 1
        self.assertFalse(validate_config(c)["valid"])

    def test_dense_twenty_mb_annotation_and_motor_limits(self):
        c = preset_config("twenty_mb")
        c["cohesin"]["count"] = 100
        template = c["genes"][0]
        c["genes"] = [dict(copy.deepcopy(template), id=f"locus_{i}",
                           name=f"Locus {i}", start=7*i, end=7*i+3)
                      for i in range(MAX_GENES)]
        c["chromatin"]["initial_marks"] = [
            {"start": 3*i, "end": 3*i+1, "mark": "H3K27ac"}
            for i in range(512)]
        c["chromatin"]["barriers"] = [
            {"bead": 3*i, "direction": "right" if i % 2 else "left"}
            for i in range(512)]
        report = validate_config(c)
        self.assertTrue(report["valid"], report["errors"])
        c["cohesin"]["count"] = MAX_COHESIN
        self.assertTrue(validate_config(c)["valid"])
        c["cohesin"]["count"] = MAX_COHESIN + 1
        self.assertFalse(validate_config(c)["valid"])
        c["cohesin"]["count"] = 100
        c["genes"].append(dict(copy.deepcopy(template), id="extra_gene"))
        self.assertFalse(validate_config(c)["valid"])

    def test_cohesin_count_also_obeys_available_chromatin_sites(self):
        c = default_config()
        c["cohesin"]["count"] = c["chromatin"]["beads"] // 2
        self.assertTrue(validate_config(c)["valid"])
        c["cohesin"]["count"] += 1
        self.assertFalse(validate_config(c)["valid"])

    def test_symmetric_pair_cannot_be_counted_twice(self):
        c = default_config()
        pair = copy.deepcopy(c["interactions"][-1])
        pair["a"], pair["b"] = pair["b"], pair["a"]
        c["interactions"].append(pair)
        self.assertFalse(validate_config(c)["valid"])

    def test_malformed_types_return_errors(self):
        for value in (None, [], "configuration", {}, {"chromatin": []}):
            with self.subTest(value=value):
                self.assertFalse(validate_config(value)["valid"])

    def test_regulatory_presets_preserve_genomic_mapping_and_explicit_coarse_units(self):
        small = preset_config("nucleosome_locus")
        large = preset_config("regulatory_20mb")
        for config, span, beads, bp, representation in (
                (small, 80_000, 400, 200, "nucleosome_unit"),
                (large, 20_000_000, 2000, 10_000, "coarse")):
            with self.subTest(representation=representation):
                report = validate_config(config)
                self.assertTrue(report["valid"], report["errors"])
                self.assertEqual(report["genomic_span_bp"], span)
                self.assertEqual(config["chromatin"]["beads"], beads)
                self.assertEqual(config["chromatin"]["bp_per_bead"], bp)
                self.assertEqual(config["chromatin"]["representation"], representation)
                genes = {gene["id"]: gene for gene in config["genes"]}
                for region in config["regulatory_regions"]:
                    self.assertGreaterEqual(region["start"], 0)
                    self.assertLess(region["end"], beads)
                    if region["kind"] == "promoter":
                        self.assertLessEqual(region["start"], genes[region["target_gene"]]["start"])
                        self.assertGreaterEqual(region["end"], genes[region["target_gene"]]["start"])
        self.assertEqual(small["genes"][0]["start"] * 200, 28_000)
        self.assertEqual((small["genes"][0]["end"] + 1) * 200, 40_000)
        self.assertTrue(any("Wrapped DNA" in warning for warning in validate_config(small)["warnings"]))

    def test_twenty_mb_at_two_hundred_bp_is_rejected_at_one_hundred_thousand_beads(self):
        config = preset_config("nucleosome_locus")
        config["chromatin"]["beads"] = 100_000
        report = validate_config(config)
        self.assertFalse(report["valid"])
        self.assertEqual(report["genomic_span_bp"], 20_000_000)
        self.assertTrue(any("chromatin.beads" in error for error in report["errors"]))
        self.assertTrue(any("4000" in error for error in report["errors"]))

    def test_regulatory_references_promoter_position_and_unique_enabled_promoter(self):
        base = preset_config("nucleosome_locus")
        for field, value in (("target_gene", "missing"), ("tf_species", "missing"),
                             ("coactivator_species", "missing"), ("kind", "silencer")):
            with self.subTest(field=field):
                config = copy.deepcopy(base)
                config["regulatory_regions"][0][field] = value
                self.assertFalse(validate_config(config)["valid"])
        config = copy.deepcopy(base)
        promoter = next(region for region in config["regulatory_regions"] if region["kind"] == "promoter")
        promoter["start"] -= 3
        promoter["end"] -= 3
        report = validate_config(config)
        self.assertTrue(any("contain the target gene start" in error for error in report["errors"]))
        config = copy.deepcopy(base)
        promoter = copy.deepcopy(next(region for region in config["regulatory_regions"] if region["kind"] == "promoter"))
        promoter["id"] = "extra_promoter"
        config["regulatory_regions"].append(promoter)
        self.assertFalse(validate_config(config)["valid"])
        promoter["enabled"] = False
        self.assertTrue(validate_config(config)["valid"])

    def test_missing_tf_coactivator_rule_warns_without_inventing_recruitment(self):
        config = preset_config("nucleosome_locus")
        config["interactions"] = [rule for rule in config["interactions"]
                                  if {rule["a"], rule["b"]} != {"activator", "p300"}]
        report = validate_config(config)
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(any("recruitment requires" in warning for warning in report["warnings"]))
        config["interactions"].append({"a": "p300", "b": "activator", "strength": 1,
                                       "on_rate": 1, "off_rate": 1, "radius": 1})
        self.assertFalse(any("recruitment requires" in warning for warning in validate_config(config)["warnings"]))

    def test_state_pair_symmetry_reserved_state_and_configurable_sequence_gap(self):
        config = preset_config("nucleosome_locus")
        rule = config["chromatin"]["state_interactions"][0]
        self.assertEqual(rule["min_separation"], 3)
        rule["min_separation"] = 1
        self.assertTrue(validate_config(config)["valid"])
        rule["b"] = "H3K27ac"
        config["chromatin"]["state_interactions"].append({**rule, "a": rule["b"], "b": rule["a"]})
        report = validate_config(config)
        self.assertTrue(any("symmetric chromatin-state pair" in error for error in report["errors"]))
        config = preset_config("nucleosome_locus")
        config["chromatin"]["state_interactions"][0]["a"] = "absent_mark"
        self.assertFalse(validate_config(config)["valid"])
        config = preset_config("nucleosome_locus")
        config["marks"].append({**config["marks"][0], "id": "unmodified", "name": "Ambiguous state"})
        self.assertFalse(validate_config(config)["valid"])

    def test_malformed_optional_regulatory_and_state_rules_return_errors(self):
        for field in ("regulatory_regions", "state_interactions"):
            for value in (None, {}, [None], ["not an object"], [{}]):
                with self.subTest(field=field, value=value):
                    config = preset_config("nucleosome_locus")
                    target = config["chromatin"] if field == "state_interactions" else config
                    target[field] = value
                    report = validate_config(config)
                    self.assertFalse(report["valid"])
                    self.assertTrue(report["errors"])
        for group, field in (("regulatory_regions", "activation_gain"),
                             ("regulatory_regions", "start"),
                             ("state_interactions", "strength"),
                             ("state_interactions", "min_separation")):
            for value in (True, "1", float("nan"), float("inf"), 2**10000):
                with self.subTest(group=group, field=field, value_type=type(value).__name__):
                    config = preset_config("nucleosome_locus")
                    target = config["chromatin"] if group == "state_interactions" else config
                    target[group][0][field] = value
                    self.assertFalse(validate_config(config)["valid"])


if __name__ == "__main__":
    unittest.main()
