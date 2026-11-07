import copy
import unittest

from chromatin_studio.schema import default_config, preset_config, validate_config


class OccupancySchemaTests(unittest.TestCase):
    def test_new_presets_have_opt_in_occupancy_and_individual_rna(self):
        for name, span in (("nucleosome_transcription", 80_000), ("transcription_20mb", 20_000_000)):
            config = preset_config(name)
            report = validate_config(config)
            self.assertTrue(report["valid"], report["errors"])
            self.assertEqual(report["genomic_span_bp"], span)
            self.assertTrue(config["chromatin"]["nucleosomes"]["enabled"])
            self.assertTrue(config["transcription"]["enabled"])
            self.assertIn("remodeler", {p["id"] for p in config["proteins"]})
        self.assertNotIn("nucleosomes", default_config()["chromatin"])
        self.assertNotIn("transcription", preset_config("regulatory_20mb"))

    def test_occupancy_regions_and_probabilities_are_validated(self):
        base = preset_config("nucleosome_transcription")
        for field, value in (("enabled", 1), ("initial_occupancy", 1.01),
                             ("initial_occupancy", True), ("eviction_rate", -1),
                             ("vacant_persistence_multiplier", float("nan")),
                             ("occupied_persistence_multiplier", 1e308),
                             ("initial_regions", [{"start": 0, "end": 400, "occupied": True}]),
                             ("initial_regions", [{"start": 0, "end": 1, "occupied": 0}])):
            with self.subTest(field=field, value=str(value)):
                config = copy.deepcopy(base)
                config["chromatin"]["nucleosomes"][field] = value
                self.assertFalse(validate_config(config)["valid"])

    def test_transcript_budget_diffusion_and_gene_controls_are_bounded(self):
        for field, value in (("max_transcripts", 0), ("max_transcripts", 5001),
                             ("max_transcripts", 2.5), ("rna_diffusion", -1),
                             ("mature_rna_diffusion", float("inf")), ("enabled", "true")):
            config = preset_config("nucleosome_transcription")
            config["transcription"][field] = value
            self.assertFalse(validate_config(config)["valid"])
        for field, value in (("pause_rate", -1), ("resume_rate", 21),
                             ("occupied_elongation_multiplier", -1),
                             ("nucleosome_eviction_rate", True)):
            config = preset_config("nucleosome_transcription")
            config["genes"][0][field] = value
            self.assertFalse(validate_config(config)["valid"])

    def test_protein_preferences_and_remodeling_allow_zeros_reject_unknown_states(self):
        config = preset_config("nucleosome_transcription")
        protein = config["proteins"][0]
        protein["nucleosome_affinity"] = {"occupied": 0, "vacant": 0}
        self.assertTrue(validate_config(config)["valid"])
        protein["nucleosome_affinity"]["missing"] = 1
        self.assertFalse(validate_config(config)["valid"])
        protein["nucleosome_affinity"] = []
        self.assertFalse(validate_config(config)["valid"])
        protein["nucleosome_affinity"] = {"occupied": 1}
        protein["nucleosome_reassembly_rate"] = -1
        self.assertFalse(validate_config(config)["valid"])

    def test_coarse_occupancy_is_explained_without_mislabeling_individual_nucleosomes(self):
        report = validate_config(preset_config("transcription_20mb"))
        self.assertTrue(any("effective binary state" in text for text in report["warnings"]))
        self.assertTrue(any("cap blocks new initiation" in text for text in report["warnings"]))


if __name__ == "__main__":
    unittest.main()
