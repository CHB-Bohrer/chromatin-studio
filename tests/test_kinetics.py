import copy
import json
import math
import unittest

from chromatin_studio.kinetics import AVOGADRO, estimate_binding
from chromatin_studio.schema import default_config


class KineticsTests(unittest.TestCase):
    def setUp(self):
        self.config = default_config()
        self.request = {"a": "srrm2", "b": "activator", "length_nm": 100,
                        "time_s": 0.01, "target_kd_um": 2, "efficiency": 0.25}

    def estimate(self, **updates):
        return estimate_binding(self.config, dict(self.request, **updates))

    def test_absolute_units_and_molar_conversion(self):
        result = self.estimate()
        # A 0.1 um length unit and 0.01 s time unit give a D scale of 1 um²/s.
        self.assertAlmostEqual(result["physical_D_a_um2_s"], 0.35)
        self.assertAlmostEqual(result["physical_D_b_um2_s"], 0.8)
        self.assertAlmostEqual(result["physical_relative_D_um2_s"], 1.15)
        self.assertAlmostEqual(result["physical_radius_um"], 0.11)
        expected_volume_kernel = 4 * math.pi * 1.15 * 0.11
        expected_molar = expected_volume_kernel * 6.02214076e8
        self.assertAlmostEqual(result["kon_limit_um3_s"], expected_volume_kernel)
        self.assertAlmostEqual(result["kon_limit_M_inv_s"] / expected_molar, 1)
        self.assertAlmostEqual(result["kon_M_inv_s"] / expected_molar, 0.25)
        self.assertEqual(result["association_stoichiometric_factor"], 1)
        self.assertEqual(result["kon_M_inv_s"], result["macroscopic_kon_M_inv_s"])
        self.assertAlmostEqual(result["nominal_koff_s_inv"], 2e-6 * expected_molar * 0.25)
        self.assertAlmostEqual(result["nominal_koff_reduced"], result["nominal_koff_s_inv"] * 0.01)
        self.assertEqual(AVOGADRO, 6.02214076e23)

    def test_identical_molecule_convention_counts_each_pair_once(self):
        homo = self.estimate(a="hp1", b="hp1", efficiency=1)
        self.assertEqual(homo["pair"]["kd_convention"], "[A]^2/[A2]")
        self.assertTrue(homo["pair"]["homodimer"])
        self.assertAlmostEqual(homo["physical_relative_D_um2_s"], 1.0)
        self.assertEqual(homo["association_stoichiometric_factor"], 0.5)
        self.assertEqual(homo["macroscopic_kon_M_inv_s"], 0.5 * homo["kon_M_inv_s"])
        self.assertAlmostEqual(homo["nominal_koff_s_inv"] / homo["kon_M_inv_s"], 1e-6)
        # The same encounter kernel with unlike labels has twice the event rate.
        twin = copy.deepcopy(self.config["proteins"][0])
        twin.update(id="hp1_other", name="Other HP1")
        self.config["proteins"].append(twin)
        self.config["interactions"].append(dict(self.config["interactions"][0], b="hp1_other"))
        hetero = self.estimate(a="hp1", b="hp1_other", efficiency=1)
        self.assertEqual(hetero["kon_M_inv_s"], homo["kon_M_inv_s"])
        self.assertEqual(hetero["nominal_koff_s_inv"], 2 * homo["nominal_koff_s_inv"])

    def test_symmetric_pair_lookup_preserves_physics(self):
        ab = self.estimate()
        ba = self.estimate(a="activator", b="srrm2")
        for key in ("physical_relative_D_um2_s", "physical_radius_um", "kon_M_inv_s",
                    "macroscopic_kon_M_inv_s", "nominal_koff_s_inv"):
            self.assertEqual(ab[key], ba[key])
        self.assertEqual(ab["physical_D_a_um2_s"], ba["physical_D_b_um2_s"])

    def test_idr_diffusion_is_molecule_com_not_domain_diffusion(self):
        result = self.estimate()
        self.config["proteins"][2]["idr_beads"] = 12
        longer = self.estimate()
        self.assertEqual(result["kon_M_inv_s"], longer["kon_M_inv_s"])
        self.assertAlmostEqual(result["physical_D_a_um2_s"], 0.35)
        self.assertTrue(any("capture radius acts between sticker/domain beads" in text
                            for text in result["limitations"]))

    def test_chromatin_uses_engine_capture_radius_and_bare_bead_mobility(self):
        result = self.estimate(a="hp1", b="chromatin")
        self.assertEqual(result["pair"]["kind"], "protein_chromatin_proxy")
        self.assertAlmostEqual(result["model"]["radius_reduced"], 0.28 + 0.85)
        self.assertAlmostEqual(result["physical_radius_um"], 0.113)
        self.assertAlmostEqual(result["physical_D_b_um2_s"], 0.08)
        self.assertEqual(result["association_stoichiometric_factor"], 1)
        self.assertEqual(result["model"]["unmodified_affinity_multiplier"], 0.04)
        self.assertTrue(any("does not infer a DNA-site Kd" in text for text in result["limitations"]))
        # Neither a histone multiplier nor bp resolution supplies physical affinity.
        self.config["proteins"][0]["chromatin_affinity"]["H3K9me3"] = 100
        self.config["chromatin"]["bp_per_bead"] = 100000
        self.assertEqual(result["kon_M_inv_s"], self.estimate(a="hp1", b="chromatin")["kon_M_inv_s"])

    def test_contact_hazards_do_not_masquerade_as_molar_kon(self):
        original = self.estimate()
        self.config["interactions"][-1].update(on_rate=0, off_rate=17)
        result = self.estimate()
        for key in ("kon_M_inv_s", "macroscopic_kon_M_inv_s", "nominal_koff_s_inv"):
            self.assertEqual(original[key], result[key])
        self.assertEqual(result["model"]["on_rate_reduced"], 0)
        self.assertEqual(result["model"]["off_rate_reduced"], 17)
        self.assertTrue(any("zero association hazard" in text for text in result["limitations"]))

    def test_length_time_efficiency_scalings(self):
        initial = self.estimate()
        length = self.estimate(length_nm=200)
        time = self.estimate(time_s=0.02)
        efficiency = self.estimate(efficiency=0.5)
        self.assertAlmostEqual(length["kon_M_inv_s"] / initial["kon_M_inv_s"], 8)
        self.assertAlmostEqual(time["kon_M_inv_s"] / initial["kon_M_inv_s"], 0.5)
        self.assertAlmostEqual(time["nominal_koff_reduced"], initial["nominal_koff_reduced"])
        self.assertAlmostEqual(efficiency["kon_M_inv_s"] / initial["kon_M_inv_s"], 2)
        self.assertEqual(efficiency["kon_limit_M_inv_s"], initial["kon_limit_M_inv_s"])

    def test_efficiency_defaults_to_absorbing_limit_only(self):
        request = dict(self.request)
        del request["efficiency"]
        result = estimate_binding(self.config, request)
        self.assertEqual(result["inputs"]["efficiency"], 1)
        self.assertEqual(result["kon_M_inv_s"], result["kon_limit_M_inv_s"])
        for required in ("length_nm", "time_s", "target_kd_um"):
            with self.subTest(required=required):
                incomplete = dict(request)
                del incomplete[required]
                with self.assertRaisesRegex(ValueError, required):
                    estimate_binding(self.config, incomplete)

    def test_finite_positive_calibration_and_efficiency_validation(self):
        for field in ("length_nm", "time_s", "target_kd_um", "efficiency"):
            for value in (0, -1, True, "1", None, [], float("nan"), float("inf"), 10**400):
                with self.subTest(field=field, value=str(value)[:24]):
                    with self.assertRaisesRegex(ValueError, field):
                        self.estimate(**{field: value})
        with self.assertRaisesRegex(ValueError, "efficiency"):
            self.estimate(efficiency=1.01)

    def test_extreme_calibrations_reject_overflow_and_underflow(self):
        for updates in ({"length_nm": 1e308}, {"time_s": 1e-320},
                        {"length_nm": 1e-320}, {"target_kd_um": 1e-320},
                        {"target_kd_um": 1e308},
                        {"length_nm": 1e-80, "efficiency": 1e-200}):
            with self.subTest(updates=updates):
                with self.assertRaisesRegex(ValueError, "unrepresentable"):
                    self.estimate(**updates)

    def test_stationary_partner_allowed_but_no_relative_diffusion_rejected(self):
        self.config["proteins"][2]["diffusion"] = 0
        result = self.estimate()
        self.assertEqual(result["physical_D_a_um2_s"], 0)
        self.assertAlmostEqual(result["physical_relative_D_um2_s"], 0.8)
        self.config["proteins"][4]["diffusion"] = 0
        with self.assertRaisesRegex(ValueError, "relative diffusion"):
            self.estimate()

    def test_invalid_request_and_pair_fail_without_guessing_a_radius(self):
        for request in (None, [], "protein"):
            with self.subTest(request=request), self.assertRaises(ValueError):
                estimate_binding(self.config, request)
        for updates in ({"a": None}, {"b": ""}, {"a": "unknown"},
                        {"b": "unknown"}, {"a": "hp1", "b": "activator"}):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                self.estimate(**updates)
        self.config["proteins"][1]["id"] = "chromatin"
        with self.assertRaisesRegex(ValueError, "conflicts"):
            self.estimate(b="chromatin")

    def test_read_only_repeatable_and_finite_json(self):
        config_copy, request_copy = copy.deepcopy(self.config), copy.deepcopy(self.request)
        first = self.estimate()
        self.assertEqual(self.config, config_copy)
        self.assertEqual(self.request, request_copy)
        self.assertEqual(first, self.estimate())
        json.dumps(first, allow_nan=False)
        # Returned metadata cannot mutate a later response through global aliases.
        first["sources"][0]["title"] = "mutated"
        first["limitations"].append("mutated")
        second = self.estimate()
        self.assertNotEqual(second["sources"][0]["title"], "mutated")
        self.assertNotIn("mutated", second["limitations"])


if __name__ == "__main__":
    unittest.main()
