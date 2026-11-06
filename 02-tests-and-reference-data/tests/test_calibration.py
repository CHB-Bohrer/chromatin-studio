"""Independent checks of measured MSD, clock recovery and saved unit provenance."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from chromatin_studio.calibration import (
    analyze_trajectory, config_fingerprint, fit_time_mapping, main, unit_summary, validate_profile,
)
from chromatin_studio.schema import default_config


def fit_request(experimental=True, exponent=.7, time_scale=2.5, length_nm=40):
    simulation_lags = np.geomspace(.01, 100, 80)
    physical_lags = np.geomspace(.1, 10, 24)
    return {"simulation": {"lag": simulation_lags.tolist(), "msd": (3*simulation_lags**exponent).tolist(),
               "dimension": 3, "observable": "locus", "drift_correction": "none", "source": "Test simulated power law"},
            "experiment": {"lag": physical_lags.tolist(),
               "msd": ((length_nm/1000)**2*3*(physical_lags/time_scale)**exponent).tolist(),
               "dimension": 3, "observable": "locus", "drift_correction": "none", "source": "Declared test reference; not biological measurements"},
            "length_nm": length_nm, "length_source": "Test length assumption",
            "fit_min_s": .1, "fit_max_s": 10, "source_kind": "experimental" if experimental else "synthetic", "notes": "Unit test fixture"}


def trajectory():
    times = np.arange(12, dtype=float)
    positions = np.zeros((12, 3, 3))
    positions[:, 0, 0] = times
    positions[:, 1, 0] = -times
    positions[:, 2, 2] = 2*times
    return {"time": times.tolist(), "positions": positions.tolist(), "beads": [0, 1],
            "dimension": 3, "observable": "locus", "drift_correction": "none", "source": "Deterministic trajectory for testing"}


class ClockFitTests(unittest.TestCase):
    def test_recovers_known_fixed_length_clock(self):
        config = default_config()
        report = fit_time_mapping(fit_request(), config)
        self.assertEqual(report["profile"]["status"], "fitted")
        self.assertAlmostEqual(report["profile"]["time_s"], 2.5, places=10)
        self.assertLess(report["diagnostics"]["log_rmse"], 1e-12)
        self.assertEqual(validate_profile(report["profile"]), report["profile"])
        json.dumps(report, allow_nan=False)

    def test_diffusive_formula_agrees_with_curve_fit(self):
        request = fit_request(exponent=1, time_scale=4, length_nm=50)
        report = fit_time_mapping(request)
        # MSDsim=6*Dsim*t=3*t, and Dphys=l0²Dsim/t0.
        d_sim = .5
        d_phys = (.05**2)*d_sim/4
        self.assertAlmostEqual(report["profile"]["time_s"], .05**2*d_sim/d_phys, places=10)

    def test_recovers_shift_across_non_power_law_interpolation_segments(self):
        request = fit_request(time_scale=1.7)
        xs = np.log(request["simulation"]["lag"])
        # Varying logarithmic slope exercises interval changes in the global
        # shift optimization, rather than only a straight log-log power law.
        ys = .7*xs + .2*np.sin(2*xs)
        request["simulation"]["msd"] = np.exp(ys).tolist()
        xe = np.log(request["experiment"]["lag"])
        request["experiment"]["msd"] = ((.04**2)*np.exp(np.interp(xe-np.log(1.7), xs, ys))).tolist()
        report = fit_time_mapping(request)
        self.assertAlmostEqual(report["profile"]["time_s"], 1.7, places=10)
        self.assertLess(report["diagnostics"]["log_rmse"], 1e-12)

    def test_synthetic_is_never_experimental_fit(self):
        report = fit_time_mapping(fit_request(experimental=False))
        self.assertEqual(report["profile"]["status"], "demo")
        self.assertTrue(any("Synthetic" in warning for warning in report["warnings"]))
        profile = copy.deepcopy(report["profile"])
        profile["status"] = "fitted"
        with self.assertRaises(ValueError):
            validate_profile(profile)

    def test_slope_mismatch_requires_review(self):
        request = fit_request()
        t = np.asarray(request["experiment"]["lag"])
        request["experiment"]["msd"] = ((.04**2)*3*(t/2.5)**.3).tolist()
        report = fit_time_mapping(request)
        self.assertEqual(report["profile"]["status"], "needs_review")
        self.assertFalse(report["diagnostics"]["shape_compatible"])
        self.assertAlmostEqual(report["diagnostics"]["slope_difference"], .4, places=10)

    def test_flat_curve_cannot_identify_clock(self):
        request = fit_request()
        request["simulation"]["msd"] = [1.0]*80
        request["experiment"]["msd"] = [.04**2]*24
        report = fit_time_mapping(request)
        self.assertEqual(report["profile"]["status"], "needs_review")
        self.assertFalse(report["diagnostics"]["time_identifiable"])

    def test_no_extrapolation_or_boundary_acceptance(self):
        request = fit_request(time_scale=1e5)
        report = fit_time_mapping(request)
        self.assertEqual(report["profile"]["status"], "needs_review")
        self.assertTrue(report["diagnostics"]["boundary_optimum"])
        request["simulation"]["lag"] = [1, 2, 3, 4]
        request["simulation"]["msd"] = [1, 2, 3, 4]
        with self.assertRaisesRegex(ValueError, "entire experimental"):
            fit_time_mapping(request)

    def test_curve_conventions_must_match(self):
        for field, value in (("dimension", 2), ("observable", "center_of_mass"), ("drift_correction", "center_of_mass")):
            request = fit_request()
            request["experiment"][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                fit_time_mapping(request)
        request = fit_request()
        for key in ("simulation", "experiment"):
            request[key]["observable"] = "relative_pair"
        self.assertEqual(fit_time_mapping(request)["profile"]["status"], "fitted")

    def test_malformed_and_unattributed_inputs_rejected(self):
        mutations = [lambda r: r.pop("source_kind"), lambda r: r.update(length_nm=True),
                     lambda r: r["simulation"].update(source=""),
                     lambda r: r["simulation"]["msd"].__setitem__(0, float("nan")),
                     lambda r: r["simulation"]["lag"].__setitem__(1, r["simulation"]["lag"][0]),
                     lambda r: r["experiment"].update(lag_unit="ms"),
                     lambda r: r.update(fit_min_s=2, fit_max_s=2.05)]
        for change in mutations:
            request = fit_request()
            change(request)
            with self.subTest(change=change), self.assertRaises(ValueError):
                fit_time_mapping(request)

    def test_saved_fit_cannot_silently_change_time_or_diagnostics(self):
        for field in ("time", "diagnostics", "status"):
            profile = fit_time_mapping(fit_request())["profile"]
            if field == "time":
                profile["time_s"] *= 2
            elif field == "diagnostics":
                profile["diagnostics"]["log_rmse"] = .2
            else:
                profile["status"] = "needs_review"
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_profile(profile)


class TrajectoryMSDTests(unittest.TestCase):
    def test_locus_com_relative_pair_are_distinct_observables(self):
        request = trajectory()
        result = analyze_trajectory(request)
        np.testing.assert_allclose(result["msd"], np.asarray(result["lag"])**2)
        request["observable"] = "center_of_mass"
        np.testing.assert_array_equal(analyze_trajectory(request)["msd"], 0)
        request["observable"] = "relative_pair"
        relative = analyze_trajectory(request)
        np.testing.assert_allclose(relative["msd"], 4*np.asarray(relative["lag"])**2)
        self.assertEqual(relative["counts"][0], 11)

    def test_projection_and_whole_chromatin_com_subtraction(self):
        request = trajectory()
        request["beads"] = [2]
        result = analyze_trajectory(request)
        np.testing.assert_allclose(result["msd"], 4*np.asarray(result["lag"])**2)
        request["dimension"] = 2
        np.testing.assert_array_equal(analyze_trajectory(request)["msd"], 0)
        request["dimension"] = 3
        request["drift_correction"] = "center_of_mass"
        corrected = analyze_trajectory(request)
        np.testing.assert_allclose(corrected["msd"], (4/3*np.asarray(corrected["lag"]))**2)

    def test_irregular_final_frame_requires_explicit_drop(self):
        request = trajectory()
        request["time"][-1] -= .5
        with self.assertRaisesRegex(ValueError, "regular frame spacing"):
            analyze_trajectory(request)
        request["drop_final_irregular"] = True
        result = analyze_trajectory(request)
        self.assertTrue(result["analysis"]["dropped_final_frame"])
        self.assertEqual(result["analysis"]["retained_frames"], 11)
        request["time"][3] += .1
        with self.assertRaises(ValueError):
            analyze_trajectory(request)

    def test_invalid_selection_or_conventions_rejected(self):
        for beads in ([99], [0, 0], [True], [{}]):
            request = trajectory()
            request["beads"] = beads
            with self.subTest(beads=beads), self.assertRaises(ValueError):
                analyze_trajectory(request)
        request = trajectory()
        request.update(observable="relative_pair", beads=[0])
        with self.assertRaises(ValueError):
            analyze_trajectory(request)
        request.update(beads=[0, 1], drift_correction="center_of_mass")
        with self.assertRaises(ValueError):
            analyze_trajectory(request)

    def test_snapshot_cli_output_is_compatible_and_preserves_provenance(self):
        request = trajectory()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            snapshots, output = root/"snapshots.jsonl", root/"msd.json"
            snapshots.write_text("\n".join(json.dumps({"time": t, "chromatin": x})
                                            for t, x in zip(request["time"], request["positions"])))
            config = default_config()
            (root/"model.json").write_text(json.dumps(config))
            self.assertEqual(main(["msd", "--snapshots", str(snapshots), "--beads", "0,1",
                                   "--observable", "relative_pair", "--output", str(output)]), 0)
            result = json.loads(output.read_text())
            self.assertEqual(result["baseline_fingerprint"], config_fingerprint(config))
            self.assertEqual(result["observable"], "relative_pair")


class UnitMappingTests(unittest.TestCase):
    def test_clock_is_retained_for_perturbations_with_explicit_warning(self):
        config = default_config()
        profile = fit_time_mapping(fit_request(), config)["profile"]
        summary = unit_summary(config, profile)
        self.assertTrue(summary["baseline_matches"])
        self.assertAlmostEqual(summary["simulation"]["dt_s"], .002*2.5)
        self.assertAlmostEqual(summary["chromatin"]["bare_bead_diffusion_um2_s"], .08*.04**2/2.5)
        self.assertAlmostEqual(summary["genes"][0]["elongation_nt_s"], 3*1000/2.5)
        self.assertAlmostEqual(summary["proteins"][0]["on_rate_s_inv"], 1.2/2.5)
        original_hash = config_fingerprint(config)
        config.update(seed=999, name="renamed", calibration=profile)
        self.assertEqual(config_fingerprint(config), original_hash)
        config["chromatin"]["diffusion"] *= 2
        changed = unit_summary(config, profile)
        self.assertFalse(changed["baseline_matches"])
        self.assertEqual(changed["time_s"], profile["time_s"])
        self.assertTrue(any("retained" in warning for warning in changed["warnings"]))
        self.assertEqual(unit_summary(config, None)["status"], "uncalibrated")

    def test_curve_with_different_baseline_is_not_accepted(self):
        config, other = default_config(), default_config()
        other["chromatin"]["diffusion"] *= 2
        request = fit_request()
        request["simulation"]["baseline_fingerprint"] = config_fingerprint(other)
        profile = fit_time_mapping(request, config)["profile"]
        self.assertEqual(profile["status"], "needs_review")
        self.assertFalse(profile["diagnostics"]["curve_baseline_matches"])
        self.assertEqual(validate_profile(profile), profile)

    def test_assumed_mapping_is_explicit_and_bounded(self):
        profile = {"version": 1, "length_nm": 10, "time_s": .1, "status": "assumed",
                   "source_kind": "assumed", "sources": {"simulation": "none", "experiment": "none", "length": "Assumption"},
                   "notes": "", "diagnostics": None, "baseline_fingerprint": None, "raw_input": None}
        self.assertEqual(validate_profile(profile)["status"], "assumed")
        profile["time_s"] = 1e-100
        with self.assertRaises(ValueError):
            validate_profile(profile)


if __name__ == "__main__":
    unittest.main()
