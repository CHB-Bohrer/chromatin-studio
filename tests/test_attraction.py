"""Attraction-mode protein interactions: schema, CPU force, scope, GPU parity."""

import copy
import importlib.util
import unittest

import numpy as np

from chromatin_studio.engine import Simulation
from chromatin_studio.openmm_native import unsupported_features
from chromatin_studio.schema import checked_config, preset_config, validate_config

HAS_OPENMM = importlib.util.find_spec("openmm") is not None


def small_condensate(scaffolds=120, clients=30, seed=5):
    config = preset_config("phase_separation")
    config["chromatin"]["beads"] = 40
    config["simulation"]["confinement_radius"] = 6.0
    config["proteins"][0]["count"] = scaffolds
    config["proteins"][1]["count"] = clients
    config["seed"] = seed
    return checked_config(config)


class SchemaTests(unittest.TestCase):
    def test_preset_is_valid_and_warns_about_interpretation(self):
        report = validate_config(preset_config("phase_separation"))
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(any("Attraction-mode" in w for w in report["warnings"]))

    def test_rejects_unknown_mode_and_rates_in_attraction_mode(self):
        for mutate, text in ((lambda r: r.update(mode="glue"), "mode must be bonds or attraction"),
                             (lambda r: r.update(on_rate=1.0), "on_rate must be 0 in attraction mode"),
                             (lambda r: r.update(radius=0.0), "radius"),
                             (lambda r: r.update(radius=0.5), "must exceed the contact distance")):
            config = preset_config("phase_separation")
            mutate(config["interactions"][0])
            report = validate_config(config)
            self.assertFalse(report["valid"])
            self.assertTrue(any(text in e for e in report["errors"]), report["errors"])

    def test_bond_mode_models_are_unchanged(self):
        config = preset_config("nuclear")
        self.assertTrue(all("mode" not in r for r in config["interactions"]))
        self.assertIsNone(Simulation(config)._attraction)


class ForceTests(unittest.TestCase):
    def test_force_is_minus_energy_gradient(self):
        sim = Simulation(small_condensate())
        sim.positions[:] *= 0.75         # crowd the proteins so many wells and cores overlap
        force, _ = sim._mechanics()
        h = 1e-7
        for index in sim._attraction["particles"][:25]:
            for d in range(3):
                sim.positions[index, d] += h
                up = sim._mechanics()[1]
                sim.positions[index, d] -= 2 * h
                down = sim._mechanics()[1]
                sim.positions[index, d] += h
                numeric = -(up - down) / (2 * h)
                self.assertLess(abs(numeric - force[index, d]), 1e-4 * max(1.0, abs(force[index, d])))

    def test_well_value_and_cutoff(self):
        sim = Simulation(small_condensate(scaffolds=2, clients=0))
        a, b = sim._attraction["particles"]
        far = np.array([50.0, 0, 0])
        sim.positions[:] = far + np.arange(len(sim.positions))[:, None] * 3.0
        depth = sim._attraction["eps"][0, 0] if sim._attraction["names"][0] == "scaffold" else sim._attraction["eps"][1, 1]
        # contact 0.6, reach 1.0: -depth at contact, zero beyond reach; inside
        # contact a WCA core adds 4*((w/r)^12-(w/r)^6)+1 with (w/r)^6 = (0.6/r)^6/2.
        core = lambda r: 4 * (((0.6 / r) ** 6 / 2) ** 2 - (0.6 / r) ** 6 / 2) + 1
        for r, expected in ((0.55, -depth + core(0.55)), (0.6, -depth), (0.8, -depth * (1 - 0.25) ** 2), (1.2, 0.0)):
            sim.positions[a] = [0, 0, 0]
            sim.positions[b] = [r, 0, 0]
            _, energy, _, _ = sim._attraction_terms()
            self.assertAlmostEqual(energy, expected, places=12)

    def test_inside_contact_only_the_core_acts_and_it_repels(self):
        sim = Simulation(small_condensate(scaffolds=2, clients=0))
        a, b = sim._attraction["particles"]
        sim.positions[:] = 50.0 + np.arange(len(sim.positions))[:, None] * 3.0
        sim.positions[a] = [0, 0, 0]
        sim.positions[b] = [0.5, 0, 0]
        force, _, _, _ = sim._attraction_terms()
        self.assertLess(force[a, 0], 0)      # pushed apart
        self.assertGreater(force[b, 0], 0)

    def test_attraction_creates_no_bonds(self):
        sim = Simulation(small_condensate())
        sim.step(300)
        self.assertEqual([b for b in sim.bindings if b["kind"] == "protein"], [])
        self.assertTrue(np.all(np.isfinite(sim.positions)))

    def test_same_molecule_beads_do_not_attract(self):
        config = small_condensate(scaffolds=1, clients=0)
        config["proteins"][0].update(representation="idr", idr_beads=4, sticker_fraction=1.0)
        sim = Simulation(checked_config(config))
        chain = sim._attraction["particles"]
        self.assertEqual(len(chain), 4)
        sim.positions[chain] = [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]]
        force, energy, _, pairs = sim._attraction_terms()
        self.assertEqual(energy, 0.0)
        self.assertEqual(len(pairs), 0)
        self.assertFalse(np.any(force))


class ScopeTests(unittest.TestCase):
    def test_attraction_only_proteins_run_natively(self):
        self.assertEqual(unsupported_features(preset_config("phase_separation")), [])

    def test_reactive_proteins_still_need_the_host(self):
        config = copy.deepcopy(preset_config("phase_separation"))
        config["proteins"][0].update(on_rate=1.0, chromatin_affinity={"unmodified": 0.5})
        self.assertIn("protein-chromatin binding", unsupported_features(config))


@unittest.skipUnless(HAS_OPENMM, "OpenMM unavailable")
class ParityTests(unittest.TestCase):
    def assert_parity(self, platform, precision, force_tol, energy_tol):
        from chromatin_studio.openmm_native import NativePolymerSimulation
        config = small_condensate()
        reference = Simulation(config)
        native = NativePolymerSimulation(config, platform=platform, precision=precision)
        for scale in (1.0, 0.6, 0.4):
            native.set_positions(reference.positions * scale)
            device_force, device_energy = native.device_forces_and_energy()
            host_force, host_energy = native._sim._mechanics()
            self.assertGreater(len(native._sim._attraction_terms()[3]), 0)
            scale_force = np.max(np.linalg.norm(host_force, axis=1))
            self.assertLess(np.max(np.linalg.norm(device_force - host_force, axis=1)) / scale_force, force_tol)
            self.assertLess(abs(device_energy - host_energy) / abs(host_energy), energy_tol)

    def test_reference_platform_matches_engine(self):
        self.assert_parity("Reference", "double", 1e-12, 1e-12)

    def test_cuda_mixed_precision_parity(self):
        import openmm as mm
        names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
        if "CUDA" not in names:
            self.skipTest("CUDA plugin not registered; no device test was run")
        try:
            self.assert_parity("CUDA", "mixed", 1e-4, 1e-5)
        except Exception as exc:
            if "could not initialize" in str(exc):
                self.skipTest(f"CUDA device could not initialize: {exc}")
            raise


if __name__ == "__main__":
    unittest.main()
