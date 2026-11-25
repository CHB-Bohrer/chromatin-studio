"""Native OpenMM bare-polymer backend: scope guard, exact force parity, dynamics."""

import copy
import importlib.util
import json
from pathlib import Path
import unittest

import numpy as np

from chromatin_studio.openmm_native import unsupported_features
from chromatin_studio.schema import checked_config

ROOT = Path(__file__).resolve().parents[1]
HAS_OPENMM = importlib.util.find_spec("openmm") is not None


def small_polymer(beads=120, radius=6.0, seed=31):
    config = json.loads((ROOT / "references/hansen/jobs/baseline-20mb-10kb-model.json").read_text())
    config["chromatin"]["beads"] = beads
    config["simulation"]["confinement_radius"] = radius
    config["seed"] = seed
    return checked_config(config)


def platform_registered(name):
    import openmm as mm
    return name in [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]


class ScopeTests(unittest.TestCase):
    def test_bare_polymer_is_supported(self):
        self.assertEqual(unsupported_features(small_polymer()), [])

    def test_every_chemistry_feature_is_reported(self):
        base = small_polymer()
        cases = {
            "protein particles": lambda c: c["proteins"].append({"id": "p", "count": 3}),
            "histone marks": lambda c: c["marks"].append({"id": "m"}),
            "state-selective chromatin attraction": lambda c: c["chromatin"].update(state_interactions=[{"a": "m"}]),
            "nuclear bodies": lambda c: c["bodies"].append({"id": "b"}),
            "cohesin extrusion": lambda c: c["cohesin"].update(enabled=True),
            "genes": lambda c: c["genes"].append({"id": "g"}),
            "explicit transcription": lambda c: c["transcription"].update(enabled=True),
            "nucleosome occupancy dynamics": lambda c: c["chromatin"]["nucleosomes"].update(enabled=True),
        }
        for reason, mutate in cases.items():
            config = copy.deepcopy(base)
            mutate(config)
            self.assertIn(reason, unsupported_features(config), reason)

    @unittest.skipUnless(HAS_OPENMM, "OpenMM unavailable")
    def test_chemistry_models_are_refused_not_silently_simplified(self):
        from chromatin_studio.openmm_native import NativePolymerSimulation
        config = small_polymer()
        config["cohesin"]["enabled"] = True
        with self.assertRaisesRegex(ValueError, "bare polymers only"):
            NativePolymerSimulation(config, platform="Reference")


@unittest.skipUnless(HAS_OPENMM, "OpenMM unavailable")
class ParityTests(unittest.TestCase):
    def states(self, sim):
        """Relaxed, wall-penetrating/stretched and overlapping conformations."""
        sim.step(400)
        relaxed = sim.positions.copy()
        return {"relaxed": relaxed, "expanded": 1.9 * relaxed, "compressed": 0.6 * relaxed}

    def assert_parity(self, platform, precision, force_tol, energy_tol):
        from chromatin_studio.engine import Simulation
        from chromatin_studio.openmm_native import NativePolymerSimulation
        config = small_polymer()
        reference = Simulation(config)
        native = NativePolymerSimulation(config, platform=platform, precision=precision)
        # Same seed -> exactly the reference engine's initial conformation.
        np.testing.assert_array_equal(native.positions, reference.positions)
        for label, positions in self.states(reference).items():
            native.set_positions(positions)
            device_force, device_energy = native.device_forces_and_energy()
            host_force, host_energy = native._sim._mechanics()
            scale = np.max(np.linalg.norm(host_force, axis=1))
            error = np.max(np.linalg.norm(device_force - host_force, axis=1))
            self.assertLess(error / scale, force_tol, label)
            self.assertLess(abs(device_energy - host_energy) / abs(host_energy), energy_tol, label)

    def test_reference_platform_reproduces_engine_forces_exactly(self):
        self.assert_parity("Reference", "double", 1e-12, 1e-12)

    def test_cuda_mixed_precision_parity(self):
        if not platform_registered("CUDA"):
            self.skipTest("CUDA plugin not registered; no device test was run")
        try:
            self.assert_parity("CUDA", "mixed", 1e-4, 1e-5)
        except Exception as exc:  # unavailable device is a skip, not a pass
            if "could not initialize" in str(exc):
                self.skipTest(f"CUDA device could not initialize: {exc}")
            raise

    def test_noise_seed_controls_only_the_noise_stream(self):
        from chromatin_studio.openmm_native import NativePolymerSimulation
        config = small_polymer()
        runs = {}
        for label, noise_seed in (("a", None), ("b", None), ("c", 999)):
            sim = NativePolymerSimulation(config, platform="Reference", noise_seed=noise_seed)
            runs[label + "0"] = sim.positions.copy()
            sim.step(20)
            runs[label] = sim.positions.copy()
        # Initial conformation always follows the model seed.
        np.testing.assert_array_equal(runs["a0"], runs["c0"])
        np.testing.assert_array_equal(runs["a"], runs["b"])  # same seed: reproducible
        self.assertGreater(np.abs(runs["a"] - runs["c"]).max(), 1e-3)  # new noise stream

    def test_bookkeeping_and_free_diffusion_scale(self):
        from chromatin_studio.openmm_native import NativePolymerSimulation
        config = small_polymer(beads=3000, radius=60.0)
        config["chromatin"].update(bond_k=0.0, persistence_length=0.0, repulsion=0.0)
        config["simulation"]["confinement_stiffness"] = 0.0
        config["simulation"]["dt"] = 0.05
        sim = NativePolymerSimulation(config, platform="Reference")
        start = sim.positions.copy()
        sim.step(4)
        self.assertEqual(sim.step_count, 4)
        self.assertAlmostEqual(sim.time, 0.2)
        msd = np.mean(np.sum((sim.positions - start) ** 2, axis=1))
        expected = 6 * float(config["chromatin"]["diffusion"]) * sim.time
        self.assertLess(abs(msd / expected - 1), 0.07)
        self.assertTrue(sim.backend_metadata["capabilities"]["native_gpu_forces"])


if __name__ == "__main__":
    unittest.main()
