"""Device selection contract plus optional real accelerator smoke checks."""

import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np

from chromatin_studio.openmm_backend import (
    BackendUnavailableError, PolychromDynamics, backend_status, validate_backend_options,
)


def particle_model(n=4):
    sim = SimpleNamespace(positions=np.zeros((n, 3)), diffusion=np.full(n, 0.2), dt=0.05)
    sim.get_forces = lambda: np.zeros_like(sim.positions)
    return sim


class BackendOptionsTests(unittest.TestCase):
    def test_existing_defaults_and_case_normalization(self):
        self.assertEqual(validate_backend_options(),
                         {"platform": "CPU", "device_index": None, "precision": "double"})
        self.assertEqual(validate_backend_options(platform="cuda", device_index=2, precision="mixed"),
                         {"platform": "CUDA", "device_index": "2", "precision": "mixed"})

    def test_bad_options_fail_before_dependency_loading(self):
        cases = [dict(platform="automatic"), dict(platform=None), dict(precision="half"),
                 dict(device_index=0), dict(platform="Reference", device_index="0"),
                 *[dict(platform="CUDA", device_index=value)
                   for value in (-1, True, 1.5, "0,1", "-1", "gpu0", "")]]
        with patch("chromatin_studio.openmm_backend._dependencies") as dependencies:
            for options in cases:
                with self.subTest(options=options), self.assertRaises(ValueError):
                    PolychromDynamics(particle_model(), **options)
            dependencies.assert_not_called()

    def _initialize_fake(self, *, actual_platform=None, reported_precision=None, **options):
        requested_platform = options.get("platform", "CPU")
        observed = {}
        mm, unit = MagicMock(), MagicMock()
        unit.MOLAR_GAS_CONSTANT_R.value_in_unit.return_value = 0.008314

        def constructor(**kwargs):
            observed["kwargs"] = kwargs
            props = {}
            if kwargs["platform"] in ("CUDA", "OpenCL") and kwargs["GPU"] != "default":
                props.update(DeviceIndex=kwargs["GPU"], Precision=kwargs["precision"])
            polymer = SimpleNamespace(properties=props)
            polymer.set_data = MagicMock()
            polymer.add_force = MagicMock()

            def apply_forces():
                observed["properties_at_context_creation"] = dict(props)
                actual_properties = dict(props)
                if requested_platform in ("CUDA", "OpenCL"):
                    actual_properties.setdefault("DeviceIndex", "0")
                    actual_properties["DeviceName"] = "Fake device for selection contract test"
                    if reported_precision is not None:
                        actual_properties["Precision"] = reported_precision
                platform = SimpleNamespace(
                    getName=lambda: actual_platform or requested_platform,
                    getPropertyNames=lambda: list(actual_properties),
                    getPropertyValue=lambda context, name: actual_properties[name],
                )
                polymer.context = SimpleNamespace(getPlatform=lambda: platform)

            polymer._apply_forces = apply_forces
            return polymer

        with patch("chromatin_studio.openmm_backend._dependencies",
                   return_value=(mm, unit, SimpleNamespace(Simulation=constructor))):
            dynamics = PolychromDynamics(particle_model(), **options)
        return dynamics, observed

    def test_unset_device_does_not_override_scheduler_and_precision_is_applied(self):
        dynamics, observed = self._initialize_fake(platform="CUDA", precision="mixed")
        self.assertEqual(observed["kwargs"]["GPU"], "default")
        self.assertEqual(observed["properties_at_context_creation"], {"Precision": "mixed"})
        self.assertIsNone(dynamics.metadata["requested_device_index"])
        self.assertEqual(dynamics.metadata["device_index"], "0")
        self.assertEqual(dynamics.metadata["precision"], "mixed")
        self.assertTrue(dynamics.metadata["gpu_initialization_succeeded"])
        self.assertFalse(dynamics.metadata["capabilities"]["native_gpu_forces"])

    def test_explicit_logical_device_is_passed_through(self):
        dynamics, observed = self._initialize_fake(platform="OpenCL", device_index=3, precision="single")
        self.assertEqual(observed["properties_at_context_creation"],
                         {"DeviceIndex": "3", "Precision": "single"})
        self.assertEqual(dynamics.metadata["device_index"], "3")
        # OpenCL can initialize on a CPU; registration/Context success alone
        # must not be advertised as evidence of GPU hardware.
        self.assertIsNone(dynamics.metadata["gpu_initialization_succeeded"])

    def test_cpu_does_not_claim_configurable_gpu_precision(self):
        dynamics, observed = self._initialize_fake(platform="CPU")
        self.assertEqual(observed["properties_at_context_creation"], {})
        self.assertEqual(observed["kwargs"]["precision"], "double")
        self.assertIsNone(dynamics.metadata["precision"])
        self.assertFalse(dynamics.metadata["precision_configurable"])

    def test_no_platform_or_precision_fallback(self):
        with self.assertRaisesRegex(BackendUnavailableError, "No fallback engine"):
            self._initialize_fake(platform="CUDA", actual_platform="CPU")
        with self.assertRaisesRegex(BackendUnavailableError, "requested precision"):
            self._initialize_fake(platform="CUDA", precision="mixed", reported_precision="single")

    def test_dependency_status_does_not_claim_gpu_initialization(self):
        status = backend_status()
        self.assertFalse(status["initialization_tested"])
        self.assertIsNone(status["gpu_initialization_succeeded"])
        self.assertEqual(status["registered_platforms"], status["platforms"])


class OptionalDeviceSmokeTests(unittest.TestCase):
    def _check_device(self, platform):
        if platform not in backend_status()["registered_platforms"]:
            self.skipTest(f"{platform} plugin is not registered; no device test was run")
        sim = particle_model(2400)
        try:
            dynamics = PolychromDynamics(sim, seed=137, platform=platform, precision="mixed")
        except BackendUnavailableError as exc:
            self.skipTest(f"{platform} device could not initialize: {exc}")
        dynamics(sim)
        self.assertEqual(dynamics.metadata["platform"], platform)
        self.assertEqual(dynamics.metadata["precision"], "mixed")
        self.assertTrue(np.isfinite(sim.positions).all())
        observed_msd = np.sum(sim.positions**2, axis=1).mean()
        self.assertLess(abs(observed_msd/(6*0.2*sim.dt) - 1), 0.08)
        frozen = sim.positions.copy()
        sim.diffusion[:] = 0
        dynamics(sim)
        np.testing.assert_allclose(sim.positions, frozen, rtol=1e-6, atol=1e-7)

    def test_cuda_brownian_smoke(self):
        self._check_device("CUDA")

    def test_opencl_brownian_smoke(self):
        self._check_device("OpenCL")


if __name__ == "__main__":
    unittest.main()
