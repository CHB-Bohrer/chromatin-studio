"""Meaningful optional-backend checks; skipped when dependencies are absent."""

import unittest
from types import SimpleNamespace

import numpy as np

from chromatin_studio.openmm_backend import PolychromDynamics, backend_status


@unittest.skipUnless(backend_status()["available"], "optional polychrom/OpenMM unavailable")
class OpenMMTests(unittest.TestCase):
    def test_diffusion_and_zero_mobility(self):
        # Ensemble statistics test both the noise scale and per-particle D.
        n = 2400
        diffusion = np.repeat([0.0, 0.05, 0.5], n)
        sim = SimpleNamespace(positions=np.zeros((3*n, 3)), diffusion=diffusion, dt=0.2)
        sim.get_forces = lambda: np.zeros_like(sim.positions)
        dynamics = PolychromDynamics(sim, seed=741, platform="Reference")
        dynamics(sim)
        np.testing.assert_array_equal(sim.positions[:n], 0)
        msd = np.sum(sim.positions**2, axis=1)
        for index, coefficient in enumerate((0.05, 0.5), start=1):
            observed = msd[index*n:(index+1)*n].mean()
            expected = 6*coefficient*sim.dt
            self.assertLess(abs(observed/expected - 1), 0.07)
        # Mobility may change after a reaction. Zero mobility must take effect
        # on the very next integrator step, without rebuilding the context.
        frozen = sim.positions.copy()
        sim.diffusion[:] = 0
        dynamics(sim)
        np.testing.assert_array_equal(sim.positions, frozen)

    def test_force_direction_and_thermal_scale(self):
        n = 5000
        sim = SimpleNamespace(positions=np.zeros((n, 3)), diffusion=np.full(n, 0.4), dt=0.5)
        sim.get_forces = lambda: np.tile([2.0, 0.0, 0.0], (n, 1))
        dynamics = PolychromDynamics(sim, seed=812, platform="Reference")
        dynamics(sim)
        # Mean dx = D * F / kBT * dt = 0.4. This catches the common
        # accidental 300-K energy factor and reversed force sign.
        self.assertAlmostEqual(sim.positions[:, 0].mean(), 0.4, delta=0.035)
        self.assertAlmostEqual(sim.positions[:, 1].mean(), 0.0, delta=0.035)
        self.assertAlmostEqual(sim.positions[:, 0].var(), 0.4, delta=0.035)
        # Reaction-dependent forces must replace the prior force immediately.
        before = sim.positions.copy()
        sim.get_forces = lambda: np.tile([-2.0, 0.0, 0.0], (n, 1))
        dynamics(sim)
        self.assertAlmostEqual((sim.positions - before)[:, 0].mean(), -0.4, delta=0.035)

    def test_finite_full_system(self):
        # Loaded lazily because schema/engine are owned by other modules.
        from chromatin_studio.openmm_backend import create_polychrom_simulation
        from chromatin_studio.schema import default_config

        sim = create_polychrom_simulation(default_config(), platform="Reference")
        sim.step(12)
        state = sim.snapshot()
        self.assertEqual(state["engine"], "polychrom")
        self.assertGreater(state["time"], 0)
        self.assertTrue(np.isfinite(sim.positions).all())
        self.assertEqual(state["step"], 12)


if __name__ == "__main__":
    unittest.main()
