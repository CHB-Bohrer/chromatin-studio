"""Optional, real polychrom/OpenMM Brownian dynamics backend.

Chemistry and mechanical force evaluation are shared with the NumPy engine.
Forces are uploaded in a single per-DOF array before *each single
Euler--Maruyama step*, never frozen across a reaction block. OpenMM performs
the stochastic integration. CPU model evaluation remains a bottleneck; this
is not a native OpenMM force implementation.

Reduced units: kBT=1, one length unit is represented internally by one nm and
one time unit by one ps. These are numerical conventions, NOT a physical
calibration. The force is supplied through an integrator variable, so the
OpenMM Context energy must not be interpreted as the model potential.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import re
from typing import Any

import numpy as np


class BackendUnavailableError(RuntimeError):
    """The explicitly requested optional dynamics backend could not start."""


def validate_backend_options(*, platform: str = "CPU", device_index=None,
                             precision: str = "double") -> dict[str, Any]:
    """Normalize an explicit single-device request without querying hardware.

    An omitted device leaves OpenMM's device selection untouched, including a
    scheduler's CUDA_VISIBLE_DEVICES mapping. An explicit index is a *logical*
    index in that visible set, never a physical host GPU identifier.
    CPU and Reference have fixed arithmetic; ``precision`` configures GPU
    platforms only and metadata distinguishes requested from actual precision.
    """
    platforms = {name.lower(): name for name in ("CPU", "Reference", "CUDA", "OpenCL")}
    if not isinstance(platform, str) or platform.lower() not in platforms:
        raise ValueError("platform must be CPU, Reference, CUDA, or OpenCL.")
    platform = platforms[platform.lower()]
    if not isinstance(precision, str) or precision not in ("single", "mixed", "double"):
        raise ValueError("precision must be single, mixed, or double.")
    if device_index is not None:
        if (isinstance(device_index, bool)
                or not isinstance(device_index, (int, str))
                or not re.fullmatch(r"[0-9]+", str(device_index))):
            raise ValueError("device_index must be one nonnegative logical device index or null.")
        if platform not in ("CUDA", "OpenCL"):
            raise ValueError("device_index is only supported for CUDA or OpenCL.")
        device_index = str(int(device_index))
    return {"platform": platform, "device_index": device_index, "precision": precision}


def _capabilities() -> dict[str, Any]:
    return {
        "integration": "OpenMM Brownian integration",
        "mechanical_forces": "CPU",
        "reactions": "CPU",
        "rna_diffusion": "CPU",
        "host_device_transfer": "positions, forces, and mobilities every substep",
        "native_gpu_forces": False,
        "multi_gpu_per_run": False,
        "physical_time_calibrated": False,
    }


def _dependencies():
    try:
        mm = importlib.import_module("openmm")
        unit = importlib.import_module("openmm.unit")
        polymer = importlib.import_module("polychrom.simulation")
    except Exception as exc:
        raise BackendUnavailableError(
            "The polychrom backend requires working OpenMM and polychrom "
            "installations in the Python environment running this app. "
            "Run the optional-backend setup and relaunch with that environment. "
            f"Import failed: {exc}"
        ) from exc
    return mm, unit, polymer


def backend_status() -> dict[str, Any]:
    """Inspect registration only; a plugin is not proof of a working device.

    Do not initialize a device merely to render the setup UI. Initialization
    and actual property checks happen when an explicit simulation is created.
    """
    try:
        mm, _, _ = _dependencies()
        try:
            polychrom_version = importlib.metadata.version("polychrom")
        except importlib.metadata.PackageNotFoundError:
            polychrom_version = "source checkout"
        platforms = [mm.Platform.getPlatform(i).getName()
                     for i in range(mm.Platform.getNumPlatforms())]
        return {
            "id": "polychrom", "available": bool(platforms),
            "detail": f"OpenMM {mm.__version__}; polychrom {polychrom_version}; platforms: "
                      + ", ".join(platforms)
                      + ". Registered platforms only; device initialization has not been tested. "
                      "Real Brownian integration; CPU forces and reactions, with transfers every substep.",
            "platforms": platforms,
            "registered_platforms": platforms,
            "registered_gpu_platforms": [name for name in platforms if name in ("CUDA", "OpenCL")],
            "initialization_tested": False,
            "gpu_initialization_succeeded": None,
            "availability_scope": "dependency imports and platform registration",
            "plugin_load_failures": list(mm.Platform.getPluginLoadFailures()),
            "capabilities": _capabilities(),
            "versions": {"openmm": mm.__version__, "polychrom": polychrom_version},
        }
    except BackendUnavailableError as exc:
        return {"id": "polychrom", "available": False, "detail": str(exc),
                "platforms": [], "registered_platforms": [], "registered_gpu_platforms": [],
                "initialization_tested": False, "gpu_initialization_succeeded": None,
                "capabilities": _capabilities()}


class PolychromDynamics:
    """Callable dynamics adapter for engine.Simulation's one-step hook.

    ``sim`` must expose positions (N,3), diffusion (N,), dt and get_forces().
    Friction changes are captured by uploading diffusion on every call.
    No reaction or time bookkeeping occurs here.
    """

    def __init__(self, sim, *, seed: int = 1, platform: str = "CPU",
                 device_index=None, precision: str = "double"):
        requested = validate_backend_options(platform=platform, device_index=device_index,
                                             precision=precision)
        platform, device_index, precision = (requested[key]
                                             for key in ("platform", "device_index", "precision"))
        mm, unit, polymer = _dependencies()
        self.mm, self.unit = mm, unit
        self.n = len(sim.positions)
        if self.n == 0:
            raise ValueError("Polychrom requires at least one particle.")
        self._validate_arrays(sim)

        integrator = mm.CustomIntegrator(float(sim.dt))
        integrator.addGlobalVariable("thermal_energy", 1.0)
        integrator.addPerDofVariable("mobility", 0.0)
        integrator.addPerDofVariable("model_force", 0.0)
        integrator.addUpdateContextState()
        integrator.addComputePerDof(
            "x", "x + dt*mobility*model_force + sqrt(2*thermal_energy*mobility*dt)*gaussian"
        )
        integrator.setRandomNumberSeed(int(seed) % (2**31 - 2) + 1)
        self.integrator = integrator

        # kBT=1 kJ/mol avoids an implicit factor of 2.494 at 300 K.
        gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(
            unit.kilojoule_per_mole / unit.kelvin
        )
        try:
            self.polymer = polymer.Simulation(
                N=self.n, platform=platform, integrator=(integrator, "brownian"),
                collision_rate=1.0, timestep=float(sim.dt) * 1000,
                temperature=1.0 / gas_constant, reporters=[], verbose=False,
                precision=precision, GPU="default" if device_index is None else device_index,
                save_decimals=False,
            )
            if platform in ("CUDA", "OpenCL"):
                # polychrom 0.1.1 skips both DeviceIndex and Precision when
                # GPU="default". Preserve default device selection but set
                # precision independently before it constructs the Context.
                self.polymer.properties["Precision"] = precision
            # Do not add polychrom's default coordinate jitter: the initial
            # conformation must be exactly the reference engine's conformation.
            self.polymer.set_data(sim.positions, random_offset=0, report=False)
            # A zero-energy placeholder initializes the polychrom Context.
            # Actual model forces are uploaded directly to the integrator,
            # avoiding thousands of individual parameter updates per step.
            force = mm.CustomExternalForce("0")
            force.name = "reference_force"
            for i in range(self.n):
                force.addParticle(i, [])
            self.force = force
            self.polymer.add_force(force)
            # Supported upstream path for custom-integrator simulations that
            # step directly, also used by polychrom's loop-extrusion examples.
            self.polymer._apply_forces()
            self.context = self.polymer.context
            actual_platform = self.context.getPlatform()
            self.platform = actual_platform.getName()
            if self.platform != platform:
                raise RuntimeError(f"requested {platform}, initialized {self.platform}")
            self.context_properties = {
                name: actual_platform.getPropertyValue(self.context, name)
                for name in actual_platform.getPropertyNames()
            }
            if platform in ("CUDA", "OpenCL"):
                if self.context_properties.get("Precision") != precision:
                    raise RuntimeError("the initialized device did not honor the requested precision")
                if (device_index is not None
                        and self.context_properties.get("DeviceIndex") != device_index):
                    raise RuntimeError("the initialized device did not honor the requested logical index")
            self.metadata = {
                "requested_platform": platform,
                "platform": self.platform,
                "requested_device_index": device_index,
                "device_index": self.context_properties.get("DeviceIndex"),
                "requested_precision": precision,
                "precision": self.context_properties.get("Precision"),
                "precision_configurable": platform in ("CUDA", "OpenCL"),
                "context_properties": dict(self.context_properties),
                "initialization_succeeded": True,
                "gpu_initialization_succeeded": True if platform == "CUDA" else None,
                "device_type": ("CUDA GPU" if platform == "CUDA" else
                                "OpenCL device (GPU/CPU type not verified)" if platform == "OpenCL"
                                else "CPU"),
                "capabilities": _capabilities(),
            }
        except Exception as exc:
            raise BackendUnavailableError(
                f"polychrom/OpenMM could not initialize platform {platform}: {exc}. "
                "No fallback engine was selected."
            ) from exc

    def _validate_arrays(self, sim):
        positions = np.asarray(sim.positions, dtype=float)
        diffusion = np.asarray(sim.diffusion, dtype=float)
        if positions.shape != (self.n, 3) or diffusion.shape != (self.n,):
            raise ValueError("Particle count or shape changed after backend initialization.")
        if not np.isfinite(positions).all():
            raise ValueError("Particle positions must be finite.")
        if not np.isfinite(diffusion).all() or np.any(diffusion < 0):
            raise ValueError("Particle diffusion must be finite and nonnegative.")
        if not np.isfinite(sim.dt) or sim.dt <= 0:
            raise ValueError("Integration timestep must be finite and positive.")

    def __call__(self, sim):
        self._validate_arrays(sim)
        substeps = int(sim.integration_substeps()) if hasattr(sim, "integration_substeps") else 1
        if substeps < 1:
            raise ValueError("Integration substeps must be positive.")
        sim.last_substeps = substeps
        for _ in range(substeps):
            self._advance_one(sim, float(sim.dt) / substeps)

    def _advance_one(self, sim, dt):
        force_values = np.asarray(sim.get_forces(), dtype=float)
        if force_values.shape != (self.n, 3) or not np.isfinite(force_values).all():
            raise ValueError("Mechanical forces must be a finite (N, 3) array.")
        # Synchronize positions as well as mobilities: shared reactions or a
        # future boundary rule may modify reference state between calls.
        self.context.setPositions(np.asarray(sim.positions) * self.unit.nanometer)
        self.integrator.setPerDofVariableByName("model_force", force_values)
        self.integrator.setStepSize(dt)
        self.integrator.setPerDofVariableByName(
            "mobility", np.repeat(np.asarray(sim.diffusion)[:, None], 3, axis=1)
        )
        self.integrator.step(1)
        positions = self.context.getState(getPositions=True).getPositions(asNumpy=True)
        updated = np.asarray(positions.value_in_unit(self.unit.nanometer), dtype=float)
        if not np.isfinite(updated).all():
            raise FloatingPointError("OpenMM integration produced nonfinite positions; reduce dt.")
        if hasattr(sim, "validate_positions"):
            sim.validate_positions(updated)
        sim.positions[:] = updated
        self.polymer.data = positions


def create_polychrom_simulation(config, *, platform: str = "CPU", device_index=None,
                               precision: str = "double"):
    """Create a reference chemistry model with genuine OpenMM integration."""
    from .engine import Simulation

    # Reference initialization defines topology and chemistry. Attach the
    # dynamics only after its particle arrays exist.
    sim = Simulation(config)
    dynamics = PolychromDynamics(sim, seed=int(config.get("seed", 1)), platform=platform,
                                device_index=device_index, precision=precision)
    sim.dynamics = dynamics
    sim.engine = "polychrom"
    sim.backend_metadata = dynamics.metadata
    sim.warnings.append(
        f"polychrom/OpenMM ({dynamics.platform}) integrates Brownian motion. "
        "Mechanical forces and reactions run on CPU; forces upload every timestep. "
        "GPU-scale force throughput is not established. Selecting OpenMM does not set a physical clock; use the separate calibration profile."
    )
    return sim


class PolychromSimulation:
    """Server-facing simulation API backed by the shared chemistry engine."""

    def __init__(self, config, *, platform: str = "CPU", device_index=None,
                 precision: str = "double"):
        self._simulation = create_polychrom_simulation(
            config, platform=platform, device_index=device_index, precision=precision,
        )

    def step(self, steps: int):
        return self._simulation.step(steps)

    def snapshot(self):
        return self._simulation.snapshot()

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_simulation"), name)
