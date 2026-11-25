"""Native OpenMM forces for bare chromatin polymers (GPU-resident dynamics).

The existing polychrom adapter integrates on the device but evaluates every
force on the host, so each substep pays a host/device round trip. For models
with *no chemistry* the whole step can instead stay on the device. This module
expresses the reference engine's mechanics exactly as OpenMM forces:

    backbone    0.5*k*(r-L)^2                          HarmonicBondForce
    bending     kappa_i*(1-cos(theta_bond)),           CustomAngleForce
                kappa_i = persistence_i/L,
                theta_bond = pi - theta_OpenMM  =>  kappa_i*(1+cos(theta))
    repulsion   eps/2*(1-r^2/sigma^2)^2, r<sigma,      CustomNonbondedForce
                sigma = r_i+r_j, bonded neighbours excluded
    confinement 0.5*k_c*(|x|-R)^2, |x|>R               CustomExternalForce

and integrates them with the reference Euler--Maruyama update
x += D*F*dt + sqrt(2*D*dt)*N(0,1) (kBT=1, reduced units: 1 length = 1 nm,
1 time = 1 ps, 1 energy = 1 kJ/mol inside OpenMM; numerical conventions only).

Scope: only models whose reaction step is a no-op (no proteins, marks,
state interactions, bodies, cohesin, genes, regulatory regions, transcription
or enabled nucleosomes) are accepted; anything else raises instead of silently
dropping physics. The initial conformation is the reference engine's for the
same seed. Trajectories are statistically equivalent to the NumPy engine, not
identical: the noise comes from OpenMM's random stream. Not a physical clock.
"""

from __future__ import annotations

import math

import numpy as np

from .openmm_backend import BackendUnavailableError, validate_backend_options


def unsupported_features(config) -> list[str]:
    """Reasons the model needs host chemistry; empty when native is exact."""
    reasons = []
    chrom = config.get("chromatin", {})
    proteins = [p for p in config.get("proteins", []) if int(p.get("count", 0))]
    # One-bead proteins whose only activity is attraction-mode pair wells are
    # pure mechanics and run on the device; anything reactive needs the host.
    if any(p.get("representation", "bead") != "bead" for p in proteins):
        reasons.append("IDR-chain proteins")
    if any(float(p.get("on_rate", 0)) > 0 and any(float(v) > 0 for v in p.get("chromatin_affinity", {}).values())
           for p in proteins):
        reasons.append("protein-chromatin binding")
    if any(p.get("writer_mark") or p.get("eraser_mark") for p in proteins):
        reasons.append("mark writers or erasers")
    if any(r.get("mode", "bonds") != "attraction" and float(r.get("on_rate", 0)) > 0
           for r in config.get("interactions", [])):
        reasons.append("bond-mode protein interactions")
    if config.get("marks"):
        reasons.append("histone marks")
    if chrom.get("initial_marks"):
        reasons.append("initial marks")
    if chrom.get("state_interactions"):
        reasons.append("state-selective chromatin attraction")
    if any(b.get("enabled", True) for b in config.get("bodies", [])):
        reasons.append("nuclear bodies")
    if config.get("cohesin", {}).get("enabled"):
        reasons.append("cohesin extrusion")
    if config.get("genes"):
        reasons.append("genes")
    if config.get("regulatory_regions"):
        reasons.append("regulatory regions")
    if config.get("transcription", {}).get("enabled"):
        reasons.append("explicit transcription")
    if chrom.get("nucleosomes", {}).get("enabled"):
        reasons.append("nucleosome occupancy dynamics")
    if config.get("nucleus", {}).get("enabled"):
        reasons.append("nuclear envelope and tethers")
    return reasons


def _openmm():
    try:
        import openmm as mm
        import openmm.unit as unit
    except Exception as exc:  # pragma: no cover - depends on environment
        raise BackendUnavailableError(f"The native backend requires OpenMM: {exc}") from exc
    return mm, unit


def build_system(sim):
    """Return (System, forces-by-name) reproducing ``sim``'s mechanics."""
    mm, _ = _openmm()
    n = len(sim.positions)
    system = mm.System()
    for _ in range(n):
        system.addParticle(1.0)  # mass is unused by the custom integrator
    chrom = sim.chrom
    nc = sim.n_chromatin
    L, k = float(sim.bond_length), float(chrom["bond_k"])
    bonds = [(i, i + 1) for i in range(nc - 1)]

    backbone = mm.HarmonicBondForce()
    for i, j in bonds:
        backbone.addBond(i, j, L, k)
    backbone.setForceGroup(0)
    system.addForce(backbone)

    forces = {"backbone": backbone}
    if nc >= 3:
        bending = mm.CustomAngleForce("kappa*(1+cos(theta))")
        bending.addPerAngleParameter("kappa")
        kappa = sim.local_persistence()[1:-1] / L
        for i in range(nc - 2):
            if kappa[i] != 0:
                bending.addAngle(i, i + 1, i + 2, [float(kappa[i])])
        bending.setForceGroup(1)
        system.addForce(bending)
        forces["bending"] = bending

    eps = float(chrom.get("repulsion", 0))
    if eps:
        repulsion = mm.CustomNonbondedForce(
            "step(sig-r)*0.5*eps*(1-(r/sig)^2)^2; sig=radius1+radius2")
        repulsion.addGlobalParameter("eps", eps)
        repulsion.addPerParticleParameter("radius")
        for radius in sim.radii:
            repulsion.addParticle([float(radius)])
        repulsion.createExclusionsFromBonds(bonds, 1)
        repulsion.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        repulsion.setCutoffDistance(2 * float(np.max(sim.radii)))
        repulsion.setForceGroup(2)
        system.addForce(repulsion)
        forces["repulsion"] = repulsion

    attraction = getattr(sim, "_attraction", None)
    if attraction is not None:
        # Same flat-bottomed well as the reference engine (see
        # Simulation._build_attraction), eps and rc^2 looked up per species
        # pair; non-participants use row 0 of an extra "none" species (depth 0).
        names = attraction["names"]
        size = len(names) + 1
        eps_table = np.zeros((size, size))
        rc2_table = np.ones((size, size))
        eps_table[1:, 1:] = attraction["eps"]
        rc2_table[1:, 1:] = attraction["rc2"]
        well = mm.CustomNonbondedForce(
            "-e*(1-u*u)^2 + select(e, step(sig-r)*(4*(x6*x6-x6)+1), 0); "
            "x6=(sig/max(r, 0.5*sig))^6/2; u=min(max((r-sig)/(sqrt(rc2)-sig), 0), 1); sig=radius1+radius2; "
            "e=eps_table(s1,s2); rc2=rc2_table(s1,s2)")
        well.addTabulatedFunction("eps_table", mm.Discrete2DFunction(size, size, eps_table.ravel(order="F").tolist()))
        well.addTabulatedFunction("rc2_table", mm.Discrete2DFunction(size, size, rc2_table.ravel(order="F").tolist()))
        well.addPerParticleParameter("s")
        well.addPerParticleParameter("radius")
        participating = np.zeros(n, dtype=bool)
        participating[attraction["particles"]] = True
        for i in range(n):
            well.addParticle([float(attraction["species"][i] + 1) if participating[i] else 0.0, float(sim.radii[i])])
        if eps:
            # CUDA requires identical exclusions on every nonbonded force; the
            # excluded backbone pairs are chromatin, which never attracts here.
            well.createExclusionsFromBonds(bonds, 1)
        well.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        well.setCutoffDistance(float(attraction["cutoff"]))
        well.setForceGroup(4)
        system.addForce(well)
        forces["attraction"] = well

    if sim.confinement_stiffness > 0:
        # The tiny offset keeps the radial derivative finite at the origin
        # (the reference chain starts with bead 0 exactly there); its effect
        # on |x| is below 1e-15 length units.
        wall = mm.CustomExternalForce(
            "step(rr-R)*0.5*kc*(rr-R)^2; rr=sqrt(x*x+y*y+z*z+1e-30)")
        wall.addGlobalParameter("R", float(sim.radius))
        wall.addGlobalParameter("kc", float(sim.confinement_stiffness))
        for i in range(n):
            wall.addParticle(i, [])
        wall.setForceGroup(3)
        system.addForce(wall)
        forces["confinement"] = wall
    return system, forces


class NativePolymerSimulation:
    """Bare-polymer simulation with forces and integration on one device.

    Exposes the attributes the batch tools use (positions, step_count, time,
    snapshot, potential_energy, get_forces ...) by delegating to a reference
    engine whose coordinates are refreshed from the device after each
    ``step`` call. Host-side force/energy evaluation is used only for
    inspection, never for dynamics.
    """

    def __init__(self, config, *, platform="CUDA", device_index=None, precision="mixed",
                 noise_seed=None):
        """``noise_seed`` (default: the model seed) seeds OpenMM's Gaussian
        stream only; the initial conformation always follows the model seed.
        Use a new noise seed when continuing a trajectory from saved positions,
        so the continuation does not replay the original noise."""
        from .engine import Simulation

        reasons = unsupported_features(config)
        if reasons:
            raise ValueError("The native OpenMM backend supports bare polymers and one-bead attraction-mode proteins only; this model has: "
                             + ", ".join(reasons) + ". Use the numpy or polychrom engine.")
        requested = validate_backend_options(platform=platform, device_index=device_index,
                                             precision=precision)
        mm, unit = _openmm()
        self._sim = sim = Simulation(config)
        self.substeps = int(sim.integration_substeps())  # static: no bonds change
        self._fine_steps, self._fine_substeps = 0, self.substeps
        if sim._attraction is not None:
            # Attraction contacts change during the run, so use a static bound
            # for the run: up to 12 simultaneous close contacts (thermal overlaps
            # into the WCA core are shallow), each with the shell curvature bound
            # plus the core curvature at contact (72/sigma^2). The initial
            # placement allows overlaps down to 0.65*sigma, so the first
            # FINE_START steps use the reference engine's contact-based count.
            a = sim._attraction
            bound = np.zeros(len(sim.positions))
            protein_radius = float(np.min(sim.radii[a["particles"]]))
            for index in a["particles"]:
                row = a["species"][index]
                sigma = sim.radii[index] + protein_radius
                width = np.maximum(np.sqrt(a["rc2"][row]) - sigma, 1e-6)
                pair = 4 * a["eps"][row] / width ** 2 + 1.54 * a["eps"][row] / (width * sigma)
                core = 72.0 / sigma ** 2
                bound[index] = 12 * float(np.max(np.where(a["eps"][row] > 0, pair + core, 0.0)))
            run_substeps = int(math.ceil(sim.dt * float(np.max(sim.diffusion * bound)) / 0.35))
            stiffness_floor = int(math.ceil(sim.dt * float(np.max(sim.diffusion * (2 * float(sim.chrom["bond_k"])
                                    + 2 * float(sim.chrom.get("repulsion", 0)) / (2 * sim.radii) ** 2
                                    + sim.confinement_stiffness))) / 0.35))
            self._fine_steps, self._fine_substeps = self.FINE_START, max(self.substeps, run_substeps)
            self.substeps = max(1, run_substeps, stiffness_floor)
        sim.last_substeps = self.substeps
        system, self.forces = build_system(sim)
        integrator = mm.CustomIntegrator(sim.dt / self.substeps)
        integrator.addGlobalVariable("thermal_energy", 1.0)
        integrator.addPerDofVariable("mobility", 0.0)
        integrator.addUpdateContextState()
        integrator.addComputePerDof(
            "x", "x + dt*mobility*f + sqrt(2*thermal_energy*mobility*dt)*gaussian")
        seed = int(config.get("seed", 1)) if noise_seed is None else int(noise_seed)
        integrator.setRandomNumberSeed(seed % (2**31 - 2) + 1)
        try:
            device = mm.Platform.getPlatformByName(requested["platform"])
            properties = {}
            if requested["platform"] in ("CUDA", "OpenCL"):
                properties["Precision"] = requested["precision"]
                if requested["device_index"] is not None:
                    properties["DeviceIndex"] = requested["device_index"]
            self.context = mm.Context(system, integrator, device, properties)
        except Exception as exc:
            raise BackendUnavailableError(
                f"OpenMM could not initialize platform {requested['platform']}: {exc}. "
                "No fallback engine was selected.") from exc
        self.integrator = integrator
        self.context.setPositions(sim.positions * unit.nanometer)
        integrator.setPerDofVariableByName(
            "mobility", [mm.Vec3(d, d, d) for d in sim.diffusion])
        actual = self.context.getPlatform()
        props = {name: actual.getPropertyValue(self.context, name)
                 for name in actual.getPropertyNames()}
        if actual.getName() != requested["platform"]:
            raise BackendUnavailableError(f"requested {requested['platform']}, got {actual.getName()}")
        sim.engine = "openmm-native"
        sim.backend_metadata = {
            "requested_platform": requested["platform"], "platform": actual.getName(),
            "requested_precision": requested["precision"], "precision": props.get("Precision"),
            "device_index": props.get("DeviceIndex"), "device_name": props.get("DeviceName"),
            "context_properties": props, "integration_substeps": self.substeps,
            "capabilities": {"integration": "OpenMM device", "mechanical_forces": "OpenMM device",
                             "reactions": "none (bare polymer only)",
                             "host_device_transfer": "positions only when step() returns",
                             "native_gpu_forces": True, "physical_time_calibrated": False},
        }
        sim.warnings.append(
            f"Native OpenMM ({actual.getName()}) computes the bare-polymer forces and integration on the device. "
            "Statistically equivalent to the NumPy engine, not trajectory-identical. Reduced units.")

    FINE_START = 500

    def step(self, steps: int):
        if isinstance(steps, bool) or not isinstance(steps, (int, np.integer)) or steps < 0:
            raise ValueError("steps must be a nonnegative integer")
        if steps == 0:
            return
        fine = min(int(steps), self._fine_steps)
        if fine:
            self.integrator.setStepSize(self._sim.dt / self._fine_substeps)
            self.integrator.step(fine * self._fine_substeps)
            self._fine_steps -= fine
            self.integrator.setStepSize(self._sim.dt / self.substeps)
        self.integrator.step((int(steps) - fine) * self.substeps)
        self._pull_positions()
        self._sim.time += float(steps) * self._sim.dt
        self._sim.step_count += int(steps)

    def _pull_positions(self):
        _, unit = _openmm()
        state = self.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        self._sim.validate_positions(positions)
        self._sim.positions[:] = positions

    def device_forces_and_energy(self):
        """Forces/energy evaluated by OpenMM at the current device state."""
        _, unit = _openmm()
        state = self.context.getState(getForces=True, getEnergy=True)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        return np.asarray(forces, dtype=float), float(energy)

    def set_positions(self, positions):
        _, unit = _openmm()
        positions = np.asarray(positions, dtype=float)
        self._sim.positions[:] = positions
        self.context.setPositions(positions * unit.nanometer)

    def snapshot(self):
        return self._sim.snapshot()

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_sim"), name)
