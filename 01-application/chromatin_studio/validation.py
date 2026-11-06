"""Portable analytical acceptance checks for CPU and requested OpenMM devices.

These are implementation checks, not biological validation or proof of
equilibration. Statistical tests compare independent ensembles, not trajectories.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform as host_platform
import time
from types import SimpleNamespace

import numpy as np

from . import __version__
from .engine import Simulation
from .schema import default_config, preset_config


def _check(name, observed, expected, tolerance, **details):
    error = abs(float(observed) - float(expected))
    return {"name": name, "passed": bool(math.isfinite(error) and error <= tolerance),
            "observed": float(observed) if math.isfinite(observed) else None,
            "expected": float(expected) if math.isfinite(expected) else None,
            "absolute_tolerance": float(tolerance),
            "absolute_error": error if math.isfinite(error) else None, **details}


def _fraction(name, count, total, expected):
    # Six binomial standard errors; the floor handles tiny expected counts.
    tolerance = max(6 * math.sqrt(expected * (1 - expected) / total), 2 / total)
    return _check(name, count / total, expected, tolerance, samples=total,
                  criterion="six binomial standard errors (minimum two counts)")


def _minimal_config(beads=16):
    config = default_config()
    config.update(seed=471, marks=[], proteins=[], interactions=[], bodies=[], genes=[])
    config["cohesin"] = {"enabled": False, "count": 0, "step_rate": 0, "off_rate": 0, "spring_k": 0}
    config["simulation"].update(dt=0.002, confinement_radius=200)
    config["chromatin"].update(beads=beads, bp_per_bead=200, initial_marks=[], barriers=[],
                               bond_length=1, bond_k=30, persistence_length=2, repulsion=0)
    return config


def _integrator_fixture(positions, diffusion, dt, seed, force):
    sim = SimpleNamespace(positions=np.asarray(positions, dtype=float),
                          diffusion=np.asarray(diffusion, dtype=float), dt=dt,
                          rng=np.random.default_rng(seed), radius=1e6)
    sim.get_forces = lambda: force(sim.positions)
    sim.integration_substeps = lambda: 1
    sim.validate_positions = lambda values: Simulation.validate_positions(sim, values)
    return sim


def _mechanics(platform, precision, device_index):
    checks, measurements, metadata = [], {}, {"actual_platform": "numpy", "integration": "CPU"}

    def adapter(sim, seed):
        nonlocal metadata
        if platform == "numpy":
            return lambda s: Simulation._brownian_move(s)
        from .openmm_backend import PolychromDynamics
        dynamics = PolychromDynamics(sim, seed=seed, platform=platform,
                                     precision=precision, device_index=device_index)
        metadata = dynamics.metadata
        return dynamics

    n, dt = 2048, .2
    sim = _integrator_fixture(np.zeros((n * 3, 3)), np.repeat([0, .05, .5], n), dt, 741,
                              lambda p: np.zeros_like(p))
    advance = adapter(sim, 741)
    advance(sim)
    checks.append(_check("zero diffusion stays immobile", np.max(np.abs(sim.positions[:n])), 0, 1e-7))
    for group, coefficient in enumerate((.05, .5), start=1):
        values = sim.positions[group*n:(group+1)*n]
        observed = float(np.mean(np.sum(values * values, axis=1)))
        expected = 6 * coefficient * dt
        checks.append(_check(f"free diffusion D={coefficient}: MSD=6Dt", observed, expected,
                             6 * expected * math.sqrt(2 / (3*n)), samples=n))
        measurements[f"diffusion_{coefficient}"] = {"value": observed,
                                                    "standard_error": expected * math.sqrt(2 / (3*n))}
    before = sim.positions.copy()
    sim.diffusion[:] = 0
    advance(sim)
    checks.append(_check("changed mobility applies next step", np.max(np.abs(sim.positions-before)), 0, 1e-7))
    del advance

    # Sudden force reversal tests sign, thermal units and CPU-to-device updates.
    n, D, dt = 4096, .4, .25
    sim = _integrator_fixture(np.zeros((n, 3)), np.full(n, D), dt, 812,
                              lambda p: np.tile([2., 0., 0.], (len(p), 1)))
    advance = adapter(sim, 812)
    for sign in (1, -1):
        sim.get_forces = lambda sign=sign: np.tile([2. * sign, 0., 0.], (n, 1))
        before = sim.positions.copy()
        advance(sim)
        increments = sim.positions - before
        checks.append(_check(f"force response {sign:+d}: mean dx=D F dt", increments[:, 0].mean(),
                             sign * D * 2 * dt, 6 * math.sqrt(2*D*dt/n), samples=n))
        checks.append(_check(f"force response {sign:+d}: thermal variance", increments[:, 0].var(),
                             2*D*dt, 6 * 2*D*dt*math.sqrt(2/(n-1)), samples=n))
    del advance

    # Known Ornstein-Uhlenbeck transition at two timesteps. Exact discrete
    # Euler moments isolate implementation errors from expected step-size bias.
    for dt in (.04, .02):
        n, D, k, duration = 4096, .4, 2., .4
        steps = round(duration / dt)
        sim = _integrator_fixture(np.ones((n, 3)), np.full(n, D), dt, 950,
                                  lambda p: -k*p)
        advance = adapter(sim, 950)
        for _ in range(steps):
            advance(sim)
        a = 1 - D*k*dt
        mean = a**steps
        variance = 2*D*dt * (1-a**(2*steps))/(1-a*a)
        samples = n*3
        checks.append(_check(f"harmonic relaxation dt={dt}: mean", sim.positions.mean(), mean,
                             6*math.sqrt(variance/samples), samples=samples,
                             continuous_time_target=math.exp(-D*k*duration)))
        checks.append(_check(f"harmonic relaxation dt={dt}: variance", sim.positions.var(), variance,
                             6*variance*math.sqrt(2/(samples-1)), samples=samples,
                             continuous_time_target=(1-math.exp(-2*D*k*duration))/k))
        measurements[f"harmonic_variance_{dt}"] = {
            "value": float(sim.positions.var()), "standard_error": variance*math.sqrt(2/(samples-1))}
        del advance
    return {"platform": platform, "metadata": metadata, "checks": checks, "measurements": measurements}


def _force_gradient():
    config = _minimal_config()
    config["chromatin"].update(repulsion=3, nucleosomes={"enabled": True, "initial_occupancy": 1,
                                "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": .4})
    config["marks"] = [{"id": "M", "decay_rate": 0, "spread_rate": 0}]
    config["chromatin"]["initial_marks"] = [{"start": 0, "end": 15, "mark": "M"}]
    sim = Simulation(config)
    sim.nucleosome_occupied[3:6] = False
    sim._clear_vacant_marks()
    # The initializer avoids overlaps. Bring a non-neighbor pair into the
    # smooth interior of the soft core so repulsion is actually exercised.
    sim.positions[-1] = sim.positions[0] + [.7, .15, .1]
    contact_pairs, _ = sim._repulsion_pairs()
    if not len(contact_pairs):
        raise AssertionError("The force-gradient fixture has no active repulsion contacts.")
    force = sim.get_forces()
    numerical, h = np.zeros_like(force), 1e-6
    for i in range(len(sim.positions)):
        for d in range(3):
            sim.positions[i, d] += h
            plus = sim.potential_energy()
            sim.positions[i, d] -= 2*h
            minus = sim.potential_energy()
            sim.positions[i, d] += h
            numerical[i, d] = -(plus-minus)/(2*h)
    return [_check("shared chromatin force equals negative energy gradient", np.max(np.abs(force-numerical)), 0, 2e-5,
                   active_repulsion_pairs=len(contact_pairs)),
            _check("internal forces sum to zero", np.max(np.abs(force.sum(axis=0))), 0, 1e-9)]


def _occupancy_checks():
    config = _minimal_config(256)
    config["marks"] = [{"id": "M", "decay_rate": 0, "spread_rate": 0}]
    config["chromatin"]["nucleosomes"] = {"enabled": True, "initial_occupancy": 1,
        "eviction_rate": .7, "reassembly_rate": .3}
    sim = Simulation(config)
    vacant = cleared = rebuilt = 0
    for _ in range(128):
        sim.nucleosome_occupied[:] = True
        sim.marks["M"][:] = True
        sim._advance_nucleosomes(.4)
        vacant += int((~sim.nucleosome_occupied).sum())
        cleared += int((~sim.marks["M"]).sum())
        sim.nucleosome_occupied[:] = False
        sim._advance_nucleosomes(.4)
        rebuilt += int(sim.nucleosome_occupied.sum())
    total = 128*256
    return [_fraction("occupancy two-state eviction endpoint", vacant, total, .7*(1-math.exp(-.4))),
            _fraction("occupancy retains memory of any mark-erasing eviction", cleared, total, 1-math.exp(-.7*.4)),
            _fraction("occupancy two-state reassembly endpoint", rebuilt, total, .3*(1-math.exp(-.4)))]


def _binding_checks():
    config = _minimal_config(64)
    protein = copy.deepcopy(default_config()["proteins"][0])
    protein.update(id="probe", count=64, representation="bead", valence=1, on_rate=.8,
                   off_rate=.7, chromatin_affinity={"unmodified": 1},
                   nucleosome_affinity={"occupied": .1, "vacant": 1})
    config["proteins"] = [protein]
    config["chromatin"]["nucleosomes"] = {"enabled": True, "initial_occupancy": 1}
    sim = Simulation(config)
    sim.positions[:64] = np.column_stack((np.arange(64)*3, np.zeros(64), np.zeros(64)))
    sim.positions[64:] = sim.positions[:64] + [0, .5, 0]
    checks = []
    for occupied, factor in ((True, .1), (False, 1)):
        sim.nucleosome_occupied[:] = occupied
        captured = 0
        for _ in range(128):
            sim.bindings = []
            sim._advance_binding(.2)
            captured += len(sim.bindings)
        checks.append(_fraction(f"isolated TF encounter, occupied={occupied}", captured, 8192,
                                1-math.exp(-.8*factor*.2)))
    sim.species["probe"]["on_rate"] = 0
    survived = 0
    for _ in range(128):
        sim.bindings = [{"a": i, "b": i+64, "kind": "chromatin", "off_rate": .7,
                         "k": 1, "rest": .5} for i in range(64)]
        sim._advance_binding(.5)
        survived += len(sim.bindings)
    checks.append(_fraction("bound-state survival exp(-koff*t), rebinding disabled", survived, 8192, math.exp(-.7*.5)))
    return checks


def _rna_checks():
    config = _minimal_config()
    gene = copy.deepcopy(default_config()["genes"][0])
    gene.update(id="probe", start=2, end=5, on_rate=0, off_rate=0, initiation_rate=0,
                elongation_rate=2, processing_rate=0, rna_decay_rate=0,
                activator_species=None, required_mark=None, repressive_mark=None,
                pause_rate=0, persistence_multiplier=1)
    config["genes"] = [gene]
    config["transcription"] = {"enabled": True, "max_transcripts": 5000,
                               "rna_diffusion": 0, "mature_rna_diffusion": 0}
    sim = Simulation(config)
    explicit = sim.explicit_transcription
    rates = {"probe": {"effective_on_rate": 0, "effective_initiation_rate": 0}}
    record = explicit._new_nascent(sim.genes[0], 0)
    explicit.advance(.75, rates)
    checks = [_check("nascent RNA length = elongation speed * elapsed time", record["length_nt"], 300, 1e-9)]
    sim.time = .75
    explicit.advance(1.5, rates)
    checks += [_check("release time = gene length / speed", record["released_time"], 2, 1e-9),
               _check("completed RNA includes last genomic bin", record["length_nt"], 800, 1e-9)]
    sim = Simulation(config)
    explicit = sim.explicit_transcription
    for _ in range(4000):
        r = explicit._new_nascent(sim.genes[0], 0)
        r.update(state="unprocessed", bead_position=None, released_time=0,
                 length_nt=800, position=np.zeros(3))
        sim.genes[0]["completed"] += 1
    sim.genes[0]["config"].update(processing_rate=1.2, rna_decay_rate=.7)
    explicit.advance(.8, rates)
    totals = explicit.snapshot()["totals"]
    survival = math.exp(-.7*.8)
    unprocessed = math.exp(-1.9*.8)
    checks += [_fraction("RNA competing processing/decay: unprocessed", totals["unprocessed"], 4000, unprocessed),
               _fraction("RNA competing processing/decay: processed", totals["mature"], 4000, survival-unprocessed),
               _check("RNA molecule accounting", totals["live"]+totals["degraded"], 4000, 0)]
    return checks


def run_validation(platform="numpy", precision="mixed", device_index=None):
    """Run device acceptance and host kinetics; a missing requested GPU fails."""
    started = time.perf_counter()
    report = {"schema_version": 1, "application_version": __version__,
              "requested": {"platform": platform, "precision": precision, "device_index": device_index},
              "host": {"python": host_platform.python_version(), "system": host_platform.system()},
              "checks": [], "mechanics": [], "errors": [],
              "scope": "Analytical component acceptance; not biological validation, GPU performance certification, or full-model timestep convergence."}
    report["versions"] = {}
    for package in ("numpy", "scipy", "openmm", "polychrom"):
        try:
            report["versions"][package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            report["versions"][package] = None
    report["source_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(Path(__file__).parent.glob("*.py"))}
    try:
        for target in (["numpy"] if platform == "numpy" else ["numpy", platform]):
            result = _mechanics(target, precision, device_index if target != "numpy" else None)
            report["mechanics"].append(result)
            report["checks"].extend({**c, "execution": target} for c in result["checks"])
        if len(report["mechanics"]) == 2:
            baseline, candidate = report["mechanics"]
            for name, reference in baseline["measurements"].items():
                measured = candidate["measurements"][name]
                tolerance = 6*math.hypot(reference["standard_error"], measured["standard_error"])
                report["checks"].append(_check(f"CPU/device ensemble agreement: {name}", measured["value"],
                    reference["value"], tolerance, execution="cross-platform", criterion="six combined standard errors; independent RNG streams"))
        for probe in (_force_gradient, _occupancy_checks, _binding_checks, _rna_checks):
            report["checks"].extend({**c, "execution": "shared CPU model"} for c in probe())
        from .batch import inspect_state, _create_simulation
        coupled = _create_simulation(preset_config("nucleosome_transcription"),
            {"platform": platform, "precision": precision, "device_index": device_index})
        # Seed controlled nascent and released RNA, rather than letting a
        # short smoke test pass with empty RNA pools and vacuous accounting.
        explicit, gene = coupled.explicit_transcription, coupled.genes[0]
        explicit._new_nascent(gene, 0)
        released = explicit._new_nascent(gene, 0)
        released.update(state="unprocessed", bead_position=None, released_time=0,
            length_nt=(gene["config"]["end"]-gene["config"]["start"]+1)*coupled.chrom["bp_per_bead"],
            position=coupled.positions[gene["config"]["end"]].copy())
        gene["completed"] += 1
        explicit._synchronize()
        coupled.step(50)
        state = inspect_state(coupled)
        report["checks"].append({"name": "coupled 80 kb model: finite state, valence and RNA accounting",
            "passed": True, "execution": platform, "steps": 50,
            "rna_totals": state["transcription"]["totals"],
            "initial_seeded_rna": {"nascent": 1, "unprocessed": 1},
            "note": "Controlled seeded RNA; short integration smoke test, not equilibration or biological accuracy."})
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
    report["passed"] = not report["errors"] and bool(report["checks"]) and all(c["passed"] for c in report["checks"])
    report["elapsed_seconds"] = time.perf_counter()-started
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=["numpy", "CPU", "Reference", "CUDA", "OpenCL"], default="numpy")
    parser.add_argument("--precision", choices=["single", "mixed", "double"], default="mixed")
    parser.add_argument("--device-index", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists; choose a new validation report path.")
    report = run_validation(args.platform, args.precision, args.device_index)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
    passed = sum(c["passed"] for c in report["checks"])
    print(f"{'PASS' if report['passed'] else 'FAIL'}: {passed}/{len(report['checks'])} analytical checks. Report: {args.output}")
    for error in report["errors"]:
        print(error)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
