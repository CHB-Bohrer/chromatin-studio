"""Reproducible stochastic benchmarks, distinct from biological calibration.

The Gaussian force fixture calls the application's actual Brownian integrator.
It does not silently substitute zero-rest springs into the production model.
The eigenmode calculation is an independent prediction of the measured MSD.
See Yuan et al., Phys. Rev. E 109, 044502 (2024), Eqs. 1, 20, 21, 26.
No external implementation or experimental dataset is embedded here.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import inspect
import json
import math
from pathlib import Path
import platform
import time
from types import SimpleNamespace

import numpy as np

from . import __version__
from .engine import Simulation

PAPER_URL = "https://doi.org/10.1103/PhysRevE.109.044502"
PREPRINT_URL = "https://arxiv.org/html/2304.11266v1"
BOOTCAMP_URL = "https://github.com/ahansenlab/polysim_bootcamp/blob/main/tutorial/calibration.ipynb"


def graph_laplacian(beads, loops=()):
    """Unit-weight, connected open chain with optional extra static springs."""
    if isinstance(beads, bool) or not isinstance(beads, (int, np.integer)) or beads < 2:
        raise ValueError("beads must be an integer of at least two")
    matrix = np.zeros((beads, beads), dtype=float)
    edges = [(i, i+1) for i in range(beads-1)]
    for edge in loops:
        if len(edge) != 2 or any(isinstance(i, bool) or not isinstance(i, (int, np.integer)) for i in edge):
            raise ValueError("Each loop must contain two integer bead indices")
        a, b = edge
        if not 0 <= a < beads or not 0 <= b < beads or a == b:
            raise ValueError("Loop endpoints must be distinct and within the chain")
        if tuple(sorted((a, b))) in {tuple(sorted(e)) for e in edges}:
            raise ValueError("Duplicate springs are not supported in this benchmark")
        edges.append((a, b))
    for a, b in edges:
        matrix[a, a] += 1; matrix[b, b] += 1
        matrix[a, b] -= 1; matrix[b, a] -= 1
    return matrix


def eigenmode_msd(matrix, times, diffusion=1., spring_k=3., *, dt=None,
                  remove_com=True):
    """3D bead-averaged MSD from equilibrium Gaussian internal modes.

    With dt=None this is the continuous-time OU result. For dt>0 it is the
    exact Euler-Maruyama prediction, initialized from the *continuous*
    equilibrium covariance, matching the stochastic fixture below.
    """
    matrix = np.asarray(matrix, dtype=float)
    times = np.asarray(times, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("The spring matrix must be square")
    if not np.all(np.isfinite(matrix)) or not np.allclose(matrix, matrix.T):
        raise ValueError("The spring matrix must be finite and symmetric")
    if np.any(~np.isfinite(times)) or np.any(times < 0):
        raise ValueError("Lag times must be finite and nonnegative")
    if not math.isfinite(diffusion) or diffusion <= 0 or not math.isfinite(spring_k) or spring_k <= 0:
        raise ValueError("Diffusion and stiffness must be positive and finite")
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues[0] < -1e-9 or np.sum(eigenvalues < 1e-9) != 1:
        raise ValueError("The spring graph must have one translational zero mode")
    values = eigenvalues[eigenvalues > 1e-9]
    variance = 1 / (spring_k * values)  # kBT=1, per Cartesian coordinate
    rates = diffusion * spring_k * values
    if dt is None:
        contributions = 2 * variance * (-np.expm1(-times[:, None]*rates))
    else:
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be finite and positive")
        a = 1-rates*dt
        if np.any(np.abs(a) >= 1):
            raise ValueError("The Euler step is outside its stability interval")
        steps = np.rint(times/dt).astype(int)
        if not np.allclose(times, steps*dt, atol=1e-10, rtol=1e-10):
            raise ValueError("Discrete lag times must be multiples of dt")
        decay = a[None, :]**steps[:, None]
        contributions = variance*(1-decay)**2 + 2*diffusion*dt*(1-decay**2)/(1-a*a)
    result = 3 * np.sum(contributions, axis=1) / len(matrix)
    if not remove_com:
        result += 6*diffusion*times/len(matrix)
    return result


def _gaussian_fixture(beads, replicas, loops, diffusion, spring_k, dt, seed):
    matrix = graph_laplacian(beads, loops)
    values, vectors = np.linalg.eigh(matrix)
    rng = np.random.default_rng(seed)
    modes = rng.normal(size=(replicas, beads, 3))
    modes[:, 0] = 0
    modes[:, 1:] /= np.sqrt(spring_k*values[1:])[None, :, None]
    positions = np.einsum('ij,rjd->rid', vectors, modes).reshape(-1, 3)
    sim = SimpleNamespace(positions=positions, diffusion=np.full(replicas*beads, diffusion),
                          dt=dt, rng=rng, radius=1e6)

    def force():
        p = sim.positions.reshape(replicas, beads, 3)
        result = np.zeros_like(p)
        differences = p[:, 1:] - p[:, :-1]
        result[:, :-1] += spring_k*differences
        result[:, 1:] -= spring_k*differences
        for a, b in loops:
            delta = spring_k*(p[:, b]-p[:, a])
            result[:, a] += delta; result[:, b] -= delta
        return result.reshape(-1, 3)

    sim.get_forces = force
    sim.integration_substeps = lambda: 1
    sim.validate_positions = lambda p: Simulation.validate_positions(sim, p)
    return sim, matrix


def _summarize(values):
    return float(np.mean(values)), float(np.std(values, ddof=1)/math.sqrt(len(values)))


def _sample_chain(beads, replicas, loops, diffusion, spring_k, dt, seed, sample_steps):
    sim, matrix = _gaussian_fixture(beads, replicas, loops, diffusion, spring_k, dt, seed)
    initial = sim.positions.reshape(replicas, beads, 3).copy()
    internal, com = [], []
    wanted = set(sample_steps)
    for step in range(1, max(sample_steps)+1):
        Simulation._brownian_move(sim)
        if step not in wanted:
            continue
        delta = sim.positions.reshape(replicas, beads, 3)-initial
        center = delta.mean(axis=1)
        internal_values = np.mean(np.sum((delta-center[:, None, :])**2, axis=2), axis=1)
        internal.append(_summarize(internal_values))
        com.append(_summarize(np.sum(center*center, axis=1)))
    return matrix, internal, com


def _sample_free(particles, diffusion, dt, seed, sample_steps):
    sim = SimpleNamespace(positions=np.zeros((particles, 3)), diffusion=np.full(particles, diffusion),
                          dt=dt, rng=np.random.default_rng(seed), radius=1e6)
    sim.get_forces = lambda: np.zeros_like(sim.positions)
    sim.integration_substeps = lambda: 1
    sim.validate_positions = lambda p: Simulation.validate_positions(sim, p)
    measurements, wanted = [], set(sample_steps)
    for step in range(1, max(sample_steps)+1):
        Simulation._brownian_move(sim)
        if step in wanted:
            measurements.append(_summarize(np.sum(sim.positions**2, axis=1)))
    return measurements


def _curve(identifier, label, kind, times, measured, exact, discrete):
    return {"id": identifier, "label": label, "kind": kind,
            "x_label": "Lag time (reduced)", "y_label": "3D MSD (reduced length²)",
            "points": [{"lag": float(t), "observed": float(m[0]), "standard_error": float(m[1]),
                        "exact": float(e), "discrete_exact": float(d)}
                       for t, m, e, d in zip(times, measured, exact, discrete)]}


def _check(identifier, label, observed, expected, tolerance, unit="dimensionless", **details):
    error = abs(observed-expected)
    return {"id": identifier, "label": label, "passed": bool(math.isfinite(error) and error <= tolerance),
            "observed": float(observed), "expected": float(expected), "tolerance": float(tolerance),
            "absolute_error": float(error), "unit": unit, "details": details}


def _curve_agreement(curve):
    # Replicate chains are the independent sampling units, not individual beads
    # and not overlapping lags. Cross-lag correlation is not used to shrink SE.
    zscores = [abs(p["observed"]-p["discrete_exact"])/max(p["standard_error"], 1e-15)
               for p in curve["points"]]
    return _check(curve["id"]+"_agreement", curve["label"]+": sampled versus exact discrete MSD",
                  max(zscores), 0., 6., "standard errors", lag_points=len(zscores),
                  criterion="Every lag is within six independent-ensemble standard errors; lags are correlated")


def run_benchmarks(*, beads=64, replicas=192, dt=.005, duration=16., seed=924731):
    """Run deterministic-seed CPU reference fixtures without touching live state."""
    for name, value, low, high in (("beads", beads, 16, 256), ("replicas", replicas, 32, 2048)):
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            raise ValueError(f"{name} must be an integer from {low} to {high}")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**32-3:
        raise ValueError("seed must be an integer from zero to 2**32-3")
    if not math.isfinite(dt) or dt <= 0 or not math.isfinite(duration) or duration < 8:
        raise ValueError("dt must be positive and duration at least eight reduced time units")
    steps = round(duration/dt)
    if not math.isclose(steps*dt, duration, rel_tol=1e-10) or not 32 <= steps <= 100_000:
        raise ValueError("duration must be an integer multiple of dt, with 32–100000 steps")
    started = time.perf_counter()
    diffusion, spring_k = 1., 3.
    loops = [(beads//8, 3*beads//8), (3*beads//8, 5*beads//8),
             (5*beads//8, 7*beads//8), (beads//8, 7*beads//8)]
    sample_steps = np.unique(np.rint(np.geomspace(1, steps, 30)).astype(int)).tolist()
    times = np.array(sample_steps)*dt
    # Check stability before spending time on stochastic integration.
    eigenmode_msd(graph_laplacian(beads, loops), times, diffusion, spring_k, dt=dt)
    curves = []
    free = _sample_free(2048, diffusion, dt, seed, sample_steps)
    curves.append(_curve("free", "Free particles", "free", times, free,
                         6*diffusion*times, 6*diffusion*times))
    for identifier, label, links, offset in (("rouse", "Rouse chain · internal", [], 1),
                                           ("looped", "Static loops · internal", loops, 2)):
        matrix, internal, com = _sample_chain(beads, replicas, links, diffusion, spring_k,
                                               dt, seed+offset, sample_steps)
        exact = eigenmode_msd(matrix, times, diffusion, spring_k)
        discrete = eigenmode_msd(matrix, times, diffusion, spring_k, dt=dt)
        curves.append(_curve(identifier, label, identifier, times, internal, exact, discrete))
        curves.append(_curve(identifier+"_com", label.split(" ·")[0]+" · center of mass", "com",
                             times, com, 6*diffusion*times/beads, 6*diffusion*times/beads))
    checks = [_curve_agreement(curve) for curve in curves]
    by_id = {c["id"]: c for c in curves}
    for identifier, expected in (("free", diffusion), ("rouse_com", diffusion/beads),
                                  ("looped_com", diffusion/beads)):
        point = by_id[identifier]["points"][-1]
        checks.append(_check(identifier+"_diffusion", identifier.replace("_", " ")+": inferred diffusion",
                             point["observed"]/(6*point["lag"]), expected,
                             6*point["standard_error"]/(6*point["lag"]), "reduced length² / time",
                             inference="D = MSD/(6 lag), final lag; no fitted intercept"))
    for identifier in ("rouse", "looped"):
        bias = max(abs(p["discrete_exact"]-p["exact"])/p["exact"]
                   for p in by_id[identifier]["points"])
        checks.append(_check(identifier+"_timestep_bias", identifier+": maximum theoretical timestep bias",
                             bias, 0., .05, "relative error", criterion="Five percent continuum-bias ceiling"))
    slope_points = [p for p in by_id["rouse"]["points"] if .5 <= p["lag"] <= 4]
    log_times = np.log([p["lag"] for p in slope_points])
    observed_slope = float(np.polyfit(log_times, np.log([p["observed"] for p in slope_points]), 1)[0])
    exact_slope = float(np.polyfit(log_times, np.log([p["exact"] for p in slope_points]), 1)[0])
    checks.append(_check("rouse_subdiffusion", "Intermediate Rouse slope follows the finite-chain prediction",
                         observed_slope, exact_slope, .10, "log–log slope", lag_window=[.5, 4.],
                         asymptotic_infinite_chain_slope=.5,
                         caveat="Finite chain and COM subtraction shift the slope; this is not a fitted anomalous diffusion model"))
    checks.append(_check("rouse_subdiffusive_regime", "The finite-chain reference lies in an intermediate subdiffusive regime",
                         exact_slope, .5, .15, "log–log slope"))
    plain, loop = by_id["rouse"]["points"][-1], by_id["looped"]["points"][-1]
    ratio = loop["observed"]/plain["observed"]
    expected_ratio = loop["discrete_exact"]/plain["discrete_exact"]
    ratio_se = ratio*math.sqrt((loop["standard_error"]/loop["observed"])**2 +
                               (plain["standard_error"]/plain["observed"])**2)
    checks.append(_check("static_loop_ratio", "Static-loop suppression agrees with the predicted MSD ratio",
                         ratio, expected_ratio, 6*ratio_se, "looped / unlooped internal MSD",
                         lag=float(times[-1]), uncertainty="Delta-method SE from independent chain ensembles"))
    checks.append({"id": "static_loop_suppression", "label": "Static loops reduce internal mobility",
                   "passed": bool(loop["observed"]+6*math.hypot(loop["standard_error"], plain["standard_error"]) < plain["observed"]),
                   "observed": ratio, "expected": "less than 1 with six-SE separation", "tolerance": None,
                   "unit": "looped / unlooped internal MSD", "details": {"lag": float(times[-1])}})
    passed = sum(c["passed"] for c in checks)
    report = {
        "schema_version": 1, "status": "passed" if passed == len(checks) else "failed",
        "title": "Gaussian-chain stochastic reference benchmarks",
        "scope": "Implementation reference, not an experimentally calibrated model or full paper replication",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "provenance": {"app_version": __version__, "python": platform.python_version(), "numpy": np.__version__,
                       "platform": "CPU / NumPy", "integrator": "Simulation._brownian_move",
                       "integrator_source_sha256": hashlib.sha256(inspect.getsource(Simulation._brownian_move).encode()).hexdigest(),
                       "benchmark_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       "elapsed_seconds": time.perf_counter()-started, "external_code_copied": False,
                       "experimental_data_used": False, "full_paper_replication": False},
        "parameters": {"beads": beads, "replicas": replicas, "free_particles": 2048, "dt": dt,
                       "duration": duration, "steps": steps, "bead_diffusion": diffusion, "spring_k": spring_k,
                       "kBT": 1., "dimensions": 3, "seeds": {"free": seed, "rouse": seed+1, "looped": seed+2},
                       "loops": [list(e) for e in loops], "sample_steps": sample_steps,
                       "initialization": "Independent exact continuous-equilibrium internal Gaussian modes; COM starts at zero",
                       "boundaries": "Free chain ends; no confinement or repulsion",
                       "force": "Zero-rest springs: F_i = -k sum_j L_ij r_j",
                       "MSD_ensemble": "Independent replicas from one equilibrium origin; internal curves subtract each chain COM",
                       "time_units_calibrated": False, "length_units_calibrated": False},
        "checks": checks, "curves": curves,
        "summary": {"passed": passed, "total": len(checks), "rouse_slope_observed": observed_slope,
                    "rouse_slope_exact": exact_slope, "loop_suppression_observed": ratio,
                    "loop_suppression_exact": loop["exact"]/plain["exact"]},
        "limitations": [
            "This controlled fixture exercises the actual NumPy Brownian integrator with custom Gaussian forces. It bypasses production forces, reaction scheduling, confinement, adaptive subdivision and GPU integration.",
            "The production app uses finite-rest springs, repulsion and optional regulatory interactions; these Gaussian reference curves are not predictions for that model.",
            "The static-loop locations and small chain size are illustrative choices. Published loop-extrusion distributions, nonuniform friction, biological parameters and experimental datasets are not reproduced.",
            "The paper removes the translational mode. Internal curves here do the same after sampling; separate COM curves explicitly test free translation with D/N.",
            "Error bars are standard errors across independent replicas. Lags are correlated. The six-SE checks are regression tolerances, not biological confidence intervals.",
            "No seconds, nanometers or experimentally fitted parameters are inferred. Biological calibration requires matched observables and measured data.",
        ],
        "sources": [{"title": "Yuan et al. (2024), Gaussian-chain and loop dynamical-matrix theory", "url": PAPER_URL},
                    {"title": "Open preprint: equations 1, 20, 21 and 26", "url": PREPRINT_URL},
                    {"title": "Hansen Lab polymer bootcamp: two-locus length and time calibration workflow", "url": BOOTCAMP_URL}],
    }
    # Fail rather than publish NaN/Infinity in an API response.
    json.dumps(report, allow_nan=False)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replicas", type=int, default=192)
    parser.add_argument("--seed", type=int, default=924731)
    parser.add_argument("--dt", type=float, default=.005)
    parser.add_argument("--duration", type=float, default=16.)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("Output already exists; choose a new path to preserve prior benchmark provenance")
    report = run_benchmarks(replicas=args.replicas, seed=args.seed, dt=args.dt, duration=args.duration)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(f"{report['status']}: {report['summary']['passed']}/{report['summary']['total']} checks; {args.output}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
