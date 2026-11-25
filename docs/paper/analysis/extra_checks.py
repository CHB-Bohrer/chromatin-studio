"""Small validation analyses for the software paper that no component study covers yet.

Each check isolates one mechanism, runs it briefly with the real engine and
compares it with a closed-form answer. Output is one JSON file of curves and
check results; docs/paper/figures/make_figures.py draws them.

    python docs/paper/analysis/extra_checks.py --output docs/paper/data/extra-checks.json

Checks
------
ideal_chain   Bare chain, no bending, no repulsion, no wall contact. Bond vectors are then
              independent, so <R^2(s)> = s <b^2> exactly, and contact probability approaches
              the Gaussian law P(s) ~ s^(-3/2).
binding       Proteins held at fixed positions next to chromatin (frozen coordinates).
              Each molecule binds with total hazard kon x (beads in reach) and releases with
              koff, so the bound fraction relaxes to K/(K+koff) at rate K+koff.
marks         One mark with background decay only (frozen coordinates): f(t) = exp(-d t).
tethers       Envelope tethering element with frozen coordinates. Beads within capture distance
              that carry the required mark attach with ka and detach with kd, so their attached
              fraction relaxes to ka/(ka+kd); unmarked beads never attach.

These are short runs on the NumPy engine. They show that each rule is implemented as
stated; they are not calibrations and say nothing about physical units.
"""

import argparse
import datetime
import json
import math
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from chromatin_studio.engine import Simulation  # noqa: E402
from chromatin_studio.components.core import (bare_config, chemistry_only_step, mean_sem, protein,  # noqa: E402
                                              validated)


RULE = "every checked time within max(0.01, 3 x SEM over replicas) of the closed form"


def _worst(stats, expected, tolerance):
    """Largest absolute error, and whether every point is within max(tolerance, 3 SEM)."""
    errors = [abs(m - e) for (m, _), e in zip(stats, expected)]
    ok = all(abs(m - e) <= max(tolerance, 3 * sem) for (m, sem), e in zip(stats, expected))
    return float(max(errors)), bool(ok)


def ideal_chain(replicas=4, beads=64, steps=200_000, burn_in=40_000, every=200, seed=701):
    r2 = np.zeros((replicas, beads - 1))
    contact = np.zeros((replicas, beads - 1))
    bond2 = np.zeros(replicas)
    cutoff = 1.0
    for rep in range(replicas):
        config = bare_config(beads, seed=seed + rep, dt=0.01, radius=60.0, persistence_length=0.0,
                             repulsion=0.0, bond_k=10.0, diffusion=1.0)
        sim = Simulation(validated(config))
        sim.step(burn_in)
        samples = 0
        for _ in range((steps - burn_in) // every):
            sim.step(every)
            x = sim.positions[:beads]
            for s in range(1, beads):
                d = x[s:] - x[:-s]
                dist2 = np.einsum("ij,ij->i", d, d)
                r2[rep, s - 1] += dist2.mean()
                contact[rep, s - 1] += (dist2 < cutoff ** 2).mean()
            b = x[1:] - x[:-1]
            bond2[rep] += np.einsum("ij,ij->i", b, b).mean()
            samples += 1
        r2[rep] /= samples
        contact[rep] /= samples
        bond2[rep] /= samples
    s = np.arange(1, beads)
    b2, b2_sem = mean_sem(bond2.tolist())
    r2_mean, r2_sem = r2.mean(0), r2.std(0, ddof=1) / math.sqrt(replicas)
    c_mean, c_sem = contact.mean(0), contact.std(0, ddof=1) / math.sqrt(replicas)
    # Gaussian contact probability for a segment of s bonds with <R^2> = s b^2:
    # P(|R| < a) for a 3D Gaussian with per-axis variance s b^2 / 3.
    def gaussian_contact(n):
        sigma = math.sqrt(n * b2 / 3.0)
        z = cutoff / sigma
        return math.erf(z / math.sqrt(2)) - math.sqrt(2 / math.pi) * z * math.exp(-z * z / 2)
    exact_r2 = s * b2
    rel = np.abs(r2_mean - exact_r2) / exact_r2
    r2_ok = bool(np.all(np.abs(r2_mean - exact_r2) <= np.maximum(0.05 * exact_r2, 3 * r2_sem)))
    window = (s >= 8) & (s <= 32)
    slope = float(np.polyfit(np.log(s[window]), np.log(np.maximum(c_mean[window], 1e-12)), 1)[0])
    return {
        "params": dict(replicas=replicas, beads=beads, steps=steps, burn_in=burn_in, sample_every=every,
                       bond_k=10.0, dt=0.01, contact_cutoff=cutoff),
        "s": s.tolist(), "r2": r2_mean.tolist(), "r2_sem": r2_sem.tolist(), "r2_exact": exact_r2.tolist(),
        "contact": c_mean.tolist(), "contact_sem": c_sem.tolist(),
        "contact_gaussian": [gaussian_contact(n) for n in s],
        "bond2": b2, "bond2_sem": b2_sem,
        "checks": [
            {"id": "r2_scaling", "label": "<R²(s)> = s<b²> (largest relative error, all s)",
             "observed": float(rel.max()), "rule": "every s within max(5%, 3 x SEM)", "passed": r2_ok},
            {"id": "contact_exponent", "label": "Contact exponent fitted over s = 8 to 32",
             "observed": slope, "expected": -1.5, "passed": bool(abs(slope + 1.5) < 0.15)},
        ],
    }


def _frozen_binding_sim(seed, kon, koff, count):
    beads = 4 * count
    config = bare_config(beads, seed=seed, dt=0.01, radius=60.0, repulsion=0.0, bp_per_bead=200)
    config["proteins"] = [protein("tf", count, chromatin_affinity={"unmodified": 1.0}, on_rate=kon,
                                  off_rate=koff, radius=0.28)]
    sim = Simulation(validated(config))
    # Lay the chain straight, spacing 1, and put each protein 0.3 above every fourth bead,
    # so the beads in reach of one protein are fixed and no two proteins share a bead.
    sim.positions[:beads] = np.column_stack([np.arange(beads) - beads / 2, np.zeros(beads), np.zeros(beads)])
    for k, particle in enumerate(sim.protein_particles):
        sim.positions[particle["index"]] = sim.positions[4 * k + 1] + np.array([0.0, 0.3, 0.0])
    sim.bindings = []
    capture = 0.28 + 1.0
    reach = []
    for k, particle in enumerate(sim.protein_particles):
        d = np.linalg.norm(sim.positions[:beads] - sim.positions[particle["index"]], axis=1)
        reach.append(int(np.count_nonzero(d < capture)))
    return sim, reach


def binding(replicas=32, count=500, kon=0.05, koff=0.2, interval=0.05, times=(1, 2, 4, 8, 16, 40), seed=801):
    runs = []
    reach = None
    for rep in range(replicas):
        sim, reach = _frozen_binding_sim(seed + rep, kon, koff, count)
        values = []
        for t in times:
            chemistry_only_step(sim, round((t - sim.time) / interval), interval)
            values.append(len({sim._particle_info[b["b"]]["molecule"] for b in sim.bindings
                               if b["kind"] == "chromatin"}) / count)
        runs.append(values)
    n_reach = reach[0]
    total_on = kon * n_reach
    def expected(t):
        f_inf = total_on / (total_on + koff)
        return f_inf * (1 - math.exp(-(total_on + koff) * t))
    stats = [mean_sem([r[i] for r in runs]) for i in range(len(times))]
    grid = np.linspace(0, max(times), 120)
    err, ok = _worst(stats, [expected(t) for t in times], 0.01)
    return {"params": dict(replicas=replicas, molecules=count, kon=kon, koff=koff, interval=interval,
                           beads_in_reach=n_reach),
            "times": list(times), "bound": [m for m, _ in stats], "bound_sem": [e for _, e in stats],
            "grid": grid.tolist(), "theory": [expected(t) for t in grid],
            "checks": [{"id": "bound_fraction", "label": "Bound fraction versus two-state solution (largest error)",
                        "observed": err, "rule": RULE, "passed": ok}]}


def marks(replicas=8, beads=2000, decay=0.05, interval=0.1, times=(2, 5, 10, 20, 40, 60), seed=901):
    runs = []
    for rep in range(replicas):
        config = bare_config(beads, seed=seed + rep, dt=0.01, radius=60.0, repulsion=0.0, bp_per_bead=200,
                             initial_marks=[{"start": 0, "end": beads - 1, "mark": "m"}])
        config["marks"] = [{"id": "m", "name": "mark", "color": "#8888ff", "decay_rate": decay,
                            "spread_rate": 0.0, "contact_radius": 1.0}]
        sim = Simulation(validated(config))
        values = []
        for t in times:
            chemistry_only_step(sim, round((t - sim.time) / interval), interval)
            values.append(float(np.mean(sim.marks["m"])))
        runs.append(values)
    stats = [mean_sem([r[i] for r in runs]) for i in range(len(times))]
    grid = np.linspace(0, max(times), 120)
    err, ok = _worst(stats, [math.exp(-decay * t) for t in times], 0.01)
    return {"params": dict(replicas=replicas, beads=beads, decay_rate=decay, interval=interval),
            "times": list(times), "marked": [m for m, _ in stats], "marked_sem": [e for _, e in stats],
            "grid": grid.tolist(), "theory": [math.exp(-decay * t) for t in grid],
            "checks": [{"id": "mark_decay", "label": "Marked fraction versus exp(-d t) (largest error)",
                        "observed": err, "rule": RULE, "passed": ok}]}


def tethers(replicas=8, beads=1600, ka=0.3, kd=0.1, interval=0.05, times=(1, 2, 4, 8, 16, 40), seed=1001):
    runs, unmarked_attached, eligible_counts = [], [], []
    for rep in range(replicas):
        config = bare_config(beads, seed=seed + rep, dt=0.01, radius=8.0, repulsion=0.0, bp_per_bead=200,
                             initial_marks=[{"start": 0, "end": beads // 2 - 1, "mark": "lad"}])
        config["marks"] = [{"id": "lad", "name": "LAD mark", "color": "#8888ff", "decay_rate": 0.0,
                            "spread_rate": 0.0, "contact_radius": 1.0}]
        config["nucleus"] = {"enabled": True, "center": [0.0, 0.0, 0.0], "radii": [12.0, 12.0, 12.0],
                             "envelope_stiffness": 30.0,
                             "tethers": [{"id": "t", "name": "tether", "color": "#7a8fa6", "start": 0,
                                          "end": beads - 1, "enabled": True, "anchor": "envelope",
                                          "attach_rate": ka, "detach_rate": kd, "capture_distance": 7.0,
                                          "spring_k": 2.0, "rest_length": 0.5, "required_mark": "lad"}]}
        sim = Simulation(validated(config))
        distance = np.abs(sim.nucleus.signed_distance(sim.positions[:beads]))
        within = distance <= 7.0
        marked = np.asarray(sim.marks["lad"], dtype=bool)
        eligible = within & marked
        eligible_counts.append(int(eligible.sum()))
        values = []
        for t in times:
            chemistry_only_step(sim, round((t - sim.time) / interval), interval)
            attached = np.zeros(beads, dtype=bool)
            attached[[tt["bead"] for tt in sim.tethers]] = True
            values.append(float(attached[eligible].mean()) if eligible.any() else 0.0)
        unmarked_attached.append(int((attached & ~marked).sum()))
        runs.append(values)
    def expected(t):
        return ka / (ka + kd) * (1 - math.exp(-(ka + kd) * t))
    stats = [mean_sem([r[i] for r in runs]) for i in range(len(times))]
    grid = np.linspace(0, max(times), 120)
    err, ok = _worst(stats, [expected(t) for t in times], 0.01)
    return {"params": dict(replicas=replicas, beads=beads, attach_rate=ka, detach_rate=kd, interval=interval,
                           capture_distance=7.0, eligible_beads=eligible_counts),
            "times": list(times), "attached": [m for m, _ in stats], "attached_sem": [e for _, e in stats],
            "grid": grid.tolist(), "theory": [expected(t) for t in grid],
            "unmarked_attached": unmarked_attached,
            "checks": [{"id": "tether_two_state", "label": "Attached fraction of eligible beads versus two-state "
                        "solution (largest error)", "observed": err, "rule": RULE, "passed": ok},
                       {"id": "tether_mark_gate", "label": "Beads without the required mark that attached",
                        "observed": int(sum(unmarked_attached)), "expected": 0,
                        "passed": sum(unmarked_attached) == 0}]}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", required=True)
    parser.add_argument("--only", nargs="*", default=None)
    args = parser.parse_args()
    out = Path(args.output)
    result = json.loads(out.read_text()) if out.exists() else {}
    for name, fn in [("ideal_chain", ideal_chain), ("binding", binding), ("marks", marks), ("tethers", tethers)]:
        if args.only and name not in args.only:
            continue
        print("running", name, flush=True)
        result[name] = fn()
        for c in result[name]["checks"]:
            print("  ", "PASS" if c["passed"] else "FAIL", c["label"], c["observed"])
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True,
                                text=True).stdout.strip()
    except OSError:
        commit = None
    result["provenance"] = {"created": datetime.date.today().isoformat(), "source_commit": commit,
                            "python": platform.python_version(), "numpy": np.__version__,
                            "engine": "NumPy reference", "script": "docs/paper/analysis/extra_checks.py"}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
