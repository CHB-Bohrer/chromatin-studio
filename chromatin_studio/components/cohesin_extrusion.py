"""Cohesin alone on a 1D lattice: growth speed, residence, barriers.

One extruder at a time, so there are no collisions. Each leg takes Poisson
steps at ``step_rate`` and the complex unloads at ``off_rate`` (and reloads
immediately, as the engine documents). Expected:

    loop growth speed      = 2 x step_rate   beads per time
    mean residence         = 1 / off_rate
    mean loop at unloading = 2 x step_rate / off_rate   (far from chain ends)

The rates are the literature-derived values of the provisional 10 kb model, so
the measured speed can be converted to kb/s. That conversion is input-matched:
it confirms the unit bookkeeping, not an independent prediction.
Coordinates are frozen (chemistry only): extrusion here is a 1D process.
"""

import numpy as np

from ..engine import Simulation
from .core import (ANALYTICAL, EXPERIMENTAL, NUMERICAL, PROVISIONAL_CLOCK, bare_config, check,
                   chemistry_only_step, mean_sem, source, validated)

STUDY = {
    "id": "cohesin_extrusion",
    "title": "Loop extrusion by a single cohesin",
    "components": ["cohesin stepping", "unloading and reloading", "directional barriers"],
    "question": "Do loop growth, residence and barrier stalling follow the configured rates?",
    "evidence": [NUMERICAL, ANALYTICAL, EXPERIMENTAL],
    "sources": [
        source("Kim, Shi, Zhang, Finkelstein & Yu (2019) Science 366:1345", "https://doi.org/10.1126/science.aaz4475",
               "Human cohesin-NIPBL extrudes DNA loops at an average of about 0.5 kb/s in vitro."),
        source("Davidson et al. (2019) Science 366:1338", "https://doi.org/10.1126/science.aaz3418",
               "Single human cohesin complexes extrude loops at up to 2.1 kb/s in vitro."),
        source("Hansen, Pustova, Cattoglio, Tjian & Darzacq (2017) eLife 6:e25776", "https://doi.org/10.7554/eLife.25776",
               "Cohesin chromatin residence of about 22 min in mouse embryonic stem cells (G1)."),
    ],
    "scales": {
        "smoke": {"beads": 600, "bp_per_bead": 10000, "step_rate": 0.0051270739, "off_rate": 0.00015536588,
                  "reaction_interval": 2.0, "residences": 40, "replicas": 2, "seed": 51},
        "standard": {"beads": 2000, "bp_per_bead": 10000, "step_rate": 0.0051270739, "off_rate": 0.00015536588,
                     "reaction_interval": 2.0, "residences": 1500, "replicas": 8, "seed": 51},
    },
}


def cases(params):
    n = params["beads"]
    return [{"id": "free", "barriers": []},
            {"id": "barriers", "barriers": [{"bead": n // 2 - 40, "direction": "left"},
                                            {"bead": n // 2 + 40, "direction": "right"}]}]


def simulate(params, case, seed, save_partial):
    n = params["beads"]
    interval = params["reaction_interval"]
    config = bare_config(n, seed=seed, radius=40.0, bp_per_bead=params["bp_per_bead"],
                         barriers=case["barriers"])
    config["cohesin"] = {"enabled": True, "count": 1, "step_rate": params["step_rate"],
                         "off_rate": params["off_rate"], "spring_k": 2.5}
    sim = Simulation(validated(config))
    left_barrier = next((b["bead"] for b in case["barriers"] if b["direction"] == "left"), None)
    right_barrier = next((b["bead"] for b in case["barriers"] if b["direction"] == "right"), None)
    current = sim.cohesin[0]
    start_time, previous = sim.time, list(current)
    durations, final_loops, censored = [], [], []
    touched_end = False
    growth, growth_time = 0.0, 0.0
    inside, violations, looped_time, total_time = None, 0, 0.0, 0.0

    def loaded_inside(pair):
        return left_barrier is not None and left_barrier <= pair[0] and pair[1] <= right_barrier

    inside = loaded_inside(current)
    while len(durations) < params["residences"]:
        chemistry_only_step(sim, 1, interval)
        pair = sim.cohesin[0]
        if pair is not current:           # unloaded and reloaded within this interval
            durations.append(sim.time - start_time)
            final_loops.append(previous[1] - previous[0])
            censored.append(touched_end)
            current, start_time, touched_end = pair, sim.time, False
            inside = loaded_inside(pair)
        else:
            at_end = pair[0] == 0 or pair[1] == n - 1 or previous[0] == 0 or previous[1] == n - 1
            touched_end = touched_end or at_end
            if not at_end and not case["barriers"]:
                growth += (pair[1] - pair[0]) - (previous[1] - previous[0])
                growth_time += interval
            if inside and (pair[0] < left_barrier or pair[1] > right_barrier):
                violations += 1
            if left_barrier is not None:
                looped_time += interval * (pair[0] == left_barrier and pair[1] == right_barrier)
                total_time += interval
        previous = list(pair)
    keep = ~np.asarray(censored)
    return {"mean_residence": float(np.mean(durations)), "residences": len(durations),
            "growth_speed": growth / growth_time if growth_time else None,
            "mean_final_loop_uncensored": float(np.mean(np.asarray(final_loops)[keep])) if keep.any() else None,
            "censored_fraction": float(1 - keep.mean()), "barrier_violations": int(violations),
            "fully_looped_fraction": looped_time / total_time if total_time else None}


def analyze(params, results):
    free = [r["observables"] for r in results if r.get("complete") and r["case"]["id"] == "free"]
    walls = [r["observables"] for r in results if r.get("complete") and r["case"]["id"] == "barriers"]
    rate, off = params["step_rate"], params["off_rate"]
    checks, findings = [], []
    if free:
        speed, speed_sem = mean_sem([o["growth_speed"] for o in free])
        checks.append(check("growth_speed", "Loop growth speed", speed, 2 * rate, evidence=ANALYTICAL,
                            tolerance=0.03, sem=speed_sem, unit="beads per reduced time"))
        residence, residence_sem = mean_sem([o["mean_residence"] for o in free])
        checks.append(check("residence", "Mean residence", residence, 1 / off, evidence=ANALYTICAL,
                            tolerance=0.05, sem=residence_sem, unit="reduced time",
                            note="Interval-endpoint timing adds up to one reaction interval."))
        loop, loop_sem = mean_sem([o["mean_final_loop_uncensored"] for o in free])
        checks.append(check("processivity", "Mean loop size at unloading", loop, 1 + 2 * rate / off,
                            evidence=ANALYTICAL, tolerance=0.10, sem=loop_sem, unit="beads",
                            note="Residences that reached a chain end are excluded, which biases this slightly low."))
        kb_per_s = speed * params["bp_per_bead"] / 1000.0 / PROVISIONAL_CLOCK["time_s_per_unit"]
        kb_sem = speed_sem * params["bp_per_bead"] / 1000.0 / PROVISIONAL_CLOCK["time_s_per_unit"] if speed_sem else None
        checks.append(check("speed_kb_per_s", "Loop growth in physical units", kb_per_s, 0.5, evidence=EXPERIMENTAL,
                            tolerance=0.05, sem=kb_sem, unit="kb/s",
                            note="Input-matched to Kim 2019 (0.5 kb/s average) through the provisional clock "
                                 "(needs_review); Davidson 2019 reports up to 2.1 kb/s. Confirms unit conversion only."))
        minutes = residence * PROVISIONAL_CLOCK["time_s_per_unit"] / 60.0
        checks.append(check("residence_min", "Residence in physical units", minutes, 22.0, evidence=EXPERIMENTAL,
                            tolerance=0.05, sem=residence_sem * PROVISIONAL_CLOCK["time_s_per_unit"] / 60.0
                            if residence_sem else None, unit="min",
                            note="Input-matched to Hansen 2017 (about 22 min) through the provisional clock."))
    if walls:
        checks.append(check("barrier_violations", "Legs crossing a blocking barrier",
                            float(sum(o["barrier_violations"] for o in walls)), 0.0, evidence=NUMERICAL,
                            low=0, high=0, unit="events"))
        looped, looped_sem = mean_sem([o["fully_looped_fraction"] for o in walls])
        findings.append(f"With one extruder and impermeable convergent barriers 800 kb apart, both anchors are held "
                        f"{looped:.3f} of the time" + (f" (± {looped_sem:.3f})" if looped_sem else "") +
                        ". Barriers never release a stalled leg; only cohesin unloading does.")
    return {"checks": checks, "findings": findings,
            "summary": {"expected_speed": 2 * rate, "expected_residence": 1 / off, "mapping": PROVISIONAL_CLOCK}}
