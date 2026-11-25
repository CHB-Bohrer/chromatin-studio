"""Nucleosome occupancy alone: two-state turnover on fixed sites.

Every site switches occupied -> vacant with hazard E and back with hazard R.
From a fully occupied start the occupied fraction must follow

    f(t) = R/(E+R) + E/(E+R) exp(-(E+R) t)

Coordinates are frozen (chemistry only): occupancy turnover does not depend on
3D motion when no protein or polymerase contributes.
"""

import math

import numpy as np

from ..engine import Simulation
from .core import (ANALYTICAL, bare_config, check, chemistry_only_step, error_bar, figure, mean_sem, series, source,
                   validated)

STUDY = {
    "id": "nucleosome_turnover",
    "title": "Background nucleosome turnover",
    "components": ["nucleosome occupancy (eviction, reassembly)"],
    "question": "Does site occupancy relax to R/(E+R) at rate E+R?",
    "evidence": [ANALYTICAL],
    "sources": [
        source("Closed-form solution of the two-state master equation (derivation in the module docstring)", None,
               "Occupancy relaxation of independent two-state sites."),
        source("Deal, Henikoff & Henikoff (2010) Science 328:1161", "https://doi.org/10.1126/science.1186777",
               "Context only: nucleosome turnover of roughly 1-1.5 h at active Drosophila genes. Not compared here.",
               used_as="context"),
    ],
    "scales": {
        "smoke": {"beads": 400, "eviction_rate": 0.02, "reassembly_rate": 0.06, "reaction_interval": 0.5,
                  "times": [10.0, 40.0, 200.0], "replicas": 2, "seed": 31},
        "standard": {"beads": 4000, "eviction_rate": 0.02, "reassembly_rate": 0.06, "reaction_interval": 0.5,
                     "times": [5.0, 10.0, 20.0, 40.0, 200.0], "replicas": 8, "seed": 31},
    },
}


def _expected(params, time):
    e, r = params["eviction_rate"], params["reassembly_rate"]
    return r / (e + r) + e / (e + r) * math.exp(-(e + r) * time)


def simulate(params, case, seed, save_partial):
    # Sites are independent, so the chain geometry is irrelevant; a tiny
    # diffusion keeps the model valid without mechanics being advanced.
    config = bare_config(params["beads"], seed=seed, radius=60.0, repulsion=0.0,
                         bp_per_bead=200,
                         nucleosomes={"enabled": True, "initial_occupancy": 1,
                                      "eviction_rate": params["eviction_rate"],
                                      "reassembly_rate": params["reassembly_rate"],
                                      "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": 1})
    sim = Simulation(validated(config))
    occupied = {}
    for time in params["times"]:
        chemistry_only_step(sim, round((time - sim.time) / params["reaction_interval"]), params["reaction_interval"])
        occupied[str(time)] = float(np.mean(sim.nucleosome_occupied))
    return {"occupied_fraction": occupied, "counters": dict(sim.nucleosome_counters)}


def analyze(params, results):
    done = [r["observables"] for r in results if r.get("complete")]
    checks, means, sems = [], [], []
    for time in params["times"]:
        if not done:
            break
        value, sem = mean_sem([o["occupied_fraction"][str(time)] for o in done])
        means.append(value)
        sems.append(sem)
        checks.append(check(f"occupancy_t{time:g}", f"Occupied fraction at t = {time:g}", value,
                            _expected(params, time), evidence=ANALYTICAL, tolerance=0.01, sem=sem))
    figures = []
    if means:
        grid = np.geomspace(min(params["times"]) / 4, max(params["times"]), 80)
        bars = [error_bar(m, s) for m, s in zip(means, sems)]
        figures.append(figure(
            "relaxation", "Occupancy relaxes to the two-state steady state",
            x_label="reduced time", y_label="occupied fraction of sites", x_scale="log",
            series=[series("Two-state solution", grid, [_expected(params, t) for t in grid], role="reference",
                           style="line"),
                    series("Engine, mean ± SEM over replicas", params["times"], means,
                           low=[b[0] for b in bars], high=[b[1] for b in bars])],
            note=f"Every site starts occupied; eviction {params['eviction_rate']:g}, "
                 f"reassembly {params['reassembly_rate']:g} per reduced time."))
    return {"checks": checks, "findings": [], "figures": figures,
            "summary": {"steady_state": params["reassembly_rate"] / (params["eviction_rate"] + params["reassembly_rate"]),
                        "relaxation_rate": params["eviction_rate"] + params["reassembly_rate"]}}
