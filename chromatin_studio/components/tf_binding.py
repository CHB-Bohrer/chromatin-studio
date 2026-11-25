"""DNA plus one DNA-binding protein, against live-cell single-molecule tracking.

Experiment (Chen et al. 2014, Sox2 in mouse embryonic stem cells): a factor
finds a specific site by trial and error. It makes 84-97 rounds of 3D diffusion
(3.3-3.7 s each) separated by brief nonspecific contacts (0.75-0.9 s), then
dwells 12.0-14.6 s on the specific site. About 23.5% of molecules are chromatin
bound (fluorescence correlation spectroscopy, an independent measurement).

Model: a coarse chromatin chain and one monovalent species, nothing else. One
species has one ``off_rate``, so the two dwell classes use the occupancy
factors as a representational device: "vacant" beads stand for accessible
specific sites and "occupied" beads for generic chromatin, with turnover off.
The two dwell times are inputs. The association hazard is scanned. For the scan
point whose free interval matches 3.3-3.7 s, the study then asks whether the
trial count and the total bound fraction come out in the measured range.

Units: the provisional 10 kb polymer clock converts seconds to reduced time.
It is needs_review, so every comparison here is conditional on it.
"""

import math

import numpy as np

from ..engine import Simulation
from .core import (EXPERIMENTAL, NUMERICAL, PROVISIONAL_CLOCK, band, bare_config, check, error_bar, figure,
                   mean_sem, protein, series,
                   source, validated)

STUDY = {
    "id": "tf_binding",
    "title": "Transcription-factor search and residence on bare chromatin",
    "components": ["chromatin mechanics", "protein diffusion", "protein-chromatin association and release"],
    "question": "With measured dwell times as inputs, do the free interval, trial count and bound fraction agree with Sox2 tracking?",
    "evidence": [NUMERICAL, EXPERIMENTAL],
    "sources": [
        source("Chen et al. (2014) Cell 156:1274", "https://doi.org/10.1016/j.cell.2014.01.062",
               "Sox2/Oct4 in mouse ES cells: 84-97 rounds of 3D diffusion (3.3-3.7 s) with nonspecific contacts "
               "(0.75-0.9 s) before a specific dwell (12.0-14.6 s); about 23.5% of Sox2 chromatin bound by FCS."),
    ],
    "scales": {
        "smoke": {"beads": 90, "specific_every": 30, "radius": 7.2, "count": 20, "on_rates": [0.165],
                  "diffusions": [0.8], "specific_dwell_s": 13.3, "nonspecific_dwell_s": 0.82,
                  "time_s_per_unit": PROVISIONAL_CLOCK["time_s_per_unit"], "dt": 0.002, "steps": 20000,
                  "burn_in_steps": 4000, "block": 50, "replicas": 2, "seed": 61},
        "standard": {"beads": 450, "specific_every": 90, "radius": 12.2, "count": 60,
                     "on_rates": [0.08, 0.165, 0.33], "diffusions": [0.8, 5.0],
                     "specific_dwell_s": 13.3, "nonspecific_dwell_s": 0.82,
                     "time_s_per_unit": PROVISIONAL_CLOCK["time_s_per_unit"], "dt": 0.002, "steps": 3000000,
                     "burn_in_steps": 200000, "block": 50, "replicas": 2, "seed": 61},
    },
}

EXPERIMENT = {"free_interval_s": (3.3, 3.7), "trials": (84.0, 97.0), "bound_fraction": 0.235,
              "specific_dwell_s": (12.0, 14.6), "nonspecific_dwell_s": (0.75, 0.9)}


def cases(params):
    return [{"id": f"on{rate:g}-D{diffusion:g}", "on_rate": rate, "diffusion": diffusion}
            for diffusion in params["diffusions"] for rate in params["on_rates"]]


def _rates(params):
    specific = params["time_s_per_unit"] / params["specific_dwell_s"]
    return specific, params["specific_dwell_s"] / params["nonspecific_dwell_s"]


def build_config(params, case, seed):
    off_rate, nonspecific_factor = _rates(params)
    n = params["beads"]
    sites = list(range(params["specific_every"] // 2, n, params["specific_every"]))
    config = bare_config(n, seed=seed, dt=params["dt"], radius=params["radius"], bp_per_bead=10000,
                         bond_length=0.85, bond_k=45.0, diffusion=0.08, persistence_length=2.0, repulsion=8.0,
                         nucleosomes={"enabled": True, "initial_occupancy": 1, "eviction_rate": 0, "reassembly_rate": 0,
                                      "initial_regions": [{"start": s, "end": s, "occupied": False} for s in sites],
                                      "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": 1})
    config["proteins"] = [protein("tf", params["count"], diffusion=case["diffusion"], valence=1,
                                  chromatin_affinity={"unmodified": 1.0}, on_rate=case["on_rate"], off_rate=off_rate,
                                  nucleosome_affinity={"occupied": 1.0, "vacant": 1.0},
                                  nucleosome_off_rate_multiplier={"occupied": nonspecific_factor, "vacant": 1.0})]
    return validated(config), sites


def simulate(params, case, seed, save_partial):
    config, sites = build_config(params, case, seed)
    sim = Simulation(config)
    start = params["burn_in_steps"] * params["dt"]
    bound_at, bound_class, released_at = {}, {}, {}
    dwell = {"specific": [], "nonspecific": []}
    binds = {"specific": 0, "nonspecific": 0}
    free_intervals = []
    occupancy = {"specific": 0.0, "nonspecific": 0.0}
    samples, last_id, dropped = 0, 0, 0
    for _ in range(params["steps"] // params["block"]):
        sim.step(params["block"])
        fresh = []
        for event in reversed(sim.genomic_events):
            if event["id"] <= last_id:
                break
            fresh.append(event)
        if fresh and fresh[-1]["id"] != last_id + 1:
            dropped += fresh[-1]["id"] - last_id - 1      # buffer overflow would bias the statistics
        for event in reversed(fresh):
            molecule, time = event.get("molecule"), event["time"]
            if event["type"] == "chromatin_bind":
                kind = "specific" if not sim.nucleosome_occupied[event["bead"]] else "nonspecific"
                bound_at[molecule], bound_class[molecule] = time, kind
                if time >= start:
                    binds[kind] += 1
                    if molecule in released_at and released_at[molecule] >= start:
                        free_intervals.append(time - released_at[molecule])
            elif event["type"] == "chromatin_unbind":
                if molecule in bound_at and bound_at[molecule] >= start:
                    dwell[bound_class[molecule]].append(time - bound_at[molecule])
                released_at[molecule] = time
                bound_at.pop(molecule, None)
        if fresh:
            last_id = fresh[0]["id"]
        if sim.time >= start:
            samples += 1
            for bond in sim.bindings:
                occupancy["specific" if not sim.nucleosome_occupied[bond["a"]] else "nonspecific"] += 1
    free = np.asarray(free_intervals)
    count = params["count"]
    return {"sites": len(sites),
            "dwell_mean": {k: float(np.mean(v)) if v else None for k, v in dwell.items()},
            "dwell_count": {k: len(v) for k, v in dwell.items()}, "binds": binds,
            "free_interval_mean": float(free.mean()) if free.size else None,
            "free_interval_median": float(np.median(free)) if free.size else None,
            "free_interval_count": int(free.size),
            "bound_fraction": {k: v / samples / count for k, v in occupancy.items()} if samples else None,
            "dropped_events": int(dropped), "max_substeps": int(sim.last_substeps)}


def _well_mixed(params, case, sites):
    """Idealized dilute estimate, used only as a diagnostic of rebinding and crowding."""
    off_rate, factor = _rates(params)
    reach = 0.28 + 0.85
    volume = 4.0 / 3.0 * math.pi * params["radius"] ** 3
    density = case["on_rate"] * (4.0 / 3.0 * math.pi * reach ** 3) / volume
    nonspecific = density * (params["beads"] - sites) / (off_rate * factor)
    specific = density * sites / off_rate
    return {"free_interval": 1.0 / (density * params["beads"]),
            "bound_fraction": (nonspecific + specific) / (1 + nonspecific + specific),
            "trials": (params["beads"] - sites) / sites}


def _figures(table):
    """Association-hazard scans, one series per protein diffusion coefficient."""
    if not table:
        return []
    diffusions = sorted({row["diffusion"] for row in table})

    def scan(key, sem_key=None, estimate_key=None):
        plotted = []
        for diffusion in diffusions:
            rows = sorted((row for row in table if row["diffusion"] == diffusion and row.get(key) is not None),
                          key=lambda row: row["on_rate"])
            if not rows:
                continue
            bars = [error_bar(row[key], row.get(sem_key) if sem_key else None) for row in rows]
            group = f"D = {diffusion:g}"
            plotted.append(series(f"Engine, protein D = {diffusion:g}", [row["on_rate"] for row in rows],
                                  [row[key] for row in rows], low=[b[0] for b in bars], high=[b[1] for b in bars],
                                  style="connected", group=group))
            if estimate_key:
                plotted.append(series(f"Dilute well-mixed estimate, D = {diffusion:g}",
                                      [row["on_rate"] for row in rows], [row[estimate_key] for row in rows],
                                      role="reference", style="line", group=group))
        return plotted

    low, high = EXPERIMENT["free_interval_s"]
    bound = EXPERIMENT["bound_fraction"]
    note = "Error bars: ±1 SEM over replicas. The 3D search uses the provisional clock (needs_review)."
    return [
        figure("free_interval", "Time between chromatin contacts", x_label="association hazard (reduced)",
               y_label="free interval (s)", x_scale="log", y_scale="log",
               series=scan("free_interval_s", "free_interval_sem_s", "well_mixed_free_interval_s"),
               bands=[band("Sox2 tracking, 3.3–3.7 s", y=(low, high))], note=note),
        figure("bound_fraction", "Chromatin-bound fraction", x_label="association hazard (reduced)",
               y_label="bound fraction", x_scale="log",
               series=scan("bound_fraction", "bound_fraction_sem", "well_mixed_bound_fraction"),
               bands=[band("Sox2 FCS 23.5% (±20% assumed)", y=(0.8 * bound, 1.2 * bound))], note=note),
        figure("trials", "Nonspecific contacts per specific binding", x_label="association hazard (reduced)",
               y_label="contacts per specific binding", x_scale="log",
               series=scan("trials_per_specific_binding", "trials_sem"),
               bands=[band("Sox2 tracking, 84–97 trials", y=EXPERIMENT["trials"])],
               note="Set mainly by the assumed share of specific sites, one bead in specific_every."),
    ]


def analyze(params, results):
    off_rate, factor = _rates(params)
    unit = params["time_s_per_unit"]
    checks, findings, table = [], [], []
    for case in cases(params):
        done = [r["observables"] for r in results if r.get("complete") and r["case"]["id"] == case["id"]]
        if not done:
            continue
        row = {"case": case["id"], "on_rate": case["on_rate"], "diffusion": case["diffusion"], "replicas": len(done)}
        for kind, expected in (("specific", 1 / off_rate), ("nonspecific", 1 / (off_rate * factor))):
            value, sem = mean_sem([o["dwell_mean"][kind] for o in done if o["dwell_mean"][kind] is not None])
            row[f"{kind}_dwell_s"] = None if value is None else value * unit
            checks.append(check(f"{case['id']}_{kind}_dwell", f"{case['id']}: {kind} dwell", value, expected,
                                evidence=NUMERICAL, tolerance=0.08, sem=sem, unit="reduced time",
                                note="Dwell times are inputs; this confirms the release rule, including the class factor."))
        free, free_sem = mean_sem([o["free_interval_mean"] for o in done if o["free_interval_mean"] is not None])
        trials, trials_sem = mean_sem([o["binds"]["nonspecific"] / o["binds"]["specific"] for o in done
                                       if o["binds"]["specific"]])
        bound, bound_sem = mean_sem([sum(o["bound_fraction"].values()) for o in done if o["bound_fraction"]])
        ideal = _well_mixed(params, case, done[0]["sites"])
        row.update(free_interval_s=None if free is None else free * unit,
                   free_interval_sem_s=None if free_sem is None else free_sem * unit,
                   median_over_mean_free=mean_sem([o["free_interval_median"] / o["free_interval_mean"]
                                                   for o in done if o["free_interval_mean"]])[0],
                   trials_per_specific_binding=trials, trials_sem=trials_sem,
                   bound_fraction=bound, bound_fraction_sem=bound_sem,
                   specific_bound_fraction=mean_sem([o["bound_fraction"]["specific"] for o in done
                                                     if o["bound_fraction"]])[0],
                   well_mixed_free_interval_s=ideal["free_interval"] * unit,
                   well_mixed_bound_fraction=ideal["bound_fraction"], well_mixed_trials=ideal["trials"],
                   dropped_events=sum(o["dropped_events"] for o in done))
        table.append(row)
    usable = [row for row in table if row["free_interval_s"] is not None]
    if usable:
        target = sum(EXPERIMENT["free_interval_s"]) / 2
        best = min(usable, key=lambda row: abs(math.log(row["free_interval_s"] / target)))
        label = best["case"]
        low, high = EXPERIMENT["free_interval_s"]
        checks.append(check("free_interval", f"{label}: free interval between contacts", best["free_interval_s"],
                            target, evidence=EXPERIMENTAL, low=low, high=high, sem=best["free_interval_sem_s"], unit="s",
                            note="Selection criterion: the scanned association hazard closest to the measured interval. "
                                 "Tuned, not predicted."))
        low, high = EXPERIMENT["trials"]
        checks.append(check("trials", f"{label}: nonspecific contacts per specific binding",
                            best["trials_per_specific_binding"], sum(EXPERIMENT["trials"]) / 2, evidence=EXPERIMENTAL,
                            low=low, high=high, sem=best["trials_sem"],
                            note="Depends on the assumed share of specific sites (one bead in specific_every)."))
        checks.append(check("bound_fraction", f"{label}: chromatin-bound fraction", best["bound_fraction"],
                            EXPERIMENT["bound_fraction"], evidence=EXPERIMENTAL, tolerance=0.20,
                            sem=best["bound_fraction_sem"],
                            note="Independent FCS measurement; its uncertainty was not retrieved, +/-20% assumed."))
        findings.append(
            f"Closest scan point {label}: free interval {best['free_interval_s']:.2f} s "
            f"(dilute well-mixed estimate {best['well_mixed_free_interval_s']:.2f} s), bound fraction "
            f"{best['bound_fraction']:.3f} (estimate {best['well_mixed_bound_fraction']:.3f}), "
            + (f"{best['trials_per_specific_binding']:.0f}" if best["trials_per_specific_binding"] else "no")
            + f" contacts per specific binding (site ratio {best['well_mixed_trials']:.0f}).")
        findings.append("Median/mean of the free interval is 0.69 for a memoryless search; lower values indicate "
                        "fast rebinding to nearby chromatin: " +
                        ", ".join(f"{row['case']} {row['median_over_mean_free']:.2f}" for row in usable
                                  if row["median_over_mean_free"]) + ".")
        if any(row["dropped_events"] for row in table):
            findings.append("Some engine events were dropped from the bounded buffer; statistics are biased. "
                            "Reduce 'block'.")
    return {"checks": checks, "findings": findings, "figures": _figures(table),
            "summary": {"cases": table, "experiment": EXPERIMENT, "mapping": PROVISIONAL_CLOCK,
                        "inputs": {"off_rate_reduced": off_rate, "nonspecific_release_factor": factor}}}
