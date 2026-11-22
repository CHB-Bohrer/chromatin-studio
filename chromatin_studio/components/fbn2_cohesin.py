"""DNA plus cohesin (plus CTCF-like barriers) against Fbn2 two-locus imaging.

Experiment (Gabriele et al. 2022): two fluorescent labels 515 kb apart at the
anchors of the Fbn2 domain in mouse ES cells. Their 3D separation shrinks from
RAD21-depleted cells (no cohesin) to CTCF-depleted cells (cohesin, no anchors)
to wild type (cohesin and anchors). The fully looped state is rare (about
3-6.5% of the time) and short-lived (median about 10-30 min). The processed
distance statistics for every condition are in references/hansen.

Model: the bare 2,000 x 10 kb polymer of the baseline calibration, plus the
literature-derived cohesin numbers of the provisional model (129 extruders per
20 Mb, 22 min residence, 0.5 kb/s) and impermeable convergent barriers. Nothing
is fitted here. Three cases:

    no_cohesin        analog of RAD21 depletion
    cohesin           analog of CTCF depletion (extrusion without barriers)
    cohesin_barriers  analog of wild type (a periodic array of 510 kb domains)

The headline comparison is the ratio of mean-square anchor separations between
cases, which does not depend on the length mapping. Times use the provisional
clock (needs_review). The barrier array is an idealization: real CTCF sites are
irregular, partly occupied and can release cohesin.

The 1D cohesin lattice is advanced alone before the 3D run so loops start near
their steady state; the polymer then relaxes during the burn-in. Warm-start
conformations (equilibrated bare polymer) are run data and are passed by path.
"""

import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

from ..engine import Simulation
from .core import EXPERIMENTAL, PROVISIONAL_CLOCK, band, check, error_bar, figure, mean_sem, series, source, validated

REPO = Path(__file__).resolve().parents[2]
BASELINE_MODEL = REPO / "references" / "hansen" / "jobs" / "baseline-20mb-10kb-model.json"
REFERENCE_REPORT = REPO / "references" / "hansen" / "analysis-report.json"

STUDY = {
    "id": "fbn2_cohesin",
    "title": "Cohesin and barriers versus Fbn2 anchor separation",
    "components": ["chromatin mechanics", "cohesin extrusion", "directional barriers"],
    "question": "Do literature cohesin numbers reproduce how the Fbn2 anchor separation changes with cohesin and CTCF?",
    "evidence": [EXPERIMENTAL],
    "sources": [
        source("Gabriele et al. (2022) Science 376:496", "https://doi.org/10.1126/science.abn6583",
               "Fbn2 anchor separations in 12 conditions (data: https://zenodo.org/records/5770531, CC-BY-4.0); "
               "looped fraction about 3-6.5%, median loop lifetime about 10-30 min."),
        source("Hansen et al. (2017) eLife 6:e25776", "https://doi.org/10.7554/eLife.25776",
               "Cohesin residence about 22 min (input)."),
        source("Cattoglio et al. (2019) eLife 8:e40164", "https://doi.org/10.7554/eLife.40164",
               "About 86,900 RAD21 per mouse ES cell (input to the extruder count)."),
        source("Kim et al. (2019) Science 366:1345", "https://doi.org/10.1126/science.aaz4475",
               "Loop extrusion at about 0.5 kb/s in vitro (input)."),
    ],
    "scales": {
        "smoke": {"beads": 240, "radius": 9.9, "domain_start": 20, "domain_period": 60, "domain_length": 25,
                  "cohesin_count": 16, "step_rate": 0.0051270739, "off_rate": 0.00015536588, "spring_k": 2.5,
                  "lattice_prerun": 4000.0, "steps": 3000, "burn_in_steps": 1000, "sample_every": 250,
                  "partial_every": 4, "loop_threshold": 3, "warm_starts": [], "replicas": 1, "seed": 71},
        "standard": {"beads": 2000, "radius": 20.0, "domain_start": 74, "domain_period": 100, "domain_length": 51,
                     "cohesin_count": 129, "step_rate": 0.0051270739, "off_rate": 0.00015536588, "spring_k": 2.5,
                     "lattice_prerun": 32000.0, "steps": 6000000, "burn_in_steps": 2000000, "sample_every": 10000,
                     "partial_every": 20, "loop_threshold": 3, "warm_starts": [], "replicas": 3, "seed": 71},
    },
}


def cases(params):
    return [{"id": "no_cohesin", "cohesin": False, "barriers": False},
            {"id": "cohesin", "cohesin": True, "barriers": False},
            {"id": "cohesin_barriers", "cohesin": True, "barriers": True}]


def _domains(params):
    starts = range(params["domain_start"], params["beads"] - params["domain_length"] - 1, params["domain_period"])
    return [(start, start + params["domain_length"]) for start in starts]


def build_config(params, case, seed):
    config = json.loads(BASELINE_MODEL.read_text(encoding="utf-8"))
    config["seed"] = int(seed)
    config["chromatin"]["beads"] = params["beads"]
    config["simulation"]["confinement_radius"] = params["radius"]
    if case["barriers"]:
        config["chromatin"]["barriers"] = [item for left, right in _domains(params)
                                           for item in ({"bead": left, "direction": "left"},
                                                        {"bead": right, "direction": "right"})]
    if case["cohesin"]:
        config["cohesin"] = {"enabled": True, "count": params["cohesin_count"], "step_rate": params["step_rate"],
                             "off_rate": params["off_rate"], "spring_k": params["spring_k"]}
    return validated(config)


def _effective_separation(sim, domains):
    """Shortest path between anchors along the backbone with cohesin shortcuts."""
    n = sim.n_chromatin
    rows = list(range(n - 1)) + [pair[0] for pair in sim.cohesin]
    cols = list(range(1, n)) + [pair[1] for pair in sim.cohesin]
    weights = [1.0] * len(rows)          # a cohesin link counts as one step, like a backbone bond
    graph = coo_matrix((weights, (rows, cols)), shape=(n, n)).tocsr()
    distance = dijkstra(graph, directed=False, indices=[left for left, _ in domains])
    return [int(distance[index, right]) for index, (_, right) in enumerate(domains)]


def simulate(params, case, seed, save_partial):
    sim = Simulation(build_config(params, case, seed))
    n = sim.n_chromatin
    warm = None
    if params["warm_starts"]:
        warm = params["warm_starts"][seed % len(params["warm_starts"])]
        positions = np.load(warm)
        if positions.shape != (n, 3):
            raise ValueError(f"Warm start {warm} has shape {positions.shape}, expected {(n, 3)}.")
        sim.positions[:n] = positions
    if case["cohesin"]:
        # 1D lattice only: loop positions do not depend on coordinates.
        for _ in range(int(params["lattice_prerun"] / 2.0)):
            sim._advance_cohesin(2.0)
    domains = _domains(params)
    left = np.array([d[0] for d in domains])
    right = np.array([d[1] for d in domains])
    observables = {"warm_start": warm, "domains": domains, "sample_time": params["sample_every"] * sim.dt,
                   "burn_in_samples": params["burn_in_steps"] // params["sample_every"],
                   "vectors": [], "effective_separation": [], "mean_loop": [], "rg": [], "steps_done": 0}
    for index in range(params["steps"] // params["sample_every"]):
        sim.step(params["sample_every"])
        x = sim.positions[:n]
        observables["vectors"].append(np.round(x[right] - x[left], 3).tolist())
        observables["rg"].append(round(float(np.sqrt(np.mean(np.sum((x - x.mean(axis=0)) ** 2, axis=1)))), 3))
        if case["cohesin"]:
            observables["effective_separation"].append(_effective_separation(sim, domains))
            observables["mean_loop"].append(round(float(np.mean([p[1] - p[0] for p in sim.cohesin])), 2))
        observables["steps_done"] = sim.step_count
        if (index + 1) % params["partial_every"] == 0:
            save_partial(observables)
    return observables


def _corrected_reference():
    """Noise-corrected mean-square separations (um^2) with crude intervals, from the repository data."""
    report = json.loads(REFERENCE_REPORT.read_text(encoding="utf-8"))
    stats = {c["condition"]: c["analysis"]["distance_statistics"] for c in report["conditions"]}
    base = stats["Rad21_2_hr"]
    offset = base["mean_squared_distance_um2"] - base["noise_corrected_mean_squared_distance_um2"]

    def corrected(name):
        s = stats[name]
        low, high = s["mean_squared_distance_95pct_interval_um2"]
        return {"value": s["mean_squared_distance_um2"] - offset, "low": low - offset, "high": high - offset,
                "raw_mean_distance_um": s["mean_distance_um"]}
    return {"noise_offset_um2": offset, "no_cohesin": corrected("Rad21_2_hr"),
            "cohesin": corrected("CTCF_4_hr"), "cohesin_barriers": corrected("C36")}


def _run_lengths(flags):
    lengths, run = [], 0
    for flag in list(flags) + [False]:
        if flag:
            run += 1
        elif run:
            lengths.append(run)
            run = 0
    return lengths


def analyze(params, results):
    reference = _corrected_reference()
    length_um = PROVISIONAL_CLOCK["length_nm_per_unit"] / 1000.0
    clock = PROVISIONAL_CLOCK["time_s_per_unit"]
    table, per_case = {}, {}
    for case in cases(params):
        r2, looped, lifetimes, loops, used = [], [], [], [], 0
        for record in results:
            if record["case"]["id"] != case["id"]:
                continue
            o = record["observables"]
            vectors = np.asarray(o["vectors"], dtype=float)[o["burn_in_samples"]:]
            if not len(vectors):
                continue
            used += 1
            r2.append(float(np.mean(np.sum(vectors ** 2, axis=2))))
            if o["effective_separation"]:
                bridged = np.asarray(o["effective_separation"])[o["burn_in_samples"]:] <= params["loop_threshold"]
                looped.append(float(bridged.mean()))
                for column in bridged.T:
                    lifetimes.extend(length * o["sample_time"] for length in _run_lengths(column))
                loops.append(float(np.mean(o["mean_loop"][o["burn_in_samples"]:])))
        mean, sem = mean_sem(r2)
        per_case[case["id"]] = (mean, sem)
        table[case["id"]] = {
            "replicas_with_samples": used, "mean_square_separation_reduced": mean, "sem": sem,
            "mean_square_separation_um2": None if mean is None else mean * length_um ** 2,
            "experiment_corrected_um2": reference[case["id"]],
            "looped_fraction": mean_sem(looped)[0], "looped_fraction_sem": mean_sem(looped)[1],
            "median_loop_lifetime_min": float(np.median(lifetimes)) * clock / 60.0 if lifetimes else None,
            "loop_episodes": len(lifetimes), "mean_cohesin_loop_kb": mean_sem(loops)[0] * 10 if loops else None}
    checks, findings = [], []
    base, base_sem = per_case.get("no_cohesin", (None, None))
    if base:
        target = reference["no_cohesin"]
        checks.append(check("bare_separation", "No cohesin: mean-square anchor separation", base * length_um ** 2,
                            target["value"], evidence=EXPERIMENTAL, low=target["low"], high=target["high"],
                            sem=base_sem * length_um ** 2 if base_sem else None, unit="um^2",
                            note="RAD21-depleted cells, noise corrected. Uses the provisional length (needs_review), "
                                 "which was fitted to this same condition: a consistency check, not a prediction."))
        for name, label in (("cohesin", "Cohesin without barriers / no cohesin (CTCF-depleted / RAD21-depleted)"),
                            ("cohesin_barriers", "Cohesin with barriers / no cohesin (wild type / RAD21-depleted)")):
            value, sem = per_case.get(name, (None, None))
            if not value:
                continue
            ratio = value / base
            ratio_sem = ratio * ((sem / value) ** 2 + (base_sem / base) ** 2) ** 0.5 if sem and base_sem else None
            target, bare = reference[name], reference["no_cohesin"]
            checks.append(check(f"{name}_ratio", label, ratio, target["value"] / bare["value"], evidence=EXPERIMENTAL,
                                low=target["low"] / bare["high"], high=target["high"] / bare["low"], sem=ratio_sem,
                                note="Ratio of mean-square separations: independent of the length mapping. The same "
                                     "localization-noise offset is assumed for every condition. Parameter-free prediction."))
    wild = table.get("cohesin_barriers", {})
    if wild.get("looped_fraction") is not None:
        checks.append(check("looped_fraction", "Fully looped fraction with barriers", wild["looped_fraction"], 0.05,
                            evidence=EXPERIMENTAL, low=0.03, high=0.065, sem=wild["looped_fraction_sem"],
                            note=f"Looped = anchors within {params['loop_threshold']} beads along backbone plus cohesin "
                                 "shortcuts. The experiment infers looping from 3D distance, a different definition."))
        if wild["median_loop_lifetime_min"] is not None:
            checks.append(check("loop_lifetime", "Median looped-state lifetime", wild["median_loop_lifetime_min"], 20.0,
                                evidence=EXPERIMENTAL, low=10.0, high=30.0, unit="min",
                                note="Provisional clock (needs_review). Episodes cut by the run length bias this low."))
    partial = sum(not r.get("complete") for r in results)
    if partial:
        findings.append(f"{partial} result(s) are partial: values will change as the runs continue.")
    if any(not r["observables"].get("warm_start") for r in results):
        findings.append("At least one run started from a freshly grown chain rather than an equilibrated warm start; "
                        "its burn-in is unlikely to be sufficient.")
    names = {"no_cohesin": "No cohesin (RAD21-depleted)", "cohesin": "Cohesin, no CTCF (CTCF-depleted)",
             "cohesin_barriers": "Cohesin and CTCF (wild type)"}
    figures = []
    plotted = [(names[key], row) for key, row in table.items() if row["mean_square_separation_um2"] is not None]
    if plotted:
        scale = length_um ** 2
        bars = [error_bar(row["mean_square_separation_um2"], row["sem"] * scale if row["sem"] else None)
                for _, row in plotted]
        figures.append(figure(
            "anchor_separation", "Fbn2 anchor separation by condition", x_label="condition",
            y_label="mean-square 3D separation (µm²)", categories=list(names.values()),
            series=[series("Engine, mean ± SEM", [name for name, _ in plotted],
                           [row["mean_square_separation_um2"] for _, row in plotted],
                           low=[b[0] for b in bars], high=[b[1] for b in bars]),
                    series("Live-cell imaging (noise corrected)", list(names.values()),
                           [reference[key]["value"] for key in names], role="experiment")],
            bands=[band("Experimental 95% interval", category=names[key],
                        y=(reference[key]["low"], reference[key]["high"])) for key in names],
            note="Provisional length mapping (needs_review), fitted to the RAD21-depleted condition; "
                 "the other two conditions are parameter-free predictions."))
    return {"checks": checks, "findings": findings, "figures": figures,
            "summary": {"cases": table, "reference": reference, "mapping": PROVISIONAL_CLOCK,
                        "domains": len(_domains(params))}}
