"""HP1α residence on chromatin arrays as a function of H3K9me3 density.

Experiment (Kilic et al. 2015, Nat Commun 6:7313; single-molecule TIRF on
12 x 601 nucleosome arrays, 1 nM HP1α): the fast dissociation phase lasts
0.13 s on an H3K9me3 mononucleosome, 0.15/0.17/0.19 s on arrays with 25/50/75%
H3K9me3 and 0.25 s on fully methylated arrays; unmodified arrays give 0.08 s
and almost no binding. The association rate falls with lower H3K9me3 density.
Adding unlabelled HP1α (50, 250, 1000 nM) shortens the residence to 0.19, 0.13
and 0.11 s. The authors attribute the density effect to rapid rebinding at
neighbouring methylated nucleosomes.

Model: one 12-bead array (one bead per nucleosome) and one tagged monovalent
reader that binds H3K9me3 beads with hazard ``on_rate`` and unmodified beads
with a small affinity. The mononucleosome is one marked bead on a chain of
inert spacers (the schema needs at least 8 beads). Bonds of the tagged reader
separated by an unbound gap shorter than ``merge_gap`` form one residence
episode, as rebinding faster than a camera frame stays within one TIRF spot;
longer gaps count as new arrivals. A time criterion, unlike a distance
criterion, does not depend on how much of the confinement the array fills.
Unmodified beads release faster through the optional
``mark_off_rate_multiplier`` (input from the unmodified-array residence).
Cases where binding is rare (mononucleosome, unmodified array) run longer.

Inputs: the per-bond release rate (sets the time unit through the
mononucleosome residence, 0.13 s) and the unmodified-bead factors. One
parameter is tuned: the association hazard, by the scan point whose
fully-methylated/mononucleosome residence ratio is closest to 0.25/0.13.
Predictions: residence and association at 25-75% density, and residence with
competitor occupancy (unlabelled readers), interpolated at the occupancy
implied by the reported apparent K_D of 1.16 µM.
"""

import math

import numpy as np

from ..engine import Simulation
from .core import (EXPERIMENTAL, NUMERICAL, EventReader, band, bare_config, check, error_bar, figure, mean_sem,
                   protein, series, source, validated)

DOI = "https://doi.org/10.1038/ncomms8313"
ARRAY = 12
# Table 1 of Kilic et al. 2015 (fast phase tau_off,1 in s; k_on in 1e6 / M / s).
RESIDENCE = {"mn": (0.13, 0.01), "f25": (0.15, 0.02), "f50": (0.17, 0.03), "f75": (0.19, 0.01),
             "f100": (0.25, 0.03), "me0": (0.08, 0.03)}
ASSOCIATION = {"f25": (1.10, 0.33), "f50": (1.73, 0.55), "f75": (1.63, 0.82), "f100": (3.64, 1.56)}
COMPETITION = {50: (0.19, 0.05), 250: (0.13, 0.02), 1000: (0.11, 0.01)}     # nM unlabelled -> tau_off,1
KD_NM = 1160.0                                                              # apparent, fast phase
DENSITY = {"mn": None, "f25": 0.25, "f50": 0.5, "f75": 0.75, "f100": 1.0, "me0": 0.0}
LABELS = {"mn": "H3K9me3 mononucleosome", "f25": "25% H3K9me3", "f50": "50% H3K9me3", "f75": "75% H3K9me3",
          "f100": "100% H3K9me3", "me0": "Unmodified array"}
NO_REBINDING = "Mononucleosome level (no neighbours to rebind)"

STUDY = {
    "id": "hp1_mark_density",
    "title": "HP1α residence versus H3K9me3 density (mark-dependent binding)",
    "components": ["protein diffusion", "mark-dependent association", "mark-dependent release",
                   "rebinding within a nucleosome array"],
    "question": "With one tuned association hazard, does rebinding among methylated nucleosomes reproduce how "
                "HP1α residence and association depend on H3K9me3 density and on competing HP1α?",
    "evidence": [NUMERICAL, EXPERIMENTAL],
    "sources": [source("Kilic et al. (2015) Nat Commun 6:7313", DOI,
                       "HP1α fast-phase residence 0.13 s (H3K9me3 mononucleosome), 0.15/0.17/0.19/0.25 s on arrays "
                       "with 25/50/75/100% H3K9me3, 0.08 s unmodified; k_on 1.10/1.73/1.63/3.64 x 10^6 M^-1 s^-1 "
                       "(25-100%); with 50/250/1000 nM unlabelled HP1α 0.19/0.13/0.11 s; apparent K_D 1.16 µM.")],
    "scales": {
        "smoke": {"radius": 5.0, "diffusion": 1.0, "off_rate": 0.5, "on_rates": [2.0], "unmodified_affinity": 0.05,
                  "unmodified_release": 1.6, "merge_gap": 1.0, "competitor_targets": [0.2],
                  "dt": 0.01, "block": 10, "time": 300.0, "time_factor": {"mn": 2, "me0": 1},
                  "replicas": 1, "seed": 91},
        "standard": {"radius": 5.0, "diffusion": 1.0, "off_rate": 0.5, "on_rates": [0.5, 1.0, 2.0, 4.0],
                     "unmodified_affinity": 0.05, "unmodified_release": 1.6, "merge_gap": 1.0,
                     "competitor_targets": [0.06, 0.25, 0.7], "dt": 0.01, "block": 10, "time": 5000.0,
                     "time_factor": {"mn": 4, "me0": 2}, "replicas": 3, "seed": 91},
    },
}


def _capture(params):
    return 0.28 + 1.0          # reader radius + bond length: the engine's capture distance


def _competitors(params, on_rate, target):
    """Rough count giving the target occupancy (reaction-limited estimate; the
    analysis uses the measured occupancy, so this only needs to bracket)."""
    per_bead = on_rate * (4 / 3 * math.pi * _capture(params) ** 3) / (4 / 3 * math.pi * params["radius"] ** 3)
    free = target / (1 - target) * params["off_rate"] / per_bead
    return max(1, int(round(free + target * ARRAY)))


def cases(params):
    rows, factor = [], params.get("time_factor") or {}
    for rate in params["on_rates"]:
        for key in DENSITY:
            rows.append({"id": f"h{rate:g}-{key}", "on_rate": rate, "density": key, "competitors": 0,
                         "time": params["time"] * factor.get(key, 1)})
        for target in params["competitor_targets"]:
            rows.append({"id": f"h{rate:g}-comp{target:g}", "on_rate": rate, "density": "f100",
                         "competitors": _competitors(params, rate, target), "target_occupancy": target,
                         "time": params["time"]})
    return rows


def merge_episodes(bonds, open_bind, merge_gap, start, end):
    """Residence episodes and arrivals of one reader between ``start`` and ``end``.

    ``bonds`` are time-ordered complete (bind, unbind) pairs and ``open_bind``
    the start of a bond still open at ``end`` (or None). Bonds separated by an
    unbound gap shorter than ``merge_gap`` join one episode. Returns the
    durations of the complete episodes, the number of arrivals (episode starts)
    and the search time at risk of an arrival: unbound time outside the first
    ``merge_gap`` after each episode, when a binding would still have merged.
    Arrivals / search time is then the Poisson rate estimate."""
    episodes = []
    for bind, unbind in bonds:
        if episodes and bind - episodes[-1][1] < merge_gap:
            episodes[-1][1] = unbind
        else:
            episodes.append([bind, unbind])
    starts, complete = [e[0] for e in episodes], episodes
    if open_bind is not None:
        if episodes and open_bind - episodes[-1][1] < merge_gap:
            complete = episodes[:-1]                      # the last episode is still running
        else:
            starts.append(open_bind)
    search, previous = 0.0, None
    for index, began in enumerate(starts):
        search += began - (start if previous is None else previous + merge_gap)
        previous = episodes[index][1] if index < len(episodes) else None
    if open_bind is None:                                 # censored final search
        search += max(0.0, end - (episodes[-1][1] + merge_gap if episodes else start))
    return [e[1] - e[0] for e in complete], len(starts), search


def _marked(case, seed):
    """Methylated bead indices. Partial densities are random per replica, like a
    statistical reconstitution; a separate generator keeps the engine RNG untouched."""
    key = case["density"]
    if key == "mn":
        return [ARRAY // 2]
    count = int(round(ARRAY * DENSITY[key]))
    return sorted(np.random.default_rng(seed + 7919).choice(ARRAY, size=count, replace=False).tolist())


def build_config(params, case, seed):
    key, marked = case["density"], _marked(case, seed)
    # binding_capacity 1: a bound reader blocks its nucleosome, so unlabelled
    # readers compete for sites (the basis of the competition experiment).
    config = bare_config(ARRAY, seed=seed, dt=params["dt"], radius=params["radius"], bp_per_bead=177,
                         representation="nucleosome_unit", bond_length=1.0, bond_k=30.0, diffusion=0.02,
                         persistence_length=0.5, repulsion=0.0, binding_capacity=1)
    config["simulation"]["confinement_stiffness"] = 10.0
    config["marks"] = [{"id": "H3K9me3", "name": "H3K9me3", "color": "#b55d7a", "decay_rate": 0, "spread_rate": 0,
                        "contact_radius": 1.0}]
    config["chromatin"]["initial_marks"] = [{"mark": "H3K9me3", "start": b, "end": b} for b in marked]
    # The schema needs at least 8 beads. For the mononucleosome the other beads
    # are inert spacers (no affinity): a single nucleosome has no neighbours.
    unmodified = 0.0 if key == "mn" else params["unmodified_affinity"]
    reader = dict(diffusion=params["diffusion"], valence=1, on_rate=case["on_rate"], off_rate=params["off_rate"],
                  chromatin_affinity={"unmodified": unmodified, "H3K9me3": 1.0},
                  mark_off_rate_multiplier={"unmodified": params["unmodified_release"], "H3K9me3": 1.0})
    config["proteins"] = [protein("hp1", 1, color="#c05c86", **reader)]
    if case["competitors"]:
        config["proteins"].append(protein("hp1_unlabelled", case["competitors"], color="#8a8f99", **reader))
    return validated(config), marked


def simulate(params, case, seed, save_partial):
    config, marked = build_config(params, case, seed)
    sim = Simulation(config)
    tagged = next(i for i, m in enumerate(sim.molecules) if m["species"] == "hp1")
    marked_set, reader = set(marked), EventReader()
    start = float(sim.time)
    blocks = int(round(case.get("time", params["time"]) / (params["block"] * params["dt"])))
    bonds, bound = [], None                    # bonds: (bind, unbind, on a marked bead)
    occupancy, samples = 0.0, 0
    for index in range(blocks):
        sim.step(params["block"])
        for event in reader.fresh(sim):
            if event.get("molecule") != tagged:
                continue
            if event["type"] == "chromatin_bind":
                bound = (event["time"], event["bead"] in marked_set)
            elif event["type"] == "chromatin_unbind" and bound is not None:
                bonds.append((bound[0], event["time"], bound[1]))
                bound = None
        if case["competitors"] and marked:
            occupancy += sum(1 for bond in sim.bindings if bond["kind"] == "chromatin" and bond["a"] in marked_set
                             and sim._particle_info[bond["b"]]["species"] == "hp1_unlabelled") / len(marked)
            samples += 1
        if index and index % 50000 == 0:
            save_partial({"bonds": len(bonds)})
    durations, arrivals, search = merge_episodes([(a, b) for a, b, _ in bonds], None if bound is None else bound[0],
                                                 params["merge_gap"], start, float(sim.time))
    marked_dwell = [b - a for a, b, on_mark in bonds if on_mark]
    plain_dwell = [b - a for a, b, on_mark in bonds if not on_mark]
    return {"marked": marked, "episodes": len(durations), "arrivals": arrivals,
            "residence": float(np.mean(durations)) if durations else None,
            "bonds_per_episode": len(bonds) / arrivals if arrivals else None,
            "association_rate": arrivals / search if search > 0 else None,
            "bond_dwell_marked": float(np.mean(marked_dwell)) if marked_dwell else None,
            "bond_count_marked": len(marked_dwell),
            "bond_dwell_unmodified": float(np.mean(plain_dwell)) if plain_dwell else None,
            "bond_count_unmodified": len(plain_dwell),
            "competitor_occupancy": occupancy / samples if samples else 0.0,
            "dropped_events": int(reader.dropped), "max_substeps": int(sim.last_substeps)}


def _ratio_interval(a, b):
    """Measured ratio a/b with first-order uncertainty from (value, error) pairs."""
    value = a[0] / b[0]
    error = value * math.sqrt((a[1] / a[0]) ** 2 + (b[1] / b[0]) ** 2)
    return value, value - error, value + error


def _stats(results, case_id, key):
    return mean_sem([r["observables"][key] for r in results if r.get("complete") and r["case"]["id"] == case_id
                     and r["observables"][key] is not None])


def analyze(params, results):
    checks, findings, scan = [], [], []
    for rate in params["on_rates"]:
        base = _stats(results, f"h{rate:g}-mn", "residence")
        full = _stats(results, f"h{rate:g}-f100", "association_rate")
        if base[0] is None:
            continue
        row = {"on_rate": rate, "mononucleosome_residence": base[0], "residence_ratio": {}, "association_ratio": {},
               "competition": [],
               "bonds_per_episode": {key: _stats(results, f"h{rate:g}-{key}", "bonds_per_episode")[0]
                                     for key in DENSITY}}
        for key in DENSITY:
            value = _stats(results, f"h{rate:g}-{key}", "residence")
            if value[0] is not None:
                row["residence_ratio"][key] = (value[0] / base[0],
                                               value[0] / base[0] * math.sqrt((value[1] or 0) ** 2 / value[0] ** 2
                                                                              + (base[1] or 0) ** 2 / base[0] ** 2))
            association = _stats(results, f"h{rate:g}-{key}", "association_rate")
            if association[0] is not None and full[0]:
                row["association_ratio"][key] = (association[0] / full[0], None)
        for target in params["competitor_targets"]:
            case_id = f"h{rate:g}-comp{target:g}"
            value = _stats(results, case_id, "residence")
            theta = _stats(results, case_id, "competitor_occupancy")
            if value[0] is not None:
                row["competition"].append((theta[0], value[0] / base[0]))
        scan.append(row)
    target = RESIDENCE["f100"][0] / RESIDENCE["mn"][0]
    usable = [row for row in scan if "f100" in row["residence_ratio"]]
    marked_dwell = mean_sem([r["observables"]["bond_dwell_marked"] for r in results if r.get("complete")
                             and r["observables"]["bond_dwell_marked"] is not None])
    if marked_dwell[0] is not None:
        checks.append(check("bond_release_marked", "Bond lifetime on H3K9me3 beads", marked_dwell[0],
                            1 / params["off_rate"], evidence=NUMERICAL, tolerance=0.05, sem=marked_dwell[1],
                            unit="reduced time", note="Input release rate; rebinding does not enter single bonds."))
    plain = [r["observables"] for r in results if r.get("complete") and r["case"]["density"] == "me0"
             and r["observables"]["bond_count_unmodified"] >= 20]
    if plain:
        value, sem = mean_sem([o["bond_dwell_unmodified"] for o in plain])
        checks.append(check("bond_release_unmodified", "Bond lifetime on unmodified beads (mark-dependent release)",
                            value, 1 / (params["off_rate"] * params["unmodified_release"]), evidence=NUMERICAL,
                            tolerance=0.08, sem=sem, unit="reduced time",
                            note="Tests the optional mark_off_rate_multiplier."))
    best = None
    if usable:
        best = min(usable, key=lambda row: abs(math.log(row["residence_ratio"]["f100"][0] / target)))
        label = f"association hazard {best['on_rate']:g}"
        ratio, low, high = _ratio_interval(RESIDENCE["f100"], RESIDENCE["mn"])
        checks.append(check("residence_f100", f"{label}: residence, 100% / mononucleosome",
                            best["residence_ratio"]["f100"][0], ratio, evidence=EXPERIMENTAL, low=low, high=high,
                            sem=best["residence_ratio"]["f100"][1],
                            note="Tuned: the scanned hazard closest to the measured ratio is selected."))
        for key in ("f75", "f50", "f25"):
            if key in best["residence_ratio"]:
                ratio, low, high = _ratio_interval(RESIDENCE[key], RESIDENCE["mn"])
                checks.append(check(f"residence_{key}", f"{label}: residence, {LABELS[key]} / mononucleosome",
                                    best["residence_ratio"][key][0], ratio, evidence=EXPERIMENTAL, low=low,
                                    high=high, sem=best["residence_ratio"][key][1],
                                    note="Prediction of rebinding among methylated neighbours (random placement, "
                                         "averaged over replicas)."))
            if key in best["association_ratio"]:
                ratio, low, high = _ratio_interval(ASSOCIATION[key], ASSOCIATION["f100"])
                checks.append(check(f"association_{key}", f"{label}: association, {LABELS[key]} / 100%",
                                    best["association_ratio"][key][0], ratio, evidence=EXPERIMENTAL, low=low,
                                    high=high,
                                    note="Prediction. Array-level episodes per unit search time; measured k_on "
                                         "uncertainties are large."))
        if "me0" in best["residence_ratio"]:
            ratio, low, high = _ratio_interval(RESIDENCE["me0"], RESIDENCE["mn"])
            checks.append(check("residence_me0", f"{label}: residence, unmodified array / mononucleosome",
                                best["residence_ratio"]["me0"][0], ratio, evidence=EXPERIMENTAL, low=low, high=high,
                                sem=best["residence_ratio"]["me0"][1],
                                note="Input-matched: the unmodified release factor was set from this ratio, "
                                     "assuming no rebinding; rebinding among unmodified beads raises it."))
        if "me0" in best["association_ratio"]:
            checks.append(check("association_me0", f"{label}: association, unmodified / 100%",
                                best["association_ratio"]["me0"][0], 0.0, evidence=EXPERIMENTAL, low=0.0, high=0.23,
                                note="Measured only as a bound (λon > 100 s vs 22.9 s). Consistency of the assumed "
                                     "unmodified affinity."))
        points = sorted([(0.0, best["residence_ratio"]["f100"][0])] + [p for p in best["competition"] if p[0] is not None])
        for concentration, measured in COMPETITION.items():
            theta = concentration / (concentration + KD_NM)
            if len(points) < 2 or theta > points[-1][0]:
                continue
            predicted = float(np.interp(theta, [p[0] for p in points], [p[1] for p in points]))
            ratio, low, high = _ratio_interval(measured, RESIDENCE["mn"])
            checks.append(check(f"competition_{concentration}", f"{label}: residence with {concentration} nM "
                                f"unlabelled HP1α / mononucleosome", predicted, ratio, evidence=EXPERIMENTAL,
                                low=low, high=high,
                                note=f"Prediction at competitor occupancy {theta:.2f} (= c/(c + 1.16 µM)), "
                                     "interpolated between simulated occupancies. HP1α dimerization at high "
                                     "concentration is not modelled."))
        findings.append("Hazard scan, residence 100% / mononucleosome (target "
                        f"{target:.2f}): " + ", ".join(f"{row['on_rate']:g} → {row['residence_ratio']['f100'][0]:.2f}"
                                                       for row in usable) + ".")
        findings.append(f"Selected association hazard {best['on_rate']:g}: residence ratios "
                        + ", ".join(f"{LABELS[k]} {v[0]:.2f}" for k, v in best["residence_ratio"].items()) + ".")
    findings.append(f"Bonds separated by unbound gaps shorter than {params['merge_gap']:g} reduced time count as "
                    "one residence episode (rebinding faster than a camera frame).")
    return {"checks": checks, "findings": findings, "figures": _figures(scan, best),
            "summary": {"scan": [{k: v for k, v in row.items()} for row in scan],
                        "selected_on_rate": best["on_rate"] if best else None,
                        "experiment": {"residence_s": RESIDENCE, "association": ASSOCIATION,
                                       "competition_s": {str(k): v for k, v in COMPETITION.items()},
                                       "kd_nM": KD_NM}}}


def _figures(scan, best):
    figures, densities = [], ["f25", "f50", "f75", "f100"]
    x = [DENSITY[k] * 100 for k in densities]
    tuned = [row for row in scan if "f100" in row["residence_ratio"]]
    if tuned:
        rates = [row["on_rate"] for row in tuned]
        values = [row["residence_ratio"]["f100"] for row in tuned]
        bars = [error_bar(*v) for v in values]
        target, low, high = _ratio_interval(RESIDENCE["f100"], RESIDENCE["mn"])
        figures.append(figure(
            "hazard_scan", "Choosing the association hazard", x_label="association hazard (reduced units)",
            y_label="residence, 100% array / mononucleosome",
            series=[series("Engine", rates, [v[0] for v in values], low=[b[0] for b in bars],
                           high=[b[1] for b in bars], style="connected"),
                    series(NO_REBINDING, [rates[0], rates[-1]], [1.0, 1.0], role="reference", style="line")],
            bands=[band("Kilic et al. 2015, ±1 propagated error", role="experiment", y=(low, high))],
            x_scale="log" if len(set(rates)) > 1 else "linear",
            note=f"The scan point closest to the measured ratio ({target:.2f}) is selected: this is the one tuned "
                 "comparison. Every other comparison uses the selected hazard."))
    if best:
        values = [best["residence_ratio"].get(k, (None, None)) for k in densities]
        bars = [error_bar(*v) for v in values]
        measured = [_ratio_interval(RESIDENCE[k], RESIDENCE["mn"]) for k in densities]
        figures.append(figure(
            "residence_density", "Residence grows with H3K9me3 density", x_label="H3K9me3 nucleosomes in the array (%)",
            y_label="residence / mononucleosome residence",
            series=[series(f"Engine, hazard {best['on_rate']:g} (selected)", x, [v[0] for v in values],
                           low=[b[0] for b in bars], high=[b[1] for b in bars], style="connected"),
                    series("Kilic et al. 2015", x, [m[0] for m in measured], low=[m[1] for m in measured],
                           high=[m[2] for m in measured], role="experiment"),
                    series(NO_REBINDING, [x[0], x[-1]], [1.0, 1.0], role="reference", style="line")],
            note="Error bars: engine ±1 SEM over replicas; measured ratios ±1 propagated error. The 100% point is "
                 "tuned; 25-75% are predictions."))
        values = [best["association_ratio"].get(k, (None, None))[0] for k in densities]
        measured = [_ratio_interval(ASSOCIATION[k], ASSOCIATION["f100"]) for k in densities]
        figures.append(figure(
            "association_density", "Association scales with methylated sites", x_label="H3K9me3 nucleosomes (%)",
            y_label="association rate / 100% array",
            series=[series("Proportional to methylated sites", [0, 100], [0, 1], role="reference", style="line"),
                    series(f"Engine, hazard {best['on_rate']:g}", x, values, style="connected"),
                    series("Kilic et al. 2015", x, [m[0] for m in measured], low=[m[1] for m in measured],
                           high=[m[2] for m in measured], role="experiment")]))
        points = sorted([(0.0, best["residence_ratio"]["f100"][0])]
                        + [p for p in best["competition"] if p[0] is not None])
        measured = [(c / (c + KD_NM), _ratio_interval(v, RESIDENCE["mn"])) for c, v in COMPETITION.items()]
        span = [0.0, max([p[0] for p in points] + [m[0] for m in measured])]
        figures.append(figure(
            "competition", "Competing HP1α shortens residence", x_label="array occupancy by unlabelled HP1α",
            y_label="residence / mononucleosome residence",
            series=[series(f"Engine, hazard {best['on_rate']:g}", [p[0] for p in points], [p[1] for p in points],
                           style="connected"),
                    series("Kilic et al. 2015 (occupancy from K_D 1.16 µM)", [m[0] for m in measured],
                           [m[1][0] for m in measured], low=[m[1][1] for m in measured],
                           high=[m[1][2] for m in measured], role="experiment"),
                    series(NO_REBINDING, span, [1.0, 1.0], role="reference", style="line")],
            note="Occupancy of the measured points uses the reported apparent K_D, itself uncertain (±0.54 µM). "
                 "With every neighbour occupied, the tagged reader can rebind only to the site it left, so site "
                 "exclusion takes residence down towards the dashed mononucleosome level; clearly below it would "
                 "need a mechanism that is not modelled, such as facilitated dissociation."))
    return figures
