"""Transcription-factor kinetics at naked versus nucleosome-embedded sites.

Experiments (single-molecule TIRF on immobilized DNA or nucleosomes, Poirier lab):

- Donovan et al. 2019 (eLife): nucleosomes slow Reb1 association 53-fold, but
  also slow its dissociation 132-fold, so Reb1 binds nucleosomal and naked
  sites with similar ensemble affinity ("dissociation rate compensation").
- Donovan et al. 2023 (Mol Cell): Cbf1 (DNA-binding domain) association is
  about 270-fold slower and dissociation about 76-fold slower on nucleosomes;
  Pho4 association is about 100-fold slower while its dissociation is unchanged,
  so nucleosomes strongly inhibit Pho4 but not Cbf1.
- Donovan et al. 2023 (NAR): Gal4 at a site 31 bp into the nucleosome binds
  about 130-fold slower than at 26 bp, with dissociation only 1.3-fold faster.

Model: one chain whose binding sites are either naked (nucleosomes disabled)
or nucleosome occupied (no turnover), and one monovalent factor that binds
only marked site beads. Two release rules are compared for each factor:

- ``association_only``: the nucleosome changes association only
  (``nucleosome_affinity``), the factor's release is unchanged;
- ``association_and_release``: the measured change in dissociation is added
  (``nucleosome_off_rate_multiplier``), standing for histone contacts.

Both kinetic ratios are inputs (input-matched). What is tested: (1) the engine
realizes them in a 3D search (numerical); (2) the ensemble affinity ratio that
follows from the simulated occupancy, against the independently titrated
apparent affinities (S1/2) of the same papers (prediction, factor-of-2 rule).
Time is mapped per factor by its measured naked-site dissociation rate.
"""

import math

import numpy as np

from ..engine import Simulation
from .core import (EXPERIMENTAL, NUMERICAL, EventReader, band, bare_config, check, error_bar, figure, mean_sem,
                   protein, series, source, survival, validated)

DOI_2019 = "https://doi.org/10.7554/eLife.43008"
DOI_MOLCELL = "https://doi.org/10.1016/j.molcel.2023.03.006"
DOI_NAR = "https://doi.org/10.1093/nar/gkac1267"

# Measured values, quoted from the papers' text, figure legends or author response.
# kon/koff ratios are nucleosomal site / reference site. affinity = S1/2 ratio
# (nucleosomal / reference) from ensemble titrations; "low_only" marks a bound.
FACTORS = {
    "reb1": {"label": "Reb1", "reference_site": "naked DNA", "nucleosome_site": "601 nucleosome, site P8",
             "kon_ratio": 1 / 53, "koff_ratio": 1 / 132, "koff_reference_per_s": 0.58,
             "affinity_ratio": (0.47, 0.65), "doi": DOI_2019,
             "measured": "kon 53-fold and koff 132-fold slower on nucleosomes; DNA koff 0.58 s^-1; "
                         "S1/2 P8 nucleosome 1.5 nM vs DNA 2.3 nM (EMSA), 2.4 vs 5.1 nM (FRET/PIFE)"},
    "cbf1": {"label": "Cbf1 (ΔN)", "reference_site": "naked DNA", "nucleosome_site": "601 nucleosome, 8 bp from edge",
             "kon_ratio": 0.0006 / 0.16, "koff_ratio": 0.0051 / 0.39, "koff_reference_per_s": 0.39,
             "affinity_ratio": (2.7, 4.1), "doi": DOI_MOLCELL,
             "measured": "kon 0.16 vs 0.0006 s^-1 nM^-1, koff 0.39 vs 0.0051 s^-1 (DNA vs nucleosome); "
                         "S1/2 ratio 3.4 ± 0.7"},
    "pho4": {"label": "Pho4", "reference_site": "naked DNA", "nucleosome_site": "601 nucleosome, 8 bp from edge",
             "kon_ratio": 0.01, "koff_ratio": 1.0, "koff_reference_per_s": None,
             "affinity_ratio": (40.0, None), "doi": DOI_MOLCELL,
             "measured": "association about 2 orders of magnitude slower; dissociation the same from nucleosomes "
                         "and DNA; S1/2 ratio more than 40-fold (nucleosome 1.08 nM, DNA below 30 pM)"},
    "gal4": {"label": "Gal4 (P31 vs P26)", "reference_site": "nucleosomal site P26 (outer)",
             "nucleosome_site": "nucleosomal site P31 (inner)",
             "kon_ratio": 1 / 130, "koff_ratio": 1.3, "koff_reference_per_s": 0.66,
             "affinity_ratio": (500 / 4.3 * 0.8, 500 / 4.3 * 1.2), "doi": DOI_NAR,
             "measured": "kon,eff 0.11 vs 0.0008 nM^-1 s^-1 (130 ± 60-fold), dissociation 0.66 vs 0.83 s^-1 "
                         "(1.3 ± 0.1-fold faster); S1/2 4.3 ± 0.2 vs 500 ± 100 nM"},
}
MODELS = {"association_only": "Occlusion only (release unchanged)",
          "association_and_release": "Occlusion + measured release change"}

STUDY = {
    "id": "nucleosome_tf_kinetics",
    "title": "Transcription-factor kinetics with and without a nucleosome",
    "components": ["protein diffusion", "protein-chromatin association and release",
                   "nucleosome occupancy (static)"],
    "question": "Do measured nucleosome effects on association and dissociation (Donovan et al.) reproduce the "
                "measured nucleosome/DNA affinity ratios, and is occlusion alone enough?",
    "evidence": [NUMERICAL, EXPERIMENTAL],
    "sources": [
        source("Donovan et al. (2019) eLife 8:e43008", DOI_2019, FACTORS["reb1"]["measured"]),
        source("Donovan et al. (2023) Mol Cell 83:1251", DOI_MOLCELL,
               FACTORS["cbf1"]["measured"] + "; Pho4: " + FACTORS["pho4"]["measured"]),
        source("Donovan et al. (2023) Nucleic Acids Res 51:1139", DOI_NAR, FACTORS["gal4"]["measured"]),
    ],
    "scales": {
        "smoke": {"beads": 32, "site_every": 4, "radius": 4.0, "count": 12, "on_rate": 0.25, "off_rate": 0.5,
                  "diffusion": 2.0, "dt": 0.005, "burn_in": 20.0, "time": 200.0, "block": 20,
                  "factors": ["reb1", "pho4"], "max_samples": 2000, "replicas": 1, "seed": 81},
        "standard": {"beads": 32, "site_every": 4, "radius": 4.0, "count": 80, "on_rate": 0.25, "off_rate": 0.5,
                     "diffusion": 2.0, "dt": 0.01, "burn_in": 500.0, "time": 10000.0, "block": 25,
                     "factors": ["reb1", "cbf1", "pho4", "gal4"], "max_samples": 5000, "replicas": 3, "seed": 81},
    },
}


def cases(params):
    rows = [{"id": "reference", "substrate": "reference"}]
    for factor in params["factors"]:
        for model in MODELS:
            if model == "association_and_release" and FACTORS[factor]["koff_ratio"] == 1.0:
                continue        # identical to association_only for this factor
            rows.append({"id": f"{factor}-{model}", "substrate": "nucleosome", "factor": factor, "model": model})
    return rows


def _sites(params):
    return list(range(params["site_every"] // 2, params["beads"], params["site_every"]))


def build_config(params, case, seed):
    sites = _sites(params)
    # binding_capacity 1: one factor per site, as in the single-site titrations.
    chromatin = dict(bp_per_bead=200, bond_length=1.0, bond_k=30.0, diffusion=0.05, persistence_length=0.5,
                     repulsion=0.0, binding_capacity=1)
    if case["substrate"] == "nucleosome":
        chromatin["nucleosomes"] = {"enabled": True, "initial_occupancy": 1, "eviction_rate": 0, "reassembly_rate": 0,
                                    "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": 1}
    config = bare_config(params["beads"], seed=seed, dt=params["dt"], radius=params["radius"], **chromatin)
    config["simulation"]["confinement_stiffness"] = 10.0
    config["marks"] = [{"id": "site", "name": "Binding site", "color": "#d6a34a", "decay_rate": 0, "spread_rate": 0,
                        "contact_radius": 1.0}]
    config["chromatin"]["initial_marks"] = [{"mark": "site", "start": s, "end": s} for s in sites]
    fields = {}
    if case["substrate"] == "nucleosome":
        factor = FACTORS[case["factor"]]
        fields["nucleosome_affinity"] = {"occupied": factor["kon_ratio"], "vacant": 1.0}
        release = factor["koff_ratio"] if case["model"] == "association_and_release" else 1.0
        fields["nucleosome_off_rate_multiplier"] = {"occupied": release, "vacant": 1.0}
    config["proteins"] = [protein("tf", params["count"], diffusion=params["diffusion"], valence=1,
                                  chromatin_affinity={"unmodified": 0.0, "site": 1.0}, on_rate=params["on_rate"],
                                  off_rate=params["off_rate"], **fields)]
    return validated(config), sites


def simulate(params, case, seed, save_partial):
    config, sites = build_config(params, case, seed)
    sim = Simulation(config)
    site_set, reader = set(sites), EventReader()
    block_time = params["block"] * params["dt"]
    bound_at, dwell, sample = {}, [], []
    binds, exposure, occupied, free_molecules, blocks = 0, 0.0, 0.0, 0.0, 0
    blocks_total = int(round((params["burn_in"] + params["time"]) / block_time))
    for index in range(blocks_total):
        sim.step(params["block"])
        measuring = sim.time > params["burn_in"] + 1e-9
        for event in reader.fresh(sim):
            if event.get("bead") not in site_set:
                continue
            molecule = event.get("molecule")
            if event["type"] == "chromatin_bind":
                bound_at[molecule] = event["time"]
                binds += measuring and event["time"] >= params["burn_in"]
            elif event["type"] == "chromatin_unbind" and molecule in bound_at:
                start = bound_at.pop(molecule)
                if start >= params["burn_in"]:
                    dwell.append(event["time"] - start)
        if measuring:
            bonds = sum(1 for bond in sim.bindings if bond["kind"] == "chromatin" and bond["a"] in site_set)
            free = params["count"] - bonds
            blocks += 1
            occupied += bonds / len(sites)
            free_molecules += free
            exposure += (len(sites) - bonds) * free * block_time
        if index and index % 20000 == 0:
            save_partial({"dwell_count": len(dwell), "binds": binds})
    dwell = np.asarray(dwell)
    keep = dwell[:params["max_samples"]]
    return {"sites": len(sites), "binds": int(binds), "dwell_count": int(dwell.size),
            "dwell_mean": float(dwell.mean()) if dwell.size else None,
            "dwell_sample": [round(float(v), 4) for v in keep],
            "association_rate": binds / exposure if exposure else None,
            "occupancy": occupied / blocks if blocks else None,
            "free_molecules": free_molecules / blocks if blocks else None,
            "measured_time": blocks * block_time, "dropped_events": int(reader.dropped),
            "max_substeps": int(sim.last_substeps)}


def _per_case(results, case_id):
    return [r["observables"] for r in results if r.get("complete") and r["case"]["id"] == case_id]


def _ratio(numerator, denominator):
    """Ratio of two (mean, sem) pairs with first-order error propagation."""
    (a, sa), (b, sb) = numerator, denominator
    if a is None or b is None or b == 0:
        return None, None
    value = a / b
    if sa is None or sb is None:
        return value, None
    return value, abs(value) * math.sqrt((sa / a) ** 2 + (sb / b) ** 2) if a else None


def _dissociation_constant(o):
    """Apparent K in units of free molecules: F (1 - theta) / theta."""
    theta = o["occupancy"]
    return o["free_molecules"] * (1 - theta) / theta if theta and 0 < theta < 1 else None


def analyze(params, results):
    reference = _per_case(results, "reference")
    ref_kon = mean_sem([o["association_rate"] for o in reference])
    ref_dwell = mean_sem([o["dwell_mean"] for o in reference])
    ref_k = mean_sem([_dissociation_constant(o) for o in reference])
    checks, findings, rows = [], [], []
    for case in cases(params)[1:]:
        done = _per_case(results, case["id"])
        if not done or ref_kon[0] is None:
            continue
        factor = FACTORS[case["factor"]]
        release = factor["koff_ratio"] if case["model"] == "association_and_release" else 1.0
        kon = _ratio(mean_sem([o["association_rate"] for o in done]), ref_kon)
        koff = _ratio(ref_dwell, mean_sem([o["dwell_mean"] for o in done]))
        affinity = _ratio(mean_sem([_dissociation_constant(o) for o in done]), ref_k)
        kinetic = (release / factor["kon_ratio"], None)
        rows.append({"case": case["id"], "factor": case["factor"], "model": case["model"],
                     "kon_ratio": kon[0], "kon_ratio_sem": kon[1], "koff_ratio": koff[0], "koff_ratio_sem": koff[1],
                     "affinity_ratio": affinity[0], "affinity_ratio_sem": affinity[1],
                     "input_kon_ratio": factor["kon_ratio"], "input_koff_ratio": release,
                     "replicas": len(done), "binding_events": sum(o["binds"] for o in done),
                     "dropped_events": sum(o["dropped_events"] for o in done)})
        label = f"{factor['label']}, {MODELS[case['model']].lower()}"
        checks.append(check(f"{case['id']}_kon", f"{label}: association ratio realized", kon[0], factor["kon_ratio"],
                            evidence=NUMERICAL, tolerance=0.10, sem=kon[1],
                            note="Input multiplier. The realized ratio may fall a few % short: a factor released at a "
                                 "naked site sometimes rebinds it at once (geminate rebinding), which inflates the "
                                 "reference rate; at the slower nucleosomal site this is negligible."))
        checks.append(check(f"{case['id']}_koff", f"{label}: dissociation ratio realized", koff[0], release,
                            evidence=NUMERICAL, tolerance=0.08, sem=koff[1], note="Input multiplier."))
        checks.append(check(f"{case['id']}_balance", f"{label}: affinity from occupancy vs from kinetics",
                            affinity[0], kinetic[0], evidence=NUMERICAL, tolerance=0.15, sem=affinity[1],
                            note="Detailed balance: the occupancy-derived affinity ratio should equal "
                                 "(dissociation ratio)/(association ratio) of the inputs."))
        low, high = factor["affinity_ratio"]
        if high is None:
            checks.append(check(f"{case['id']}_affinity", f"{label}: nucleosome/reference affinity ratio",
                                affinity[0], low, evidence=EXPERIMENTAL, low=low, sem=affinity[1],
                                note="Prediction from the two kinetic inputs; the measured titration gives only a "
                                     "lower bound."))
        else:
            centre = math.sqrt(low * high)
            checks.append(check(f"{case['id']}_affinity", f"{label}: nucleosome/reference affinity ratio",
                                affinity[0], centre, evidence=EXPERIMENTAL, low=low / 2, high=high * 2,
                                sem=affinity[1],
                                note=f"Prediction from the two kinetic inputs versus independent titrations "
                                     f"(S1/2 ratio {low:.2g}-{high:.2g}). Different assays: agreement within a "
                                     "factor of 2 of that range is the criterion."))
    dropped = sum(o["dropped_events"] for o in (r["observables"] for r in results))
    if dropped:
        findings.append(f"{dropped} engine events were lost from the bounded buffer; kinetics are biased. Reduce 'block'.")
    by = {(row["factor"], row["model"]): row for row in rows}
    for key, factor in FACTORS.items():
        only, both = by.get((key, "association_only")), by.get((key, "association_and_release"))
        if only and only.get("affinity_ratio"):
            text = (f"{factor['label']}: occlusion alone gives an affinity ratio of {only['affinity_ratio']:.3g}")
            if both and both.get("affinity_ratio"):
                text += f"; with the measured release change, {both['affinity_ratio']:.3g}"
            low, high = factor["affinity_ratio"]
            text += f" (measured {low:.2g}" + (f"-{high:.2g})." if high else "+).")
            findings.append(text)
    return {"checks": checks, "findings": findings, "figures": _figures(params, results, rows, ref_dwell),
            "summary": {"cases": rows, "factors": FACTORS, "models": MODELS,
                        "reference": {"association_rate": ref_kon[0], "dwell": ref_dwell[0], "affinity": ref_k[0]}}}


def _figures(params, results, rows, ref_dwell):
    if not rows:
        return []
    labels = {key: FACTORS[key]["label"] for key in params["factors"]}
    names = [labels[key] for key in params["factors"]]
    figures = []
    measured_kon = series("Association, measured", names, [FACTORS[k]["kon_ratio"] for k in params["factors"]],
                          role="experiment", group="association")
    measured_koff = series("Dissociation, measured", names, [FACTORS[k]["koff_ratio"] for k in params["factors"]],
                           role="experiment", group="dissociation")
    engine = []
    for quantity, title, group in (("kon", "Association, engine", "association"),
                                   ("koff", "Dissociation, engine", "dissociation")):
        chosen = []
        for key in params["factors"]:
            pick = next((r for r in rows if r["factor"] == key and r["model"] == "association_and_release"), None) \
                or next((r for r in rows if r["factor"] == key), None)
            chosen.append(pick)
        bars = [error_bar(r[f"{quantity}_ratio"], r[f"{quantity}_ratio_sem"]) if r else (None, None) for r in chosen]
        engine.append(series(title, names, [r[f"{quantity}_ratio"] if r else None for r in chosen],
                             low=[b[0] for b in bars], high=[b[1] for b in bars], group=group))
    figures.append(figure(
        "kinetic_ratios", "Nucleosome effect on association and dissociation", x_label="factor",
        y_label="nucleosomal / reference rate", y_scale="log", categories=names,
        series=[measured_kon, engine[0], measured_koff, engine[1]],
        note="Both measured ratios are model inputs; the engine points show that the 3D search realizes them. "
             "Ratios below 1 mean slower on the nucleosome."))
    affinity_series = []
    for model, title in MODELS.items():
        chosen = [next((r for r in rows if r["factor"] == key and r["model"] == model), None)
                  or (next((r for r in rows if r["factor"] == key and r["model"] == "association_only"), None)
                      if model == "association_and_release" else None) for key in params["factors"]]
        bars = [error_bar(r["affinity_ratio"], r["affinity_ratio_sem"]) if r else (None, None) for r in chosen]
        affinity_series.append(series(f"Engine: {title}", names, [r["affinity_ratio"] if r else None for r in chosen],
                                      low=[b[0] for b in bars], high=[b[1] for b in bars]))
    bands = []
    for key in params["factors"]:
        low, high = FACTORS[key]["affinity_ratio"]
        bands.append(band("Measured S1/2 ratio", category=labels[key], y=(low, high)))   # high None: open-ended
    figures.append(figure(
        "affinity_ratio", "Does nucleosome binding cost affinity?", x_label="factor",
        y_label="K nucleosomal / K reference (S1/2 ratio)", y_scale="log", categories=names,
        series=affinity_series, bands=bands,
        note="Above 1: the nucleosome weakens binding. Pioneer factors (Reb1, Cbf1) need slower release from "
             "nucleosomes to match the titrations; occlusion alone suffices for Pho4 and Gal4."))
    for key in ("reb1", "cbf1"):
        if key not in params["factors"] or not FACTORS[key]["koff_reference_per_s"] or ref_dwell[0] is None:
            continue
        # s per reduced time unit: the simulated naked-site dwell maps to 1 / measured k_off(DNA).
        seconds = 1.0 / (FACTORS[key]["koff_reference_per_s"] * ref_dwell[0])
        plotted, times = [], np.geomspace(0.05, 2000, 40)
        for case_id, title, group in (("reference", "naked DNA", "naked DNA"),
                                      (f"{key}-association_only", "nucleosome, occlusion only", "occlusion only"),
                                      (f"{key}-association_and_release", "nucleosome, + release change",
                                       "nucleosome")):
            samples = [v * seconds for o in _per_case(results, case_id) for v in o["dwell_sample"]]
            if samples:
                plotted.append(series(f"Engine, {title}", times, survival(samples, times), style="line",
                                      group=group))
        rate_dna = FACTORS[key]["koff_reference_per_s"]
        rate_nuc = rate_dna * FACTORS[key]["koff_ratio"]

        def measured(rate):     # stop where a dwell sample of this size could no longer resolve it
            return [value if value >= 1e-3 else None for value in np.exp(-rate * times)]
        plotted += [series("Measured, naked DNA", times, measured(rate_dna), role="experiment", style="line",
                           group="naked DNA"),
                    series("Measured, nucleosome", times, measured(rate_nuc), role="experiment", style="line",
                           group="nucleosome")]
        figures.append(figure(
            f"dwell_{key}", f"{FACTORS[key]['label']} residence: dissociation compensation",
            x_label="time bound (s)", y_label="fraction still bound", x_scale="log", y_scale="log", series=plotted,
            note="Time mapped by the measured naked-DNA dissociation rate (input-matched). Measured curves are "
                 "single exponentials with the reported primary rates."))
    return figures
