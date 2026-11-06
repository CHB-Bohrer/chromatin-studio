"""Isolated mechanistic studies, not fits to the accompanying manuscript.

The reference is an exact continuous-time Markov model of four promoters and
nascent-RNA birth/removal. The production polymer pilot is a different model;
its independent-replica null is explicitly not an experimental trans control.
All times and spatial distances remain in reduced units.
"""

import copy
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .schema import checked_config, validate_config


_CONDITIONS = [
    ("baseline", "All configured mechanisms", "Retain the selected model and its stated assumptions.", ("reference", "polymer")),
    ("uncoupled", "Independent promoters", "Remove cis competition, cis environmental modulation and cell-wide extrinsic variation.", ("reference",)),
    ("no_competition", "No cis competition", "Set mutual promoter inhibition to zero; retain shared environments and extrinsic variation.", ("reference",)),
    ("no_shared_environment", "No cis environmental drive", "Remove homolog-specific environmental modulation; retain competition and cell-wide variation.", ("reference",)),
    ("no_extrinsic_noise", "No cell-wide variation", "Use the same baseline activation factor in every cell; retain cis mechanisms.", ("reference",)),
    ("competition_only", "Cis competition only", "Remove environmental and cell-wide modulation; retain mutual cis inhibition.", ("reference",)),
    ("environment_only", "Shared cis environment only", "Remove competition and cell-wide variation; retain independent homolog-specific environmental switches.", ("reference",)),
    ("extrinsic_only", "Cell-wide variation only", "Remove both cis mechanisms; the shared cell factor can correlate cis and trans alike.", ("reference",)),
    ("no_regulatory_gain", "No regulatory activation gains", "Remove regional activation/initiation gains and legacy activator gain; retain binding and mark gates.", ("polymer",)),
    ("no_cohesin", "No cohesin", "Disable loop extrusion and its mechanical springs.", ("polymer",)),
    ("no_protein_binding", "No protein–protein binding", "Set protein-pair association hazards to zero; retain protein–chromatin binding.", ("polymer",)),
    ("no_mark_spreading", "No mark spreading", "Set contact-driven spreading to zero; retain bound-enzyme catalysis and decay.", ("polymer",)),
    ("no_chromatin_attraction", "No state-selective attraction", "Remove effective mark-selective chromatin attraction rules.", ("polymer",)),
    ("no_body_affinity", "No landmark affinity", "Remove mechanical attraction/exclusion; retain body geometry, optional gene body responses and processing effects.", ("polymer",)),
    ("no_speckle_processing", "No speckle processing gain", "Remove RNA processing enhancement; this is not a nascent-initiation perturbation.", ("polymer",)),
    ("no_transcription_stiffness", "No transcription stiffness feedback", "Set gene persistence multipliers to one.", ("polymer",)),
    ("no_body_response", "No direct body response", "Remove optional gene body-response hypotheses; retain body mechanics and speckle RNA processing.", ("polymer",)),
    ("no_confinement", "No nuclear wall", "Set soft-wall stiffness to zero while keeping the same initialization volume.", ("polymer",)),
    ("no_nucleosome_access", "No occupancy-dependent TF access", "Set occupied/vacant capture and release factors to one; retain occupancy, mark loss and mechanical effects.", ("polymer",)),
    ("no_nucleosome_turnover", "Frozen nucleosome occupancy", "Remove spontaneous, protein-driven and transcription-driven turnover and sliding; retain the initial occupancy pattern.", ("polymer",)),
]
_BY_ID = {item[0]: item for item in _CONDITIONS}
DEFAULTS = {
    "reference": {"engine": "reference", "conditions": ["baseline", "uncoupled", "no_competition", "no_shared_environment"],
                  "cells": 1024, "duration": 35.0, "seed": 42, "detection_efficiency": 1.0,
                  "on_rate": .5, "off_rate": 1.0, "initiation_rate": 2.0, "clearance_rate": 1.0,
                  "extrinsic_cv": .5, "competition_strength": 10.0, "environment_gain": 8.0,
                  "environment_on_rate": .08, "environment_off_rate": .08, "bootstrap_samples": 200},
    "polymer": {"engine": "polymer", "conditions": ["baseline", "no_regulatory_gain"],
                "gene_a": None, "gene_b": None, "replicas": 4, "warmup_steps": 100,
                "samples_per_replica": 40, "sample_every": 10, "seed": 42,
                "contact_radius": 1.5, "detection_efficiency": 1.0, "bootstrap_samples": 200},
}
LIMITS = {
    "reference": {"cells": {"min": 64, "max": 4096}, "duration": {"min": 1, "max": 200},
                  "bootstrap_samples": {"min": 0, "max": 1000}, "max_work": 200000},
    "polymer": {"replicas": {"min": 2, "max": 12}, "warmup_steps": {"min": 0, "max": 10000},
                "samples_per_replica": {"min": 4, "max": 500}, "sample_every": {"min": 1, "max": 1000},
                "bootstrap_samples": {"min": 0, "max": 1000}, "max_particle_steps": 20000000},
}
SCOPE = [
    "These are mechanistic sufficiency and ablation experiments, not a reproduction or fit of the manuscript's figures or 97% parameter scan.",
    "No live simulation is used or modified. A study starts independent seeded simulations from copied configurations.",
    "Reference nascent activity means RNA count > 0, separately from promoter ON. Exponential RNA residence is a reduced observation model, not explicit elongation or an intron-FISH forward model.",
    "Reference cis and trans compare two homolog pairs within the same simulated cell, sharing its extrinsic factor; there is no spatial polymer in that reference.",
    "The polymer engine has one chain. Its independent-replica null removes shared-cell effects and must not be called trans or subtracted as the manuscript's cis–trans contrast.",
    "A sampled minimum distance does not establish that a pair never contacted. Physical nanometer thresholds require independent length calibration and an imaging model.",
    "Bootstrap fractions are support measures, not p-values. Short pilot runs do not establish stationarity, converged correlations, or biological causality.",
]


def study_catalog():
    return {"engines": [{"id": "reference", "name": "Four-promoter reference",
                         "description": "Paired homologs, cis competition/environment, shared cell variation and nascent RNA."},
                        {"id": "polymer", "name": "Single-chain spatial pilot",
                         "description": "Independent copies of the selected polymer; within-chain association and contact analysis."}],
            "conditions": [{"id": i, "name": n, "description": d, "engines": list(e)}
                           for i, n, d, e in _CONDITIONS],
            "defaults": copy.deepcopy(DEFAULTS), "limits": copy.deepcopy(LIMITS), "scope": list(SCOPE)}


def _number(value, name, low, high, integer=False):
    if (isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value)
            or not low <= value <= high or (integer and value != int(value))):
        raise ValueError(f"{name} must be {'an integer ' if integer else ''}between {low} and {high}.")
    return int(value) if integer else float(value)


def _options(config, supplied):
    if supplied is None:
        supplied = {}
    if not isinstance(supplied, dict):
        raise ValueError("Study options must be an object.")
    engine = supplied.get("engine", "reference")
    if engine not in DEFAULTS:
        raise ValueError("Study engine must be reference or polymer.")
    unknown = set(supplied) - set(DEFAULTS[engine])
    if unknown:
        raise ValueError("Unknown study options: " + ", ".join(sorted(unknown)))
    result = copy.deepcopy(DEFAULTS[engine])
    result.update(copy.deepcopy(supplied))
    conditions = result["conditions"]
    if (not isinstance(conditions, list) or not conditions or len(conditions) > len(_CONDITIONS)
            or any(not isinstance(c, str) or c not in _BY_ID or engine not in _BY_ID[c][3] for c in conditions)):
        raise ValueError("Choose a nonempty list of conditions supported by this study engine.")
    if len(set(conditions)) != len(conditions):
        raise ValueError("Study conditions must not be duplicated.")
    for key, bounds in LIMITS[engine].items():
        if isinstance(bounds, dict):
            result[key] = _number(result[key], key, bounds["min"], bounds["max"], key != "duration")
    result["seed"] = _number(result["seed"], "seed", 0, 2**32 - 1, True)
    result["detection_efficiency"] = _number(result["detection_efficiency"], "detection_efficiency", 0, 1)
    if engine == "reference":
        for key in ("on_rate", "off_rate", "initiation_rate", "competition_strength", "environment_gain",
                    "environment_on_rate", "environment_off_rate"):
            result[key] = _number(result[key], key, 0, 20)
        result["clearance_rate"] = _number(result["clearance_rate"], "clearance_rate", .01, 20)
        result["extrinsic_cv"] = _number(result["extrinsic_cv"], "extrinsic_cv", 0, 2)
    else:
        result["contact_radius"] = _number(result["contact_radius"], "contact_radius", .01, 100)
        genes = [g["id"] for g in config["genes"]]
        if len(genes) < 2:
            raise ValueError("A spatial gene-pair study requires at least two configured genes.")
        result["gene_a"] = result["gene_a"] if result["gene_a"] is not None else genes[0]
        result["gene_b"] = result["gene_b"] if result["gene_b"] is not None else genes[1]
        if result["gene_a"] not in genes or result["gene_b"] not in genes or result["gene_a"] == result["gene_b"]:
            raise ValueError("Choose two different configured gene IDs.")
    return result


def _ablate(config, condition):
    c = copy.deepcopy(config)
    if condition == "no_regulatory_gain":
        for region in c.get("regulatory_regions", []):
            region["activation_gain"] = region["initiation_gain"] = 0.0
        for gene in c["genes"]:
            gene["activator_species"] = None
    elif condition == "no_cohesin":
        c["cohesin"]["enabled"] = False
    elif condition == "no_protein_binding":
        for interaction in c["interactions"]:
            interaction["on_rate"] = 0.0
    elif condition == "no_mark_spreading":
        for mark in c["marks"]:
            mark["spread_rate"] = 0.0
    elif condition == "no_chromatin_attraction":
        c["chromatin"]["state_interactions"] = []
    elif condition == "no_body_affinity":
        for body in c["bodies"]:
            body["strength"] = 0.0
    elif condition == "no_speckle_processing":
        for gene in c["genes"]:
            gene["speckle_gain"] = 0.0
    elif condition == "no_transcription_stiffness":
        for gene in c["genes"]:
            gene["persistence_multiplier"] = 1.0
    elif condition == "no_body_response":
        for gene in c["genes"]:
            gene["body_responses"] = {}
    elif condition == "no_confinement":
        c["simulation"]["confinement_stiffness"] = 0.0
    elif condition == "no_nucleosome_access":
        for protein in c["proteins"]:
            protein["nucleosome_affinity"] = {"occupied": 1.0, "vacant": 1.0}
            protein["nucleosome_off_rate_multiplier"] = {"occupied": 1.0, "vacant": 1.0}
    elif condition == "no_nucleosome_turnover":
        if "nucleosomes" in c["chromatin"]:
            c["chromatin"]["nucleosomes"]["eviction_rate"] = 0.0
            c["chromatin"]["nucleosomes"]["reassembly_rate"] = 0.0
        for protein in c["proteins"]:
            protein["nucleosome_eviction_rate"] = protein["nucleosome_reassembly_rate"] = 0.0
            protein["nucleosome_sliding_rate"] = 0.0
        for gene in c["genes"]:
            gene["nucleosome_eviction_rate"] = 0.0
    return c


def _reference_condition(options, condition):
    result = copy.deepcopy(options)
    if condition in ("uncoupled", "no_competition", "environment_only", "extrinsic_only"):
        result["competition_strength"] = 0.0
    if condition in ("uncoupled", "no_shared_environment", "competition_only", "extrinsic_only"):
        result["environment_gain"] = 0.0
    if condition in ("uncoupled", "no_extrinsic_noise", "competition_only", "environment_only"):
        result["extrinsic_cv"] = 0.0
    return result


def build_study(config, options=None):
    """Build a downloadable plan; copied polymer conditions are valid model JSON."""
    engine = options.get("engine", "reference") if isinstance(options, dict) else "reference"
    base = checked_config(config) if engine == "polymer" else None
    normalized = _options(base, options)
    conditions = []
    for identifier in normalized["conditions"]:
        _, name, description, _ = _BY_ID[identifier]
        row = {"id": identifier, "name": name, "description": description}
        if engine == "polymer":
            row["config"] = checked_config(_ablate(base, identifier))
        else:
            row["options"] = _reference_condition(normalized, identifier)
        conditions.append(row)
    warnings = ["Pilot duration and sample count are numerical choices, not established equilibration or power requirements."]
    if normalized["detection_efficiency"] < 1:
        warnings.append("Detection independently misses a fraction of nascent-positive alleles. No universal multiplicative correlation correction is applied.")
    if engine == "polymer":
        warnings += ["Polymer confidence intervals resample independent replicas, retaining all frames in each replica.",
                     "A common seed across ablations does not keep reaction histories paired after mechanisms diverge."]
    else:
        warnings += ["The default toy coupling strengths deliberately make contrasting mechanisms visible; they were not fitted to the manuscript.",
                     "Cells start with promoters OFF and no RNA. Duration is a transient observation time; check longer durations before assuming stationarity."]
    return {"engine": engine, "options": normalized, "conditions": conditions, "scope": list(SCOPE),
            "warnings": warnings, "version": 1}


def correlation(x, y):
    """Pearson correlation, equal to phi for binary inputs; constants are undefined."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("Correlation inputs must be equal-length vectors.")
    if len(x) < 2 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return None
    x, y = x - x.mean(), y - y.mean()
    variance = float(np.dot(x, x) * np.dot(y, y))
    if variance <= 0:
        return None
    return float(np.clip(np.dot(x, y) / math.sqrt(variance), -1, 1))


def detection_phi(latent_phi, p_a, p_b, q_a, q_b=None):
    """Exact independent false-negative attenuation for binary allele detection.

    Observed covariance=q_a*q_b*latent covariance. The attenuation depends on
    both marginal probabilities, so a common detection efficiency alone does
    not justify one universal correction factor.
    """
    q_b = q_a if q_b is None else q_b
    for name, value in (("p_a", p_a), ("p_b", p_b), ("q_a", q_a), ("q_b", q_b)):
        _number(value, name, 0, 1)
    if latent_phi is None or not (0 < p_a < 1 and 0 < p_b < 1 and q_a > 0 and q_b > 0):
        return None
    return float(latent_phi * math.sqrt(q_a * q_b * (1-p_a) * (1-p_b)
                                       / ((1-q_a*p_a) * (1-q_b*p_b))))


def _four_correlations(values):
    values = np.asarray(values)
    if values.ndim != 2 or values.shape[1] != 4:
        raise ValueError("Four-promoter observations must have columns A, B, A_prime, B_prime.")
    pairs = [correlation(values[:, a], values[:, b]) for a, b in ((0, 1), (2, 3), (0, 3), (2, 1))]
    cis = sum(pairs[:2]) / 2 if all(v is not None for v in pairs[:2]) else None
    trans = sum(pairs[2:]) / 2 if all(v is not None for v in pairs[2:]) else None
    delta = cis - trans if cis is not None and trans is not None else None
    return {"phi_cis": cis, "phi_trans": trans, "delta_phi": delta, "pair_correlations": pairs,
            "reason": None if delta is not None else "At least one allele has no variation; the matched contrast is undefined."}


def four_promoter_statistics(values, *, bootstrap_samples=200, seed=1):
    """Matched correlations with a bootstrap over complete independent cells."""
    values = np.asarray(values)
    result = _four_correlations(values)
    _number(bootstrap_samples, "bootstrap_samples", 0, 1000, True)
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(bootstrap_samples):
        delta = _four_correlations(values[rng.integers(len(values), size=len(values))])["delta_phi"]
        if delta is not None:
            boot.append(delta)
    result.update(bootstrap_fraction=float(np.mean(np.asarray(boot) >= 0)) if boot else None,
                  ci95=np.quantile(boot, [.025, .975]).tolist() if boot else [None, None],
                  valid_bootstraps=len(boot), bootstrap_unit="independent cell",
                  bootstrap_target="delta_phi", n_cells=len(values))
    return result


def simulate_reference(options, *, seed=None, max_events=2_000_000, deadline=None):
    """Exact Gillespie paths, vectorized over independent cells, observed once.

    Every cell starts with promoters OFF and no RNA at t=0. Its lognormal
    extrinsic factor has mean one and is constant across all four alleles for
    that cell. Independent homolog environment switches modulate A/B and A'/B'.
    Event clocks are resampled only after an event, as in the direct SSA.
    """
    n, horizon = int(options["cells"]), float(options["duration"])
    rng = np.random.default_rng(options["seed"] if seed is None else seed)
    promoters = np.zeros((n, 4), dtype=bool)
    rna = np.zeros((n, 4), dtype=np.int64)
    env_on, env_off = options["environment_on_rate"], options["environment_off_rate"]
    environment = rng.random((n, 2)) < (env_on / (env_on + env_off) if env_on + env_off else 0)
    sigma = math.sqrt(math.log1p(options["extrinsic_cv"] ** 2))
    extrinsic = np.exp(rng.normal(-sigma*sigma/2, sigma, n))
    clocks = np.zeros(n)
    active = np.arange(n)
    total_events, rounds = 0, 0
    while active.size:
        if rounds % 100 == 0 and deadline is not None and time.monotonic() > deadline:
            raise ValueError("Study pilot exceeded 30 seconds. Reduce cells/duration/conditions and use the saved plan for larger runs.")
        p, counts, env = promoters[active], rna[active], environment[active]
        activation = options["on_rate"] * extrinsic[active, None]
        activation = activation * (1 + options["environment_gain"] * np.repeat(env, 2, axis=1))
        activation = activation / (1 + options["competition_strength"] * p[:, [1, 0, 3, 2]])
        rates = np.concatenate((np.where(p, options["off_rate"], activation),
                                options["initiation_rate"] * p,
                                options["clearance_rate"] * counts,
                                np.where(env, env_off, env_on)), axis=1)
        totals = rates.sum(axis=1)
        mobile = totals > 0
        active, rates, totals = active[mobile], rates[mobile], totals[mobile]
        if not active.size:
            break
        next_times = clocks[active] + rng.exponential(1 / totals)
        event = next_times <= horizon
        active, rates, totals = active[event], rates[event], totals[event]
        if not active.size:
            break
        clocks[active] = next_times[event]
        choices = np.argmax(np.cumsum(rates, axis=1) > rng.random(len(active))[:, None] * totals[:, None], axis=1)
        toggle = choices < 4
        promoters[active[toggle], choices[toggle]] ^= True
        initiation = (choices >= 4) & (choices < 8)
        rna[active[initiation], choices[initiation] - 4] += 1
        removal = (choices >= 8) & (choices < 12)
        rna[active[removal], choices[removal] - 8] -= 1
        environmental = choices >= 12
        environment[active[environmental], choices[environmental] - 12] ^= True
        total_events += len(active)
        rounds += 1
        if total_events > max_events:
            raise ValueError("Study pilot exceeded the two-million-event safety limit; reduce cells, rates or duration.")
    return {"promoter_on": promoters, "nascent_counts": rna, "nascent_positive": rna > 0,
            "environment_on": environment, "extrinsic_factor": extrinsic, "events": total_events}


def _condition_seed(seed, identifier, replica=0, stream=0):
    code = int.from_bytes(hashlib.sha256(identifier.encode()).digest()[:4], "little")
    return int(np.random.SeedSequence([seed, code, replica, stream]).generate_state(1)[0])


def _run_reference(plan, progress, deadline):
    settings = plan["options"]
    if settings["cells"] * settings["duration"] * len(plan["conditions"]) > LIMITS["reference"]["max_work"]:
        raise ValueError("Reference pilot exceeds the cells × duration × conditions budget; reduce the pilot or export its plan.")
    rows, events = [], 0
    for index, condition in enumerate(plan["conditions"]):
        options = condition["options"]
        seed = _condition_seed(options["seed"], condition["id"])
        result = simulate_reference(options, seed=seed, max_events=2_000_000-events, deadline=deadline)
        events += result["events"]
        rng = np.random.default_rng(_condition_seed(options["seed"], condition["id"], stream=1))
        observed = result["nascent_positive"] & (rng.random((options["cells"], 4)) < options["detection_efficiency"])
        statistics = {key: four_promoter_statistics(values, bootstrap_samples=options["bootstrap_samples"],
                                                   seed=_condition_seed(options["seed"], condition["id"], stream=2))
                      for key, values in (("latent", result["nascent_positive"]), ("observed", observed),
                                          ("promoter", result["promoter_on"]))}
        statistics.update(mean_nascent_probability=result["nascent_positive"].mean(axis=0).tolist(),
                          observed_nascent_probability=observed.mean(axis=0).tolist(),
                          mean_nascent_counts=result["nascent_counts"].mean(axis=0).tolist(),
                          mean_promoter_on=result["promoter_on"].mean(axis=0).tolist())
        warnings = []
        if any(p < .05 or p > .95 for p in statistics["mean_nascent_probability"]):
            warnings.append("At least one nascent allele is rarely positive or nearly always positive; correlation is poorly resolved.")
        rows.append({"id": condition["id"], "name": condition["name"], "metrics": statistics,
                     "options": options, "seed": seed, "events": result["events"], "warnings": warnings,
                     "samples": {"columns": ["A", "B", "A_prime", "B_prime"],
                                 "nascent_counts": result["nascent_counts"].tolist(),
                                 "promoter_on": result["promoter_on"].astype(int).tolist(),
                                 "observed_positive": observed.astype(int).tolist(),
                                 "extrinsic_factor": result["extrinsic_factor"].tolist()}})
        if progress:
            progress({"condition": condition["id"], "completed": index+1, "total": len(plan["conditions"])})
    return rows


def _pair_vectors(samples, gene_a, gene_b):
    if not isinstance(samples, list) or not samples:
        raise ValueError("Gene-pair analysis requires a nonempty sample list.")
    try:
        a = [sample["genes"][gene_a] for sample in samples]
        b = [sample["genes"][gene_b] for sample in samples]
        counts_a, counts_b = np.array([x["nascent"] for x in a]), np.array([x["nascent"] for x in b])
        distances = np.linalg.norm(np.array([x["position"] for x in a]) - np.array([x["position"] for x in b]), axis=1)
        if not (np.all(np.isfinite(counts_a)) and np.all(np.isfinite(counts_b)) and np.all(np.isfinite(distances))):
            raise ValueError("Samples must be finite.")
        return a, b, counts_a, counts_b, distances
    except (KeyError, TypeError, IndexError) as exc:
        raise ValueError("Each sample needs both genes, nascent counts, promoter ON and 3D promoter positions.") from exc


def analyze_gene_pair(samples, gene_a, gene_b, *, contact_radius=1.5, bootstrap_seed=1, bootstrap_samples=200):
    """Analyze within-chain samples, with replica-clustered uncertainty.

    Each sample has replica/time and genes[id]={nascent,on,position,...}.
    Distances are between gene-start bead centers. Rows are observations, not
    independent cells; the bootstrap resamples complete replicas.
    """
    if gene_a == gene_b:
        raise ValueError("Choose distinct genes.")
    _number(contact_radius, "contact_radius", .01, 100)
    _number(bootstrap_samples, "bootstrap_samples", 0, 1000, True)
    a, b, counts_a, counts_b, distances = _pair_vectors(samples, gene_a, gene_b)
    latent_a, latent_b = counts_a > 0, counts_b > 0
    observed_a = np.array([x.get("observed_positive", x["nascent"] > 0) for x in a], dtype=bool)
    observed_b = np.array([x.get("observed_positive", x["nascent"] > 0) for x in b], dtype=bool)
    groups = {}
    for index, sample in enumerate(samples):
        groups.setdefault(sample["replica"], []).append(index)
    ids = sorted(groups)
    rng, boot = np.random.default_rng(bootstrap_seed), []
    if len(ids) >= 2:
        for _ in range(bootstrap_samples):
            indices = np.concatenate([groups[ids[int(i)]] for i in rng.integers(len(ids), size=len(ids))])
            value = correlation(latent_a[indices], latent_b[indices])
            if value is not None:
                boot.append(value)
    null_a, null_b = [], []
    if len(ids) >= 2:
        for index, replica in enumerate(ids):
            left = sorted(groups[replica], key=lambda i: samples[i]["time"])
            right = sorted(groups[ids[(index+1) % len(ids)]], key=lambda i: samples[i]["time"])
            for i, j in zip(left, right):
                if not math.isclose(samples[i]["time"], samples[j]["time"], rel_tol=1e-9, abs_tol=1e-12):
                    continue
                null_a.append(latent_a[i]); null_b.append(latent_b[j])
    result = {"gene_a": gene_a, "gene_b": gene_b,
              "nascent_phi": correlation(latent_a, latent_b), "observed_nascent_phi": correlation(observed_a, observed_b),
              "promoter_phi": correlation([x["on"] for x in a], [x["on"] for x in b]),
              "count_correlation": correlation(counts_a, counts_b),
              "mean_distance": float(distances.mean()), "median_distance": float(np.median(distances)),
              "min_distance": float(distances.min()), "contact_fraction": float(np.mean(distances <= contact_radius)),
              "contact_radius": contact_radius, "distance_definition": "gene-start bead center separation in reduced length units",
              "independent_replica_null_phi": correlation(null_a, null_b),
              "null_definition": "Gene A from replica r versus gene B at matching saved time in replica r+1; not trans",
              "genomic_separation_bp": abs(a[0].get("start_bp", 0) - b[0].get("start_bp", 0)),
              "bootstrap_ci95": np.quantile(boot, [.025, .975]).tolist() if boot else [None, None],
              "valid_bootstraps": len(boot), "bootstrap_unit": "independent replica", "bootstrap_target": "nascent_phi",
              "n_replicas": len(ids), "n_frames": len(samples),
              "mean_nascent_probability": [float(latent_a.mean()), float(latent_b.mean())],
              "mean_nascent_counts": [float(counts_a.mean()), float(counts_b.mean())],
              "reason": None if correlation(latent_a, latent_b) is not None else "One or both genes have constant nascent-positive state."}
    marks = sorted(set(a[0].get("mark_fractions", {})) & set(b[0].get("mark_fractions", {})))
    result["mark_cofluctuations"] = {mark: correlation([x["mark_fractions"][mark] for x in a],
                                                      [x["mark_fractions"][mark] for x in b]) for mark in marks}
    return result


def _polymer_sample(sim, replica, selected, detection_rng, efficiency):
    genes = {}
    for gene in sim.genes:
        spec = gene["config"]
        if spec["id"] not in selected:
            continue
        start, end = int(spec["start"]), int(spec["end"])
        nascent = len(gene["polymerases"])
        genes[spec["id"]] = {"nascent": nascent, "on": bool(gene["on"]),
                              "observed_positive": bool(nascent > 0 and detection_rng.random() < efficiency),
                              "position": sim.positions[start].tolist(),
                              "start_bp": start * sim.chrom["bp_per_bead"],
                              "mark_fractions": {mark: float(np.mean(flags[start:end+1])) for mark, flags in sim.marks.items()},
                              "speckle_contact": sim._gene_speckle_contact(spec)}
    return {"replica": replica, "time": float(sim.time), "genes": genes}


def _run_polymer(plan, progress, deadline):
    from .engine import Simulation
    options = plan["options"]
    steps = options["warmup_steps"] + options["samples_per_replica"] * options["sample_every"]
    particles = max(validate_config(c["config"])["particle_count"] for c in plan["conditions"])
    if particles * steps * options["replicas"] * len(plan["conditions"]) > LIMITS["polymer"]["max_particle_steps"]:
        raise ValueError("Spatial pilot exceeds the particle-step budget. Download its plan for a larger batch experiment or reduce pilot size.")
    def advance(sim, count):
        while count:
            if time.monotonic() > deadline:
                raise ValueError("Spatial pilot exceeded 30 seconds; reduce steps/replicas or run exported condition models as batch experiments.")
            chunk = min(count, 25)
            sim.step(chunk)
            count -= chunk
    rows = []
    for condition in plan["conditions"]:
        samples, seeds, blocked = [], [], 0
        for replica in range(options["replicas"]):
            seed = _condition_seed(options["seed"], "polymer", replica)
            config = copy.deepcopy(condition["config"])
            config["seed"] = seed
            sim = Simulation(config)
            detection_rng = np.random.default_rng(_condition_seed(options["seed"], condition["id"], replica, 1))
            advance(sim, options["warmup_steps"])
            for _ in range(options["samples_per_replica"]):
                if time.monotonic() > deadline:
                    raise ValueError("Spatial pilot exceeded 30 seconds; reduce steps/replicas or run exported condition models as batch experiments.")
                advance(sim, options["sample_every"])
                samples.append(_polymer_sample(sim, replica, {options["gene_a"], options["gene_b"]},
                                               detection_rng, options["detection_efficiency"]))
            if sim.config.get("transcription", {}).get("enabled"):
                blocked += int(sim.snapshot()["transcription"]["totals"]["blocked_initiations"])
            seeds.append(seed)
            if progress:
                progress({"condition": condition["id"], "replica": replica+1, "total_replicas": options["replicas"]})
        metrics = analyze_gene_pair(samples, options["gene_a"], options["gene_b"], contact_radius=options["contact_radius"],
                                    bootstrap_seed=options["seed"], bootstrap_samples=options["bootstrap_samples"])
        warnings = ["Saved frames are autocorrelated; the interval and duration have not been selected from a convergence study."]
        if metrics["nascent_phi"] is None:
            warnings.append("No interpretable nascent correlation: one or both genes stayed constant in this short pilot.")
        if options["replicas"] < 8:
            warnings.append("Fewer than eight independent replicas: bootstrap intervals are exploratory and may be unstable.")
        if blocked:
            warnings.append("RNA capacity blocked initiation; these expression statistics are not quantitatively interpretable.")
        config_hash = hashlib.sha256(json.dumps(condition["config"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        rows.append({"id": condition["id"], "name": condition["name"], "metrics": metrics, "seeds": seeds,
                     "config_sha256": config_hash,
                     "samples": samples, "blocked_initiations": blocked, "warnings": warnings})
    return rows


def run_study(config, options=None, *, progress=None):
    """Run bounded independent pilots; no reference to the live engine is accepted."""
    plan = build_study(config, options)
    start = time.monotonic()
    rows = (_run_reference(plan, progress, start+30) if plan["engine"] == "reference"
            else _run_polymer(plan, progress, start+30))
    return {"engine": plan["engine"], "options": plan["options"], "conditions": rows,
            "scope": plan["scope"], "warnings": plan["warnings"], "elapsed_seconds": time.monotonic()-start,
            "version": 1, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def reference_verification():
    """Reproducible numerical controls, deliberately separate from a paper fit."""
    from datetime import datetime, timezone
    import platform
    options = {"conditions": ["uncoupled", "competition_only", "environment_only", "extrinsic_only"],
               "cells": 2048, "duration": 24, "seed": 1741, "bootstrap_samples": 300,
               "competition_strength": 20, "environment_gain": 12, "extrinsic_cv": 1,
               "environment_on_rate": .06, "environment_off_rate": .06, "detection_efficiency": .3}
    result = run_study(None, options)
    controls = {condition["id"]: condition for condition in result["conditions"]}
    checks = []
    def check(name, observed, expected, tolerance, passed, explanation):
        checks.append({"name": name, "observed": observed, "expected": expected, "tolerance": tolerance,
                       "passed": bool(passed), "explanation": explanation})
    uncoupled = controls["uncoupled"]
    m = uncoupled["metrics"]["latent"]
    for label in ("phi_cis", "phi_trans", "delta_phi"):
        check("uncoupled_"+label, m[label], 0.0, .07, abs(m[label]) < .07,
              "Independent promoters have zero population correlation; finite-cell estimates fluctuate.")
    on, off, initiation, clearance = (result["options"][k] for k in ("on_rate", "off_rate", "initiation_rate", "clearance_rate"))
    expected_mean = initiation * on / ((on+off)*clearance)
    expected_variance = expected_mean + initiation**2 * on * off / (clearance*(on+off)**2*(clearance+on+off))
    counts = np.asarray(uncoupled["samples"]["nascent_counts"]).ravel()
    tolerance = 6 * math.sqrt(expected_variance/len(counts))
    check("uncoupled_telegraph_mean", float(counts.mean()), expected_mean, tolerance,
          abs(counts.mean()-expected_mean) < tolerance,
          "Stationary birth/clearance telegraph mean; six analytical standard errors. A duration sweep remains necessary for other parameters.")
    ratio = on/off
    double = ratio**2 / (1+options["competition_strength"])
    z = 1+2*ratio+double
    p, p11 = (ratio+double)/z, double/z
    expected_phi = (p11-p*p)/(p*(1-p))
    observed = controls["competition_only"]["metrics"]["promoter"]["phi_cis"]
    check("competition_stationary_promoter_phi", observed, expected_phi, .06,
          abs(observed-expected_phi) < .06,
          "Exact four-state stationary promoter distribution with weights 1,r,r,r²/(1+c). This target is promoter ON, not nascent RNA.")
    negative = controls["competition_only"]["metrics"]["latent"]
    check("competition_nascent_sign", negative["delta_phi"], "negative", "upper 95% cell-bootstrap bound < 0",
          negative["ci95"][1] < 0, "Mutual inhibition is sufficient for cis-reduced nascent coupling under these toy parameters.")
    positive = controls["environment_only"]["metrics"]["latent"]
    check("shared_environment_nascent_sign", positive["delta_phi"], "positive", "lower 95% cell-bootstrap bound > 0",
          positive["ci95"][0] > 0, "Independent homolog environments positively coordinate two cis promoters without direct competition.")
    extrinsic = controls["extrinsic_only"]["metrics"]["latent"]
    check("cell_wide_noise_preserves_cis_trans_symmetry", extrinsic["delta_phi"], 0.0, .07,
          abs(extrinsic["delta_phi"]) < .07 and extrinsic["phi_cis"] > .03 and extrinsic["phi_trans"] > .03,
          "Exchangeable alleles share a cell factor; both correlations rise while their population difference remains zero.")
    target_pairs = []
    probabilities = controls["environment_only"]["metrics"]["mean_nascent_probability"]
    for value, (a, b) in zip(positive["pair_correlations"], ((0, 1), (2, 3), (0, 3), (2, 1))):
        target_pairs.append(detection_phi(value, probabilities[a], probabilities[b], options["detection_efficiency"]))
    predicted_delta = float(np.mean(target_pairs[:2])-np.mean(target_pairs[2:]))
    measured_delta = controls["environment_only"]["metrics"]["observed"]["delta_phi"]
    check("false_negative_detection_forward_model", measured_delta, predicted_delta, .08,
          abs(measured_delta-predicted_delta) < .08,
          "Independent allele-level detection attenuates covariance by q²; correlation attenuation also depends on both latent marginal probabilities.")
    for condition in result["conditions"]:
        samples = condition.pop("samples")
        condition["samples_sha256"] = hashlib.sha256(json.dumps(samples, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"status": "passed" if all(c["passed"] for c in checks) else "needs_review",
            "created_utc": datetime.now(timezone.utc).isoformat(), "checks": checks,
            "passed_checks": sum(c["passed"] for c in checks), "total_checks": len(checks),
            "reproduce": "python -m chromatin_studio.studies --verify --output references/manuscript/new-verification.json",
            "environment": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()},
            "reference_result": result,
            "sources": ["https://doi.org/10.1073/pnas.0803850105", "https://doi.org/10.1073/pnas.162041399"],
            "limitations": ["These are reduced-model controls, not a manuscript figure or parameter-scan reproduction.",
                            "The main manuscript's supporting kinetic methods and original analysis data were not supplied.",
                            "The lognormal cell factor is a positive illustrative alternative to the manuscript's stated normal modulation.",
                            "Bootstrap intervals use independent simulated cells; they do not incorporate biological model or parameter uncertainty."]}


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Run the independent reference-model controls.")
    parser.add_argument("--output", required=True, help="New JSON output path; existing files are not overwritten.")
    args = parser.parse_args(argv)
    if not args.verify:
        parser.error("Choose --verify.")
    output = Path(args.output)
    if output.exists():
        parser.error("Output already exists; choose a new path.")
    report = reference_verification()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"{report['status']}: {report['passed_checks']}/{report['total_checks']} controls; saved {output}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
