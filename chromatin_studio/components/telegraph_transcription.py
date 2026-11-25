"""One gene alone: two-state promoter, initiation, release and decay.

With constant rates the released-RNA copy number follows the telegraph model
(Peccoud & Ycart 1995). Elongation only delays release and processing only
relabels RNA, so neither changes the stationary statistics of total RNA:

    ON fraction  p   = k_on / (k_on + k_off)
    mean RNA         = p k_ini / delta
    Fano factor      = 1 + k_ini k_off / ((k_on + k_off)(k_on + k_off + delta))

Both transcription representations are tested: the count-based default and the
individual-polymerase/RNA mode. Coordinates are frozen (chemistry only).
"""

import numpy as np

from ..engine import Simulation
from .core import ANALYTICAL, bare_config, check, chemistry_only_step, mean_sem, source, validated

STUDY = {
    "id": "telegraph_transcription",
    "title": "Single-gene bursting versus the telegraph model",
    "components": ["promoter switching", "initiation", "elongation", "RNA processing and decay"],
    "question": "Do the ON fraction, mean and Fano factor of RNA copy number match the telegraph model?",
    "evidence": [ANALYTICAL],
    "sources": [
        source("Peccoud & Ycart (1995) Theor. Popul. Biol. 48:222", "https://doi.org/10.1006/tpbi.1995.1027",
               "Stationary moments of the two-state (telegraph) gene expression model."),
        source("Raj, Peskin, Tranchina, Vargas & Tyagi (2006) PLoS Biol. 4:e309",
               "https://doi.org/10.1371/journal.pbio.0040309",
               "Context only: mammalian mRNA is made in bursts consistent with a two-state promoter. "
               "No copy-number distribution is compared here.", used_as="context"),
    ],
    "scales": {
        "smoke": {"on_rate": 0.05, "off_rate": 0.2, "initiation_rate": 2.0, "decay_rate": 0.02,
                  "processing_rate": 0.5, "elongation_rate": 20.0, "reaction_interval": 1.0, "burn_in": 400.0,
                  "sample_every": 200.0, "samples": 60, "replicas": 2, "seed": 41},
        "standard": {"on_rate": 0.05, "off_rate": 0.2, "initiation_rate": 2.0, "decay_rate": 0.02,
                     "processing_rate": 0.5, "elongation_rate": 20.0, "reaction_interval": 1.0, "burn_in": 1000.0,
                     "sample_every": 200.0, "samples": 3000, "replicas": 8, "seed": 41},
    },
}


def cases(params):
    return [{"id": "count", "explicit": False}, {"id": "explicit", "explicit": True}]


def _theory(params):
    on, off, ini, decay = (params[k] for k in ("on_rate", "off_rate", "initiation_rate", "decay_rate"))
    p = on / (on + off)
    return {"on_fraction": p, "mean": p * ini / decay,
            "fano": 1 + ini * off / ((on + off) * (on + off + decay))}


def simulate(params, case, seed, save_partial):
    config = bare_config(8, seed=seed, radius=30.0, repulsion=0.0)
    interval = params["reaction_interval"]
    config["genes"] = [{"id": "gene", "name": "gene", "start": 2, "end": 5,
                        "on_rate": params["on_rate"], "off_rate": params["off_rate"],
                        "initiation_rate": params["initiation_rate"], "elongation_rate": params["elongation_rate"],
                        "processing_rate": params["processing_rate"], "rna_decay_rate": params["decay_rate"],
                        "activator_species": None, "required_mark": None, "repressive_mark": None,
                        "speckle_gain": 0.0, "persistence_multiplier": 1.0}]
    if case["explicit"]:
        # Zero RNA diffusion: copy-number statistics do not depend on RNA position.
        config["transcription"] = {"enabled": True, "max_transcripts": 5000,
                                   "rna_diffusion": 0.0, "mature_rna_diffusion": 0.0}
    sim = Simulation(validated(config))
    gene = sim.genes[0]
    chemistry_only_step(sim, round(params["burn_in"] / interval), interval)
    counts, on_time, total_time = [], 0.0, 0.0
    block = round(params["sample_every"] / interval)
    for _ in range(params["samples"]):
        for _ in range(block):
            chemistry_only_step(sim, 1, interval)
            on_time += gene["last_transcription_interval"]["on_time"]
            total_time += interval
        counts.append(int(gene["unprocessed"]) + int(gene["mature"]))
    counts = np.asarray(counts, dtype=float)
    blocked = sim.explicit_transcription.snapshot()["totals"]["blocked_initiations"] if case["explicit"] else 0
    return {"mean": float(counts.mean()), "variance": float(counts.var(ddof=1)),
            "on_fraction": on_time / total_time, "samples": len(counts),
            "initiated": int(gene["initiated"]), "completed": int(gene["completed"]), "blocked": int(blocked)}


def analyze(params, results):
    theory = _theory(params)
    checks = []
    for case in cases(params):
        done = [r["observables"] for r in results if r.get("complete") and r["case"]["id"] == case["id"]]
        if not done:
            continue
        for key, label, tolerance in (("on_fraction", "promoter ON fraction", 0.03), ("mean", "mean RNA copy number", 0.05)):
            value, sem = mean_sem([o[key] for o in done])
            checks.append(check(f"{case['id']}_{key}", f"{case['id']}: {label}", value, theory[key],
                                evidence=ANALYTICAL, tolerance=tolerance, sem=sem))
        value, sem = mean_sem([o["variance"] / o["mean"] for o in done])
        checks.append(check(f"{case['id']}_fano", f"{case['id']}: Fano factor", value, theory["fano"],
                            evidence=ANALYTICAL, tolerance=0.10, sem=sem,
                            note="Samples are spaced four RNA lifetimes apart."))
    return {"checks": checks, "findings": [], "summary": {"theory": theory}}
