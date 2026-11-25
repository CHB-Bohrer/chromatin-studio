"""One non-binding protein species: is its diffusion what the input says?

Proteins with zero association hazard diffuse in a large sphere next to a short
inert chain. The centre-of-mass mean-square displacement must be 6 D t for a
single bead and, by the engine's stated convention, also for an IDR chain whose
beads each carry n x D.

There is no experimental pass/fail here. The report converts the input D to
physical units under the provisional polymer mapping so the gap to measured
transcription-factor mobility is visible.
"""

import numpy as np

from ..engine import Simulation
from .core import (ANALYTICAL, PROVISIONAL_CLOCK, bare_config, check, error_bar, figure, mean_sem, protein, series,
                   source, validated)

STUDY = {
    "id": "protein_diffusion",
    "title": "Free protein diffusion",
    "components": ["protein Brownian motion", "IDR chain drag convention"],
    "question": "Does a free molecule's centre of mass diffuse with the configured coefficient?",
    "evidence": [ANALYTICAL],
    "sources": [
        source("Einstein relation for free Brownian motion in 3D", "https://doi.org/10.1002/andp.19053220806",
               "MSD = 6 D t."),
    ],
    "scales": {
        "smoke": {"count": 12, "diffusion": 0.5, "dt": 0.002, "steps": 3000, "record_every": 50,
                  "lags": [1, 4], "radius": 40.0, "replicas": 2, "seed": 21},
        "standard": {"count": 40, "diffusion": 0.5, "dt": 0.002, "steps": 400000, "record_every": 50,
                     "lags": [1, 4, 16], "radius": 40.0, "replicas": 4, "seed": 21},
    },
}


def cases(params):
    return [{"id": "bead", "representation": "bead"}, {"id": "idr", "representation": "idr"}]


def simulate(params, case, seed, save_partial):
    config = bare_config(8, seed=seed, dt=params["dt"], radius=params["radius"], repulsion=8.0)
    config["proteins"] = [protein("probe", params["count"], diffusion=params["diffusion"],
                                  representation=case["representation"], idr_beads=4)]
    sim = Simulation(validated(config))
    molecules = [np.asarray(m["indices"]) for m in sim.molecules]

    def centres():
        return np.array([sim.positions[index].mean(axis=0) for index in molecules])

    frames = [centres()]
    for _ in range(params["steps"] // params["record_every"]):
        sim.step(params["record_every"])
        frames.append(centres())
    frames = np.asarray(frames)
    msd = {}
    for lag in params["lags"]:
        delta = frames[lag:] - frames[:-lag]
        msd[str(lag)] = float(np.mean(np.einsum("tmi,tmi->tm", delta, delta)))
    return {"msd": msd, "frame_time": params["record_every"] * params["dt"],
            "max_radius": float(np.max(np.linalg.norm(frames, axis=2))), "max_substeps": int(sim.last_substeps)}


def analyze(params, results):
    checks, plotted = [], []
    frame_time = params["record_every"] * params["dt"]
    for case in cases(params):
        done = [r["observables"] for r in results if r.get("complete") and r["case"]["id"] == case["id"]]
        points = []
        for lag in params["lags"]:
            if not done:
                continue
            value, sem = mean_sem([o["msd"][str(lag)] for o in done])
            time = lag * frame_time
            points.append((time, value, sem))
            checks.append(check(f"{case['id']}_lag{lag}", f"{case['id']} MSD at t = {time:g}", value,
                                6 * params["diffusion"] * time, evidence=ANALYTICAL, tolerance=0.05, sem=sem,
                                unit="reduced length^2",
                                note="Dilute probes; the confining wall and soft collisions are negligible here."))
        if points:
            bars = [error_bar(v, s) for _, v, s in points]
            plotted.append(series(f"Engine, {case['id']} probe (mean ± SEM)", [t for t, _, _ in points],
                                  [v for _, v, _ in points], low=[b[0] for b in bars], high=[b[1] for b in bars]))
    figures = []
    if plotted:
        times = [lag * frame_time for lag in params["lags"]]
        grid = np.geomspace(min(times) / 2, max(times) * 2, 40)
        figures.append(figure(
            "msd", "Free-protein mean-square displacement", x_label="lag (reduced time)",
            y_label="mean-square displacement (reduced length²)", x_scale="log", y_scale="log",
            series=[series("Free diffusion, 6Dt", grid, 6 * params["diffusion"] * grid, role="reference",
                           style="line")] + plotted,
            note="IDR probes are tracked by their centre of mass."))
    length_um = PROVISIONAL_CLOCK["length_nm_per_unit"] / 1000.0
    scale = length_um ** 2 / PROVISIONAL_CLOCK["time_s_per_unit"]
    findings = [
        f"Under the provisional 10 kb polymer mapping (needs_review), D = {params['diffusion']} reduced is "
        f"{params['diffusion'] * scale:.4f} um^2/s and the schema maximum of 5 is {5 * scale:.3f} um^2/s. "
        "Free nuclear transcription factors are typically reported at roughly 1-10 um^2/s by single-molecule "
        "tracking (order of magnitude; not re-verified in this study), so protein mobility in the coarse "
        "presets is not on the polymer clock and must be calibrated separately."]
    return {"checks": checks, "findings": findings, "figures": figures,
            "summary": {"reduced_to_um2_per_s": scale, "mapping": PROVISIONAL_CLOCK}}
