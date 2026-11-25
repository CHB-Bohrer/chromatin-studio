"""DNA only: does the bending term give a worm-like chain with the intended stiffness?

Model: one free chain, no excluded volume, no confinement contact. The engine's
bending energy is kappa (1 - cos theta) per joint with kappa = persistence_length
/ bond_length, and successive joints are independent. For rigid bonds this is
exactly solvable:

    <cos theta> = coth(kappa) - 1/kappa                 (Langevin function)
    <t_i . t_j>  = <cos theta>^|i-j|
    effective persistence length = -bond / ln <cos theta>

so the configured ``persistence_length`` is a bending modulus, and the realized
persistence length is somewhat shorter on a discrete chain. The experimental
anchor is the persistence length of double-stranded DNA, about 50 nm.
"""

import math

import numpy as np

from ..engine import Simulation
from .core import ANALYTICAL, EXPERIMENTAL, bare_config, check, error_bar, figure, mean_sem, series, source, validated

NM_PER_BP = 0.34

STUDY = {
    "id": "polymer_wlc",
    "title": "Bare DNA as a worm-like chain",
    "components": ["chromatin mechanics (bonds, bending)"],
    "question": "Does a chain configured with the 50 nm persistence length of DNA show that persistence length?",
    "evidence": [ANALYTICAL, EXPERIMENTAL],
    "sources": [
        source("Bustamante, Marko, Siggia & Smith (1994) Science 265:1599", "https://doi.org/10.1126/science.8079175",
               "Force-extension of single lambda-DNA fits a worm-like chain with persistence length ~53 nm."),
        source("Baumann, Smith, Bloomfield & Bustamante (1997) PNAS 94:6185", "https://doi.org/10.1073/pnas.94.12.6185",
               "DNA persistence length 45-50 nm at physiological ionic strength."),
    ],
    "scales": {
        "smoke": {"beads": 24, "bp_per_bead": 25, "persistence_nm": 50.0, "bond_k": 100.0, "diffusion": 0.5,
                  "dt": 0.002, "steps": 4000, "burn_in_steps": 1000, "sample_every": 20, "replicas": 2, "seed": 11},
        "standard": {"beads": 48, "bp_per_bead": 25, "persistence_nm": 50.0, "bond_k": 100.0, "diffusion": 0.5,
                     "dt": 0.002, "steps": 6000000, "burn_in_steps": 1000000, "sample_every": 200,
                     "replicas": 8, "seed": 11},
    },
}


def _kappa(params):
    """Configured persistence length in bonds (bond length is 1 reduced unit)."""
    return params["persistence_nm"] / (params["bp_per_bead"] * NM_PER_BP)


def simulate(params, case, seed, save_partial):
    config = validated(bare_config(
        params["beads"], seed=seed, dt=params["dt"], radius=100.0,
        bp_per_bead=params["bp_per_bead"], bond_length=1.0, bond_k=params["bond_k"],
        diffusion=params["diffusion"], persistence_length=_kappa(params), repulsion=0.0))
    sim = Simulation(config)
    n = params["beads"]
    lags = n - 2
    cos_sum, r2_sum, bond_sum, bond2_sum, samples = np.zeros(lags), np.zeros(n - 1), 0.0, 0.0, 0
    sim.step(params["burn_in_steps"])
    for _ in range((params["steps"] - params["burn_in_steps"]) // params["sample_every"]):
        sim.step(params["sample_every"])
        x = sim.positions[:n]
        bonds = x[1:] - x[:-1]
        length = np.linalg.norm(bonds, axis=1)
        tangent = bonds / length[:, None]
        for lag in range(1, lags + 1):
            cos_sum[lag - 1] += float(np.mean(np.einsum("ij,ij->i", tangent[lag:], tangent[:-lag])))
        for separation in range(1, n):
            delta = x[separation:] - x[:-separation]
            r2_sum[separation - 1] += float(np.mean(np.einsum("ij,ij->i", delta, delta)))
        bond_sum += float(np.mean(length))
        bond2_sum += float(np.mean(length ** 2))
        samples += 1
    return {"samples": samples, "tangent_correlation": (cos_sum / samples).tolist(),
            "mean_square_distance": (r2_sum / samples).tolist(),
            "mean_bond": bond_sum / samples, "mean_bond_squared": bond2_sum / samples,
            "max_substeps": int(sim.last_substeps)}


def analyze(params, results):
    kappa = _kappa(params)
    cosine = 1.0 / math.tanh(kappa) - 1.0 / kappa
    bond_nm = params["bp_per_bead"] * NM_PER_BP
    predicted_nm = -bond_nm / math.log(cosine)
    done = [r["observables"] for r in results if r.get("complete")]
    checks, findings = [], []
    if done:
        first, first_sem = mean_sem([o["tangent_correlation"][0] for o in done])
        checks.append(check("joint_cosine", "Mean cosine between neighboring bonds", first, cosine,
                            evidence=ANALYTICAL, tolerance=0.01, sem=first_sem,
                            note="Langevin function of kappa = persistence_length / bond_length."))
        # Persistence length from the decay over the first few joints (the tail is noisy).
        fitted = []
        for o in done:
            c = np.asarray(o["tangent_correlation"][:6])
            lag = np.arange(1, len(c) + 1)
            good = c > 0.05
            if np.count_nonzero(good) >= 2:
                fitted.append(-bond_nm / np.polyfit(lag[good], np.log(c[good]), 1)[0])
        lp, lp_sem = mean_sem(fitted)
        checks.append(check("persistence_discrete", "Persistence length from tangent decay", lp, predicted_nm,
                            evidence=ANALYTICAL, tolerance=0.03, sem=lp_sem, unit="nm",
                            note="Expected value is the exact discrete-chain result, not the configured input."))
        checks.append(check("persistence_dna", "Persistence length versus double-stranded DNA", lp, 50.0,
                            evidence=EXPERIMENTAL, low=45.0, high=53.0, sem=lp_sem, unit="nm",
                            note=f"Length mapping: {params['bp_per_bead']} bp per bead at 0.34 nm per bp."))
        # End-to-end distance of the whole chain against the same discrete model.
        m = params["beads"] - 1
        b2, _ = mean_sem([o["mean_bond_squared"] for o in done])
        exact = b2 * (m * (1 + cosine) / (1 - cosine) - 2 * cosine * (1 - cosine ** m) / (1 - cosine) ** 2)
        r2, r2_sem = mean_sem([o["mean_square_distance"][-1] for o in done])
        checks.append(check("end_to_end", "Mean-square end-to-end distance", r2, exact, evidence=ANALYTICAL,
                            tolerance=0.05, sem=r2_sem, unit="reduced length^2",
                            note="Slowest observable; needs many relaxation times. Soft bonds enter through <b^2>."))
        needed = 1.0
        for _ in range(60):   # solve coth(k) - 1/k = exp(-bond/50 nm) for the modulus that realizes 50 nm
            needed = 1.0 / (1.0 / math.tanh(needed) - math.exp(-bond_nm / 50.0))
        findings.append(
            f"persistence_length is a bending modulus: {kappa:.2f} bonds configured gives {predicted_nm:.1f} nm "
            f"on this discrete chain; {needed:.2f} bonds would realize 50 nm at this resolution.")
    figures = []
    if done:
        lags = min(len(o["tangent_correlation"]) for o in done)
        contour = [bond_nm * s for s in range(1, lags + 1)]
        stats = [mean_sem([o["tangent_correlation"][s] for o in done]) for s in range(lags)]
        bars = [error_bar(m, e) for m, e in stats]
        figures.append(figure(
            "tangent_decay", "Bond orientation memory along the chain",
            x_label="contour separation (nm)", y_label="mean cosine between bonds", y_scale="log",
            series=[series("Discrete worm-like chain", contour, [cosine ** s for s in range(1, lags + 1)],
                           role="reference", style="line"),
                    series("DNA, 45 nm", contour, [math.exp(-c / 45.0) for c in contour], role="experiment",
                           style="line"),
                    series("DNA, 53 nm", contour, [math.exp(-c / 53.0) for c in contour], role="experiment",
                           style="line"),
                    series("Engine, mean ± SEM", contour, [m for m, _ in stats],
                           low=[b[0] for b in bars], high=[b[1] for b in bars])],
            note=f"{params['bp_per_bead']} bp per bond at 0.34 nm/bp. Points at large separation are noisy; "
                 "the persistence-length check fits the first six joints only."))
        sizes = min(len(o["mean_square_distance"]) for o in done)
        stats = [mean_sem([o["mean_square_distance"][s] for o in done]) for s in range(sizes)]
        bars = [error_bar(m, e) for m, e in stats]
        exact = [b2 * (m * (1 + cosine) / (1 - cosine) - 2 * cosine * (1 - cosine ** m) / (1 - cosine) ** 2)
                 for m in range(1, sizes + 1)]
        figures.append(figure(
            "internal_distance", "Mean-square distance between beads",
            x_label="separation (bonds)", y_label="mean-square distance (reduced length²)",
            x_scale="log", y_scale="log",
            series=[series("Discrete worm-like chain", list(range(1, sizes + 1)), exact, role="reference",
                           style="line"),
                    series("Engine, mean ± SEM", list(range(1, sizes + 1)), [m for m, _ in stats],
                           low=[b[0] for b in bars], high=[b[1] for b in bars])]))
    return {"checks": checks, "findings": findings, "figures": figures,
            "summary": {"configured_persistence_bonds": kappa, "bond_nm": bond_nm,
                        "predicted_persistence_nm": predicted_nm, "replicas": len(done)}}
