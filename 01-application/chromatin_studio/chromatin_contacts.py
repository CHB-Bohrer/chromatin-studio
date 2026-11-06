"""Effective mark-dependent chromatin attraction in reduced energy units.

The deepest matching well is used per physical pair; simultaneous marks do not
add duplicate interactions. Changing marks changes the energy landscape, an
active phenomenological coupling rather than an equilibrium mark mechanism.
"""
import numpy as np

from .neighbors import neighbor_pairs


def compute_chromatin_contacts(positions, marks, rules):
    positions = np.asarray(positions, dtype=float)
    n = len(positions)
    force = np.zeros_like(positions)
    stiffness = np.zeros(n)
    reports = [dict(a=r["a"], b=r["b"], strength=r["strength"],
                    radius=r["radius"], min_separation=r.get("min_separation", 3),
                    contacts=0, energy=0.0) for r in rules]
    result = dict(forces=force, stiffness=stiffness, energy=0.0,
                  report=dict(rules=reports, total_contacts=0, energy=0.0))
    positive = [r for r in rules if r["strength"] > 0]
    if not positive:
        return result
    pairs = neighbor_pairs(positions, max(r["radius"] for r in positive))
    if not len(pairs):
        return result
    i, j = pairs.T
    delta = positions[i] - positions[j]
    r2 = np.einsum("ij,ij->i", delta, delta)
    states = {key: np.asarray(value, dtype=bool) for key, value in marks.items()}
    states["unmodified"] = ~np.logical_or.reduce(list(states.values())) if states else np.ones(n, dtype=bool)
    energy = np.zeros(len(pairs))
    winner = np.full(len(pairs), -1, dtype=int)
    coefficient = np.zeros(len(pairs))
    curvature = np.zeros(len(pairs))
    for k, rule in enumerate(rules):
        epsilon, radius = float(rule["strength"]), float(rule["radius"])
        if epsilon <= 0:
            continue
        a, b = states[rule["a"]], states[rule["b"]]
        eligible = ((j - i >= rule.get("min_separation", 3)) & (r2 < radius**2)
                    & ((a[i] & b[j]) | (b[i] & a[j])))
        shape = 1 - r2 / radius**2
        well = -epsilon * shape**2
        choose = eligible & (well < energy)
        winner[choose] = k
        energy[choose] = well[choose]
        coefficient[choose] = -4 * epsilon / radius**2 * shape[choose]
        curvature[eligible] = np.maximum(curvature[eligible], 16 * epsilon / radius**2)
    pair_force = coefficient[:, None] * delta
    np.add.at(force, i, pair_force)
    np.add.at(force, j, -pair_force)
    np.add.at(stiffness, i, curvature)
    np.add.at(stiffness, j, curvature)
    for k, report in enumerate(reports):
        selected = winner == k
        report.update(contacts=int(selected.sum()), energy=float(energy[selected].sum()))
    total = float(energy.sum())
    result["energy"] = total
    result["report"].update(total_contacts=int((winner >= 0).sum()), energy=total)
    return result
