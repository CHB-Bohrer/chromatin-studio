<!--
Section: How time is advanced (numerics)
Update when: step ordering, operator splitting, sub-stepping, RNG/seed-replay rules, or backends change.
Code: chromatin_studio/engine.py (step, advance_reactions, integration_substeps, _probability),
      openmm_backend.py, openmm_native.py, batch.py, cluster.py.
Sources: AGENTS.md "Scientific distinctions", README.md "Two available integration backends", docs/PROJECT_STATUS.md "GPU and scale status".
-->

# How time is advanced

::: {.box}
**Short version.** Each step first moves all particles, then updates all the reactions using the new positions, in a fixed order. Inside each reaction family the random events are drawn exactly for rates held fixed over the step. Coupling between families is approximate (operator splitting), so results should be checked at a smaller time step.
:::

## One step

![**One simulation step.** Mechanics (possibly in several sub-steps) moves every particle. Then the reaction families run in this order, each seeing the state left by the previous one. The rates of each family are frozen for the step.](figures/fig-step.svg){#fig:step width=100%}

The order inside `advance_reactions` is: protein binding and release, nucleosome turnover, nucleosome sliding, marks, cohesin, tethers, and finally genes and RNA. A spatial search tree of the chromatin positions is built once per step and shared by binding and mark spreading.

## What is exact and what is not

- **Exact, given frozen rates.** Two-state flips (marks, nucleosomes) use the exact two-state probability. Promoter paths and initiations are drawn as exact event sequences inside the step. Sliding events and competing binding encounters are ordered by exponential clocks.
- **Approximate.** Rates depend on positions and on other states, but they are held fixed over a step, and the families run one after another. This is called **operator splitting**. It is accurate when $\Delta t$ is small compared with the fastest reaction and the fastest motion that matters. The whole simulation is therefore **not** an exact global Gillespie simulation. The safe practice is to repeat key results at a smaller $\Delta t$.
- **Mechanical sub-steps** protect stability only. They do not refine the reaction updates.

## Randomness and replay

One random number generator, seeded from the model's `seed`, drives the whole run. Two promises follow from it:

1. **Observers do not disturb.** Snapshots, rate tables, event logs and charts never draw a random number.
2. **New features are invisible when off.** An optional mechanism that is disabled draws no extra numbers, so an old model with the same seed gives the identical trajectory. For example, a `binding_capacity` that is never reached, or a model without a nucleus, replays exactly.

Changes to the core rules (such as the move to exact within-step promoter paths in v0.9) do change replay, and are recorded as such. Saved runs and cluster bundles are fresh starts with recorded seeds, not exact restart checkpoints.

## Backends and speed

| Backend | What runs where | Status |
|---|---|---|
| NumPy (reference) | Everything on the CPU. | Default; all features. About 2.1 ms per step for the 2,000-bead baseline. |
| OpenMM adapter | OpenMM integrates motion; the CPU still computes forces and reactions. | Works, but is *slower* on a GPU (2.8–3.6 ms per step) because data crosses to the device every sub-step. |
| OpenMM native | All forces and integration on the GPU. | Bare polymers, plus one-bead proteins whose only activity is attraction-mode wells (Section 3); refuses any other chemistry. About 60 µs per step for the bare baseline, about 35 times faster, and statistically equivalent to NumPy within 0.3%. Forces match the CPU engine to 10⁻¹² (Reference platform) and 10⁻⁴ (CUDA mixed precision). |

With attraction wells, the GPU uses a fixed sub-step count that allows for up to 12 close contacts per bead, and finer sub-steps for the first 500 steps, because the initial placement allows overlaps down to 0.65 of contact. The steep core is capped at half the contact distance (about 4,000 $k_BT$) so that it cannot overflow in single precision. Measured on an RTX 3090 that was shared with other runs, so the numbers are pessimistic: 0.4–0.7 ms per step for bare chains of 5,000–50,000 beads, and 9–16 ms per step for 10,000–50,000 beads plus 5,000–20,000 droplet proteins.

Models with chromatin binding, marks, cohesin or genes still run on the CPU. A native GPU path for that chemistry is future work (Section 13).

## Running at scale

From the interface you can export any draft model as a cluster bundle: the validated model, job settings, Slurm scripts and hashes. On the cluster, `python -m chromatin_studio.batch` checks the hashes and runs one replica per job, streaming results to disk with the actual backend recorded. The app prepares jobs but does not submit them.
