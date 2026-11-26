<!--
Section: Introduction and design
Update when: the scope, the architecture, the set of mechanisms, or the size limits change.
Code: chromatin_studio/engine.py (Simulation), schema.py (limits), server.py, static/app.js.
Sources: docs/PROJECT_STATUS.md "Intended scientific scope", docs/DEVELOPMENT.md "Architecture map".
-->

# Introduction: what the software is for

## The question

Genes are switched on and off by molecules that have to find them. A transcription factor diffuses through the nucleus, bumps into chromatin, sticks for a while, and leaves. Whether it sticks depends on the DNA, on whether a nucleosome is in the way, and on chemical marks on the histones. Enzymes it brings along can change those marks. Far-away enhancers can loop over to a promoter, helped or blocked by cohesin and by where the locus sits in the nucleus. The result, in the end, is a gene that fires in bursts.

Each of these steps has been studied on its own. What is hard is to see how they act **together**, because they happen on different scales of space and time and feed back on each other. Chromatin Studio is built to put them in one simulation where every step is explicit, every step can be turned off, and every number can be traced back to its source.

## The design in one picture

![**How the pieces connect.** Each box is a mechanism with its own section. Arrows show which quantity affects which. Most arrows can be switched off individually, which is how ablation studies work.](figures/fig-feedback.svg){#fig:feedback width=100%}

The simulation has two kinds of state:

1. **Positions** of beads in 3D: the chromatin chain, the protein molecules, and released RNA. These move by Brownian motion under forces (Section 2).
2. **Discrete states** on the genome and on molecules: which protein is bonded to which bead, which histone marks are present, whether a nucleosome sits on a bead, where cohesin legs are, whether a promoter is ON, where each polymerase is. These change by random events with set rates (Sections 3 to 7).

Each time step first moves the particles, then updates all the discrete states using the new positions. Section 8 explains why this split is reasonable and where it can fail.

## Reduced units

All inputs are in **reduced units**. Energy is measured in units of thermal energy ($k_BT = 1$). Length is measured in units of the chromatin bond length. Time is whatever unit makes the diffusion coefficients come out as entered. A bead can stand for 200 bp (one coarse nucleosome unit) or 10 kb (a chromatin blob); this is the setting `bp_per_bead`. Converting to nanometres and seconds needs a calibration against data (Section 11), which is not finished. Until then, compare runs with each other, not with seconds on a stopwatch.

## Size limits

The default model and the four small presets (Nuclear neighborhood, Multivalent assembly, Reader–writer memory, Loop extrusion) use one coarse nucleosome unit (200 bp) per bead; the 96 kb default is 480 beads. The 20 Mb presets stay at 10 kb per bead.

Model size is no longer capped. Above the sizes below, validation only warns that the CPU engine and the live view will be slow:

| Quantity | Warning above |
|---|---|
| Chromatin beads | 4,000 |
| Chromatin plus protein particles | 6,000 |
| Cohesin complexes | 256 |

Hard limits remain for individually tracked RNA (5,000), genes (256) and regulatory or tethering elements (128 each); up to 1,024 protein species and 256 marks are accepted. Measured cost per step on the CPU engine with full chemistry is about 3 ms at 480 beads, 35 ms at 10,000 and 220 ms at 50,000. A 20 Mb region at nucleosome resolution (100,000 beads) now validates, but is practical only on the GPU engine and for models it supports (Section 8).

## How the code is organised

The program is a small Python package with a local web server and a browser interface written without frameworks.

| Part | What it does |
|---|---|
| `schema.py` | Defines every setting, its default and allowed range, and rejects invalid models. |
| `engine.py` | Holds the state, computes forces, moves particles, runs all reactions, and writes read-only snapshots. |
| `transcription.py` | Promoter switching, polymerases and individual RNA. |
| `nucleus.py` | The ellipsoidal nuclear envelope. |
| `components/` | Single-mechanism studies with checks against references (Section 10). |
| `calibration.py`, `kinetics.py` | Tools to map reduced units to physical ones. |
| `openmm_backend.py`, `openmm_native.py` | Optional OpenMM integrators. |
| `server.py`, `static/` | The local web app. |
| `cluster.py`, `batch.py` | Export of a model as a job for a computing cluster. |

Three rules hold across the code. First, anything that only **reads** the state (snapshots, charts, rate tables) never draws a random number or changes the physics, so watching a run never alters it. Second, new mechanisms are **optional** and off by default, and a model that does not use them gives the exact same trajectory for the same seed as before. Third, the program never silently substitutes something: if you ask for a GPU backend that is not available, it stops with an error.
