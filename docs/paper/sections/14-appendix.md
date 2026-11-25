<!--
Section: Appendix (code map, how this document is maintained)
Update when: a section is added, renamed, or covers different code.
Sources: docs/paper/README.md.
-->

# Appendix: how this document is kept current {.unnumbered}

This document lives in `docs/paper/` of the repository. Each section is its own Markdown file in `docs/paper/sections/`, and every file starts with a comment that names the code it describes. When that code changes, only that section and its figures need updating.

| Section | File | Code it describes |
|---|---|---|
| 1 Introduction | `01-introduction.md` | `engine.py`, `schema.py` (limits), overall layout |
| 2 Chromatin | `02-chromatin.md` | `engine.py` mechanics and integrator, `chromatin_contacts.py` |
| 3 Proteins | `03-proteins.md` | `engine.py` binding and release, `kinetics.py` |
| 4 Marks and nucleosomes | `04-marks-nucleosomes.md` | `engine.py` marks, nucleosomes, sliding |
| 5 Cohesin | `05-cohesin.md` | `engine.py` cohesin |
| 6 Nucleus | `06-nucleus.md` | `nucleus.py`, `engine.py` tethers and bodies |
| 7 Transcription | `07-transcription.md` | `engine.py` regulation, `transcription.py` |
| 8 Numerics | `08-numerics.md` | `engine.py` stepping, OpenMM backends, `batch.py` |
| 9 Interface | `09-interface.md` | `static/`, `server.py`, `studies.py` |
| 10 Validation | `10-validation.md` | `components/`, `benchmarks.py`, `analysis/extra_checks.py` |
| 11 Calibration | `11-calibration.md` | `calibration.py` |
| 12–13 | `12-limitations.md`, `13-roadmap.md` | `docs/PROJECT_STATUS.md` |

**Figures and data are cached.** Simulation results used by the figures are stored as small JSON files in `docs/paper/data/`. The figure script redraws a figure only when its data or its drawing code has changed, and the build never reruns a simulation. A validation run is repeated only when the code of that mechanism changes, and only for that mechanism. See `docs/paper/README.md` for the commands.
