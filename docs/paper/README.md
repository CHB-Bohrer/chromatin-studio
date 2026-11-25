# The software paper

`chromatin-studio-paper.pdf` is a plain-language technical description of how Chromatin Studio works: the model, each mechanism, the numerics, the interface, validation against theory and experiments, limitations and roadmap. It is built from small, independent source files so that one change touches one section.

## Layout

| Path | What it holds |
|---|---|
| `sections/NN-name.md` | One section each, in file-name order. The comment at the top of each file says which code it describes and when it needs updating. |
| `data/*.json` | Small results the figures draw: collected component-study reports (checks and figure series only), the benchmark curves, the extra checks, and `recorded-results.json` (numbers quoted from `docs/PROJECT_STATUS.md` for runs whose data lives outside the repository). |
| `data/checks-table.md` | Generated summary table, included in Section 10 with `@include`. |
| `analysis/extra_checks.py` | Short validation analyses no component study covers (ideal chain, binding, marks, tethers). |
| `analysis/collect.py` | Copies the small, figure-relevant parts of run outputs into `data/` and regenerates the table. |
| `figures/make_figures.py` | Draws every figure from `data/`. Redraws only figures whose data or drawing code changed (`figures/.figure-cache.json`). |
| `figures/*.png` | Interface screenshots (real captures). |
| `build.py`, `style.css` | Assemble the sections and print the PDF. |

## Rebuild

Needs pandoc, matplotlib and Playwright with Chromium (`pip install matplotlib playwright && playwright install chromium`, or set `CHROMIUM_PATH` to an existing Chrome).

```sh
python docs/paper/build.py
```

This redraws stale figures and prints the PDF in a few seconds. It never runs a simulation.

## Updating after a code change (keep it cheap)

1. Find the section(s) whose header comment names the changed file. Edit only those.
2. If the change could alter a mechanism's numbers, rerun **only that mechanism's** check, outside the repository, then collect it:

   ```sh
   python -m chromatin_studio.components run cohesin_extrusion --scale standard --jobs 8 --output /data/val/cohesin_extrusion
   python docs/paper/analysis/extra_checks.py --output /data/val/extra.json --only tethers
   python -m chromatin_studio.benchmarks --output /data/val/reference-results.json
   python docs/paper/analysis/collect.py --runs /data/val
   ```

   `collect.py` only replaces the data files it finds, so the others keep their cached values.
3. Run `python docs/paper/build.py`. Only figures whose data changed are redrawn.
4. Commit the changed section, data, figures and PDF together with the code.

Documentation-only or interface-only changes need no reruns. Never commit raw run data.

## Cost of the validation runs (4 CPU cores, NumPy engine)

| Check | Time |
|---|---|
| `benchmarks`, `nucleosome_turnover`, extra checks except `ideal_chain` | under 2 min each |
| `telegraph_transcription`, `ideal_chain` | about 5 min |
| `cohesin_extrusion`, `protein_diffusion` | about 10–20 min |
| `polymer_wlc` (paper uses 1.5 M steps per replica instead of the standard 6 M) | about 1 h |
| `tf_binding`, `nucleosome_tf_kinetics`, `hp1_mark_density`, `fbn2_cohesin` | hours; results quoted from `PROJECT_STATUS.md` |

## Conventions

- Plain language, short sentences, a "Short version" box at the top of each section.
- Mark modelling choices without a measurement as **[Assumption]** and statements not checked against code or a run as **[Unverified]**.
- Every figure caption says where its numbers come from.
- Figures use at most three colours (blue, orange, aqua), dashed grey for closed-form references, hollow diamonds for measurements, and one y-axis per panel.
