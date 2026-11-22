# Component studies

A component study switches on **one mechanism at a time** and compares the engine with a reference. The purpose is to learn which parts of the model can be trusted, and against what, before mechanisms are combined. Results and their dates are in [PROJECT_STATUS.md](PROJECT_STATUS.md).

## Three levels of evidence

Every check in a report carries one of these labels. They are never merged.

| Label | Meaning | A pass shows |
| --- | --- | --- |
| `numerical_check` | The engine follows its own stated rule. | The implementation is consistent. |
| `analytical_reference` | The engine matches a closed-form result. | The mechanism is computed correctly. |
| `experimental_reference` | An output is compared with a published measurement. | The value lies in the stated experimental range **under the stated unit mapping**. |

An experimental pass is not a calibration. Several experimental checks are *input-matched* (the rate was set from the measurement, so the check only confirms unit conversion) or *tuned* (a scan point was selected to match). Each such check says so in its `note`. Only checks described as predictions test the model.

Studies that convert to seconds or nanometres use the provisional polymer mapping (`needs_review`, 2,000 beads × 10 kb only). Their comparisons are conditional on it.

## The studies

| Study | Components switched on | Reference | Cost at `standard` |
| --- | --- | --- | --- |
| `polymer_wlc` | Bare chain: bonds and bending | Exact discrete worm-like chain; DNA persistence length of about 50 nm | 8 replicas, about 30 min |
| `protein_diffusion` | One non-binding protein species, bead and IDR forms | MSD = 6Dt | minutes |
| `nucleosome_turnover` | Background eviction and reassembly | Two-state relaxation | seconds |
| `telegraph_transcription` | One gene, count-based and individual-RNA modes | Telegraph-model ON fraction, mean and Fano factor | about 1 min |
| `cohesin_extrusion` | One cohesin on a 1D lattice, with and without barriers | Growth speed, residence, loop size; in vitro extrusion speed (input-matched) | about 10 min |
| `tf_binding` | Chain plus one monovalent DNA-binding protein | Sox2 single-molecule tracking: free interval, trial count, bound fraction | 12 runs, hours |
| `nucleosome_tf_kinetics` | Chain with static nucleosomes at binding sites, one monovalent factor | Donovan et al. single-molecule kinetics for Reb1, Cbf1, Pho4 and Gal4: kinetic ratios (input-matched) and titrated affinity ratios (prediction) | 24 runs, several hours |
| `hp1_mark_density` | 12-nucleosome array with H3K9me3 marks, one monovalent reader | Kilic et al. 2015 HP1α residence and association versus H3K9me3 density and competing HP1α (one tuned hazard, the rest predictions) | 108 runs, about 65 min on 12 workers |
| `fbn2_cohesin` | Chain, cohesin, convergent barriers | Fbn2 anchor separation with and without cohesin and CTCF (data in `references/hansen`) | 9 runs, overnight |

Both nucleosome studies limit each binding-site bead to one bound factor with the optional `chromatin.binding_capacity` (see [CONTRACT.md](../CONTRACT.md)); without it a bead accepts any number of bonds, which would make site competition meaningless.

The three purely 1D studies (`nucleosome_turnover`, `telegraph_transcription`, `cohesin_extrusion`) advance reactions with frozen coordinates. Their reaction interval is a study parameter and may exceed the mechanical time step, because no mechanics is integrated.

## Commands

Run from the repository root with the application's Python environment. Output is run data: write it outside the repository (see `LOCAL_DATA.md`).

```sh
python -m chromatin_studio.components list
python -m chromatin_studio.components run polymer_wlc --scale standard --jobs 8 --output /data/component-studies/polymer_wlc
python -m chromatin_studio.components analyze /data/component-studies/polymer_wlc
python -m chromatin_studio.components show /data/component-studies/polymer_wlc
```

- `--scale smoke` runs in seconds and only proves the pipeline works. Its numbers are not evidence.
- `--set key=json` overrides one study parameter, for example `--set steps=2000000`.
- `--case ID` restricts a run to named cases.
- Rerunning the same command **resumes**: finished `(case, replica)` files are kept. A changed parameter or changed source code is refused; use a new output directory.
- `analyze` can be repeated while a long run is in progress. Partial results are counted and flagged.
- `fbn2_cohesin` needs equilibrated bare-polymer conformations: pass `--set 'warm_starts=["/path/a.npy", ...]'`. Without them it starts from a freshly grown chain and says so in the report.

On Windows with a conda environment, put both the environment root and its `Library\bin` on `PATH`; otherwise SciPy's sparse-graph import fails without a traceback.

Each output directory holds `plan.json` (parameters, cases, seeds, source hashes), `raw/` (one file per case and replica) and `report.json`. These are fresh runs with recorded seeds, not restart checkpoints.

Reports use format 2. `source_matches_current` says whether the source that produced the runs still matches the code that analyzed them, `changed_since_run` lists the files that differ, and `analysis_source_sha256` records the analyzing code. Analyzing again with different code moves the previous report to `report-history/`, so each report stays attributable to the code that produced it.

For runs that take hours, freeze the source: copy the `chromatin_studio` package (and `references/` when a study reads it) to a folder outside the checkout and run from that copy, so edits to the checkout cannot change a run in progress. Analyze afterwards with the current checkout; the report then lists what changed.

## Figures and the Evidence view

`analyze()` may also return `"figures"`, chart descriptions built with `core.series`, `core.band` and `core.figure`. They are data only. The interface draws every study with one generic renderer (`static/evidence-view.js`), so a study changes its plots in its own file.

- `series(label, x, y, role=, style=, low=, high=, group=)`. The `role` is `model` (engine output: filled marks in a hue), `reference` (closed form or model input: dashed gray) or `experiment` (published measurement: hollow diamonds, or a wide pale line). The `style` is `points`, `connected`, `line` or `step`. `low`/`high` are ±1 SEM ends unless the note says otherwise. Series that share a `group` share a hue, for example an engine curve and the measurement it is compared with.
- `band(label, role=, x=(low, high) or y=(low, high), category=)` is a shaded range; `None` leaves an end open.
- `figure(id, title, x_label=, y_label=, series=, bands=, x_scale=, y_scale=, categories=, note=)`. Axes are `linear` or `log`; `categories` makes the x axis categorical.
- One figure needs at most three hues (`core.MAX_COLORED_GROUPS`, counted by `core.colored_groups`): only the first three interface hues stay distinguishable for every pair under color-vision deficiency. Split a larger comparison into several figures. The end-to-end test checks this at tiny scale; if a fourth group slips through, the renderer draws its markers as squares.

Start the server with `--evidence <folder>` (or set `CHROMATIN_STUDIO_EVIDENCE`) and open the **Evidence** tab. It lists every `report.json` under that folder (up to four levels deep; older versions in `report-history/` are selectable from each report), draws the figures with hover values, keyboard stepping and a table view, and shows every check with its acceptance window, grouped by evidence level and tagged with the basis its note states (input-matched, tuned or prediction). The view only reads saved reports: it runs nothing and leaves the live model untouched. Absolute paths are shortened before they reach the browser.

## Changing or adding one study

Each study is one file in `chromatin_studio/components/`. A study imports only from `core.py` and the engine, never from another study. To change a study, edit its file; nothing else needs to change. To add one:

1. Create `chromatin_studio/components/<id>.py` with
   - `STUDY`: a dictionary with `id`, `title`, `components`, `question`, `evidence`, `sources` and `scales`. `scales` must contain `smoke` and `standard`; every scale needs `replicas` and `seed`.
   - `simulate(params, case, seed, save_partial)`: build a model, run it, and return JSON-serializable observables. Call `save_partial(observables)` occasionally in long runs.
   - `analyze(params, results)`: return `{"checks": [...], "findings": [...], "summary": {...}}`, and optionally `"figures": [...]` (see above). Build each check with `core.check(...)` and label its evidence level.
   - optionally `cases(params)`: a list of `{"id": ...}` dictionaries when the study compares several conditions.
2. Add one line to `REGISTRY` in `chromatin_studio/components/__init__.py`.
3. Add its shrunk parameters to `TINY` in `tests/test_components.py`. The existing tests then cover the module contract and an end-to-end run.

Rules that keep studies honest:

- Build models with `core.bare_config` and pass them through `core.validated`, so a study never runs a configuration the schema would reject.
- Read engine state; do not draw from `sim.rng` or change engine state from analysis code.
- Put only literature values you have checked into `sources`, with a DOI or URL. Mark context-only sources with `used_as="context"`.
- State in the check's `note` whether a comparison is input-matched, tuned or a prediction, and which unit mapping it uses.
- Do not edit an old report. Rerun into a new directory and record both.

`runner.py` and `__main__.py` know nothing about individual studies, and `core.py` holds only shared helpers. Changing those three files changes the recorded source hash of every study, so existing output directories can no longer be resumed.
