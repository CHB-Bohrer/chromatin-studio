# Chromatin Studio development guide

Start with [AGENTS.md](../AGENTS.md) and [PROJECT_STATUS.md](PROJECT_STATUS.md). The status file is the concise development handoff; [CONTRACT.md](../CONTRACT.md) defines implementation semantics and [README.md](../README.md) explains use and scientific interpretation.

## Architecture map

| Path | Responsibility |
| --- | --- |
| `chromatin_studio/schema.py` | Configuration defaults, presets, validation, limits and compatibility. |
| `chromatin_studio/engine.py` | Simulation state, forces, binding, marks, nucleosomes/remodelers, cohesin, stepping and read-only snapshots. |
| `chromatin_studio/transcription.py` | Promoter paths, initiation, elongation, individual RNA and associated summaries. |
| `chromatin_studio/neighbors.py`, `chromatin_contacts.py` | Spatial candidates and chromatin contact/attraction calculations. |
| `chromatin_studio/openmm_backend.py` | Optional polychrom/OpenMM integration and actual platform metadata. The CPU engine still evaluates the model. |
| `chromatin_studio/openmm_native.py` | Bare-polymer models only: the engine's mechanics as native OpenMM forces, with forces and integration on the device. Refuses models with chemistry. |
| `chromatin_studio/server.py`, `run.py` | Local HTTP API, applied workspace, static/tutorial serving and launch entry point. |
| `chromatin_studio/calibration.py`, `kinetics.py` | Physical mapping/MSD tools and separate idealized binding estimates. |
| `chromatin_studio/studies.py`, `benchmarks.py`, `validation.py` | Independent ablations, analytical reference fixtures and numerical backend checks. |
| `chromatin_studio/cluster.py`, `batch.py` | Portable job plans/bundles and headless per-replica execution. |
| `static/app.js` | Draft model editing, API coordination and integration of the browser views. |
| `static/genome-view.js`, `genome-live.css` | 1D genomic tracks, binding/occupancy detail and recorded-event display. |
| `static/interaction-view.js`, `regulation-view.js`, `transcription-view.js` | Protein context, regulatory dependencies and individual RNA inspection. |
| `static/calibration-view.js`, `cluster-view.js`, `study-view.js` | Calibration, cluster planning and isolated study interfaces. |
| `static/learn-view.js`, `model-explain.js` | Explanatory illustrations and teaching controls, separate from measured run state. |
| `examples/` | Complete editable model JSON specifications. |
| `tests/`, `benchmarks/`, `references/` | Implementation tests and recorded scientific evidence; provenance and scope matter. |
| `tutorial/storyboard.json`, `tutorial/README.md` | Editable narration/capture plan and the regeneration workflow. |

The browser edits a draft. Apply/restart validates it and creates an applied simulation. Sequential step requests advance that run; snapshots feed the 3D, 1D, interaction and RNA views. Studies create separate simulations. Cluster export copies the draft into a fresh-run bundle. Avoid coupling diagnostic/UI calculations back into the stochastic state.

## Start from a source checkout

From the project root on macOS/Linux, with Python 3.10 or later:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run.py --config examples/nucleosome_access.json --port 8766 --open
```

On Windows, use the equivalent `.venv\Scripts\python.exe` paths. Choose an unused port; do not replace another running experiment. The server binds to loopback by default. This is a local application, not a hardened public multiuser service.

The optional backend has a separate dependency file:

```sh
.venv/bin/python -m pip install -r requirements-polychrom.txt
```

No account, browser framework build or voice model is needed to run the simulator. Narration regeneration has a separate environment. Preserve the source directory layout until installed-package assets are fixed.

## Verification commands

Use the environment that runs the application. For the full Python suite, the tutorial-renderer tests also need Pillow and imageio-ffmpeg. Install them in the development environment, not as a requirement for ordinary simulation users:

```sh
.venv/bin/python -m pip install -r requirements-dev.txt
.venv/bin/python -m unittest discover -s tests -v
node --test tests/*ui.cjs
```

The Node built-in test runner requires an installed compatible Node.js runtime. Optional-backend checks can skip without their dependencies or compatible hardware. Report the actual pass/skip counts, reasons and platform; compare with the dated evidence in [PROJECT_STATUS.md](PROJECT_STATUS.md), not a promised fixed count. HTTP tests use loopback sockets.

For a focused nucleosome change:

```sh
.venv/bin/python -m unittest discover -s tests -p 'test_live_nucleosomes.py' -v
node --test tests/test_genome_ui.cjs
```

Other changes require their corresponding tests; a numerical integration change also needs distribution/convergence checks. These checks produce new reports without overwriting saved historical results:

```sh
.venv/bin/python -m chromatin_studio.validation --platform numpy --output /tmp/chromatin-validation-new.json
.venv/bin/python -m chromatin_studio.benchmarks --output /tmp/chromatin-reference-new.json
```

Use unique output paths for retained evidence. A CPU pass is not a substitute for the actual target-platform validation. A documentation-only edit normally needs link/factual checks, not these simulation workloads.

## Tutorial changes

Follow [the tutorial guide](../tutorial/README.md). Edit the storyboard and transcript together, recapture affected views from real simulations, regenerate changed narration, then render and check the assembled result. Verify image dimensions before setting focus rectangles. Preserve model/run provenance, the natural 1× voice pace and the distinction between draft, applied model, illustration and measurement.

Do not regenerate the entire video for a development-handoff-only change. After a behavior/UI change, identify which chapters are affected rather than silently leaving contradictory instructions in the video.

## Maintain a useful handoff

For each substantive change, record:

1. The concrete behavior completed, including changed defaults or semantics.
2. What ran, its outcome, source/version identity and the limits of the evidence.
3. Work still incomplete, with the next action and its completion criteria.
4. Changes needed in the contract, examples, UI and tutorial.

Keep process IDs, absolute personal paths, private manuscript details and temporary environment instructions in a local note outside the release directory. Run data (trajectories, analysis reports and the scripts that produced them) is also kept outside the repository. Each machine records its data location in a gitignored `LOCAL_DATA.md` at the repository root; read it before working with run data. Keep durable scientific/engineering state in the project. Never depend on another AI having access to prior chat messages.

A reusable opening instruction for another collaborator is:

> Read AGENTS.md, docs/PROJECT_STATUS.md, docs/DEVELOPMENT.md and CONTRACT.md. Verify the current source and evidence before claiming completion. Continue with the task I specify, preserving any active simulation. Distinguish implementation, numerical verification, empirical calibration and biological validation. Update the handoff when finished.

## Release maintenance

The user authorized public release of the tool and examples while excluding the full manuscript. Preserve that distinction and the data attribution recorded in [RELEASE.md](RELEASE.md). Citation metadata, ignore rules and automated CPU/interface checks are included. Keep the licensing status explicit until a software license is selected. Test a fresh checkout when changing distribution files. The current standalone wheel omits required UI/reference assets; source-checkout execution is the supported initial path. See [the milestone table](PROJECT_STATUS.md#next-milestones-and-completion-criteria).
