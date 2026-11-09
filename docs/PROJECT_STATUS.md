# Chromatin Studio project status

Last reviewed: **2026-11-06**. Application: **0.10.0**. Configuration format: **0.1**.

The application is a functioning local research prototype with a visual model builder, coupled simulation, live inspection, tests and an editable tutorial. The latest completed feature work is live 1D protein binding, nucleosome accessibility and chromatin remodeling. The current development task is preparing the project for sharing and making its state understandable to future contributors.

The production presets remain in **reduced units**. No accepted empirical length/time fit, full-model GPU performance validation, 20 Mb nucleosome-resolution run or complete published-figure reproduction has been completed. These are unfinished milestones, not hidden capabilities.

## Intended scientific scope

Support reusable models across diverse loci, from a small nucleosome-unit region to a coarse 20 Mb domain. The intended system connects protein diffusion and multivalent binding, chromatin state and mechanics, enhancer/promoter regulation, cohesin, nuclear environments and stochastic transcription. Users should be able to inspect mechanisms, disable them for controlled comparisons, export reproducible cluster experiments and calibrate against matching experimental observables.

The user selected individual nascent transcripts and diffusing released RNA as the first RNA representation. RNA bead chains, residue-level IDRs and fully emergent nuclear bodies remain future extensions. The interface and tutorial are part of the scientific tool, not optional presentation work.

## Implemented capabilities and their boundaries

| Area | Current implementation | Boundary that matters |
| --- | --- | --- |
| Chromatin mechanics | Bead–spring polymer, bending stiffness, soft repulsion, state-selective attraction and adjustable soft confinement | Soft potentials do not guarantee preserved topology. Resolution and density require independent justification. |
| Proteins and assemblies | Many species, molecular diffusivities, finite valence, reversible chromatin and protein links, optional IDR domain chains | Larger assemblies can form; this does not establish a phase diagram or sequence-specific protein chemistry. |
| Marks and accessibility | Binary marks, writers/erasers, optional contact spreading, occupied/vacant nucleosome sites, occupancy-dependent persistence | No histone-copy identity, wrapped DNA, tails or explicit histone reservoir. |
| Remodelers | Directly anchored or TF-recruited eviction, reassembly and adjacent-site sliding | Sliding moves occupancy and marks, conserves nucleosome number and leaves genomic coordinates/protein anchors fixed. No ATP cycle. |
| Regulatory elements | Enhancer/promoter regions, TF preferences, explicit coactivator recruitment and configured contact-dependent rate gains | Mechanisms and parameters are hypotheses. No DNA sequence search or universally fitted enhancer model. |
| Transcription | Stochastic promoter switching and initiation; optional individual polymerases, elongation, pausing, RNA release/processing/decay | Rates are frozen within each reaction interval. No polymerase traffic/exclusion, finite polymerase pool or RNA-chain forces. |
| Chromatin feedback | Occupancy and engaged transcription can alter local persistence; polymerases can affect occupancy | These are configurable phenomenological couplings. |
| Cohesin | Stochastic loop extrusion with directional barriers, collisions and unloading/reloading | No independently parameterized loading rate or completed published contact-map reproduction. |
| Nuclear environments | Separately configurable speckle and nucleolar spatial affinity fields and transcription responses | Fixed imposed fields; no nucleolar subphases, rRNA biogenesis or self-assembled organelles. |
| Inspection | 3D view, live 1D binding/occupancy, full genomic tracks, interaction context, regulatory rate decomposition and individual RNA views | Read-only snapshots; free proteins retain 3D positions instead of invented genomic positions. |
| Studies | Independent control/ablation runs and reference calculations | Quantitative manuscript conclusions and full spatial matched-allele comparisons remain unfinished. |
| Cluster execution | Validated draft export, Slurm/job-array scripts, headless execution, streamed results and provenance | Fresh runs only; the app does not submit jobs or provide exact checkpoint continuation. |

Detailed definitions and coordinate conventions are in [CONTRACT.md](../CONTRACT.md) and [README.md](../README.md).

## Latest completed increment

The v0.10 increment is complete:

- **1D genome → Live binding → Binding detail** displays actual direct-bound molecules and occupied/vacant sites at suitable resolution. A selected site exposes association and release hazards for each species.
- Each protein has independent occupied/vacant factors for capture and chromatin-bond release. Gaps need not favor every species; mark-dependent readers can lose binding opportunities.
- The **Live nucleosome access · 25.6 kb** preset uses 128 × 200 bp sites with initial gaps, distinct TF/reader preferences and an explicit remodeler.
- Remodelers can evict, reassemble and slide occupancy. Eviction clears marks; reassembly begins unmarked; sliding carries marks to an adjacent vacant site.
- Structured genomic events include binding, release, sliding and occupancy observations. The engine retains up to 2,048 events; the browser caches up to 240 and displays a recent filtered subset. Turnover observations resolve reaction-interval endpoints; event IDs follow emission order across operator-split updates.
- Study ablations neutralize both occupancy-sensitive capture and release, and disable sliding when turnover is disabled.
- The tutorial was rebuilt with 16 chapters, close-up views, natural-paced narration and current event provenance. It is a scripted walkthrough of genuine captured states, not continuous simulation footage.

## Evidence from the latest completed checks

These are recorded results from the feature implementation, not a claim that every new checkout or future edit has already passed them. Compare report hashes with the source before reusing results.

| Evidence | Recorded result | Scope and location |
| --- | --- | --- |
| Python suite | 237 tests run: 235 passed, 2 skipped | CUDA/OpenCL hardware checks skipped. [v0.10 report](../benchmarks/nucleosome-access-v0.10.json) |
| Interface helpers | 75 tests passed | Node tests; distinct from visual browser inspection. Same [report](../benchmarks/nucleosome-access-v0.10.json). |
| OpenMM Reference smoke run | 200 steps; 11 checks passed; 8 bindings and 2 releases | CPU execution with finite coordinates/rates. No sliding occurred in this short run. Same [report](../benchmarks/nucleosome-access-v0.10.json). |
| Recorded live example | 4,900 NumPy steps; 171 bindings, 146 releases, 21 slides, 15 evictions, 7 reassemblies | Occupancy balances as 117 − 15 + 7 = 109; vacant sites are unmarked. Same [report](../benchmarks/nucleosome-access-v0.10.json). |
| Analytical polymer fixture | 14 of 14 saved checks passed | Gaussian-chain/loop reference fixture; not a production-model biological fit. [Reference results](../benchmarks/reference-results.json) |
| Study reference controls | 9 of 9 saved checks passed | Limited stochastic controls, not full manuscript reproduction. [Verification](../references/manuscript/verification.json) |
| Tutorial | 610.4 seconds; 18,312 decoded frames; 16 chapters; audio in every chapter | [Media validation](../tutorial/render-validation.json) and [source validation](../tutorial/render-source-validation.json). No human listening review is claimed. |

## Physical calibration is unfinished

Completed: MSD import/analysis/fitting tools; synthetic conversion tests; public Fbn2 reference processing across 12 conditions and 2,894 trajectories; a fitting/holdout split; and a proposed baseline simulation protocol. The references and conditions must remain distinct. See [reference inventory](../references/hansen/inventory.json), [processing verification](../references/hansen/verification.json) and [calibration verification](../benchmarks/calibration-verification.json).

Not completed: equilibrated baseline trajectories, defensible spatial mapping, an accepted experimental clock, held-out agreement and timestep/replica convergence. The [protocol](../references/hansen/jobs/protocol.json) explicitly records `prepared_not_run_not_calibrated`, `simulation_submitted: false` and `physical_mapping_assigned: false`. Four replicas of eight million steps are a proposed workload; their burn-in is not established.

Even after that fit, protein diffusion, residence times, binding probabilities, catalytic rates and transcription parameters need their own observations. One shared clock provides consistent units; it does not infer all those parameters. The [kinetics estimator](../chromatin_studio/kinetics.py) does not turn the interacting model's encounter hazards into measured macroscopic binding constants.

## GPU and scale status

The [OpenMM adapter](../chromatin_studio/openmm_backend.py) performs stochastic integration while host code calculates mechanical forces, reaction searches and RNA motion. Coordinates, forces and mobilities cross the host/device boundary every integration substep. Native GPU forces and full-model speedups are unimplemented/unmeasured; local Reference/CPU checks do not establish CUDA performance.

The [schema](../chromatin_studio/schema.py) permits 4,000 chromatin beads, 6,000 chromatin/protein particles and 5,000 individual RNA records. Therefore:

- 20 Mb at 10 kb per bead uses 2,000 chromatin beads and is supported.
- 20 Mb at 5 kb per bead uses 4,000 chromatin beads and reaches the chromatin limit.
- At 200 bp per bead, the maximum chromatin span is 0.8 Mb, subject to the total particle limit.
- 20 Mb at 200 bp per bead needs 100,000 chromatin beads and is unsupported. Increasing a schema constant alone is not a scaling implementation.

## Next milestones and completion criteria

This order is a proposed continuation plan. Follow the user's current priority rather than automatically starting all items.

| Priority | Work remaining | Completion evidence |
| --- | --- | --- |
| 1 | Complete the GitHub source release | The source is on the public remote. The first web upload committed the numbered upload folders; the 2026-11-07 layout fix restores the intended flat tree (see below). Remaining: merge the fix, confirm GitHub Actions passes, and select a license. |
| 2 | Complete baseline physical calibration | Run the matching model, establish stationarity, constrain length, fit MSD shape/time, evaluate held-out data, and record uncertainty plus timestep dependence. Attach a mapping only after the declared acceptance criteria pass. |
| 3 | Calibrate and isolate biological mechanisms | Fit protein mobility/residence and transcription independently; compare controlled binding, mark, extrusion and transcription limits before combining all mechanisms. Record which published observables have actually been reproduced. |
| 4 | Validate and optimize cluster/GPU execution | Run acceptance checks and end-to-end throughput measurements on real hardware. Move bottleneck forces/searches to a suitable backend with parity/convergence checks before claiming acceleration. |
| 5 | Support larger nucleosome-resolution systems | Demonstrate memory, neighbor-search, reaction, rendering and output scaling on increasing workloads; validate physics before raising bounds toward 100,000 beads. |
| 6 | Extend missing biology when it answers a specified question | Examples include emergent bodies, RNA chains/interactions, polymerase traffic, mechanochemical remodeling and paired spatial alleles. Define the representation and observables before implementation. |

Exact restart checkpoints and portable installed-package assets are also unfinished engineering work. The current reliable distribution method is a source checkout. `pyproject.toml` currently packages Python code while the server expects sibling `static`, `tutorial` and reference directories.

## Publication and private context

On 2026-11-06 the user authorized public release of the complete tool, examples, tutorial and handoff, excluding the full manuscript. The project is now a Git repository linked to the existing public `CHB-Bohrer/chromatin-studio` remote. Citation metadata, ignore rules, source-checkout instructions and Python/Node CI have been added. A fresh source copy passed 235 Python tests (2 GPU hardware skips), 75 interface tests and an independent app/tutorial launch check; see [release verification](../benchmarks/public-release-verification.json). GitHub Actions has not run yet. A code license awaits the maintainer's selection.

**Upload and layout fix (2026-11-07).** The release was uploaded through the GitHub website in four batches. Those commits kept the numbered batch folders (`01-application/` … `04-tutorial-narration/`) and the local upload helpers, and omitted the hidden `.gitignore` and `.github/workflows/tests.yml`. That broke relative paths in the code, tests and documentation. The fix moves every file to its release path with history-preserving renames, removes the upload helpers and adds `.gitattributes` to keep LF line endings on Windows. Afterwards 274 of the 276 release files are byte-identical to the v0.10.0 upload manifest at their intended paths. The original ignore and CI files were not available, so they were rewritten (Python 3.10/3.12 unit tests, NumPy acceptance validation and Node interface tests) and are not byte-identical to the release.

Manuscript hypothesis examples, short source attribution/questions and tutorial scenes are intentionally included under that authorization. The source PDF and any full-text extraction are excluded. A file/archive audit found no complete manuscript in the tool tree. Personal absolute narration paths in two generated metadata files were normalized, with original/published hashes recorded in [publication metadata](../tutorial/publication-metadata.json). Public reference datasets keep their own attribution and license records. See [release scope](RELEASE.md).

This handoff intentionally omits manuscript content, personal paths and machine-specific process IDs. Local session notes are maintained outside this project directory. They are optional context, not a prerequisite for using the public code or proof that a process is still running.

## Updating this status

After each substantive increment, update the date, feature status, evidence and next action together. Preserve explicit distinctions between implemented, tested, fitted and planned. Record the actual environment and source identity with new results. Do not carry forward a success count from an earlier version as fresh verification.

The handoff and source-release preparation were added on 2026-11-06 after the v0.10 feature and tutorial checks. Release packaging does not establish a new empirical calibration, GPU result or biological validation.
