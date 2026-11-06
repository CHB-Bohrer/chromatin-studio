# Chromatin Studio

A local visual builder and exploratory simulator for chromatin, reacting proteins, multivalent assemblies, histone marks, cohesin and transcription. The model is configurable across loci; the included biological names are illustrative species classes, not fitted molecular parameterizations.

**This is a working research prototype in reduced units.** It now includes an explicit MSD calibration workflow and reproducible polymer reference benchmarks. No experimental clock is attached to the supplied presets. Numerical agreement with a reference model is separate from calibration and biological validation.

**Continuing development with an AI or a new collaborator?** Start with [AGENTS.md](AGENTS.md), the [current project status and unfinished milestones](docs/PROJECT_STATUS.md), and the [development guide](docs/DEVELOPMENT.md). They record what is implemented, what the evidence supports, and how to continue without relying on earlier chats.

## Open the application

Use Python 3.10 or later. On macOS or Linux, download this repository with Git and start the local app:

```sh
git clone https://github.com/CHB-Bohrer/chromatin-studio.git
cd chromatin-studio
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run.py --open
```

On Windows, use `.venv\Scripts\python.exe` instead of `.venv/bin/python`. Downloading and extracting the repository ZIP also works; run the environment/install commands inside the extracted project directory. The interface opens at [http://127.0.0.1:8765](http://127.0.0.1:8765). Keep the terminal open while using it; Control-C stops the server. On macOS, after setup, **Open Chromatin Studio.command** provides a double-click launcher.

The server listens only on this computer. It does not upload models or simulation results. The repository distributes the application; GitHub Pages alone cannot run its Python simulation server.

To enable the genuine polychrom/OpenMM integration backend:

```sh
.venv/bin/python -m pip install -r requirements-polychrom.txt
.venv/bin/python run.py --open
```

The optional dependency file pins the upstream polychrom commit used here. Source checkout execution is supported. If using `pip install -e .`, keep the source directory layout intact; standalone wheel assets are not yet packaged. To use a different port, pass `--port 8766`. For a compact binding/remodeler demonstration, add `--config examples/nucleosome_access.json`.

See [release scope and attribution](docs/RELEASE.md), [citation metadata](CITATION.cff), and [development/test instructions](docs/DEVELOPMENT.md). The tool, hypothesis examples and tutorial are public; the complete source manuscript is not distributed.

## Tutorial and experimental references

Use **Tutorial** in the application header, or open [the offline player](tutorial/index.html). The captioned walkthrough uses natural narration at a comfortable speaking pace. Chapter links let you revisit individual topics, and the player offers adjustable playback speed; start at 1× and choose the pace that suits you. It uses actual application captures: the 20 Mb overview is coarse; the 80 kb manuscript-hypothesis example demonstrates regulation and individual RNA; and the 25.6 kb nucleosome-access example shows live 1D binding, occupied/vacant sites, species-specific capture/release and remodelers. Preset changes are explicit. The [editable storyboard](tutorial/storyboard.json), [renderer instructions](tutorial/README.md), transcript, captions and original captures let you regenerate only changed scenes. This is a scripted keyframe walkthrough, not an uninterrupted recording of simulation dynamics.

**Calibrate → Import measured MSD** now imports public Hansen lab Fbn2 measurements, keeping all twelve conditions separate. The 17 curve variants include raw and noise-corrected RAD21 2h references and a movie-level fitting/holdout split. The fitting subset contains 424 trajectories from 74 movies; 113 trajectories from 18 different movies remain held out. Import retains units, sources, processing, counts and distance statistics. It neither assigns a physical clock nor changes the live run.

The [study inventory](references/hansen/inventory.json) distinguishes retrieved/processed data from studies only inventoried. All twelve Fbn2 conditions (2,894 trajectories) have been processed; newer MINFLUX, coupled enhancer–RNA and protein datasets have **not all been reanalyzed**. Publicly available conditions cannot be pooled as if they were replicates of one biological state. Data redistribution and processing provenance are recorded in [licenses](references/hansen/licenses.json), [source methods](references/hansen/source-methods.json) and [verification](references/hansen/verification.json).

A historical [20 Mb calibration cluster bundle](references/hansen/jobs/baseline-20mb-10kb-cuda.zip) and [acceptance protocol](references/hansen/jobs/protocol.json) preserve the saved polymer mechanics and turn off active biological mechanisms for a candidate RAD21-depleted baseline. That archived bundle contains **v0.8.0 code**, not the current v0.10 implementation; it is retained for provenance. Export a new bundle from the current application for current-code work, and recheck the protocol/model match. Four seeded replicas, eight million steps each, are a **proposed workload**, not an established equilibration time or a completed run. The 10 kb resolution represents the 515 kb label separation with separate 510/520 kb brackets. A 5 kb alternative resolves it exactly but needs its own parameterization and clock. The production model remains uncalibrated until stationarity, length matching, MSD shape, held-out predictions and timestep convergence pass.

## Using the workbench

1. Select a preset to initialize a new run. **80 kb nucleosome-unit locus** demonstrates TF recruitment and transcription at 200 bp per bead; **20 Mb regulatory landscape** demonstrates the same couplings at 10 kb per bead. The additional `nucleosome_transcription` and `transcription_20mb` presets enable nucleosome turnover and individual RNA tracking at those respective scales. Earlier presets remain available.
2. Explore the component editors. Add protein species, choose bead or flexible IDR-chain representation, set molecular abundance and valence, then specify pairwise interactions and chromatin affinities.
3. Define chromatin marks, initial marked regions, genes, directional barriers and nuclear landmarks. Use the left **Regulation** editor to add enhancers and promoters, and **Chromatin** to add state-selective chromatin attraction or enable nucleosome occupancy. **Proteins** exposes occupied/vacant binding preferences and remodeling rates; **Genes** enables individual transcripts, elongation and pausing. Transcription's persistence multiplier is a hypothesis switch: **1 means no mechanical effect**.
4. Apply the configuration to reset the system and initialize a reproducible run. Edits are not applied silently to an existing trajectory.
5. Run, pause, inspect particles, and watch the live measurements and reaction events. Rotate and zoom the scene to understand the component relationships.
6. Export the model JSON to reuse it, or export the experiment including the current snapshot and recent measurement history. Import accepts a model JSON or the model portion of an experiment export.
7. Use **Cluster job** (or **Cluster** in the builder) to prepare an independent batch experiment from the **draft model**. Preview and download the ZIP without changing the local run. Configure the target platform and resource requests before transferring the bundle to your cluster.

The **Live nucleosome access · 25.6 kb** example provides a compact **25.6 kb, 128 × 200 bp** experiment with TFs, occupied/vacant binding preferences and explicit remodeler species. Use it to inspect binding and occupancy changes before adding a full regulatory landscape. Its parameters are illustrative, not an empirical fit.

The example JSON files in `examples/` are the full model specifications, including `regulatory_20mb.json`, `nucleosome_locus.json`, `nucleosome_transcription.json`, `transcription_20mb.json` and `nucleosome_access.json`. They can also be loaded with `python run.py --config examples/nuclear.json`.

## Understanding the model visually

Use **How it works** in the header, **Learn** in the builder, or the contextual explanation above each editor. The illustrated guide links molecular motion, encounters, chromatin accessibility, regulation and RNA output, with separate lessons for architecture and validation. Its controls are local teaching examples: they do not edit the draft, apply changes, advance the engine, or contribute samples to the live charts. **Configure** links take you to the relevant editor; **Return to simulation** restores the workbench.

Examples distinguish reduced units, microscopic encounter hazards and hypothetical physical-rate estimates. An ideal diffusion-limited association calculation does not calibrate the interacting polymer model's Kd. Likewise, the nucleus backdrop is decorative, while the central 3D view, genome tracks and observables use the applied model and actual run state.

The workbench's dark nucleus imagery and serif headings are inspired by [the Bohrer Lab website](https://chb-bohrer.github.io/). All illustrations are local SVG or canvas drawings; the interface requires no external font, image or visualization service.

## Physical length and time calibration

Open **Calibrate** in the builder. The workflow fixes physical length first, then estimates a single clock from matching simulated and measured MSD curves:

```text
MSD_physical(t) = (length_nm / 1000)² × MSD_simulation(t / time_s)
physical diffusion [µm²/s] = reduced diffusion × (length_nm / 1000)² / time_s
physical hazard [1/s] = reduced hazard / time_s
elongation [nt/s] = reduced elongation [beads/time] × bp_per_bead / time_s
```

`time_s` means seconds per **reduced time unit**, not seconds per integration step or saved frame. A step lasts `simulation.dt × time_s`. Genomic resolution does not establish a nanometer scale. Length and clock are degenerate in a single power-law MSD amplitude, so provide independent spatial evidence for length (for example distance statistics at a matching genomic separation). A plateau can constrain length only when the same stationary observable and model geometry justify comparing it; a nearly flat curve cannot identify a clock.

1. Choose an experimental condition and matching model baseline. Disable unsupported mechanisms for the initial comparison, equilibrate, and record independent replicas. A short animation is not an equilibrium trajectory.
2. Match the measurement: single-locus MSD, selected-locus center-of-mass MSD, or **two-locus separation-vector MSD**. The latter measures changes in `r1 − r2`, not changes in its scalar length and not the average of two locus MSDs. Match 2D versus 3D, genomic separation, drift convention, localization noise and exposure blur. Whole-polymer COM subtraction removes genuine collective motion; it is not automatically a stage-drift correction.
3. Import two-column `lag,msd` CSV or analysis JSON. Simulation lag/MSD are reduced; experimental lag is seconds and MSD is µm². Supply source identifiers, length evidence and a fit window. The fitter uses log-space interpolation, requires complete window coverage and does **not** extrapolate.
4. Inspect the overlaid curves, log residuals and slopes. A fit with poor shape agreement, a boundary optimum or an unidentifiable clock is marked `needs_review`; the UI does not apply it. Thresholds are explicit diagnostic choices, not a statistical goodness-of-fit test. No confidence interval is invented from correlated lag samples. Use independent replicas and sensitivity to length, window and imaging corrections.
5. Attach an acceptable mapping to the **draft**, inspect the implied rates, and export or apply it to a new run. Mapping metadata does not change reduced dynamics. Synthetic examples remain demonstrations and cannot be attached through the UI as empirical calibration.

The saved profile includes both input curves, their conventions, sources, fitting window, diagnostics and baseline model hash. Hashes exclude the replica seed, display name and calibration itself. The comparison is conservative: different spellings of optional defaults can flag a mismatch. For clean provenance, generate trajectories from the same exported configuration used for fitting. Changed physical parameters retain the original baseline clock and generate a review warning; perturbations must not be automatically refitted to erase their predicted mobility differences.

The conversion summary covers chromatin and free-protein diffusion, RNA diffusion, persistence and bond lengths, reaction hazards, cohesin movement and transcription elongation. A clock fit does **not** fit these parameters independently. Protein binding hazards remain per-encounter inverse-time rates; they do not become molar `kon` or measured `Kd`. The interaction estimator can explicitly load the applied model's fitted scales while retaining its independent-sphere assumptions.

### Measure MSD from cluster trajectories

Cluster jobs must save snapshots at a suitable regular interval. Metrics alone, and a **Save run** export with only one coordinate snapshot, cannot supply a trajectory MSD. From an extracted job bundle or this project directory:

```sh
python -m chromatin_studio.calibration msd \
  --snapshots results/replica-0/snapshots.jsonl \
  --beads 10,20 --dimension 2 --observable relative_pair \
  --discard-before 10 --output simulation-msd.json
```

These bead indices and discard time are command examples, not calibrated choices. Use `--observable locus` or `center_of_mass` when appropriate; `--drift-correction center_of_mass` subtracts the full supplied polymer COM. The analyzer records selected beads, sample spacing, exclusions and overlapping pair counts. It finds the corresponding `model.json` beside the trajectory unless `--model` specifies another file. It rejects irregular spacing; a last partial interval can be excluded explicitly with `--drop-final-irregular`. Select an equilibration cutoff and sampling interval appropriate to the requested lag window. The initial analyzer is bounded to 256 MiB of snapshot input and 12 million coordinate entries; larger trajectory analysis needs a selected window or a future streaming analyzer.

The fitted profile is carried in `model.json`. Portable bundles and batch results also contain `physical-units.json`, and the run manifest records the same mapping, status and baseline comparison. Raw output coordinates and times stay reduced so conversions remain reproducible. No profile means **uncalibrated**, rather than an implicit one-second time unit.

### Reference models and study reproduction

The saved report at `benchmarks/reference-results.json` contains actual stochastic runs, analytic predictions, seeds, parameters, tolerances and source hashes. Reproduce it with a new output filename:

```sh
python -m chromatin_studio.benchmarks --output benchmarks/new-reference-results.json
```

The first suite exercises the application's NumPy Brownian integrator under controlled Gaussian-chain forces: free diffusion, whole-chain `D/N`, finite-chain internal MSD, intermediate subdiffusion, timestep bias and static-loop suppression. It uses independent equilibrium chain replicas and compares sampled MSD with both continuous and exact finite-Euler predictions. The saved run passed **14/14** checks: internal Rouse slope **0.489** versus **0.503** predicted, and looped/unlooped internal MSD **0.716** versus **0.720** predicted. This fixture bypasses production forces, reactions, adaptive subdivision and GPU integration; those have separate implementation checks.

This is a limited reference-model reproduction grounded in the Gaussian/loop theory of [Yuan et al. (2024)](https://doi.org/10.1103/PhysRevE.109.044502). Its small chain and fixed loops are specified test choices, not that paper's full loop-extrusion model or a recreated figure. The production app's finite-rest springs and optional interactions are not the Gaussian fixture. The UI labels full published-figure reproduction as **not yet done**.

The [Hansen Lab calibration notebook](https://github.com/ahansenlab/polysim_bootcamp/blob/main/tutorial/calibration.ipynb) is a useful workflow precedent: it compares two-locus relative motion, obtains a spatial scale from plateau statistics and adjusts time scaling. Its experimental input is referenced through a local path and is not bundled. Neither its clock nor its noise correction can be transferred to this model without the matching data and measurement conventions.

The next reproduction targets should be isolated before combining all mechanisms:

| Target | Observables to match | Remaining work |
| --- | --- | --- |
| Yuan et al. loop-dependent dynamics | MSD curves and local slopes with specified loop statistics | Reproduce the published distributions, parameters and selected figure; current static-loop fixture is only the starting test. |
| [Fudenberg et al. (2016), loop extrusion](https://pmc.ncbi.nlm.nih.gov/articles/PMC4889513/) | Contact maps, contact probability versus genomic separation, boundary orientation/deletion and cohesin removal | Match extrusion rules, density, processivity, polymer forces and the relative 1D/3D clock before comparison. |
| [Michieletto et al. (2016), epigenetic recoloring](https://doi.org/10.1103/PhysRevX.6.041047) | Compaction, mark ordering, domain stability and hysteresis | Implement a separately specified compatible recoloring benchmark. Current independent mark flags and catalytic hazards are not the paper's recoloring model. |
| Binding and transcription limits | Dilute binding curves/residence times, telegraph-model RNA moments, elongation delays | Calibrate microscopic reaction rules and measured kinetics after the polymer clock, then test coupled responses on held-out conditions. |

For a quantitative biological model, the next required evidence is a chosen experimental tracking/spatial dataset and a matched equilibrated simulation. The general workflow is ready; the presets have not been fitted to such a dataset.

## Inspecting chromatin in 1D

Use **1D genome → Live binding** to inspect the same applied model along its genomic coordinate axis. Navigate to a gene or select a bead, then choose **Binding detail** for a close view of roughly two dozen sites. At the 200 bp mapping, sufficiently zoomed views show occupied nucleosome symbols, vacant sites and actual bound protein molecules. Wider views aggregate occupancy and counts; at 10 kb per bead, each occupancy state is an effective bin rather than an individual nucleosome. **All tracks** restores genes, regulatory elements, loops, marks, persistence, polymerases, affinities and occupancy together. Zoom, pan or enter a coordinate interval to move along the locus. Select a feature or bead to connect the genomic and 3D views.

The view updates from actual simulation snapshots while running. Pause to inspect a consistent frame. Protein markers show direct molecule–chromatin attachments; spatial proximity and an indirect path through another protein are not direct occupancy. Multiple domains of an IDR molecule attached at one site still count as one molecule. These stable molecule identities describe proteins within the current run; nucleosome glyphs describe occupancy at genomic sites, not tracked histone identities.

Select a site to see per-species **Association / t** and **Release / t**. Association is the hazard for an eligible encounter at that bead, not an arrival rate including concentration or diffusion. Release is the per-bond hazard at the site’s current occupancy. Both use reduced inverse time unless an empirical physical mapping is attached; the interface does not silently call reduced time seconds.

The tracks distinguish several different quantities:

- **Genes:** configured genomic intervals with current promoter state. Gene coordinates use the bead-to-base-pair mapping described below.
- **Regulatory elements:** configured enhancers and promoters linked to their target genes. Selecting an element connects the genome track, 3D highlight and Regulation view.
- **Cohesin loops:** current pairs of loop anchors. Configured directional barriers are shown separately; they are potential stopping sites, not automatically occupied anchors. Loop-list indices are not persistent motor identities after unloading and reloading.
- **Histone marks:** current bead-level flags, including changes from writers, erasers and spreading.
- **Nucleosomes:** occupied/vacant bead states when occupancy is enabled. These are individual coarse nucleosome units only at an appropriate genomic mapping; a 10 kb bin is not one physical nucleosome.
- **Persistence:** current local persistence values, including occupancy and transcription multipliers, in reduced length units.
- **Polymerases:** current contour positions; individual tracking supplies stable identities, pause states and growing RNA lengths.
- **Protein affinity:** the exact dimensionless encounter-rate multiplier used by the engine at each bead. These values are not equilibrium dissociation constants or bound fractions. For each species, the highest applicable mapped mark affinity is used; its unmodified value is the fallback only when no mapped mark is present. An explicit zero overrides that fallback, although another positive mapped mark can still permit association. Enabled TF elements multiply this result by their region preference; overlapping elements for the same TF use the largest region multiplier. If nucleosome occupancy is enabled, the species' occupied/vacant affinity factor multiplies the result again.
- **Protein occupancy:** distinct molecules of each species directly bound to a bead. Multiple attached IDR domains from the same molecule count once. Protein–protein links alone do not count as chromatin occupancy. A molecule can remain bound after the bead's current association affinity has changed.

**Recorded binding & remodeling events** reports engine transitions with the affected bead, species/molecule when available, event identity and a sampled time or explicitly labeled interval endpoint. Binding, release and sliding record their sampled conditional clocks; eviction, reassembly and reset records summarize the turnover interval. Event IDs order recorded observations, not a globally interleaved chemical timeline: reaction families remain operator split. Binding and release are not inferred by comparing two display frames. Occupancy-turnover entries report endpoint changes or a mark-reset interval; turnover can occur and return to the same occupancy within a step. The bounded recent-event buffer is not a complete trajectory, and its retained/dropped counts make truncation visible.

The viewer uses **applied configuration and live simulation state**. Unapplied editor changes do not silently alter the plotted experiment. Affinity and occupancy are supplied by the engine rather than inferred from spatial proximity. Numeric details remain available even when an overview compresses many beads into a pixel.

Exported snapshots include `genome_tracks`, with `bp_per_bead`, `span_bp`, per-species `affinity_runs`, `occupancy_runs`, `on_rate_runs` and `off_rate_runs`, `persistence_runs`, and `nucleosome_occupancy_runs`. `binding_sites` retains each direct domain bond’s bead, particle, molecule, species and effective off-rate; the UI deduplicates molecule markers at each site. Snapshots also include `genomic_events` with run identity, retained events, total/dropped counts and the latest sequence ID. Each run is `[start_bead, end_bead_exclusive, value]`; available tracks cover the full chromosome, including zeros. The nucleosome track is empty when the feature is disabled, rather than claiming measured occupancy. This compact representation uses half-open bead intervals even though configured gene and initial-mark intervals have inclusive end-bead indices.

## Building and inspecting regulatory pathways

In the left **Regulation** editor, add an enhancer or promoter and choose its bead interval, target gene, TF species and optional coactivator. **TF region affinity ×** multiplies that TF's current mark-dependent encounter affinity only inside the enabled element. A value of 1 adds no preference; 0 suppresses association there unless another overlapping element supplies a larger multiplier. A region multiplier cannot rescue a zero mark-dependent affinity. These are specified site preferences, not DNA sequence scanning or motif-energy calculations.

Configure a TF–coactivator pair in **Interactions** so the coactivator can actually bind the TF. In **Proteins**, select the coactivator's writer or eraser mark and catalytic rate if it should modify chromatin. Recruitment requires a direct protein–protein bond to a TF that is directly bound within the element; proximity or an indirect path through other proteins is insufficient. Counts represent distinct molecules, not the number of IDR domain links. Recruited catalysis acts at the TF's anchor beads, with each coactivator counted once per bead, including when it also binds that bead directly. This is a tethered-complex approximation without separate catalytic-domain reach geometry.

An enhancer contributes transcriptional gains when its closest bead is within its **contact radius** of the target promoter interval. Without an enabled promoter element, the target is the gene's start bead. A promoter element acts locally and must contain its target gene's start bead; at most one promoter per gene can be enabled. Enhancer-local mark writing can proceed without enhancer–promoter contact. It does not automatically transfer the mark to the promoter.

For each enabled element, the dimensionless signal is `TF count + coactivator gain × recruited count`. When its contact gate passes, its activation contribution is `activation gain × signal`, and its initiation contribution is `initiation gain × signal`. Contributions from elements targeting the same gene add:

```text
promoter ON rate = base ON rate × mark gate × (1 + sum of activation contributions) × body activation factor
initiation rate when ON = base initiation rate × (1 + sum of initiation contributions)
```

Required/repressive marks at the **gene's start bead** gate the ON rate, not the initiation rate of an already ON promoter. Losing enhancer contact removes its rate gains without forcing the promoter OFF or erasing marks. Promoter switching and mark turnover therefore retain memory. Set activation and initiation gains independently to zero to isolate direct regulatory stimulation from modification-mediated effects. A gene with enabled elements uses these regional gains instead of its legacy gene-wide activator-count gain; a gene without enabled elements retains the legacy behavior.

Use the top **Regulation** view to select a gene and inspect current TF occupancy, recruited coactivators, mark fractions, enhancer contacts, effective rates, nascent polymerases and processed/unprocessed RNA. **Locate in 1D** and **Highlight in 3D** connect the selected element to the spatial views. The arrows organize configured pathways; they are not a measured causal sequence. Rates are per reduced time, distances are reduced lengths, gains are dimensionless, and mark percentages are fractions of beads. The view uses the applied model, including while editor changes are pending. Snapshots export the same measurements in `regulation` and per-gene `regulation` records.

These couplings are hypotheses with experimental precedents. Targeted p300 recruitment can acetylate and activate promoters/enhancers ([Hilton et al., 2015](https://doi.org/10.1038/nbt.3199)), but contact–transcription relationships depend on the system: sustained proximity supported activation in an engineered Drosophila reporter ([Chen et al., 2018](https://doi.org/10.1038/s41588-018-0175-z)), whereas contemporaneous Sox2 transcription did not correlate with enhancer proximity in the measured mouse stem-cell trajectories ([Alexander et al., 2019](https://doi.org/10.7554/eLife.41769)). Neither result establishes a universal contact gate or excludes unresolved transient encounters.

### Which transcription factors and stochastic processes are present?

The example species named **Activator** is an explicit diffusing transcription factor, with its own abundance, diffusion, chromatin binding and protein–protein links. Its ID remains `activator`; the name does not imply measured parameters for a particular biological TF. The builder and live protein inspector derive its **Transcription factor** role from the configured genes or enabled regulatory elements. Several species can play this role, and a TF can also have catalytic or remodeling activity.

| Preset ID | Span and resolution | Explicit Activator molecules | Regulatory elements | RNA representation |
|---|---|---:|---:|---|
| `twenty_mb` | 20 Mb; 2,000 × 10 kb beads | 50 | None; gene-wide activator gain | Aggregate counts |
| `regulatory_20mb` | 20 Mb; 2,000 × 10 kb beads | 50 | 12 | Aggregate counts |
| `transcription_20mb` | 20 Mb; 2,000 × 10 kb beads | 50 | 12 | Individual polymerases and diffusing RNA |
| `nucleosome_transcription` | 80 kb; 400 × 200 bp beads | 10 | 4 | Individual polymerases and diffusing RNA |

These are different model presets. Enabling a regulatory element does not automatically enable individual RNA, and none of these 20 Mb presets resolves individual nucleosomes. The **Experiment** and **Genes** panels show the current draft's regulator assignments and RNA mode; the scene inspector describes the applied model.

For effective promoter ON rate `a`, OFF rate `b`, and initiation rate `r` frozen during a reaction interval, the engine now samples the complete two-state switching path with exponential waiting times. If the realized ON duration is `T_on`, the number of initiations is `Poisson(r × T_on)`; their times are sampled within the actual ON intervals. New polymerases elongate only for the part of the interval after their birth. This fixes the earlier endpoint-state approximation, which assigned a full interval of initiation to a promoter found ON at its endpoint.

This promoter/initiation process is exact conditional on frozen rates. The coupling to moving chromatin, occupancy, protein binding and marks remains operator split. Mechanical adaptive substeps do not refine those reaction updates. Reduce the global timestep and check convergence of burst statistics, contacts and RNA output. Event-budget failures stop with an explicit diagnostic; live-RNA capacity blocks attempts with counters and a warning. Quantitative runs should have zero capacity-blocked initiations.

In individual-RNA mode, pause/resume and released-RNA processing/degradation use exponential waiting times, while elongation is deterministic between events at the current occupancy-dependent speed. Aggregate mode also uses within-interval promoter and initiation histories, with sequential processing/decay over each molecule's actual remaining age; it still lacks individual RNA positions and identities.

The implemented feedback path is: **marks and nucleosome occupancy → TF capture → bound TF and recruited coactivator counts + spatial enhancer contact → promoter ON/initiation rates → polymerase movement → local nucleosome eviction and gene-wide stiffness changes → altered chromatin contacts**. Mark writing can also change subsequent protein capture and chromatin–chromatin attraction. Important boundaries are:

- Chromatin affinity changes the local capture hazard, not an existing bond's off-rate. Multiple mapped marks use the greatest affinity; overlapping TF-region preferences also use their greatest multiplier. Binding is not a sequence-resolved motif model, and each coarse chromatin bead can bind multiple protein molecules.
- Coactivator recruitment for gene regulation and tethered catalysis requires one actual protein link to a directly DNA-bound TF. Higher-order assemblies exist mechanically, but arbitrary multi-hop recruitment does not contribute to this signal.
- Mark gates are evaluated at the gene's start bead. Losing a mark or contact removes future activation/gain without instantly stopping an ON promoter or engaged polymerase. At least one engaged polymerase, even if paused, applies the stiffness multiplier across the whole gene; eviction acts at its current bead.
- Cohesin changes contacts mechanically; the model has no independent cohesin-to-transcription rate multiplier. The original `speckle_gain` accelerates RNA processing only. Optional `body_responses` independently alter activation and initiation as explicit hypotheses. The imposed nucleolus has no explicit rRNA synthesis pathway.
- There is no free Pol II pool, polymerase traffic, nascent-RNA degradation, MS2 fluorescence observation model, or released RNA–protein binding. Individual RNA is a diffusing point. Keep blocked initiations at zero for quantitative interpretation.

## Manuscript experiments and confinement (v0.9)

Choose **Manuscript hypotheses · 80 kb** for a small nucleosome-unit model, or **Manuscript hypotheses · 20 Mb** for a coarse landscape. Both retain TFs, actual TF–coactivator bonds, nucleosome turnover, individual polymerases and RNA. They add illustrative enhancer, active-body and Polycomb mark patterns, an HDAC-like eraser and PRC1/PRC2-like writers. The species labels, mark choices and rates are hypotheses, not fitted parameters from the manuscript's SIF associations. The 20 Mb model is still 10 kb per bead. The supplied manuscript and the scope of each comparison are recorded in `references/manuscript/manifest.json`; no raw cell data or Supporting Material were provided.

**Confinement is already part of the mechanical model.** System now exposes both the radius `R` and optional stiffness `k` (default 30; zero disables the wall). For each mechanical particle and released RNA point, the soft harmonic wall is

```text
U(r) = (k/2) max(r − R, 0)²
F(r) = −k max(r − R, 0) r̂
```

Small penetrations are possible: positions are not clipped to a hard sphere. The panel reports particle and RNA excursions and wall energy. `R` is the modeled confining domain, not necessarily the whole nucleus around an isolated 80 kb or 20 Mb locus. Radius, stiffness, excluded volume and packing alter contact frequencies; compare wall strengths and domain sizes as controls. This custom harmonic wall is not identical to the built-in polychrom confinement potential.

**Environmental responses** are optional per gene and body:

```json
"body_responses": {
  "speckle": {"activation_multiplier": 1.5, "initiation_multiplier": 1},
  "nucleolus": {"activation_multiplier": 0.7, "initiation_multiplier": 1}
}
```

At promoter-start distance `d` from an enabled body's center and body radius `s`, the weight is `w = exp(−d²/(2s²))`; a configured multiplier `m` contributes `exp(w log(m))`. Effects multiply across bodies. A multiplier of 1 is neutral, values below 1 suppress the rate, and positive activation gains do not bypass the required/repressive-mark gate. These Gaussian prescribed responses are distinct from mechanical body affinity and released-RNA speckle processing. They do not infer causal effects from imaging correlations. The live Regulation view decomposes baseline, recruitment, mark-gate and environment factors and shows nucleosome access and feedback.

Open **Studies** to compare mechanisms without replacing or stepping the live simulation:

- **Four-promoter reference:** two homolog pairs in each simulated cell, exact continuous-time promoter competition/environment switching and nascent-RNA birth/removal. All four promoters share a positive cell-level extrinsic factor. The two cis pairs have independent local environments. Calculate Φcis from A–B and A′–B′, Φtrans from A–B′ and A′–B, and ΔΦ = Φcis − Φtrans. Separate latent nascent-positive activity, detection-thinned activity and promoter ON state. RNA presence is not identical to the promoter state. Detection efficiency is an independent allele-call dropout model, not a calibrated intron-imaging model. The positive lognormal extrinsic distribution and exponential nascent residence are explicit approximations, not an exact reconstruction of the manuscript's missing Supporting Material.
- **Single-chain spatial pilot:** copy the selected model, choose a gene pair and remove selected mechanisms across independent seeded replicas. The nucleosome-access ablation neutralizes both occupied/vacant capture and release contrasts; the nucleosome-turnover ablation disables background, anchored and polymerase turnover plus sliding, retaining the initialized occupancy state. Measure within-chain nascent association, contact/noncontact association and sampled distances. A shuffled independent-replica null is **not Φtrans**: it removes shared-cell effects and is labeled separately. Local pilots are deliberately bounded; full condition configurations can be exported for longer cluster runs.

A sampled minimum distance does not prove absence of historical contact. Nanometer thresholds require independent length calibration and an imaging/localization model. Short trajectories do not establish stationarity. Independent-cell bootstraps apply to the reference; whole-independent-replica resampling applies to the spatial pilot. Correlation uncertainty or bootstrap support is not a causal conclusion. `references/manuscript/verification.json` records actual qualitative controls and code hashes; this is not reproduction of the manuscript's 97% parameter scan or machine-learning results.

Next evidence needed for manuscript-scale inference: paired-homolog spatial simulations, experimental allele-level nascent/SIF/position tables and Supporting Material, observation/noise and cell-state models, condition-specific physical calibration, converged replica and timestep studies, and quantitative target-GPU profiling. Existing ideal Rouse and loop-control benchmarks remain a simpler mechanical baseline. Compare contact-mediated activation, shared environmental modulation and competition separately before combining them.

## Nucleosome occupancy, TF access and mechanics

In **Chromatin → Nucleosome occupancy**, select **Enable nucleosome occupancy**. **Initial occupied fraction** is a probability independently sampled for each bead, not a promise of an exactly matching count. Values 0 and 1 initialize uniformly without sampling. Inclusive initial intervals override that sampled state in list order; **Initially occupied** selects their state. The setting is optional and disabled in earlier models.

Occupancy is a binary state on the existing chromatin bead. Eviction does **not** delete DNA, change the genomic axis, alter bead size or remove backbone links. At 200 bp per bead it is a coarse nucleosome-unit occupancy hypothesis. At 10 kb per bead it is an effective occupied/vacant **bin** variable covering many possible nucleosomes. A 20 Mb run at this coarse resolution cannot report individual nucleosome eviction events.

In **Proteins → Nucleosome access & release**, set **Occupied-site capture ×** and **Vacant-site capture ×** independently for each species. **Occupied-site release ×** and **Vacant-site release ×** independently modify dissociation. Both release factors default to 1, preserving existing models; zero prevents release by that chemical hazard while the corresponding state persists. Protein–protein dissociation is unaffected. The eligible protein–chromatin encounter hazard is:

```text
association hazard = protein on_rate × mark affinity × TF-region multiplier
                     × species occupied/vacant capture factor
release hazard = protein off_rate × species current occupied/vacant release factor
```

A positive occupied-site factor permits association with occupied chromatin; zero blocks new association there. For example, an ordinary-access hypothesis could use occupied/vacant factors `0.05 / 1`, and a more permissive, pioneer-like hypothesis `0.5 / 1`. These examples are arbitrary controls, not measured constants or a universal classification. TF–nucleosome interactions vary by factor and binding arrangement ([Zhu et al., 2018](https://www.nature.com/articles/s41586-018-0549-5)). A factor multiplying an already zero mark/site affinity still yields zero.

Eviction clears **all** histone-mark flags at the affected bead. Vacant beads cannot receive writer-deposited marks or serve as marked donors/recipients in spreading. Reassembled occupancy starts unmarked; it does not restore the previous marks. Writers can subsequently deposit new marks. Because mark loss and accessibility factors act together, eviction need not increase every protein's total association hazard: losing an attractive mark may offset a larger vacant-site factor. Existing protein–chromatin bonds use the occupied/vacant release factor for their current anchor state at subsequent reaction updates; eviction does not automatically unbind them. Setting both release factors to 1 retains the base off-rate regardless of occupancy. Effective **unmodified** chromatin-attraction rules still mean no mark flags and therefore also include vacant DNA.

Define an initially vacant interval in **Chromatin → Nucleosome occupancy** to create a nucleosome-free starting region. It is an initial condition, not a permanently excluded zone: reassembly or sliding can fill it later. A gap affects each species through that species’ capture and release factors, while mark loss may alter its separate mark preference.

**Background eviction rate** and **Background reassembly rate** act at every bead. Protein remodeling rates act only at actual direct chromatin attachments or a recruited coactivator's TF anchors. A protein molecule contributes once per bead across duplicate IDR links, overlapping elements and its own direct attachment. The optional **Polymerase eviction rate** in **Genes** adds a hazard at each engaged polymerase's current bead; paused polymerases remain engaged. In inverse reduced time:

```text
E_i = background eviction + sum of anchored-protein eviction rates
      + sum of engaged-polymerase eviction rates at bead i
R_i = background reassembly + sum of anchored-protein reassembly rates
P(occupied → vacant) = E_i/(E_i+R_i) × [1 - exp(-(E_i+R_i) dt)]
P(vacant → occupied) = R_i/(E_i+R_i) × [1 - exp(-(E_i+R_i) dt)]
```

Zero total rate leaves occupancy unchanged. Reassembly draws from an implicit reservoir; there are no explicit free histones, ATP reactions or histone recycling trajectories. Reactions are ordered as protein binding/release, occupancy turnover, local remodeler sliding, mark chemistry, cohesin, then gene/RNA updates. Within the occupancy substep, rates are held fixed. For initially occupied sites, `P(any eviction) = 1 - exp(-E_i dt)` is sampled jointly with the final state, so eviction followed by reassembly in the same interval still erases old marks. `endpoint_evictions` and `endpoint_reassemblies` count changes between interval endpoints. `mark_reset_intervals` counts initially occupied bead intervals containing an eviction, even if the site ends occupied; none of these counters counts every microscopic turnover event.

### Chromatin remodelers

**Proteins → Add chromatin remodeler** creates an editable species and enables occupancy in the draft. The species still needs to reach and bind chromatin, or be recruited by an actual protein bond to a TF in a configured regulatory element. Its biological name does not silently supply interactions, measured rates or ATP dependence. Applying the draft starts a fresh simulation.

Each species can specify **Nucleosome eviction rate**, **Nucleosome reassembly rate** and **Nucleosome sliding rate**. These activities act at actual direct anchors or configured TF-recruited anchors, with one contribution per unique molecule and anchor despite duplicate IDR links or overlapping elements. Eviction removes occupancy and all marks; reassembly restores occupancy without old marks. Sliding moves an occupied anchor or neighboring nucleosome into an adjacent vacant bead. The configured sliding hazard is shared uniformly over valid local moves for that molecule and anchor; no available vacancy means no move. Sliding carries the nucleosome’s current mark flags to the destination, clears the source, and conserves total occupied sites. DNA beads, genes and existing protein tethers stay at their genomic coordinates.

These are coarse kinetic rules. An adjacent 200 bp bead move is not a measured base-pair remodeling step. There is no explicit ATP hydrolysis, wrapped DNA, histone identity, sliding directionality, remodeling motor cycle or steric exclusion between a moved nucleosome and a bound TF. Compare remodeling-disabled controls, equal occupied/vacant capture and release factors, and different reaction timesteps before interpreting a coupled effect. A nucleosome-free region need not increase every species’ occupancy or transcription.

Occupancy also multiplies the local persistence value:

```text
local persistence_i = baseline persistence
                      × product of persistence multipliers for overlapping engaged genes
                      × occupied/vacant persistence multiplier
```

Set **Occupied persistence ×** and **Vacant persistence ×** to 1 to remove this mechanical effect while retaining accessibility changes. A smaller vacant multiplier makes the local bending potential softer; it is an independent modeling hypothesis, not a deduction from improved TF access. Occupancy changes the force and energy on the next mechanical evaluation. Transcription's multiplier applies across an entire gene while at least one polymerase remains engaged, including during pauses.

## Individual elongation and RNA molecules

In **Genes → Individual RNA molecules**, enable **Track individual transcripts**. This replaces the legacy aggregate RNA update with individually identified polymerases and transcripts. The `nucleosome_transcription` and `transcription_20mb` presets enable it together with occupancy and a recruited remodeling species. Existing count-mode models remain available when the switch is off.

Every accepted initiation creates a stable `pol-N` identity and matching `rna-N` identity within the run. Nascent RNA follows its polymerase's interpolated chromatin position. All genes currently transcribe in the increasing bead-index direction, from `start` to `end + 1`. With genomic mapping `b` base pairs per bead:

```text
nascent RNA length = (polymerase bead-position - gene start) × b nucleotides
completed RNA length = (gene end - gene start + 1) × b nucleotides
elongation speed = configured beads/time × occupied elongation multiplier, if occupied
elongation speed = configured beads/time, if vacant or occupancy is disabled
```

In **Elongation & nucleosome coupling**, **Occupied elongation ×** controls speed at each occupied bead. The engine integrates movement across each bead boundary at its corresponding speed. A value of 0 stalls movement at occupied beads until occupancy changes. **Polymerase pause rate** and **Polymerase resume rate** govern additional stochastic pauses; paused polymerases have zero elongation speed. These speed and pause controls apply to individual-transcript mode. The polymerase-eviction coupling also works with the positions retained by the legacy model.

Nucleosomes can alter Pol II pausing and elongation, as shown by single-molecule measurements ([Hodges et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2775800/)). The implemented multipliers and pause process are configurable approximations; they do not resolve DNA unwrapping, histone transfer, backtracking or a measured nucleosomal energy landscape. There is no finite polymerase pool, polymerase exclusion/traffic, strand selection, nucleotide sequence, intron structure or alternative splicing. The displayed nucleotide length is the full modeled genomic interval, not an annotated mature transcript length.

On completion, RNA is released at the final gene bead as an **unprocessed passive point particle**. **Released RNA diffusion D** and **Processed RNA diffusion D** specify the two stages' diffusivities in reduced length²/time. The particle undergoes Brownian motion within the same soft nuclear confinement; it is not a resolved RNA chain and has no RNA–protein binding, excluded-volume feedback, RNA–RNA attraction or force on chromatin. The 3D **Transcripts** toggle displays actual simulated particle positions. The Regulation view's **Elongation & individual RNA** panel and 1D **Polymerases** track expose identities, progress, pause state and transcript length.

In individual mode, processing depends on the **released RNA particle's** position. Inside any configured speckle, `processing rate = base processing rate × (1 + speckle_gain)`; overlapping speckles do not multiply the gain further. Outside speckles it uses the base rate. The legacy count mode instead uses gene proximity as its processing proxy. Neither mode changes promoter activation directly through this gain. There is no nucleolus-specific RNA processing or rRNA biogenesis.

The configured `rna_decay_rate` acts on both released unprocessed and processed RNA, but not nascent RNA. Processing and degradation compete, and a newly processed transcript accrues processed-state diffusion and decay only over the remaining interval. Local processing hazards are frozen within each spatial substep, so this remains an operator-split approximation of diffusion across a speckle boundary. Accepted initiations occur at their sampled times within actual promoter ON intervals; new polymerases accrue only their post-birth elongation time. The following accounting identities hold within a run:

```text
initiated = nascent + unprocessed + mature + degraded
completed = unprocessed + mature + degraded
```

**Maximum tracked transcripts** limits live nascent and released RNA together: default **2,000**, allowed up to **5,000**, separately from the 6,000 chromatin/protein-particle limit. Existing transcripts are never removed to make room. When capacity is exhausted, a lottery allocates available slots among initiation attempts and increments **blocked initiations** for rejected attempts. This changes expression kinetics: blocked initiations must remain **zero** for quantitative interpretation. Increase capacity within the supported limit or reduce the simulated scope, then rerun. A full cap with no RNA degradation can block initiation indefinitely.

Snapshots include `nucleosomes` with occupancy, endpoint/reset counters and accepted `slides`, and `transcription` with current transcript records, identities, stages, positions, lengths, times, per-gene counts and cap diagnostics. These are current-state exports, not a saved full trajectory of every RNA particle.

## Protein preferences and chromatin attraction

The workbench exposes two distinct mechanisms. **Proteins → Chromatin encounter affinity** changes how quickly an eligible protein associates with a nearby chromatin bead. **Chromatin → State-selective chromatin attraction** directly adds a mechanical attraction between selected chromatin states, independently of protein occupancy. The Regulation view displays both the protein preference table and the live chromatin-contact table.

For a chromatin rule, select two marks (or **unmodified**), a well depth `epsilon` in kBT, an attraction radius `R` in reduced length, and a minimum contour separation in beads. The symmetric pair potential is:

```text
U(r) = -epsilon × (1 - r²/R²)²   for r < R
U(r) = 0                       for r ≥ R
```

The default minimum separation is **3**, excluding pairs fewer than three indices apart along the chain. **1 is allowed** and includes backbone neighbors; doing so changes their local mechanics as well as long-range compaction. Zero strength disables a rule. When multiple current marks make several rules eligible for a physical pair, only the **most negative current well energy** is used. Wells are not added, and the winner need not have the largest configured depth because their radii can differ. Each individual well has zero force at its cutoff; the combined minimum-energy potential can have a force kink where differently shaped wells cross. Equal-energy ties select the first rule.

Here **unmodified means no configured mark flag is present on the bead**. This differs from protein affinity's unmodified *fallback*, which applies when no present mark has an entry in that protein's mapping. Live contact counts count each physical pair once under its winning rule; the reported negative energy is the attraction contribution in kBT, not the system's total potential energy or a binding Kd. Marks are read afresh when computing forces and reports.

Mark writing can therefore alter subsequent folding, which can alter encounters and further writing. This is active phenomenological feedback, not an equilibrium mark model. State-dependent attractive polymers provide a modeling precedent ([Michieletto et al., 2016](https://doi.org/10.1103/PhysRevX.6.041047)); they do not imply that a histone mark itself exerts a universal pair force. Explicit reader-mediated bridges and an effective state attraction may count overlapping mechanisms. Compare zero-attraction and zero-reader-bridging controls separately before combining them. Strong compaction or a connected cluster alone does not establish liquid behavior.

## Protein interactions and DNA attachment

Use **Interactions** beside the 3D and 1D views to explore which species can associate and which connections currently exist. The view uses the applied configuration and live simulation state. A configured interaction is a possible reaction, not evidence that an assembly has formed.

For a selected protein pair, three diagrams separate connections where **neither**, **one**, or **both** partners are directly bound to chromatin. These categories describe direct attachment. A protein with no direct chromatin bond can still be tethered indirectly through a connected protein assembly; the population counts distinguish this from proteins in completely detached components. Sharing a chromatin bead does not itself create a protein–protein edge.

The view distinguishes physical domain links from unique molecular pairs: two IDR molecules connected through several domain contacts remain one molecular pair. Components describe connected molecules in the current frame; their identifiers are not persistent assembly lineages. The same configured protein–protein reaction rule currently applies in all three attachment contexts. Attachment changes mechanics and encounter opportunities; it does not silently introduce a different chemical affinity.

### Estimating physical binding constants

The **physical binding estimate** requires a length mapping (nanometers per reduced length), a time mapping (seconds per reduced time), and a target dissociation constant in micromolar. These are explicit modeling assumptions, not inferred from base pairs per bead. The estimate does not edit the simulation's reaction rates.

For isolated spherical partners in dilute three-dimensional solution, the Smoluchowski pair-encounter kernel is

```text
D_physical = D_model × length_um² / time_s
R_physical = R_model × length_um
k_pair,limit = 4π (D_A + D_B) R × N_A × 10⁻¹⁵     [M⁻¹ s⁻¹]
k_pair = capture_efficiency × k_pair,limit
```

Here diffusivities are in µm²/s, the encounter radius is in µm, and `10⁻¹⁵` converts a cubic micrometer to liters. Capture efficiency **1** assumes every encounter reacts. Smaller values provide a phenomenological correction for orientation or other reaction barriers, not a calibrated model of those barriers.

For distinct partners `A + B ⇌ AB`, `Kd = [A][B]/[AB]`, so the nominal dissociation rate is `koff = Kd × k_pair`. Identical partners require an explicit counting convention: for `2 A ⇌ A₂`, association-event production is `½ k_pair [A]²`, and with `Kd = [A]²/[A₂]`, the nominal `koff = Kd × ½ k_pair`. Both the pair kernel and the macroscopic association coefficient are reported. Concentrations refer to free monomers and complexes, not total protein.

**These are target-rate estimates, not measured Kd values of this simulation.** The engine's existing association parameter has units of inverse reduced time *while an eligible pair is within its capture radius*. It is not a bimolecular molar association constant. Simply substituting the calculated rate into that field would count diffusion incorrectly. Immediate rebinding after release, timestep, capture geometry, valence, bound-state mechanics and competing reactions affect the apparent equilibrium and dwell times.

The protein–chromatin estimate treats one capture region as an isolated spherical target and uses the configured bare chromatin-bead mobility. A connected polymer, a 10 kb coarse bead, and a biological DNA-binding site are different targets; mark multipliers are not converted to experimental Kd values. Likewise, the free-molecule diffusion coefficient is not the measured diffusion of a DNA-tethered assembly. The current estimate does not infer three-dimensional, sliding or confined association rates for bound partners. Multivalent complexes can exhibit avidity and apparent affinities that cannot be represented by a single pair Kd.

To make physical affinities control validated simulations, first calibrate free-particle diffusion and length/time mappings, then fit the finite-radius reaction and release scheme using dilute two-body association, dissociation and equilibrium measurements. Check timestep/capture-radius convergence and rebinding, before adding DNA attachment, multiple binding sites, marks or condensates. Preserve intrinsic pair constants separately from effective assembly affinities.

The distinction between finite-radius hazards and macroscopic association follows [Donev, Yang & Kim (2018)](https://doi.org/10.1063/1.5009464) and [Vijaykumar, Bolhuis & ten Wolde (2016)](https://doi.org/10.1039/C6FD00104A). For the consequences of rebinding beyond a simple rate substitution, see [Lawley & Keener (2016)](https://doi.org/10.1016/j.bpj.2016.10.008).

## Choosing genomic coverage and nucleosome-scale resolution

Choose **20 Mb chromatin** in the Experiment menu. This preset contains **2,000 chromatin beads at 10,000 bp per bead**, plus 314 protein beads representing 194 molecules across six species. It includes six illustrative gene regions, histone-mark domains, directional barriers, 24 cohesins, a speckle field and a nucleolar field. The full configuration is in `examples/twenty_mb.json`.

**20 Mb regulatory landscape** adds enhancer/promoter elements, activator–p300 recruitment and a weak H3K9me3-selective chromatin attraction to that coarse landscape. **80 kb nucleosome-unit locus** uses **400 beads × 200 bp**, two genes with enhancer/promoter pairs, six cohesins and the same illustrative recruitment/attraction mechanisms. At 10 kb per bead, a regulatory element represents a coarse interval; promoter and motif structure within that interval is unresolved.

The workbench displays the genomic span of both the applied model and your draft. The Chromatin editor computes:

```text
genomic span = chromatin bead count × base pairs per bead
20,000,000 bp = 2,000 × 10,000 bp
```

The limits are **4,000 chromatin beads and 6,000 total particles**, allowing a 20 Mb model at 5 kb/bead as well. Gene annotations and cohesins each support up to 256 entries; regulatory elements support 128; initial mark regions and barriers each support 512. A 1 kb resolution model would need 20,000 chromatin beads and is beyond this version's supported limit. At **200 bp per bead, the current maximum is 0.8 Mb**; covering **20 Mb requires 100,000 chromatin beads before adding proteins**, and is not supported by this workbench.

Changing resolution manually retains annotation **bead indices**. Update gene and regulatory-element intervals, initial mark regions and barriers if you want to preserve their genomic coordinates. Internal intervals use zero-based inclusive bead indices: bead `i` covers the half-open bp interval `[i × b, (i + 1) × b)`, with `b` bp per bead. A gene from bead `start` through `end` therefore covers `[start × b, (end + 1) × b)`; the UI labels its last included base explicitly. Gene elongation and cohesin stepping are configured in beads per reduced time, so changing bp/bead changes the displayed genomic speed unless you adjust those rates.

**Genomic coverage is specified; physical scale is not yet calibrated.** Assigning 10 kb to a bead does not determine its nanometer size, diffusivity, persistence length or seconds per timestep. Histone marks in this setup describe coarse genomic regions, not individual nucleosomes. Resolution comparisons require mechanical and kinetic calibration rather than changing bp labels alone.

The **nucleosome-unit** setting labels an isotropic bead at a chosen repeat length; it does not add a histone core, wrapped DNA, orientation, tails, variable linker geometry or unwrapping. The separate, optional occupancy model adds coarse eviction/reassembly and remodeling rates on that bead, without resolving these structures. A nucleosome core wraps roughly 146–147 bp, with linker DNA completing a variable repeat; 200 bp is an approximation ([Luger et al., 1997](https://doi.org/10.1038/38444)). A binary mark on such a unit is still a simplified state, not separate modifications on both copies of each histone. Explicit geometry costs more particles: a published minimal model represents 300 kb/1,500 nucleosomes with approximately 50,000 particles, retaining 5 bp DNA beads and linker mechanics ([Farr et al., 2021](https://doi.org/10.1038/s41467-021-23090-3)). Another model sampled 6,000 oriented nucleosomes plus linker DNA over more than 1 Mb using Monte Carlo, demonstrating a structural modeling scale rather than real-time Brownian throughput ([Zülske et al., 2024](https://doi.org/10.1016/j.bpj.2024.02.024)).

A future 100,000-unit model needs a benchmarked backend with native GPU forces, scalable reaction updates, streamed trajectories and reduced-detail browser rendering. The current polychrom adapter still evaluates forces and chemistry on CPU and exchanges arrays each integration substep; selecting it does not supply that architecture. A practical present workflow is to explore a coarse 20 Mb context and a separate 80 kb locus, with calibrated comparisons rather than assuming their physical scales are interchangeable.

Spatial neighbor queries replace stored all-pairs arrays for repulsion, capture searches and mark spreading. The polychrom adapter uploads force arrays in bulk. A development baseline of the 2,314-particle model took **2.93 seconds for 500 polychrom/OpenMM CPU steps** after warm-up (about 5.9 ms/step), saved in `benchmark-20Mb-polychrom.json`. That advances 1.0 reduced time unit, not one biological second.

The bounded NumPy measurements below are saved with their exact configurations, software versions and engine-source hash in `benchmark-nucleosome-baseline.json`. They were taken **before regulatory recruitment and state-selective attraction were added**, using existing models relabeled at 200 bp/bead; they are not timings of the new presets.

| Baseline case | Total particles | Initialization | First 100 steps | Next 100 steps | Snapshot + JSON, median |
|---|---:|---:|---:|---:|---:|
| 400 chromatin beads, 80 kb | 492 | 0.022 s | 0.169 s | 0.198 s | 1.36 ms |
| 4,000 chromatin beads, 0.8 Mb | 4,314 | 0.173 s | 1.242 s | 1.351 s | 8.83 ms |

Both runs completed 200 steps (0.4 reduced time) with finite coordinates and potential energy. They differ in protein counts and confinement, so their ratio is not a controlled scaling law. The measurements exclude browser/network cost and establish short-run execution, not equilibration or biological validity. **No 100,000-bead run or throughput extrapolation is claimed.**

The regulatory presets were also tested for **100 steps each on both NumPy and OpenMM Reference**, including recruitment logic and state-dependent attraction, before adding nucleosome occupancy and individual RNA. The 80 kb preset took about 2.3–2.5 ms per step; the 20 Mb preset about 14–16 ms per step. All four short runs retained finite coordinates, forces and energies. See [the regulatory benchmark](benchmark-regulation.json) for environment, configurations and measurements.

The newer [nucleosome and RNA benchmark](benchmark-nucleosome-rna.json) exercises the 80 kb preset for **4,000 NumPy steps**: 20 transcripts initiated, seven completed and one degraded, with zero blocked initiations and conserved RNA accounting. It also checks the 20 Mb preset for 100 steps on each backend. Coordinates, forces and energies remained finite. These local execution checks exclude browser/API overhead and do not establish long-run equilibration, physical time calibration or biological validity.

## What this version represents

| Component | Implemented representation | Interpretation and boundary |
|---|---|---|
| Chromatin | Bead–spring chain, bending stiffness, soft repulsion, optional state-selective attraction and confinement | Coarse mechanics; nucleosome-unit labels do not add wrapped-DNA geometry. Soft interactions can permit chain crossing; no topological conservation guarantee. |
| Protein species | Independent counts, free molecular diffusion, radii, binding valence and affinities | Up to 32 species and 6,000 chromatin/protein particles in the local workbench; individual RNA has a separate cap. |
| Higher-order association | Reversible molecule–molecule connections, bounded by valence | Assemblies can exceed dimers. Cluster size counts connected molecules, not proof of liquid–liquid phase separation. |
| Explicit IDRs | Flexible domain-bead chains with selectable sticker fraction | Spatial chain conformations and stickers are explicit. These are not residue-level sequences, electrostatics or CALVADOS parameters. |
| Chromatin binding | Reversible mechanical tethers with reciprocal forces | Reduced mobility arises from attachment and collective motion. Free diffusivity and measured effective mobility are distinct. |
| Histone modifications | Independent binary flags on chromatin beads; direct and TF-recruited writer/eraser reactions, turnover and optional contact spreading | Marks alter protein association and configured chromatin attraction. Enabled occupancy clears marks on eviction and prevents writing/spreading onto vacant beads. No automatic mutual exclusivity or explicit histone-copy states. |
| Nucleosome occupancy | Optional binary occupancy, background/anchored eviction and reassembly, mark-carrying adjacent-site sliding, polymerase-associated eviction, species-specific capture/release and local persistence factors | A nucleosome-unit or coarse-bin hypothesis with an implicit free-histone reservoir. No wrapped-DNA geometry, tracked histone identities, histone recycling or ATP-dependent motor cycle. |
| Regulatory elements | Enhancer/promoter intervals, TF site multipliers, explicit one-hop coactivator recruitment and contact-gated rate gains | Up to 128 elements. No DNA sequence scanning, resolved catalytic reach or universal enhancer mechanism is claimed. |
| Cohesin | Stochastic unloading/repositioning and outward-moving harmonic links, barriers and motor collisions | Initial placement honors the requested count. Unbound motors reload when an adjacent pair of sites is free; there is no independent loading-rate parameter yet. |
| Genes | OFF/ON promoters, initiation and contour-based elongation; optional individual polymerases with stochastic pausing/resumption and occupancy-dependent speed | Promoter and RNA kinetics are distinct. No finite free-polymerase pool, polymerase exclusion/traffic, backtracking, introns or sequence-dependent rates. |
| RNA | Legacy counts, or optional growing nascent records and passive diffusing released/processed points with stable IDs | Individual RNA positions govern speckle processing in this mode. No RNA chain mechanics, RNA–protein binding or RNA-driven forces. |
| Gene mechanics | Local persistence multiplier across a gene while any polymerase remains engaged, including paused polymerases | A configurable phenomenological coupling that multiplies the optional occupancy factor; each coupling can be independently disabled. |
| Speckle and nucleolus | Separately configurable, stationary spatial affinity fields | These imposed landmarks are not self-assembled organelles. Explicit protein networks evolve separately. Nucleolar subphases, rDNA and rRNA biogenesis are not implemented. |

## Dynamics and reaction semantics

The reference integrator uses overdamped Euler–Maruyama motion at kBT = 1:

```text
x(t+dt) = x(t) + D F dt + sqrt(2 D dt) Normal(0,1)
```

Each protein's configured diffusion is its dilute **molecular center-of-mass diffusivity**. An IDR represented by n equal-drag beads uses n times that diffusivity per bead, so its isolated center of mass has the requested diffusivity. Binding, internal constraints, crowding and confinement then change measured trajectories. Independent drag omits hydrodynamic coupling.

Protein chromatin-affinity values are dimensionless **encounter-rate multipliers**, not binding energies or experimental dissociation constants. When a bead has mapped marks, the maximum applicable marked affinity is used; the unmodified affinity is used when no mapped mark is present. Thus a single mapped repressive mark can suppress binding, but conflicting multiple marks follow the documented maximum rule. Enabled TF-region preferences multiply this value, with the largest multiplier used for overlapping regions of the same TF. Enabled nucleosome occupancy supplies a further species-specific occupied/vacant multiplier. Chromatin–chromatin state attraction instead contributes energy and reciprocal mechanical forces as described above.

Protein association respects a molecule's total valence across protein and chromatin connections. Finite-distance candidate encounters and stochastic reaction hazards govern association. Reaction parameters plus mechanical potentials are phenomenological; they have not been jointly fitted to detailed balance or macroscopic on/off rates. Energy-consuming writing, transcription and extrusion intentionally drive the model.

Contact spreading is a separate, optional effective rule representing local mark propagation. Set `spread_rate` to zero to study only explicit writer-mediated deposition and erasure. RNA processing can depend on proximity to an imposed speckle through `speckle_gain`: individual mode uses each released RNA's position, while legacy count mode uses gene proximity. This is distinct from promoter activation and can be disabled with zero gain.

Required and repressive marks gate promoter activation at the gene's start bead. **For a gene with no enabled regulatory elements**, distinct activator molecules bound anywhere in that gene multiply its ON rate by `1 + bound activator count`. Enabled elements targeting the gene replace this legacy gain with the regional activation/initiation rules above. While at least one polymerase remains engaged, including a paused polymerase, the persistence multiplier applies across that gene's entire region. Multipliers combine for overlapping engaged genes. These explicit choices allow controlled hypotheses; they are not inferred from the illustrative species names.

### Two available integration backends

- **NumPy reference:** computes forces, advances Brownian motion and executes the reactions locally.
- **Polychrom/OpenMM:** uses a real polychrom simulation and OpenMM custom Brownian integrator with per-particle mobility. The same CPU model computes forces and reactions; force arrays upload in bulk at every integration substep. It is an integration adapter, **not yet an optimized native OpenMM force implementation**. It should not be expected to provide GPU-scale acceleration.

Both engines use the same configuration and chemistry. Their random streams differ, so compare distributions and mechanisms, not frame-for-frame identity. Requesting an unavailable backend returns an error; the app does not silently substitute a different engine.

## GPU and cluster execution

**CUDA execution is supported for integration, but this version is not a GPU-accelerated implementation of the complete model.** The browser's Polychrom option retains its existing CPU default. Cluster jobs explicitly select NumPy, OpenMM CPU, Reference, CUDA, or OpenCL. CUDA runs the Brownian integrator on an NVIDIA GPU; OpenCL runs on the selected OpenCL device, which can also be a CPU. Force calculation, spatial reaction searches, mark chemistry, cohesin and individual RNA updates remain on the host. Every integration substep transfers coordinates, forces and mobilities. These transfers can outweigh GPU integration gains. Benchmark whole-run throughput before purchasing or reserving substantial GPU time.

### Prepare and transfer a job

1. Build the model in the UI, then open **Cluster job**. The exported model is the current editor draft, including unapplied changes. It starts from the configuration and seed; it does not continue the displayed trajectory.
2. Choose the target platform, precision, independent replica count, integration steps and output intervals. Leave the device index blank to use the allocated device. A populated index is a **logical index within the visible set**, not a physical node GPU number. Precision controls apply to CUDA/OpenCL; CPU and Reference use their platform's fixed precision.
3. Enter CPU cores, memory, walltime and optional Slurm partition/account. These are scheduler requests, not measured requirements. The template requests one GPU per CUDA/OpenCL replica and one node per task. Replicas have distinct seeds `model.seed + replica_index`; they do not distribute one polymer across GPUs or nodes.
4. Preview the model size, reduced simulation duration, replica seeds and scripts, then download the ZIP. Transfer it with your cluster's approved file-transfer tool. The app does not connect to the cluster, store credentials or submit jobs.
5. Extract the ZIP, prepare a Python environment following its `RUN.md`, and use an allocated compute node for the preflight. Inspect the Slurm directives. For site-specific changes, copy `submit.slurm` to `submit-local.slurm` and edit the copy, keeping the hashed originals intact. From the extracted directory, submit with `sbatch submit-local.slurm`, or run `sh run.sh 0` under another scheduler's allocation for replica 0. Activate the prepared environment first.

OpenMM's [platform documentation](https://docs.openmm.org/latest/userguide/library/04_platform_specifics.html) describes device and precision settings. Slurm [GPU management](https://slurm.schedmd.com/gres.html#GPU_Management) explains scheduler-controlled visibility, and [job arrays](https://slurm.schedmd.com/job_array.html) provide the replica index. Do not overwrite `CUDA_VISIBLE_DEVICES` in the generated scripts. Registered OpenMM plugins are not proof that a usable device exists; actual Context initialization is required. A requested device that fails initialization stops the job, without CPU fallback. No CUDA hardware was available for local development acceptance; run the bundled checks on your allocated GPU before production.

The ZIP contains the validated model, job settings, headless source, dependency requirements, scripts and SHA-256 manifest. Optional dependencies pin OpenMM and a polychrom commit; NumPy/SciPy minimum requirements are not an exact environment lock. Each run records actual package versions, source/model hashes, requested/actual backend properties, scheduler identifiers, seed, timing and completion status. Bundle contents are checked before validation and batch execution. This is provenance and accidental-change detection, not a signed software distribution.

Results stream to per-replica `metrics.jsonl` and optional `snapshots.jsonl`, plus an atomic `final-state.json` and `run-manifest.json`. Setting snapshot interval to zero disables periodic snapshots. Initial and final summary metrics and a final full state are retained even when periodic snapshots are disabled; output volume still grows with particle count and duration. Reusing an existing result directory is refused. Numerical failure records a failed manifest; abrupt scheduler termination can leave a `running` manifest and partial final JSONL line. There is no exact restart checkpoint or automatic walltime continuation yet, so pilot durations should fit within the requested walltime. Download the results from the cluster with the same approved transfer workflow.

GPU export preserves the current **4,000 chromatin bead / 6,000 mechanical particle / 5,000 live RNA** limits. A 20 Mb model at 200 bp per nucleosome would require 100,000 chromatin beads before adding proteins and remains unsupported. The 20 Mb presets run at coarse genomic resolution. A validated native force backend and new scaling measurements are prerequisites for raising these limits.

### Scientific acceptance and experimental controls

Every generated `run.sh` first runs the analytical acceptance suite on its allocated node, then runs the experiment only if it passes. Run it independently with:

```sh
python -m chromatin_studio.validation --platform numpy --output validation-cpu.json
python -m chromatin_studio.validation --platform CUDA --precision mixed --output validation-cuda.json
```

Reports give observed values, analytic targets, numerical/statistical tolerances, source hashes and actual platform details. A requested GPU that is unavailable is a **failure**, rather than a skipped success. Stochastic acceptance uses six standard errors from independent-particle or binomial ensembles; these deliberately broad regression thresholds detect gross implementation errors, not small biological effects. OpenMM checks are compared with a NumPy ensemble; identical trajectories across devices are not expected. The harmonic test checks two timesteps against exact *discrete Euler* moments and reports the continuous-time targets, separating known discretization bias from coding errors.

| Test | Expected result | What it can reveal |
|---|---|---|
| Free diffusion and frozen particles | MSD = 6Dt; D = 0 gives no displacement; changed D applies immediately | Noise amplitude, thermal units, stale mobility uploads |
| Constant force followed by reversal | Mean displacement = D F dt at kBT = 1; thermal variance = 2Ddt per coordinate | Force sign, drift/noise scaling, stale force uploads |
| Harmonic relaxation at dt and dt/2 | Ensemble mean/variance match the known discrete Ornstein–Uhlenbeck solution | Integrator mistakes and explicit timestep bias |
| Chromatin force/energy consistency | F = −∇U and internal forces sum to zero | Bending/repulsion sign or gradient errors with heterogeneous occupancy |
| Isolated TF capture and dwell | Capture = 1 − exp(−k×access×t); survival = exp(−koff×t) without rebinding | Incorrect accessibility multipliers, dissociation hazards |
| Nucleosome turnover | Exact two-state endpoint probabilities; any eviction erases old marks | Turnover bias or marks surviving hidden eviction/reassembly |
| Individual RNA | Length = speed×time; release at full gene length; analytical processing/decay fractions; conserved molecule accounting | Off-by-one genes, incorrect growth, lost/double-counted molecules |
| Short coupled model on selected device | Finite state, valid protein valence and RNA accounting | Integration/reaction compatibility, not equilibration |

The existing unit/integration suite additionally covers coactivator recruitment, histone spreading, cohesin barriers, TF/mark preferences, RNA spatial processing, pausing, capacity and UI/model consistency. The batch runner checks finite positions/forces/energies, valence, RNA identities and accounting at output boundaries, and records blocked RNA initiation as a quantitative-interpretation warning.

Before interpreting any locus experiment, use the following protocol. These full-model convergence and biological calibration studies are **planned experiments**, not completed acceptance checks:

1. **Timestep and output convergence:** compare dt, dt/2 and dt/4 at the same total reduced duration and sampling times across multiple seeds. Compare Rg, contact maps, protein occupancy/dwell distributions, cluster-size distributions, mark fractions and RNA production with uncertainty intervals. Brownian displacements must resolve capture radii and narrow fields; bounded forces or adaptive mechanical substeps do not by themselves guarantee converged reaction encounters. Thin saved frames separately to check observables are not an output-frequency artifact.
2. **Binding and assembly controls:** measure free/attached/complex center-of-mass MSD versus lag; tethered motion should cross over with tether stiffness and chromatin mobility. Sweep concentration and valence with association disabled as a negative control. Calibrate effective Kd from dilute occupancy/dwell/rebinding experiments before claiming that a diffusion-limit estimator produces the requested macroscopic affinity. For condensates, require concentration coexistence, size dependence and exchange/FRAP behavior; a visible cluster alone is insufficient.
3. **Regulatory interventions:** disable one coupling at a time—TF nucleosome access contrast, remodeler activity, coactivator recruitment/catalysis, histone spreading, chromatin attraction, cohesin, pausing, elongation slowdown and transcription-dependent stiffness. Measure expected downstream effects and confirm unrelated mechanisms stay unchanged. Allow mark loss to oppose improved access. Compare body affinity and speckle processing gain separately; the original speckle gain affects processing; optional body responses change promoter rates independently.
4. **Reproducibility and statistics:** separate equilibration from measurement, estimate autocorrelation, and use independent replicas rather than treating every frame as independent. Start with a small pilot array, then choose the number of replicas from observed uncertainty. GPU versus CPU comparisons should include the same observables and uncertainty, beyond these simple integration fixtures. Require zero blocked RNA initiations for quantitative expression results.
5. **Experimental calibration:** map reduced length/time using chromatin distances and free-protein diffusion, fit binding and remodeling rates to suitable measurements, and compare contact/RNA response distributions. Fit a subset of loci, then predict held-out loci and perturbations. This is required before claiming biological accuracy.

### Route to useful GPU throughput

First profile representative 80 kb and 20 Mb runs on the target node, separating initialization, host forces, reactions/RNA, transfers, integration and output; report walltime per simulated reduced time and scaling with particles/protein counts. Then port backbone, bending, confinement and short-range forces to native OpenMM kernels, followed by dynamic tethers, cohesin and mark-dependent attraction. Validate each force against the CPU energy gradients and frozen-state force arrays before testing ensemble dynamics. Preserve reaction timing during this port; batching more mechanical steps between chemistry updates changes the model and needs its own convergence study. Reduce host/device transfers and optimize reaction neighbor searches only after profiling identifies the dominant cost. Exact stochastic checkpoints and GPU benchmarks at increasing system sizes come before supporting a full 20 Mb nucleosome-resolution model.

## Validation and saved data

Run the checks from this directory with the Python environment used to start the app:

```sh
python -m unittest discover -s tests -v
node --test tests/test_genome_ui.cjs
node --test tests/test_interaction_ui.cjs
node --test tests/test_regulation_ui.cjs
node --test tests/test_transcription_ui.cjs
node --test tests/test_cluster_ui.cjs
```

The Python tests cover configuration integrity and resource bounds, mechanical force gradients and state-attraction selection, reaction/occupancy invariants, TF recruitment and catalytic deduplication, regulatory rate/contact semantics, nucleosome access and mark-memory loss, occupancy-dependent mechanics, individual RNA accounting/pausing/capacity/spatial processing, authoritative read-only summaries, physical-estimate units/conventions, and optional OpenMM diffusion/force response. The Node checks cover genomic coordinates, interval boundaries, navigation, interaction catalog/count semantics, regulatory interval/mark/unique-molecule summaries, and transcript identity/state displays. Optional-engine tests skip when dependencies are missing. Numerical correctness does not establish that the illustrative parameters represent a real cell.

The API rejects invalid regions, dangling species/mark references, duplicate pair rules, nonfinite values and configurations above the prototype particle limit. Numerical failures stop stepping rather than silently replacing the trajectory.

Model exports contain all configuration and seed information. Experiment exports contain the model, current spatial snapshot, recent events, gene/regulatory state, nucleosome occupancy, current individual RNA/polymerase records when enabled, chromatin-contact counts and energies, engine/version metadata, and up to **2,000 measurement frames**. The current RNA records retain molecule identities within that run; the export does not contain complete individual RNA trajectories. Exports are not full restart checkpoints and do not preserve every integration step. A new run can be initialized from the exported configuration; continuing the exact saved stochastic state is a later feature.

## Protein and transcription evidence after polymer calibration

A polymer clock supplies **one shared physical time conversion** for the coupled system. It does not determine protein diffusivities, binding lifetimes, catalytic rates or promoter switching rates. To match a measured free-protein diffusivity, set `D_reduced = D_measured × time_s / (length_nm/1000)²`. For an independently measured unimolecular rate, set `rate_reduced = rate_per_second × time_s`. These inputs then run on the same clock as chromatin; assigning a different physical clock to each species would make their coupled interactions inconsistent.

In a stationary, ideal two-state binding model, the bound fraction constrains `k_on_effective / (k_on_effective + k_off)`. A residence-time measurement can supply a separate off-rate, but only after accounting for photobleaching, missed events and the definition of the bound state. Neither an occupancy fraction nor chromatin MSD alone identifies both rates. A diffusion-limited capture estimate is an ideal upper-bound model; reaction probability, binding geometry, rebinding and local chromatin motion still determine the measured association kinetics. The simulator's conditional encounter hazard is not a molar association constant.

The following primary evidence provides distinct calibration or validation targets. These references are an **evidence inventory**, not a claim that the example parameters have been fitted to every dataset.

| Evidence | What it can constrain | What it does not establish |
|---|---|---|
| [Hansen, Woringer et al. 2018, Spot-On](https://elifesciences.org/articles/33125); [open trajectory data](https://zenodo.org/records/834781), including Sox2 in mESCs and CTCF in U2OS cells | Free/bound diffusion and population fractions, with defocalization and acquisition corrections | The displacement-mixture model ignores state transitions; these fits do not independently identify association/dissociation rates or locus-specific affinity. Cell types and proteins must remain matched. |
| [Hansen et al. 2017, CTCF/cohesin dynamics](https://elifesciences.org/articles/25776) | Bound-state survival and residence-time analysis with photobleaching controls; a strategy for constraining off-rates | A genomic average is not the off-rate of every binding site. `k_off = 1/mean lifetime` applies to a specified exponential state, not an arbitrary mixture or uncorrected FRAP recovery. |
| [Hansen & Zechner 2021, promoter context](https://pmc.ncbi.nlm.nih.gov/articles/PMC7888307/); [data and code](https://doi.org/10.5281/zenodo.2755026) | Dynamic TF input and gene-output constraints on promoter switching and memory | Yeast Msn2 reporter experiments do not supply mammalian default rates or direct molecular TF binding kinetics. Context-dependent responses may require a promoter network beyond the current two-state model. |
| [Yang, Pinholt et al. 2026, E–P dynamics and transcription](https://pmc.ncbi.nlm.nih.gov/articles/PMC13308063/); [trajectory data](https://doi.org/10.5281/zenodo.20293143) | Joint enhancer–promoter distance and nascent MS2 output; a test of contact-dependent promoter activation after matching imaging and RNA observation models | **Preprint, not peer reviewed** (version 1, June 21, 2026). One engineered mESC E–P system does not establish a universal contact radius or rate law. The app does not yet reproduce this fit or its MS2 measurement model. |
| [Waterbury et al. 2024, LSD1 IDR switch](https://doi.org/10.1016/j.molcel.2024.05.017) | Supplementary phase diagrams, co-condensation, recruitment and FRAP measurements provide coregulator/IDR material-property targets | This corepressor study shows that an IDR can inhibit association; more IDR attraction is not a universal rule. Co-condensation or recovery alone does not identify a unique microscopic pair Kd or the app's encounter hazard. |

Fit free protein motion first, then binding fractions and corrected dwell-time distributions, followed by TF/coactivator recruitment and catalytic changes. Test promoter/elongation parameters against nascent output, and processing/decay against released RNA. Match the observation process: the current RNA records are molecule states, not directly an MS2 intensity trace. Use held-out cells, loci or perturbations to validate a fitted model, and compare basal-only, no-recruitment, no-contact-gain and no-chromatin-feedback controls. A good MSD overlay establishes a motion mapping; it cannot by itself validate transcription regulation.

## Scientific next steps

| Stage | Implementation target | Evidence required before using it scientifically |
|---|---|---|
| 1. Quantitative motion and reactions | Native OpenMM forces, recorded trajectories, restart checkpoints and physical-unit mapping; spatial neighbor queries are already implemented | Timestep convergence; dilute diffusion and binding/rebinding benchmarks; equilibrium controls with consistent reaction energetics |
| 2. Protein material properties | Interchangeable coarse sticker/domain and sequence-based IDR force fields; separate folded domains and binding sites | Single-chain sizes, diffusion and binding measurements; concentration/valence sweeps; coexistence and exchange measurements rather than cluster size alone |
| 3. Dynamic nuclear bodies | Emergent scaffold/client mixtures, interacting RNA chains beyond the current passive RNA points, rDNA tethering and separate nucleolar compartments | Partition coefficients, exchange/FRAP and spatial composition; distinguish fixed-body effects from self-assembly |
| 4. Regulatory mechanisms | Calibrate implemented TF/coactivator recruitment, occupancy/access, stochastic pausing and state feedback; extend promoter networks, backtracking, traffic, finite pools and mark exclusion | Burst and dwell-time distributions, nascent versus mature RNA, modification propagation, contact-memory alternatives and perturbation controls |
| 5. Transfer across loci | Load locus annotations and optional chromatin-contact constraints; separate shared molecular parameters from locus-specific inputs | Fit a subset of loci and predict held-out loci/perturbations; propagate parameter uncertainty |

Keep body proximity effects on initiation, elongation, processing and RNA stability as separate candidate mechanisms. The original speckle_gain acts on processing; optional body_responses independently affect promoter activation/initiation. Likewise, transcription-induced stiffness should compete against a no-effect model and, later, alternatives such as active forces or changed compaction.

Before drawing conclusions, compare controls that disable protein–protein binding, recruited catalytic activity, regulatory activation/initiation gains, effective chromatin attraction, mark spreading, cohesin, body affinity and speckle processing gain one at a time. Separately test equal occupied/vacant capture and release factors, zero remodeling rates, equal occupancy persistence factors, zero pause rate, an occupied elongation multiplier of one, and a transcription persistence multiplier of one. Require zero blocked initiations when interpreting individual-RNA output quantitatively. Match molecule/locus motion, occupancy, exchange and gene-output distributions rather than visual appearance alone.

Useful primary references:

- [Ancona & Brackley: chromatin-mediated protein phase separation](https://doi.org/10.1016/j.bpj.2022.05.039).
- [Eeftens et al.: binding, modification and chromatin compaction](https://doi.org/10.1038/s41467-021-26147-5).
- [Owen et al.: polychrom and epigenetic memory](https://pmc.ncbi.nlm.nih.gov/articles/PMC11075759/).
- [Salari et al.: transcription and chromatin dynamics](https://pmc.ncbi.nlm.nih.gov/articles/PMC11199603/).
- [Bhat et al.: speckle proximity and splicing](https://pmc.ncbi.nlm.nih.gov/articles/PMC11164319/).
- [Feric et al.: nucleolar subcompartments and multivalent domain models](https://pmc.ncbi.nlm.nih.gov/articles/PMC5127388/).
- [CALVADOS: sequence-based IDR and multidomain protein models](https://github.com/KULL-Centre/CALVADOS).

## Source layout

`chromatin_studio/schema.py` defines the model, presets and validation. `engine.py` owns reference mechanics, reactions, nucleosome occupancy, regulatory recruitment and read-only genome/interaction/regulation summaries. `transcription.py` advances individual polymerases, growing nascent RNA and passive released RNA with processing/degradation. `chromatin_contacts.py` computes state-selective chromatin energies, forces and contact reports. `kinetics.py` supplies the independent physical binding estimator. `openmm_backend.py` supplies optional polychrom/OpenMM integration. `cluster.py` validates and packages jobs, `batch.py` executes streamed replicas, and `validation.py` runs analytical acceptance checks. `server.py` serves a loopback API and `static/` contains the browser interface, including `genome-view.js`, `interaction-view.js`, `regulation-view.js`, `transcription-view.js` and `cluster-view.js`. `tests/` contains numerical/integration checks and `examples/` contains editable model presets. The `benchmark-*.json` files preserve measured development baselines and their stated limitations.
