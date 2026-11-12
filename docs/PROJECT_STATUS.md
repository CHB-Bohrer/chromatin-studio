# Chromatin Studio project status

Last reviewed: **2026-11-09**. Application: **0.10.0**. Configuration format: **0.1**.

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
| Whole nucleus and tethers | Optional soft ellipsoidal envelope around the confined territory; genomic tethering elements (e.g. LADs) whose beads attach stochastically to an envelope or body-surface anchor and feel a harmonic spring; attach/detach hazards can depend on a mark and on nucleosome occupancy | Phenomenological anchors: no lamina meshwork, pores, membrane mechanics or anchor diffusion. NumPy engine only; the OpenMM backends refuse these models. |
| Regulatory elements | Enhancer/promoter regions, TF preferences, explicit coactivator recruitment and configured contact-dependent rate gains | Mechanisms and parameters are hypotheses. No DNA sequence search or universally fitted enhancer model. |
| Transcription | Stochastic promoter switching and initiation; optional individual polymerases, elongation, pausing, RNA release/processing/decay | Rates are frozen within each reaction interval. No polymerase traffic/exclusion, finite polymerase pool or RNA-chain forces. |
| Chromatin feedback | Occupancy and engaged transcription can alter local persistence; polymerases can affect occupancy | These are configurable phenomenological couplings. |
| Cohesin | Stochastic loop extrusion with directional barriers, collisions and unloading/reloading | No independently parameterized loading rate or completed published contact-map reproduction. |
| Nuclear environments | Separately configurable speckle and nucleolar spatial affinity fields and transcription responses | Fixed imposed fields; no nucleolar subphases, rRNA biogenesis or self-assembled organelles. |
| Inspection | 3D view, live 1D binding/occupancy, full genomic tracks, interaction context, regulatory rate decomposition and individual RNA views | Read-only snapshots; free proteins retain 3D positions instead of invented genomic positions. |
| Studies | Independent control/ablation runs and reference calculations | Quantitative manuscript conclusions and full spatial matched-allele comparisons remain unfinished. |
| Cluster execution | Validated draft export, Slurm/job-array scripts, headless execution, streamed results and provenance | Fresh runs only; the app does not submit jobs or provide exact checkpoint continuation. |

Detailed definitions and coordinate conventions are in [CONTRACT.md](../CONTRACT.md) and [README.md](../README.md).

## Tutorial scale-wording revision (2026-11-12)

A viewer read "each bead is 10 kb" as contradicting "bead = nucleosome". Both are valid: bead size is a setting (`bp_per_bead`); the 20 Mb overview uses 10 kb bins, and the detailed 80 kb and 25.6 kb examples use 200 bp coarse nucleosome units. The tutorial opened with the 10 kb example without saying so.

- Rewrote chapters 01 and 02 (narration, transcript and captions) to state bead size as a setting and to tie 200 bp to one coarse nucleosome unit.
- Chapter 12 now says only bare-polymer models run fully on the GPU (native backend); models with proteins still evaluate forces and reactions on the CPU.
- Audited the other chapters against the contract and status; no other scale or limit claims needed changes.
- Added a future-update checklist to [tutorial/README.md](../tutorial/README.md).
- **Regenerated later the same day:** installed Python 3.12 and Kokoro on the Windows machine (venvs outside the repo), narrated scenes `01-overview`, `02-build-and-learn`, `12-cluster` and `17-nucleus-tethers` (the other 13 reused), and rendered all 17 chapters (691.8 s, natural 1× pace). The old video lacked chapter 17.
- `render_tutorial.py` gained a Windows Arial fallback. On Windows run with `PYTHONUTF8=1`, or the renderer fails on the "≈" in a caption. The renderer wrote absolute local paths into `render-manifest.json` and `render-timeline.json`; they were rewritten to repo-relative `narration-neural/...` before committing.
- **Not done:** captures were not retaken (chapter 01 still shows the 20 Mb capture, which matches the narration). The video was rendered but not watched at 1×, so audio, key-point readability and transitions are unreviewed.
- Next concrete step: watch the video at 1× and fix any problems.

## Latest completed increment

**Whole nucleus, tethering elements and occupancy views (2026-11-09, uncommitted on `feat/nucleus-tethering`).**

- [`nucleus.py`](../chromatin_studio/nucleus.py): soft ellipsoidal envelope `U = k/2 |u|²(1 − 1/ρ)²` outside the envelope, with an exact gradient. For a sphere it reduces to the confinement-wall law. RNA points also feel the envelope.
- Tethering elements (`nucleus.tethers`, up to 128): genomic intervals whose beads attach one at a time, with these rules:
  - **Anchor:** the radially projected point on the envelope, or on an enabled body surface.
  - **Attach:** a bead attaches when it is within `capture_distance`. The hazard is `attach_rate × occupancy_attach_multiplier[state]`, gated by an optional `required_mark`.
  - **Detach:** the hazard is `detach_rate × occupancy_detach_multiplier[state]`.
  - **Events:** each transition emits `tether_attach` or `tether_detach`.
  - **Replay:** a model without a nucleus draws no extra random numbers, so legacy seeds replay unchanged.
- Release multipliers may now reach 1000. Occupied-site release can therefore follow measured nucleosome-induced dissociation ([Luo et al. 2014](https://doi.org/10.1093/nar/gkt1319); [Donovan et al. 2019](https://doi.org/10.15252/embj.2018100809)).
- Study conditions `no_tethering` and `no_nucleus` were added.
- Interface changes:
  - **3D view:** envelope rings, a **Whole nucleus** framing toggle, amber tether links with anchor markers, and hollow vacant-nucleosome beads (**Occupancy** toggle).
  - **1D Live binding:** Cohesin, mark and **Tethers** rows (element intervals, attached beads, attach/detach glyphs). The Nucleosomes row now keeps every vacant site visible at any zoom, with an occupied-fraction count.
  - **Bodies editor:** a **Whole nucleus & tethering** section for these settings.
- The tutorial gained chapter 17, "Add the whole nucleus and tethering elements". It uses four captures of real states from an isolated UI-check instance with accelerated tethers. The storyboard validates in captions-only mode (17 scenes, 760 s). Narration and the rendered video were not regenerated, because the Kokoro voice model is not installed on this machine.
- Evidence: 11 new tests in [`tests/test_nucleus.py`](../tests/test_nucleus.py) cover the sphere law, the ellipsoid gradient, force/energy consistency, events, the read-only snapshot, the capture range, occupancy/mark hazards, legacy replay, ablations, native refusal and the schema. Full suite: 250 passed. The 5 `test_tutorial_server.py` tests were deselected because they need the Windows symlink privilege (WinError 1314), which is unrelated to this change.

**Should the polymer have constraint forces?** Not as rigid bonds. In the equilibrated Fbn2 baseline the bond length is 0.90 ± 0.15 units (16%). Only 0.5% of bonds exceed 1.5× rest length, and no bead pair overlaps deeply. Holonomic constraints would therefore change little, and they would alter the dynamics the clock was fitted to. At 10 kb per bead, the bead is a fluctuating blob, so a fixed length is not physically motivated. The soft repulsion (barrier about 4 kT) does allow occasional chain crossing. Stronger excluded volume is the option to consider only if topology or entanglement matters to a question.

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

**First baseline calibration attempt (2026-11-07), status `needs_review`, not accepted.** The protocol's baseline model (fingerprint `a2d9b040…`) was run with the native CUDA backend (mixed precision, RTX 3090): 8 replicas of 24 million steps each (48,000 reduced units, seeds 20261102–20261109), plus 4 replicas at dt = 0.001 (16,000 units) for the timestep check. Two replicas were extended from 16,000 units by Markov continuation with a fresh noise seed. Raw outputs are local and not committed, and no mapping was attached to the model. The analysis used pairs 974–1025/1026 (510/520 kb), the RAD21-2h fit subset and a held-out acquisition, with no refitting:

- **Length is consistent: about 50 nm per reduced unit.** Pair 974–1025 gives 49.7 nm with a replica-bootstrap 95% interval of 44.9–56.8 nm. The value stays within 49–51 nm across the 510/520 kb pairs, burn-ins of 4,000, 8,000, 12,000 and 24,000 units, and the pooled interior-pair estimator.
- **Time fits the training curve but is not robust.** With an 8,000-unit burn-in and a 60–1,200 s window, the fit gives 0.217 s per unit (bootstrap 0.15–0.36). Fit quality: log RMSE 0.016; MSD slope 0.54 (simulation) vs 0.54 (experiment); the optimum is not at the search boundary. The 40–600 s and 100–2,000 s windows agree. The estimate depends on burn-in, though: 0.24, 0.21 and 0.18 s per unit at burn-ins of 4,000, 12,000 and 24,000. At 24,000 the simulated slope drops to 0.45.
- **Not stationary.** Interior 51-bead R² falls steadily from 178 to 156 over 48,000 units; the post-burn-in drift has z = −3.1 (−1.6 after a 24,000-unit burn-in). Rg levels off after about 18,000 units. Single-pair R² is dominated by slow fluctuations.
- **Timestep not converged.** Over matched seeds and window, dt = 0.001 changes interior R² by −5% and the fitted time at fixed length by +18–21%.
- **The held-out test fails, but the failure reflects the data.** The holdout gives log RMSE 0.355 against a threshold of 0.25, with bias −0.35. The two experimental acquisitions themselves differ by a mean log ratio of +0.353 over 60–1,200 s, with slopes 0.54 vs 0.49. No single mapping can pass this criterion as written. Their noise-corrected second moments agree more closely (predicted 0.402 vs observed 0.423 µm²).

The step proposed after this attempt was to equilibrate longer and repeat the dt check, adding dt = 0.0005. Separately, decide with the data owner whether the holdout criterion should account for between-acquisition variability.

**Extended baseline run (2026-11-07, same day), still `needs_review`.** All 8 replicas were extended by Markov continuation to 96,000 reduced units. The analysis discards the first 48,000 units. No mapping was attached to the model.

- **Stationarity is approximately established after about 48,000 units.** Rg levels off at 13.3–13.5, and bond lengths are flat. Interior 51-bead R² fell about 13% over the first half (177 → 154). In the second half, a small downward trend remains that is not significant (z = −1.4). Results barely move across burn-ins of 24,000, 36,000, 48,000 and 64,000 units: length 46.9–48.0 nm and time 0.19–0.21 s per unit. Block autocorrelation analysis is not done; uncertainties are replica bootstrap.
- **Length is about 47 nm per unit, but depends on the estimator.** The protocol's pair estimator gives 47.2 nm for 974–1025 (95% interval 40.6–53.9) and 46.7 nm for 974–1026. Pooling all interior pairs at the same separation gives 51.3 nm, about 9% higher. This means 8 replicas do not fully average the slow fluctuations of one specific pair.
- **Time is about 0.20 s per unit and fits the training curve.** Pair 974–1025 gives 0.205 s (bootstrap 0.13–0.32) and pair 974–1026 gives 0.198 s, both `fitted`. Fit quality: log RMSE 0.034; MSD slope 0.57 (simulation) vs 0.54 (experiment); the optimum is not at the search boundary. The 40–600 s and 100–2,000 s windows agree within 0.193–0.212 s. With the pooled-interior length, the time is 0.278 s per unit.
- **Timestep is converged for well-resolved observables.** Four replicas each restarted from their 96,000-unit state and ran 8,000 units at dt 0.002, 0.001 and 0.0005, with distinct noise. Bead MSD and interior relative-pair MSD agree within 0.5% at lags of 4–400 units, and within 2% at 2,000 units, inside the noise. Bond variance shifts by −0.4% to −0.6%. The single-pair time fit on these short runs differs by 22–30%, but the dt = 0.002 fit hits its search boundary there, so that comparison is too noisy to test the timestep. The earlier 18–21% dt sensitivity was confounded with relaxation.
- **The held-out test fails, unchanged: log RMSE 0.357, bias −0.35.** The two experimental acquisitions differ from each other by the same amount.

**Remaining protocol checks (2026-11-08), using the same 96,000-unit data with no new simulation.** Status is still `needs_review`.

- **The time estimate is limited by the noise model.** Length barely changes with the noise correction. In the table, s scales the localization-noise offset, where s = 1 is the reference's correction (relative-vector error 0.203 µm, taken from the polysim bootcamp tutorial and not measured for this dataset):

  | Noise offset scale s | Length (nm/unit) | Time (s/unit) | Fit log RMSE | MSD slope, sim/exp | Status |
  | --- | --- | --- | --- | --- | --- |
  | 0 (uncorrected) | 49.6 | 0.12 | 0.158 | 0.55/0.34 | shape incompatible, search boundary |
  | 0.8 | 48.1 | 0.162 | 0.092 | 0.56/0.44 | fitted |
  | 1.0 | 47.2 | 0.205 | 0.034 | 0.57/0.54 | fitted |
  | 1.2 | 46.1 | 0.304 | 0.120 | 0.60/0.74 | fitted |

  Values are for pair 974–1025; pair 974–1026 agrees within 3%.

  A ±20% change in the noise parameter therefore moves the clock by −21% to +48%. The shape agrees best at s = 1. An independent measurement of the noise for this dataset would constrain the clock more than additional simulation.
- **Distance distribution.** The comparison passes the simulated pair vectors (scaled by the fitted length) through the same isotropic Gaussian noise model and compares them with the observed |R| histograms. The second moment matches by construction. Mean |R| is 0.621 µm in the simulation, against 0.582 µm for the fit set (95% interval 0.532–0.628) and 0.617 µm for the holdout (0.573–0.665). The simulated distribution is narrower than the fit set's (10th/50th/90th percentiles 0.32/0.61/0.94 vs 0.23/0.52/1.02 µm; KS distance 0.15). It is closer to the holdout's (0.28/0.57/1.01 µm; KS 0.06). This is consistent with cell-to-cell heterogeneity that a single homogeneous chain does not represent. The frames are autocorrelated, so no formal test p-value is claimed.
- **Holdout band.** The fixed-mapping MSD prediction lies inside the fit set's movie-cluster bootstrap band at every lag from 60 to 1,200 s. It lies inside the holdout's band at none of them (prediction/holdout median 0.72). The fit-set experimental curve itself lies inside the holdout band at 0% of those lags. So the two acquisitions differ beyond their own movie-bootstrap uncertainty; this is not just sampling noise.
- **Autocorrelation.**
  - The integrated autocorrelation time (Sokal window) of the pair R² is about 3,300 reduced units (median; 1,600–5,100 across replicas), and about 2,900 units for interior R².
  - The retained 48,000 units give about 8 effective samples per replica and about 70 in total. The 48,000-unit burn-in is about 15 τ.
  - The replica SE of pair R² (26.8 on a mean of 180) is about twice the block and τ-based SEs (12–13). Replica-level variation slower than the window therefore remains, and the replica bootstrap, not overlapping-lag counts, remains the right uncertainty.
  - The pair-vs-pooled-interior length difference (R² 180 vs 152) is about 1 replica SE. It is consistent with sampling noise of the single pair. The pooled estimator is about 13× more precise in this homogeneous model.

Remaining for acceptance:
- a holdout tolerance or acquisition model that accounts for the demonstrated between-acquisition difference;
- an independent estimate of localization noise for this dataset;
- a decision whether to adopt the pooled-interior estimator, valid here because the chain is homogeneous away from its ends, or to add replicas (about 4× for half the pair SE);
- the 5 kb alternative resolution, if 5 kb bins are adopted.

Even after that fit, protein diffusion, residence times, binding probabilities, catalytic rates and transcription parameters need their own observations. One shared clock provides consistent units; it does not infer all those parameters. The [kinetics estimator](../chromatin_studio/kinetics.py) does not turn the interacting model's encounter hazards into measured macroscopic binding constants.

## Literature-timescale model and first ablation set (2026-11-08/09, exploratory)

This is a provisional exploration, not a calibration. It converts literature timescales to reduced hazards with the `needs_review` clock (0.205 s/unit, range 0.13–0.32), so every converted rate scales with that unaccepted clock. Mechanics match `manuscript_20mb` and the Fbn2 baseline (2,000 × 10 kb), so the clock's fingerprint scope applies. Models, source tables, scripts and outputs are local run data (see `LOCAL_DATA.md`), not committed.

- **v1, binding and cohesin.**
  - Residences: activator 13.3 s (Sox2/Oct4 single-molecule tracking), p300 26 s, PRC2 257 s. HP1 and the remaining species are unchanged; each has a recorded reason.
  - Cohesin: 129 bound extruders per 20 Mb, 22 min residence, 0.5 kb/s loop growth.
  - First version scaled on_rate with off_rate. A 41 s pilot showed association became far too slow.
- **v2, all feasible processes.**
  - Association: activator and PRC2 on_rates set from measured effective association, targeting 23.5% and 20% bound.
  - Genes: ON 5 min, base OFF 30 min, one initiation per 10 s while ON, processing 5 min, mRNA half-life 7.1 h.
  - Elongation: clamped at the schema floor of 29 kb/min, about 10× faster than in vivo.
  - Marks: H3K27ac half-life 1.5 h; methyl marks 1 day as a lower bound for "days"; H2AK119ub1 1.5 h (assumption, no source).
  - Nucleosomes: turnover about 1 h at active genes.
  - Rates with no literature value (spreading, catalysis, background nucleosome hazards) were scaled by the same factor as their process. This preserves hypothesis balances, which remain hypotheses.
- **Ablation set.** Baseline plus 11 study ablations (all polymer conditions except `no_confinement`). Each ran 731,350 NumPy steps, about 5 min biological time, as 1 replica with a shared seed, 16–25 ms/step with 8 concurrent. Findings, using the mean over the last 60 s compared with baseline:
  - **Cohesin.** Loops are still growing at 5 min (mean 86 kb, max 215 kb, +10 kb/min).
  - **Removing cohesin** gave the clearest effect: Rg +0.63 units (baseline drift −0.11/min); 10-bead contact probability 0.040 → 0.015.
  - **Removing body affinity** halved protein–protein bonds and expanded the chain by about 0.24 units.
  - **Removing protein binding or chromatin attraction** expanded the chain by about 0.2–0.3 units.
  - **Bound fractions overshot the targets.** The activator settled near 0.55 bound against a 0.235 target; the linear on_rate assumption fails. HP1 is about 0.98 bound and p300 about 0.57. PRC2 (4 molecules) bound 0–1 molecule. Bound fractions therefore still need iterative tuning.
  - **Slow processes barely moved.** Marks and nucleosome occupancy changed by a few beads.
  - **Transcription was not exercised.** Two genes turned ON, giving 26 initiations; no transcript completed, because transit takes about 4 min.
  - **No measurable effect.** `no_regulatory_gain`, `no_body_response` and `no_speckle_processing` produced no measurable change. The last two gave identical trajectories.
  - **Noise.** All other differences in gene, mark and bound-fraction values are within single-replica trajectory divergence. With 1 replica there is no noise estimate, and polymer observables are still relaxing (τ ≈ 11 min biological time).

- **v3 model: nucleus, LAD tethers and occupancy-dependent kinetics (2026-11-09).**
  - **Nucleus and LADs:** a 5 µm nucleus (105.9 units) whose envelope cuts the 20-unit territory at x = 15. The three H3K9me3 domains (beads 80–249, 720–899, 1450–1619) serve as LAD tethers. Each tether attaches at 1/min within 71 nm of the envelope and detaches at 1/h. Attachment needs H3K9me3 and an occupied nucleosome; vacant sites detach 10× faster. These hazards are assumptions.
  - **Activator:** occupied-site release ×300. Basis: Gal4 dwell is reduced "by orders of magnitude" on nucleosomes (Donovan et al. 2019 EMBO J), and nucleosomes accelerate dissociation up to ~1000× (Luo et al. 2014).
  - **Nucleosome readers and enzymes:** need an occupied site (vacant capture 0) and release 10× faster from vacant sites; this is an assumption.
  - **Pioneer kinetics not used:** Reb1/Cbf1 dissociation-rate compensation (Donovan et al. 2019 eLife) is recorded but not applied, because the model has no pioneer species.
- **Retuning bound fractions.** Two pilot rounds were run, 60 biological seconds each, warm-started from equilibrated Fbn2 conformations.
  - **Round 1:** 5 on_rate brackets. The occupancy kinetics collapsed the activator's bound fraction to 0.01–0.06, because the ×300 release makes most encounters transient.
  - **Fitting:** HP1 and PRC2 were fitted on logit(f) against log on_rate. PRC2 used the equilibrium prediction k_on,eff/(k_on,eff + k_off), because it relaxes over about 190 s. The activator fit was noise-limited, so its on_rate came from the linear scaling of effective association.
  - **Round 2 (tuned model):**

    | Species | Bound fraction | Target | Note |
    |---|---|---|---|
    | Activator | 0.19 / 0.24 | 0.235 | on_rate 1.4; mean completed dwell 0.7–0.8 s, not the 13.3 s SMT residence |
    | HP1 | 0.48 / 0.61 | 0.5 | target is an assumption |
    | p300 | 0.43 at 2.5× on_rate | ~0.54 | final on_rate 3.9×, chosen by log-odds slope |

  - **Activator dwell mismatch:** with occupied sites dominating, matching both the bound fraction and the SMT residence would need far more accessible sites, or a pioneer-like occupied-site release. This is an open modeling question, not a tuning failure.
  - **Outputs:** the tuned model and tuning report are in `literature-kinetics/` (`model-literature-v3-tuned.json`).
- **The 60 s pilots underestimated the 5-minute state.** Bound fractions kept rising through the run. Two full 5 bio-min baselines at the round-2 rates ended (240–300 s window, 2 replicas) at:

  | Species | Bound fraction | Target |
  |---|---|---|
  | Activator | 0.52 / 0.66 | 0.235 |
  | HP1 | 0.86 / 0.78 | 0.5 |
  | p300 | 0.95 / 0.97 | ~0.54 |
  | Writer | 0.97–1.0 | no target |
  | PRC2 | 0 | 0.2 |

  They are in `ablation-v3-5min/baseline-r*`. The first paired ablation launch was stopped after 1.3 h to avoid 11 h of runs at mis-tuned rates.
- **Tuning round 3** (6 variants × 5 bio-min, `literature-kinetics/tuning-r3`). Each fit used the 240–300 s means:
  - **HP1:** well determined (logit slope 2.2, residual 0.23); fitted on_rate 0.26.
  - **p300:** 0.59 at on_rate 0.028, so it was set to 0.025.
  - **PRC2:** with 4 molecules the fraction ranges 0–0.75 per run, mean about 0.16 over 8 runs; left unchanged.
  - **Activator:** stayed at 0.5–0.8 for every on_rate from 0.7 to 1.4. It appears to saturate the ~28 vacant sites, where release is not accelerated; 0.6 × 50 molecules ≈ 30. That is the expected consequence of Donovan-style occupied-site release: long dwells happen only at nucleosome-free sites, so site supply rather than association sets the bound fraction. Remodeler eviction may also add vacant sites over time; nucleosome occupancy fell from 0.988 to 0.986.
- **Tuning round 4.** 8 × 5 bio-min runs (`literature-kinetics/tuning-r4`), with HP1 at 0.26 and p300 at 0.025, scanning the activator on_rate.
  - **Activator response is switch-like:** 4–5 min bound fraction 0.50 at on_rate 0.7, 0.28 at 0.35, 0.02–0.05 at 0.175, and below 0.01 at 0.044 or less. This is consistent with positive feedback: TF-recruited remodeler eviction, plus long dwells only at vacant sites.
  - **Other species across the scan:** HP1 0.47–0.59 (target 0.5); p300 0.56–0.69 (target ~0.54); PRC2 averages about 0.14 (target 0.2; 4 molecules).
  - **Final rates:** `model-literature-v3-tuned2.json` uses activator 0.32 (log-odds interpolation toward 0.235), HP1 0.26, p300 0.025 and PRC2 0.00286 (with the tuning report).
  - **Caveats:** the targets are matched at 4–5 bio-min after a cold protein start, not at a proven steady state. Because the activator's response is steep, expect replica spread.
- **Ablation set v3b (complete 2026-11-10; 16/16 runs, 31 ms/step with 8 concurrent).** It uses the tuned-2 model; output is in `ablation-v3b-5min/` (`summary.json`, `summary.csv`). All values are means over 240–300 s and 2 paired replicas.
  - **Baseline bound fractions:** activator 0.19 (target 0.235), HP1 0.55 (0.5), p300 0.74 (0.54; still high), writer 0.97, PRC2 0 (4 molecules).
  - **Cohesin.** Removing cohesin gave the clearest effect again: 10-bead contact probability fell from 0.040 to 0.013, and Rg rose from 12.7 to 13.1.
  - **Occupancy-dependent access.** Removing it raised the activator from 0.19 to 0.97 bound. In this model, occupied-site release is what limits activator binding.
  - **Tethering is conformation-dependent.**
    - Replica 0 held about 13 LAD1 beads at the envelope (LAD1 within 1 unit of the envelope: 3.5%; about 0 without tethers). Replica 1 had no LAD near the envelope, so it had no attachments.
    - LAD2 and LAD3 never reached the periphery in 5 minutes: the territory touches the envelope only in a small cap.
  - **Not established with 2 replicas:**
    - **Activator under `no_tethering` and `no_nucleus`:** it reached 0.34 and 0.28. Replica 1 had no attachments in either arm, yet differed by 0.06. That difference is stochastic divergence, which shows the noise scale of this switch-like species.
    - **H3K27ac:** +8 to +12 marked beads in nearly every condition. This points to a low baseline rather than many separate effects.
    - **PRC2:** 0–0.375 with 4 molecules.
    - **Transcription:** 0–25 initiations.
  - **What would make the tethering and periphery questions answerable:**
    - Place the territory so it overlaps the envelope more.
    - Start LAD domains near the periphery, as LADs are positioned in vivo.
    - Run several more replicas and tens of minutes.
    - Use a larger p300 molecule pool, or lower its on_rate, if its bound fraction matters. Conditions: baseline, `no_tethering`, `no_nucleus`, `no_nucleosome_access`, `no_nucleosome_turnover`, `no_cohesin`, `no_protein_binding` and `no_body_affinity`. Each condition has 2 replicas of 731,350 steps (about 5 biological minutes), in a paired design: replica k uses seed 1000+k and Fbn2 warm start r{k} in every condition. `ablation-v3-5min/` keeps only the two round-2-rate baselines. `analyze_ablation_v3.py` reports replica half-ranges and flags deltas that exceed them.

Next concrete step: analyze `ablation-v3b-5min` once it completes (`analyze_ablation_v3.py`). For transcription and mark questions, start runs from a pre-equilibrated polymer and cohesin state, and use replicas and longer windows (tens of minutes to hours), or GPU-native chemistry.

## GPU and scale status

The [OpenMM adapter](../chromatin_studio/openmm_backend.py) performs stochastic integration while host code calculates mechanical forces, reaction searches and RNA motion. Coordinates, forces and mobilities cross the host/device boundary every integration substep. Measured on an RTX 3090 (2026-11-07), this makes CUDA *slower* than the NumPy engine: 2.8–3.6 vs 2.1 ms/step for the 2,000-bead baseline.

**Native bare-polymer backend (2026-11-07).** [`openmm_native.py`](../chromatin_studio/openmm_native.py) expresses the engine's backbone, bending, soft-repulsion and soft-wall terms as OpenMM forces, and integrates them with the same Euler–Maruyama update entirely on the device. It accepts only models whose reaction step is a no-op. Proteins, marks, state interactions, bodies, cohesin, genes, regulatory regions, transcription and nucleosome dynamics are refused with an explicit error. Evidence on the RTX 3090 with the 2,000-bead Fbn2 baseline:

- Forces and energies match the engine to about 1e-15 relative on OpenMM Reference (double precision), and to about 1e-6 on CUDA mixed precision. This holds for relaxed, stretched/wall-penetrating and overlapping conformations ([tests](../tests/test_openmm_native.py)).
- Comparing 8 seeds per backend (10k burn-in + 40k steps), bond length/variance, bending, Rg, single-bead MSD and the all-pairs 510 kb separation MSD at lags of 200/2,000/10,000 steps all agree within 0.3% (|z| ≤ 1.1).
- Throughput is 60 µs/step in mixed precision, about 35× the NumPy engine. One replica saturates the GPU, so concurrent replicas do not raise aggregate throughput.

Trajectories are statistically equivalent to the NumPy engine, not identical: OpenMM draws its own noise. Native GPU paths for models with chemistry remain unimplemented.

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
| 2 | Complete baseline physical calibration | First attempt and 96,000-unit extension run 2026-11-07 (`needs_review`, see above; stationarity and dt now largely addressed, holdout fails on between-acquisition variability). Remaining: run the matching model, establish stationarity, constrain length, fit MSD shape/time, evaluate held-out data, and record uncertainty plus timestep dependence. Attach a mapping only after the declared acceptance criteria pass. |
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
