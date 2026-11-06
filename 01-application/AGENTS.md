# Guidance for AI contributors

Chromatin Studio is a local research prototype for coupled chromatin mechanics, protein reactions and transcription. Continue from the implemented system; do not assume that a requested scientific capability has already been calibrated or validated.

## Read first

1. [Current project status](docs/PROJECT_STATUS.md): completed work, evidence, limitations and next milestones.
2. [Development guide](docs/DEVELOPMENT.md): architecture, commands and how to maintain this handoff.
3. [Implementation contract](CONTRACT.md): configuration, API and reaction semantics.
4. The relevant section of [README.md](README.md), then the implementation and tests for the requested change.

Use the current user request to determine scope. The roadmap describes unfinished work; it is not an instruction to launch long simulations, submit cluster jobs, publish data or redesign unrelated components.

These files can also be supplied explicitly to an AI that does not automatically read `AGENTS.md`. No particular AI service or prior chat history is required.

## Scientific distinctions to preserve

- The supplied presets have **no accepted empirical physical clock**. Calibration tools, processed reference data and synthetic tests exist; the production length/time fit has not been completed. A polymer clock converts units but does not independently fit protein or transcription kinetics.
- Protein association inputs are conditional encounter hazards and affinity multipliers. Do not relabel them as measured molar `kon` or `Kd`. The diffusion-limited estimator is a separate idealized calculation.
- OpenMM/polychrom currently integrates motion while the host CPU evaluates forces and reactions. CUDA support is not evidence of full-model GPU acceleration or completed GPU validation.
- Current bounds are 4,000 chromatin beads and 6,000 chromatin/protein particles; individual RNA has a separate 5,000-record limit. A 20 Mb model at 200 bp per bead needs 100,000 chromatin beads and is unsupported.
- Nucleosomes are occupied/vacant genomic sites, with coarse turnover and sliding. IDRs are domain-and-sticker chains. Speckles and nucleoli are imposed spatial fields. Do not imply wrapped DNA, residue-level proteins, self-assembled organelles or validated liquid–liquid phase separation.
- Promoter paths and initiation are stochastic within frozen-rate reaction intervals; mechanics and reaction families are operator split. Do not call the complete coupled simulation an exact global Gillespie process.
- Save-run exports and cluster bundles are not exact stochastic restart checkpoints. Sampling summaries and event buffers are bounded; they are not complete trajectories.
- Numerical checks, an analytical reference-model comparison, a published-figure reproduction and empirical biological validation are different levels of evidence.

## Working rules

- Preserve the user's active model and trajectory. Use isolated simulations or a separate local instance for tests and recording. Inspect current processes before starting or stopping a server; a port mentioned in a handoff is not proof that a server is running.
- Keep draft edits separate from the applied run. Apply/restart is an explicit state change. Studies and cluster export operate on independent configurations.
- Keep snapshots, rate diagnostics and UI summaries read-only; they must not draw from the model RNG or mutate physics. Demonstrate events from real engine transitions.
- Preserve optional defaults and backwards compatibility unless the user requests a behavior change. Explain changes that invalidate seed replay or physical calibration.
- Run checks that match the change. Record skipped platforms, the actual backend and the scope of any scientific claim. Never turn missing GPU hardware into a passing GPU result.
- Use existing tutorial sources and capture real UI states. Preserve natural narration at 1×; do not speed speech to fit old timings. See [tutorial instructions](tutorial/README.md).
- The user authorized public release of the tool, hypothesis examples, tutorial and development handoff on 2026-11-06, excluding the complete manuscript. Keep the full manuscript, full-text extractions, credentials and machine-specific local notes out of the repository. Included simulated examples/reference summaries are intentional public tool content. Do not treat this as permission to publish future experimental data automatically.
- Task ownership is temporary. Coordinate concurrent edits, but do not infer permanent ownership from names in an earlier chat.

## At the end of substantive work

Update [PROJECT_STATUS.md](docs/PROJECT_STATUS.md) with the date, completed behavior, evidence, remaining work and the next concrete step. Update the contract when interfaces or semantics change. Keep old reports attributable to the source that produced them; do not alter an old result to make it look current.

For documentation-only changes, check links and factual consistency. Do not rerun expensive scientific workloads merely to refresh a date. Report what actually ran in the current session.
