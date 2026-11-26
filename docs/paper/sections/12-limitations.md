<!--
Section: Limitations
Update when: a limitation is removed, or a new one is discovered.
Sources: AGENTS.md "Scientific distinctions to preserve", README.md "What this version represents", docs/PROJECT_STATUS.md table.
-->

# Limitations

::: {.box}
**Short version.** The model is coarse, phenomenological and in reduced units. Read every result as "what these rules imply", not as a measured fact about cells.
:::

**Units and calibration.**

- No preset has an accepted physical clock (Section 11).
- Association inputs are hazards per encounter, not molar $k_\text{on}$ or $K_d$.

**Physics and chemistry.**

- Soft repulsion lets chains pass through each other; topology is not preserved.
- Nucleosomes are on/off states on fixed beads, with no wrapped DNA, histone identity, tails or histone pool.
- IDRs are coarse bead chains, not sequences. Speckles and nucleoli are imposed fields, not self-assembled condensates. Large clusters are not evidence of phase separation.
- Cohesin barriers are perfect, unloaded cohesin reloads instantly, and there is no loading rate.
- Transcription has no free polymerase pool, no polymerase traffic or backtracking, no introns, and RNA is a point particle with no forces or protein binding.
- Tethers are anchor springs, with no lamina, pores or membrane mechanics.
- Many couplings (enhancer contact gates, mark gates, stiffness changes during transcription, body responses) are configurable hypotheses, not fitted laws.

**Numerics.**

- Reaction families and mechanics are operator split; the full model is not an exact stochastic simulation. Check key results at smaller $\Delta t$.
- Saved runs are not exact restart checkpoints; event logs and summaries are bounded.

**Scale and speed.**

- At most 4,000 chromatin beads; 20 Mb at nucleosome resolution is out of reach.
- The GPU only speeds up bare polymers. Models with chemistry run on the CPU.

**Validation.**

- Each mechanism has been checked on its own. The combined model has been explored with ablations, but its biological predictions have not been validated, and no published figure has been fully reproduced.
- The cohesin study against Fbn2 data has not been run.
