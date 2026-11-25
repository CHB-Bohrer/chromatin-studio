<!--
Section: Summary (front page)
Update when: the headline capabilities, the version, or the main validation results change.
Sources: docs/PROJECT_STATUS.md (top), status.json, pyproject.toml (version).
-->

# Summary {.unnumbered}

::: {.box}
**Chromatin Studio in one paragraph.** Chromatin Studio is a program that simulates a stretch of chromatin inside a cell nucleus, the proteins that diffuse around it and stick to it, and the genes on it that switch on and make RNA. You build a model in a web page, press run, and watch it in 3D and along the genome. Every piece of the model (bending of the fibre, protein binding, histone marks, nucleosomes, cohesin loops, the nuclear envelope, transcription) can be switched on or off, so you can ask what each piece does. It is a **research prototype in reduced units**: the code works and is tested, several single mechanisms have been checked against theory and against published experiments, but the full model has not yet been calibrated to real seconds and nanometres.
:::

**What is in this document.** Section 1 explains the goal and the overall design. Sections 2 to 7 describe each mechanism: what it represents, the rule the code follows, and what it deliberately leaves out. Section 8 explains how time is advanced. Section 9 shows the interface. Sections 10 and 11 report what has been checked against theory, against experiments (Donovan et al. on nucleosomes and transcription factors; Kilic et al. on HP1α and H3K9me3) and against live-cell polymer data. Sections 12 and 13 list the limits and the plan. The appendix maps each section to the code it describes.

**Key results so far.**

- The engine matches closed-form answers for nucleosome turnover, promoter bursting (the telegraph model), cohesin stepping, free diffusion and a worm-like chain (all within about 1%).
- Using measured on- and off-rate changes as inputs, the model predicts the measured nucleosome-versus-DNA affinity ratios for four yeast transcription factors, **but only if** nucleosomes are allowed to slow the factor's release. Blocking access alone is not enough for the pioneer-like factors Reb1 and Cbf1.
- Rebinding among neighbouring marked nucleosomes, with one tuned hazard, reproduces how HP1α residence grows with H3K9me3 density up to 50%. It does **not** reproduce how a competing HP1α shortens residence, which points to a missing mechanism.
- A first fit to live-cell Fbn2 locus tracking gives about 47 nm and 0.2 s per reduced unit, but this mapping is still marked *needs review* and is not attached to any preset.

**Reading guide.** Boxes like the one above give the short version of a section. Text marked **[Assumption]** is a modelling choice without a measurement behind it. Text marked **[Unverified]** is something this document states that has not been checked against a run or the code in detail. Numbers come from `docs/PROJECT_STATUS.md` unless another source is named.
