<!--
Section: Validation against theory and experiments
Update when: any component study, extra check or experimental comparison is added or rerun.
Code: chromatin_studio/components/*.py, chromatin_studio/benchmarks.py, docs/paper/analysis/extra_checks.py.
Data: docs/paper/data/*.json (collected by docs/paper/analysis/collect.py).
Sources: docs/COMPONENT_STUDIES.md, docs/PROJECT_STATUS.md "Experiment-matched component studies".
-->

# Validation: does each part work?

::: {.box}
**Short version.** Every mechanism has been run on its own and compared with a closed-form answer; all of those checks pass (figures in Sections 2 to 7). Against experiments, the nucleosome study shows that pioneer-like factors need slower release from nucleosomes, and the HP1α study reproduces the effect of mark density but not of competition.
:::

## Three kinds of evidence

Every check carries one of three labels, and they are never mixed:

| Label | Question it answers | What a pass shows |
|---|---|---|
| Numerical | Does the engine follow its own stated rule? | The code does what it says. |
| Analytical | Does it match a textbook closed-form result? | The mechanism is computed correctly. |
| Experimental | Does an output fall inside a published measurement? | Agreement **under the stated unit mapping**. Not a calibration. |

Experimental checks are further tagged. **Input-matched** means the measurement was used to set the rate, so the check only confirms units. **Tuned** means one parameter was scanned to hit the measurement. Only a **prediction** actually tests the model.

## Summary of mechanism checks

@include data/checks-table.md

The analytical checks were rerun for this document on the current code (NumPy engine, Linux, 4 processes), except where the table says otherwise. Standard-scale runs from 2026-11-15 on the project's Windows machine gave the same outcomes (`docs/PROJECT_STATUS.md`).

## Experiment 1: nucleosomes and transcription-factor kinetics (Donovan et al.)

**The question.** Single-molecule experiments measured how much a nucleosome slows a factor's binding and its release, for four yeast factors. Separately, titrations measured the overall affinity for nucleosomal versus naked DNA. If the two kinds of measurement are consistent, plugging the kinetic changes into a 3D search with rebinding should reproduce the affinity ratios.

**The model.** A short chain whose binding-site beads either have a static nucleosome or not, and one factor that binds only the site beads (one per bead). The measured binding slow-down is an input. The release rule is run two ways: release unchanged (the nucleosome only **blocks** access), or release slowed by the measured factor.

![**Nucleosome/DNA affinity ratio, simulated against measured.** Values above 1 mean the nucleosome weakens binding. Grey bars: measured ranges from titrations; a model value passes if it lies within a factor of 2 of that range, which is why the Reb1 square (0.37) passes. Blocking access alone (circles) overshoots by about 100-fold for Reb1 and Cbf1. Adding the measured slower release (squares) brings both into range. Pho4 and Gal4 pass with blocking alone. Values from `docs/PROJECT_STATUS.md` (24 runs, 2026-11-22).](figures/exp-donovan.svg){#fig:donovan width=82%}

**What it means.** Steric blocking alone cannot explain why pioneer-like factors bind nucleosomes and naked DNA about equally well; the nucleosome must also slow their release, as the experiments report. Because both kinetic ratios were inputs, this shows the kinetic and titration measurements are mutually consistent once rebinding is included. It does not explain *why* release is slower.

## Experiment 2: HP1α and H3K9me3 density (Kilic et al. 2015)

**The question.** HP1α stays longer on nucleosome arrays with more H3K9me3, which the authors attribute to fast rebinding at neighbouring marked nucleosomes. Can rebinding alone reproduce this, and the effect of adding unlabelled competitor HP1α?

**The model.** A 12-bead array (one bead per nucleosome) with 0 to 100% of beads marked, and one tagged reader. Release from marked and unmarked beads is set from the measured residences. One association hazard is tuned to the 100%-array residence ratio; everything else is a prediction.

![**HP1α residence on arrays, relative to a single marked nucleosome.** Left: against H3K9me3 density. Right: against competitor occupancy. Diamonds with bars: Kilic et al. 2015. Circles: model predictions as recorded in the project status. **[Unverified]** The status text does not say clearly whether these predictions are at the tuned hazard 1.0 or an interpolated hazard near 0.8; the 100% point (2.21) is at hazard 1.0. The model follows the density trend up to 50% but overshoots slightly at 75–100%. With competitor the model can only fall to the no-rebinding level (dashed line, ratio 1), while the measurements fall below it. Model values from `docs/PROJECT_STATUS.md` (108 runs, 2026-11-22); model error bars were not recorded there and are omitted.](figures/exp-kilic.svg){#fig:kilic width=100%}

**What it means.** 8 of 14 checks pass. Rebinding plus mark-dependent release explains the density effect. It cannot explain competition: in the model a competitor only removes sites to rebind to, which at most brings residence down to the single-nucleosome level. The measured residence goes lower. That points to a mechanism the model lacks, such as competitor-facilitated dissociation or HP1α–HP1α interactions. Interpolating the hazard scan suggests a hazard near 0.8 would fit the tuned point better; this has not yet been run.

## Experiment 3: Sox2 search on bare chromatin

The `tf_binding` study takes the measured Sox2 bound residence times as inputs and scans the association hazard. 13 of 15 checks pass. The free interval between bindings scales as one over the hazard (9.9, 4.9 and 2.5 time units for hazards 0.08, 0.165 and 0.33). The measured 3.5 s lies between scan points, so the chosen point fails the interval and the trial count while its bound fraction (0.27 against 0.235 ± 0.047) passes. Interpolation suggests a hazard near 0.235 would meet all three. Lower diffusion barely changes the averages but makes rebinding cluster in time.

## Experiment 4: cohesin and the Fbn2 locus

The `fbn2_cohesin` study compares anchor separation with and without cohesin and CTCF against the Fbn2 live-cell data. It has not been run yet.
