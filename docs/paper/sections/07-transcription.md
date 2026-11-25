<!--
Section: Regulatory elements and transcription
Update when: enhancer/promoter rules, mark gates, recruitment, promoter switching, initiation, elongation,
             pausing, RNA processing/decay, or gene-mechanics couplings change.
Code: chromatin_studio/engine.py (_regulatory_region_records, _region_recruitment, _gene_regulatory_rates,
      _legacy_gene_on_rate, _advance_genes, local_persistence), chromatin_studio/transcription.py.
Sources: README.md "Building and inspecting regulatory pathways", "Individual elongation and RNA molecules".
-->

# Regulation and transcription

::: {.box}
**Short version.** A gene has a promoter that flips between OFF and ON. While ON it starts polymerases at a set rate. Bound transcription factors, the coactivators they recruit, and enhancer–promoter contact in 3D speed up switching ON and starting. Marks at the promoter can allow or block activation. Polymerases travel along the gene, can pause, slow down at nucleosomes and evict them, and release an RNA that diffuses, is processed and decays.
:::

## Enhancers and promoters

A **regulatory element** is a genomic interval with a kind (enhancer or promoter), a target gene, a transcription factor (TF) species and, optionally, a coactivator species. For each element the engine counts, from the actual bonds at that moment:

- the number of distinct TF molecules bound directly to chromatin inside the element, and
- the number of distinct coactivator molecules linked directly to one of those TFs (one protein–protein hop; longer chains do not count).

These give a signal
$$s = n_\text{TF} + g_\text{coact}\, n_\text{recruited}.$$
An enhancer contributes only while its nearest bead is within `contact_radius` of the target promoter in 3D; a promoter element always contributes. Contributions to the same gene add up:
$$k_\text{ON} = k_\text{ON}^0 \times G_\text{mark} \times \Big(1 + \textstyle\sum g_\text{act}\, s\Big) \times F_\text{body},
\qquad
k_\text{init} = k_\text{init}^0 \times \Big(1 + \textstyle\sum g_\text{init}\, s\Big) \times F_\text{body}' .$$
$G_\text{mark}$ is 0 or 1: it is 0 if the gene's start bead lacks a `required_mark` or carries a `repressive_mark`. $F_\text{body}$ is the optional nuclear-body response: for each body, a factor $\exp(w \ln m)$ with $w = \exp(-\tfrac12 d^2/r^2)$ falling off with the promoter's distance $d$ from the body. A gene with no regulatory elements uses a simpler older rule: $k_\text{ON} = k_\text{ON}^0 (1 + \text{number of bound activators in the gene})$.

Losing contact or a mark removes future gain, but does not instantly switch off an ON promoter or stop a polymerase. The promoter therefore has memory. All these rules are **hypotheses** with experimental precedents; the literature does not support one universal contact rule.

![**A regulatory element at work** (Regulation view, illustrative 80 kb preset). The panel shows bound TFs, recruited coactivators, enhancer contact and the resulting rate terms for one gene, all read from the current state.](figures/screen-regulation.png){#fig:regulation width=72%}

## Promoter bursting

With rates frozen over one step of length $\Delta t$, the engine draws the **exact** ON/OFF path inside that step from exponential waiting times. Initiations are a Poisson process during the ON periods only, at times drawn inside those periods. This is the classic telegraph model of transcriptional bursting, and the study `telegraph_transcription` (Section 10) checks the ON fraction, the mean RNA and the Fano factor against its closed-form answers.

## Polymerases and RNA

There are two RNA modes.

- **Count mode** keeps per-gene counts of nascent, unprocessed and mature RNA.
- **Individual mode** gives every initiation a polymerase `pol-N` and transcript `rna-N`. The polymerase moves along the gene at `elongation_rate` beads per unit time, multiplied by `occupied_elongation_multiplier` on occupied nucleosome sites. It can pause and resume with exponential waiting times. While any polymerase is engaged, the gene's bending stiffness is scaled by `persistence_multiplier`, and the polymerase can evict the nucleosome under it.

At the end of the gene the RNA is released as a point particle that diffuses inside the nucleus. It is processed at `processing_rate`, increased by $(1 + \text{speckle gain})$ while it is inside a speckle, and decays at `rna_decay_rate`. The bookkeeping identity
$$\text{initiated} = \text{nascent} + \text{unprocessed} + \text{mature} + \text{degraded}$$
holds in every run. If the cap on tracked RNA (default 2,000, up to 5,000) is reached, new initiations are blocked and counted. Quantitative runs should have zero blocked initiations.

**Not represented:** a limited pool of free polymerase, polymerases blocking each other, backtracking, introns and splicing, DNA sequence, RNA chains with their own forces, and RNA binding to proteins.

## Does it work?

![**Bursting matches the telegraph model in both RNA modes.** Engine values divided by the closed-form answer for the ON fraction, mean RNA count and Fano factor (variance over mean), with count-based and individual-RNA bookkeeping. Grey boxes: accepted window. Study `telegraph_transcription`, rerun for this document (8 replicas × 3,000 samples).](figures/val-telegraph.svg){#fig:telegraph width=88%}

With switching rates 0.05 (ON) and 0.2 (OFF), initiation 2 and decay 0.02, the telegraph model gives an ON fraction of 0.2, a mean of 20 RNAs and a Fano factor of 6.93. The engine gives 0.199, 19.9 and 6.84 in count mode and 0.200, 20.2 and 7.07 in individual-RNA mode; all 6 checks pass within their windows (3% for the ON fraction, 5% for the mean and 10% for the Fano factor, or 3 standard errors). The Fano factor well above 1 is the signature of bursting.

Regulatory gains, mark gates and body responses are rules that multiply rates; they are covered by unit tests (`tests/test_regulation.py`) and by the rate breakdown in the interface, but not yet by an analysis with a figure.

## The full feedback loop

Putting Sections 3 to 7 together, the implemented chain of cause and effect is:

> marks and nucleosomes → TF capture → bound TFs and recruited coactivators, plus enhancer contact → promoter ON and initiation rates → polymerase movement → nucleosome eviction and gene stiffness → chromatin folding and contacts → (back to) encounters and marks.

Each link can be switched off on its own, which is the basis of the ablation studies in Section 10.
