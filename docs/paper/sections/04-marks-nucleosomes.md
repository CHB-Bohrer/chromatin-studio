<!--
Section: Histone marks, nucleosomes and remodelers
Update when: mark writing/erasing/spreading, nucleosome turnover, sliding, or remodeler recruitment changes.
Code: chromatin_studio/engine.py (_advance_marks, _advance_nucleosomes, _advance_nucleosome_sliding,
      _nucleosome_rates, _remodeling_targets, _clear_vacant_marks).
Sources: CONTRACT.md (marks, nucleosomes, remodelers), README.md "Nucleosome occupancy, TF access and mechanics".
-->

# Histone marks, nucleosomes and remodelers

::: {.box}
**Short version.** Each chromatin bead carries on/off flags for histone marks and, optionally, an occupied/vacant flag for its nucleosome. Bound proteins write or erase marks; marks can also spread to nearby beads. Nucleosomes can be evicted and reassembled, and remodelers can slide them to a neighbouring empty site. Losing a nucleosome wipes its marks.
:::

## Marks

A **mark** is a yes/no flag on each bead (for example H3K27ac or H3K9me3). Marks are independent; nothing makes two marks mutually exclusive unless you build that in. Each mark can change by three routes:

- **Writing** by a bound protein whose species has `writer_mark` set, at its `catalytic_rate`. A coactivator recruited by a bound transcription factor writes at that factor's bead (Section 7).
- **Erasing** by a bound protein with `eraser_mark` set, plus a background `decay_rate`.
- **Spreading** (optional): an unmarked bead gains the mark at rate `spread_rate` times the number of marked beads within `contact_radius` in 3D. This is a phenomenological rule, set to zero to study only enzyme-driven marking.

For each bead, the total "on" hazard $w$ and "off" hazard $e$ are added up, and the flag is updated with the exact two-state probability over the step:
$$P(\text{flip}) = \frac{h}{w + e}\left(1 - e^{-(w+e)\Delta t}\right),$$
where $h = e$ if the bead is marked and $h = w$ if it is not.

## Nucleosome occupancy

When enabled, each bead is a **nucleosome site** that is either occupied or vacant. This is a coarse state on a fixed piece of DNA: there is no wrapped DNA geometry, no histone identity, and an implicit unlimited pool of free histones. Occupancy matters in four ways:

1. Proteins can prefer or avoid occupied sites, both for binding and for release (Section 3).
2. Occupied and vacant sites can have different bending stiffness.
3. Marks can only exist on occupied sites. Eviction clears them, and a reassembled nucleosome starts unmarked.
4. Polymerases can move more slowly through occupied sites and can evict them (Section 7).

Each site switches with an eviction hazard $E$ (background plus any bound remodeler or passing polymerase) and a reassembly hazard $A$. The engine samples the exact two-state endpoint over each step. It also tracks whether an occupied site was evicted **and** refilled within the same step, so the marks are still wiped in that case.

![**Background nucleosome turnover against the exact two-state answer.** Starting from all sites occupied, the occupied fraction relaxes to $A/(E+A)$ with rate $E + A$. Points: engine (8 replicas × 4,000 sites, eviction 0.02, reassembly 0.06). Dashed: closed form. All 5 checked times agree within 0.4%. Study `nucleosome_turnover`, rerun for this document.](figures/val-turnover.svg){#fig:turnover width=78%}

## Does it work?

Nucleosome turnover is checked in the figure above: the occupied fraction follows the exact two-state curve at every checked time. Marks are checked the same way: with only background decay at rate $d$, the marked fraction must fall as $e^{-dt}$, and it does.

![**Mark decay follows the exact exponential.** 8 replicas × 2,000 beads, decay rate 0.05, frozen coordinates (`extra_checks.py`).](figures/val-marks.svg){#fig:marks width=62%}

Writing by bound enzymes, spreading and sliding are covered by unit tests (`tests/test_engine.py`, `tests/test_live_nucleosomes.py`) but not yet by a quantitative analysis with a figure.

## Remodelers

A protein species can carry three remodeling rates, all zero by default:

- `nucleosome_eviction_rate`: adds to $E$ at the bead it is bound to.
- `nucleosome_reassembly_rate`: adds to $A$ there.
- `nucleosome_sliding_rate`: moves an occupied nucleosome at or next to its bead one site into an adjacent vacant site. The marks move with the nucleosome, and the total number of nucleosomes stays the same.

A remodeler acts either because it is bound to chromatin directly or because it is recruited by a bound transcription factor. Each molecule counts once per bead. Sliding events are drawn one by one with exponential waiting times within each step. There is no ATP cycle and no motor mechanics; these are effective rates.
