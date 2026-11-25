<!--
Section: Cohesin loop extrusion
Update when: cohesin loading, stepping, barriers, collisions or unloading change.
Code: chromatin_studio/engine.py (_advance_cohesin, _load_cohesin), components/cohesin_extrusion.py, components/fbn2_cohesin.py.
Sources: CONTRACT.md "cohesin", docs/PROJECT_STATUS.md logic review notes.
-->

# Cohesin and loop extrusion

::: {.box}
**Short version.** Each cohesin is a pair of "legs" on the chain, joined by a spring. Both legs step outward, growing a loop, until they hit a barrier (such as a CTCF site) or another cohesin. Cohesin falls off at a set rate and immediately reloads somewhere random.
:::

A cohesin starts on two adjacent free beads. In each step, each leg takes a Poisson-distributed number of steps outward with mean `step_rate` $\times\,\Delta t$, so a loop with both legs free grows at $2 \times$ `step_rate` beads per unit time. A leg stops when:

- the bead it is on has a **barrier** facing its direction (`left`, `right` or `both`), as for a convergent CTCF site;
- the next bead is held by another cohesin leg;
- it reaches the end of the chain.

The two legs are connected by a spring (`spring_k`) that pulls the two loop anchors together in 3D, so extrusion compacts the chain mechanically. Each cohesin unloads with hazard `off_rate`; the mean residence is $1/k_\text{off}$.

**Choices worth knowing.**

- Barriers are perfect walls: they never let a leg through and never release a stalled leg. Real CTCF sites are leaky.
- An unloaded cohesin reloads at once, so the number of bound cohesins is constant. There is no separate loading rate, because the model has none to set it from.
- Cohesin affects transcription only through 3D contacts. There is no direct "cohesin makes the gene fire faster" rule.


## Does it work?

![**Cohesin follows its configured rates.** Each point is the engine's value divided by the closed-form expectation (1 = exact); grey boxes show the accepted window. The two physical-unit rows are input-matched: the rates were set from measured values (0.5 kb/s, 22 min), so they confirm the unit conversion only. Study `cohesin_extrusion`, rerun for this document (8 replicas × 1,500 residences on 2,000 beads of 10 kb).](figures/val-cohesin.svg){#fig:cohesin width=88%}

The loop at unloading is slightly smaller than predicted (65.3 against 67 beads) because residences that hit the chain ends are excluded; the check allows for this. In 8 replicas no leg ever crossed a blocking barrier.
