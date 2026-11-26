<!--
Section: Proteins, binding and assemblies
Update when: binding rules, affinity/release multipliers, valence, IDR representation, or binding_capacity change.
Code: chromatin_studio/engine.py (_advance_binding, _chromatin_affinity, _bond_off_rate, _mark_release_factors),
      chromatin_studio/kinetics.py, schema.py (proteins, interactions).
Sources: CONTRACT.md "Config" (proteins, interactions), README.md "Protein interactions and DNA attachment".
-->

# Proteins: diffusion, binding and assemblies

::: {.box}
**Short version.** Each protein molecule is a diffusing bead (or a short flexible chain). When it comes within reach of a chromatin bead, it can bind with a set hazard; the bond is a spring and breaks with another hazard. Marks, nucleosomes and regulatory regions scale these hazards up or down. Molecules can also bind each other, up to a limit on how many partners each one can have.
:::

## Species and molecules

You define protein **species** (up to 32), each with a copy number, a diffusion coefficient, a radius and a **valence** (the maximum number of bonds one molecule can make, counting both chromatin and other proteins). A species can be one bead or an **IDR chain**: a few domain beads joined by springs, of which a chosen fraction are "stickers" that can bind. The IDR is a coarse cartoon of a disordered region, not a residue-level sequence. The diffusion coefficient you enter is the molecule's free centre-of-mass value; for an $n$-bead chain each bead gets $n$ times that value so the whole chain diffuses at the requested rate.

## Binding to chromatin

At each step, for every protein site that still has free valence, the engine finds chromatin beads within the capture distance (protein radius plus bond length). For each such pair it draws a waiting time from an exponential distribution with rate
$$k_\text{on}^\text{eff} = k_\text{on} \times a_\text{mark} \times a_\text{region} \times a_\text{nuc}.$$

- $k_\text{on}$ is the species' `on_rate`, a **hazard per encounter**, not a molar association constant.
- $a_\text{mark}$ is the species' affinity for the bead's marks. If the bead carries several mapped marks, the largest value is used; with no mapped mark the `unmodified` value is used.
- $a_\text{region}$ is an optional extra preference inside enhancers or promoters (Section 7).
- $a_\text{nuc}$ is an optional factor for occupied versus vacant nucleosome sites (Section 4).

All candidate bindings in a step are sorted by their waiting times and accepted in that order, as long as valence and site limits still allow. This ordering by exponential clocks is what stops two molecules from both taking the last free slot. An accepted bond is a spring of stiffness 12 between protein and bead.

**Optional bead capacity.** By default a chromatin bead can hold any number of bonds. Setting `binding_capacity` (1 to 64) caps the total across species. The cap is checked only when a capture is accepted, so it uses no extra random numbers, and a cap that is never reached gives exactly the same run as no cap. The nucleosome studies in Section 10 use a cap of 1 so that factors compete for a site as they would in a single-site experiment.

## Release

Each existing chromatin bond breaks with probability $1 - e^{-k_\text{off}^\text{eff}\Delta t}$ per step, where
$$k_\text{off}^\text{eff} = k_\text{off} \times r_\text{nuc} \times r_\text{mark}.$$
$r_\text{nuc}$ lets a nucleosome speed up or slow down release (allowed range 0 to 1000, because measured nucleosome effects on transcription-factor release reach about 1000-fold). $r_\text{mark}$ lets a mark change release: if several listed marks are present the smallest factor is used, mirroring the "strongest affinity wins" rule for association. Both are off (equal to 1) unless set.

## Protein–protein links

Pairs of species can be given an interaction rule in one of two modes.

**Bonds** (the default): a capture radius, an on-hazard, an off-hazard and a spring strength. Links obey valence, so molecules can form chains and networks larger than dimers. The interface reports the largest connected cluster.

**Attraction**: no discrete bonds. Sticker beads of the two species feel a smooth pair well. Its depth $\varepsilon$ (in $k_BT$) applies at contact, $r \le \sigma = r_a + r_b$, and it fades to zero at the rule's radius $r_c$:

$$U(r) = -\varepsilon\,(1 - s^2)^2, \qquad s = \frac{r - \sigma}{r_c - \sigma} \text{ clipped to } [0, 1].$$

Inside contact the well is flat, so it exerts no force there. Attracting pairs also get a steep core inside contact, $4[(w/r)^{12} - (w/r)^6] + 1$ with $w = \sigma/2^{1/6}$. Without it, the bounded soft repulsion let many attracting beads collapse onto one another, giving clusters far denser than beads can pack **[Assumption: the core shape is the standard WCA form, not fitted]**. How many partners fit is then set by geometry, not by a valence number. This mode is meant for coarse phase separation with one bead per protein. It runs on the GPU when it is the only protein activity (Section 8).

A large cluster alone is **not** evidence of liquid–liquid phase separation. The tests are: a dense phase coexisting with a dilute one, a threshold concentration below which no droplet forms, and molecules exchanging between the two.

**Does the attraction mode phase-separate?** The **Phase separation · 80 kb** preset (400 × 200 bp chain, 600 one-bead scaffolds, 150 clients) was run on the GPU, with one seed per case **[Unverified: single seeds, not converged to one droplet at high concentration]**. At 600 scaffolds, a 1.0 $k_BT$ contact depth gave no droplet. 1.5 $k_BT$ gave one droplet (density about 4 per unit volume, below close packing of about 6.5), with about 13% of its members replaced every 20,000 steps. 2.0 $k_BT$ gave two droplets that did not merge, and 3.0 $k_BT$ several small arrested clusters. At 1.5 $k_BT$, no droplet formed at scaffold densities of 0.014 to 0.055. Droplets formed at 0.083 and 0.138, where the dilute phase stayed at 0.050 and 0.038, so the dilute concentration did not rise with the total amount. The preset uses 1.5 $k_BT$. The well depth is a model setting, not a measured affinity.

## What the numbers mean

Because $k_\text{on}$ is a hazard per close encounter, the bound fraction you see depends on diffusion, crowding, confinement and rebinding, not only on $k_\text{on}/k_\text{off}$. A separate tool (`kinetics.py`) computes an idealised diffusion-limited association constant for comparison, but it does not turn the simulation's hazards into measured $K_d$ values. Section 10 shows that matching a measured bound fraction usually needs a scan of $k_\text{on}$.

## Does it work?

**Free diffusion.** One non-binding species, as single beads and as 4-bead IDR chains, should spread with $\text{MSD} = 6Dt$ at the configured diffusion coefficient. All 6 checks pass. At the longest lag the MSD is about 1% below $6Dt$ for beads (4.75 against 4.80) and 1.5% below for IDR chains (4.73), more than the standard error. This small shortfall was also seen on the Windows machine and is probably soft collisions and the wall; it is not yet explained **[Unverified]**.

![**Free proteins diffuse at the configured rate.** Mean-square displacement against lag time for bead and IDR forms, against $6Dt$. Study `protein_diffusion`.](figures/val-diffusion.svg){#fig:diffusion width=100%}

**Binding and release.** To test the binding rule by itself, proteins are held at fixed positions above a straight chain, each within reach of exactly 3 beads, and only the reactions are advanced. Each molecule then binds with total hazard $K = 3 k_\text{on}$ and releases with $k_\text{off}$, so the bound fraction must rise as $\frac{K}{K + k_\text{off}}\left(1 - e^{-(K + k_\text{off})t}\right)$. It does, within the rule of 0.01 or 3 standard errors at every checked time. The two earliest points sit about 2 standard errors low (about 4%); this has not been explained and is worth a closer look **[Unverified]**.

![**Bound fraction follows the two-state solution.** 32 replicas × 500 molecules, $k_\text{on} = 0.05$ per bead, $k_\text{off} = 0.2$, frozen coordinates (`extra_checks.py`).](figures/val-binding.svg){#fig:binding width=62%}

How binding behaves when molecules really search in 3D, with rebinding, is tested against experiments in Section 10.
