<!--
Section: Chromatin as a polymer, and how particles move
Update when: any force term, the integrator, or the confinement changes.
Code: chromatin_studio/engine.py (_mechanics, _brownian_move, integration_substeps, local_persistence),
      chromatin_studio/chromatin_contacts.py, chromatin_studio/neighbors.py.
Sources: README.md "Dynamics and reaction semantics", "Protein preferences and chromatin attraction".
-->

# Chromatin as a polymer

::: {.box}
**Short version.** Chromatin is a chain of beads joined by springs. The chain resists bending, beads push each other away softly, beads with chosen marks can attract each other, and everything is kept inside a soft wall. All particles move by overdamped Brownian motion.
:::

## The chain

The chromatin is $N$ beads in a row. Neighbours are joined by a harmonic spring with stiffness $k_b$ (`bond_k`) and rest length $b$ (`bond_length`):
$$U_\text{bond} = \tfrac{1}{2} k_b \,(r_{i,i+1} - b)^2 .$$
The springs are soft on purpose. At 10 kb per bead a bead is a fluctuating blob, not a rigid unit, so a fixed bond length would not be physically motivated. In the equilibrated Fbn2 model, bonds vary by about 16% and only 0.5% stretch past 1.5 times their rest length.

## Bending stiffness

Each set of three neighbouring beads has a bending energy
$$U_\text{bend} = \kappa_i \,(1 - \cos\theta_i), \qquad \kappa_i = \frac{\ell_p \, m_i}{b},$$
where $\theta_i$ is the angle between consecutive bonds and $\ell_p$ is the setting `persistence_length`. The local factor $m_i$ lets other mechanisms change stiffness: a transcribed gene can be made stiffer or softer while a polymerase is on it, and occupied and vacant nucleosome sites can have different stiffness. With no such couplings $m_i = 1$.

One subtlety: on a chain of discrete beads, the persistence length you actually **measure** is shorter than the setting. The exact relation is
$$\ell_p^\text{realised} = \frac{-b}{\ln\!\big(\coth k - 1/k\big)}, \qquad k = \ell_p / b .$$
The polymer check in Section 10 confirms this.

## Soft repulsion

Beads that overlap push apart with a smooth, bounded energy
$$U_\text{rep} = \tfrac{\varepsilon}{2}\left(1 - \frac{r^2}{\sigma^2}\right)^2 \quad (r < \sigma),$$
where $\sigma$ is the sum of the two radii. Because the energy is finite at $r = 0$ (about 4 $k_BT$ in the presets), two chain segments can occasionally pass through each other. This is a known simplification: the model does not preserve knots or entanglements. Stronger repulsion would fix that if a question needs it.

## Mark-dependent attraction

You can add rules of the form "beads with mark A attract beads with mark B". Each rule is a smooth well of depth $\epsilon$ and range $R$:
$$U_\text{attr}(r) = -\epsilon \left(1 - \frac{r^2}{R^2}\right)^2 \quad (r < R).$$
Pairs closer than a minimum distance along the chain (default 3 beads) are skipped so that the rule acts between distant parts of the chain. If several rules apply to one pair, only the deepest well counts; wells do not add up. Because marks change during a run, the energy landscape changes too. This is a deliberate, non-equilibrium coupling and not a claim that a mark is itself a force.

## Confinement

All particles are kept inside a sphere of radius $R_c$ by a soft harmonic wall:
$$U_\text{wall} = \tfrac{1}{2} k_w \,(|\mathbf{x}| - R_c)^2 \quad (|\mathbf{x}| > R_c).$$
This is the chromosome territory. An optional larger nuclear envelope is described in Section 6.

![**The pair energies used by the model**, in reduced units with typical preset values. Soft repulsion is finite at contact, so chains can cross. The attraction well is optional and only acts between beads with chosen marks. The wall is zero inside and rises outside.](figures/fig-potentials.svg){#fig:potentials width=100%}

## How particles move

Every particle (chromatin, protein and released RNA) moves by the overdamped Langevin equation, integrated with the Euler–Maruyama method:
$$\mathbf{x}(t + \Delta t) = \mathbf{x}(t) + D\,\mathbf{F}\,\Delta t + \sqrt{2 D \Delta t}\;\boldsymbol{\xi}, \qquad \boldsymbol{\xi} \sim \mathcal{N}(0, 1).$$
Here $D$ is the particle's diffusion coefficient and $\mathbf{F}$ the total force. There is no inertia and no hydrodynamic coupling between particles.

To keep this stable, the engine looks at the stiffest force acting on each particle and, if needed, splits one step into several smaller sub-steps so that $D \cdot \text{stiffness} \cdot \Delta t_\text{sub} \le 0.35$. This is a safety net against blow-ups, not a guarantee of accuracy. If more than 1,000 sub-steps would be needed the engine refuses and asks for a smaller time step. If coordinates ever become non-finite or fly far outside the wall, the engine stops with an error instead of clipping them.

## Does it work?

Three analyses test the polymer and the integrator against exact answers.

**Brownian motion.** Free particles and a Rouse chain (beads joined by springs, nothing else) are started in equilibrium and their mean-square displacement (MSD) is compared with the exact discrete-time solution. All 14 checks pass, including the slow $t^{1/2}$ growth of a chain bead's MSD at intermediate times and the reduced motion of looped segments.

![**Particle and chain motion match the exact solution.** A free particle's MSD grows as $6Dt$; a bead inside a chain grows more slowly (sub-diffusion) until the whole chain moves; the chain's centre of mass diffuses $N$ times slower. Points: engine (192 chains of 64 beads). Dashed: exact. Source: `python -m chromatin_studio.benchmarks`, rerun for this document.](figures/val-rouse.svg){#fig:rouse width=70%}

**Ideal-chain statistics.** With bending and repulsion switched off, bond directions are independent, so the mean-square distance between two beads must grow exactly as $\langle R^2(s) \rangle = s \langle b^2 \rangle$, and the probability that two beads touch must fall as $s^{-3/2}$, the classic result for an ideal polymer. Both hold: the fitted contact exponent is −1.50 against the expected −3/2. Short separations sit about 1–2% below the line; this is the expected small time-step bias of the Euler–Maruyama method, which slightly inflates the shortest bond fluctuations, and it is inside the acceptance rule (5% or 3 standard errors).

![**An ideal chain behaves like an ideal chain.** Left: mean-square distance against separation along the chain. Right: contact probability, with the Gaussian prediction and the $s^{-3/2}$ slope. 4 replicas of a 64-bead chain, 200,000 steps each (`docs/paper/analysis/extra_checks.py`).](figures/val-ideal-chain.svg){#fig:ideal width=100%}

**Bending stiffness.** A chain of 48 beads of 25 bp, with the bending modulus set for 50 nm, should have bond directions that lose memory along the chain exactly as the discrete worm-like chain predicts. All 4 checks pass. The neighbouring-bond cosine is 0.827 against the exact 0.830 (0.3% low, probably because the bonds are slightly soft). The persistence length measured from the decay is 47.9 ± 3.6 nm, consistent with the discrete-chain prediction of 45.6 nm and inside DNA's 45–53 nm range, but too noisy to tell those two apart. The whole-chain mean-square end-to-end distance (541 ± 50 against 472) passes only through its large error, because it is the slowest quantity to average. This run used 1.5 million steps per replica instead of the study's standard 6 million; the longer run on the Windows machine gave the same picture (47.2 ± 1.7 nm). Note the setting is a bending modulus: 5.88 bonds gives 45.6 nm at this resolution, and 6.40 bonds would give 50 nm.

![**Bending stiffness gives the expected worm-like chain.** Left: average cosine between bonds against their distance along the chain, with the exact discrete prediction (dashed) and DNA's 45–53 nm persistence length (grey). Right: mean-square distance between beads against the exact discrete worm-like chain. Study `polymer_wlc`.](figures/val-wlc.svg){#fig:wlc width=100%}
