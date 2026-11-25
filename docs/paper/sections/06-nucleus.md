<!--
Section: Nuclear envelope, tethering elements and nuclear bodies
Update when: the envelope law, tethers (attach/detach rules, anchors), or bodies (speckles, nucleoli) change.
Code: chromatin_studio/nucleus.py, chromatin_studio/engine.py (_advance_tethers, _tether_hazard, _tether_geometry,
      body wells in _mechanics, _gene_regulatory_rates body_responses).
Sources: CONTRACT.md "nucleus", "bodies"; docs/PROJECT_STATUS.md "Whole nucleus, tethering elements".
-->

# The nucleus: envelope, tethers and bodies

::: {.box}
**Short version.** The chromosome territory can sit inside a larger soft ellipsoid, the nuclear envelope. Chosen stretches of the genome (for example lamina-associated domains) can stick to the envelope or to a nuclear body when they come close, held by a spring, and let go later. Speckles and nucleoli are fixed attractive regions in space.
:::

## The envelope

The optional envelope is an ellipsoid with centre $\mathbf{c}$ and semi-axes $(a, b, c)$. For a point $\mathbf{x}$, let $\mathbf{u} = \mathbf{x} - \mathbf{c}$ and $\rho = |(\mathbf{x}-\mathbf{c})/(a,b,c)|$, so $\rho = 1$ on the surface. Outside, the energy is
$$U_\text{env} = \tfrac{k}{2}\,|\mathbf{u}|^2\left(1 - \frac{1}{\rho}\right)^2 \quad (\rho > 1),$$
which, for a sphere, reduces exactly to the territory wall law of Section 2. The force is its exact gradient. The envelope acts on chromatin, proteins and released RNA, in addition to the territory wall.

## Tethering elements

A tethering element is a genomic interval, say beads 80 to 249. Each bead in it can attach to an anchor surface, either the envelope or a chosen nuclear body:

- **Attach.** A bead that is within `capture_distance` of the surface attaches with hazard `attach_rate`, multiplied by a factor for its nucleosome state. If `required_mark` is set, beads without that mark cannot attach.
- **Hold.** On attachment the anchor point is fixed at the nearest surface point (radial projection). The bead is then pulled toward it by a harmonic spring (`spring_k`, `rest_length`).
- **Detach.** An attached bead lets go with hazard `detach_rate`, again multiplied by a nucleosome-state factor.

Each transition is logged as an event, so you can see attachments appear along the genome track. These are phenomenological anchors: there is no lamina meshwork, no nuclear pores and no membrane mechanics. Models with an enabled nucleus run only on the NumPy engine; the OpenMM backends refuse them.

## Does it work?

The envelope force is checked in unit tests against its energy (exact gradient) and against the sphere law (`tests/test_nucleus.py`). For the tether rules, a chain is placed inside a spherical envelope, coordinates are frozen, and only the reactions run. Half of the beads carry the required mark. Marked beads within capture distance should attach following the two-state curve $\frac{k_a}{k_a + k_d}\left(1 - e^{-(k_a + k_d)t}\right)$, and unmarked beads should never attach. Both hold.

![**Tethers attach and detach at the configured rates.** Attached fraction of eligible beads (marked, within capture distance), 8 replicas with 243–298 eligible beads each, $k_a = 0.3$, $k_d = 0.1$. No unmarked bead attached in any replicate (`extra_checks.py`).](figures/val-tethers.svg){#fig:tethers width=62%}

## Nuclear bodies

**Speckles** and **nucleoli** are fixed Gaussian wells in space, with a centre, a radius and separate strengths for chromatin and each protein species:
$$U_\text{body}(\mathbf{x}) = -\,\alpha\, \exp\!\left(-\frac{|\mathbf{x} - \mathbf{c}|^2}{2 r^2}\right).$$
A negative $\alpha$ makes a body repel instead. They are **imposed landmarks**, not condensates that form by themselves. Bodies can also change transcription (Section 7): speckles can speed up RNA processing, and a gene can be given optional activation and initiation multipliers that fade with its distance from a body.
