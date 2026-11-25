<!--
Section: Interface and workflow
Update when: views, tabs, the build/apply/run workflow, studies, Evidence tab, cluster export or the tutorial change.
Code: static/*.js, chromatin_studio/server.py, chromatin_studio/studies.py, chromatin_studio/components/evidence.py.
Sources: README.md "Using the workbench", docs/COMPONENT_STUDIES.md "Figures and the Evidence view", tutorial/README.md.
-->

# Using the software

::: {.box}
**Short version.** Start the local server, open the page, pick a preset or build a model, press Apply, then Run. Watch it in 3D and along the genome. Run ablation studies, export a cluster job, or open the Evidence tab to see every saved validation result.
:::

## The workflow

1. **Start.** `python run.py --open` starts a server on your own computer (nothing is uploaded) and opens the app at `http://127.0.0.1:8765`.
2. **Build.** The left-hand builder has one editor per part of the model: System, Proteins, Interactions, Chromatin, Genes, Regulation, Bodies (including the nucleus and tethers), Cohesin, Studies, Calibrate and Cluster. Every editor works on a **draft**.
3. **Apply.** Apply validates the draft and starts a fresh, seeded run. Edits never change a running simulation silently.
4. **Run and inspect.** The 3D view shows the chain, proteins, bonds, envelope and tethers. The 1D genome browser shows marks, nucleosomes, cohesin, tethers and every bound molecule along the genome, with a log of real binding, release and remodeling events. The Regulation view breaks each gene's rates into their terms.
5. **Compare.** The Studies editor runs a control and chosen ablations (for example `no_cohesin`, `no_tethering`, `no_nucleosome_access`) as independent copies, without touching the live run.
6. **Scale up.** Cluster export packages the draft as a job bundle for many replicas on a cluster.
7. **Save.** Export the model JSON to reuse it, or the experiment with its current snapshot.

![**The live 1D genome browser** (25.6 kb nucleosome-access preset). Filled circles are occupied nucleosome sites, dashed circles vacant ones; coloured dots below are individual bound molecules, labelled with their simulation IDs.](figures/screen-live-genome.png){#fig:genome width=100%}

## The Evidence tab

Component studies (Section 10) save a `report.json` with every check and the data for its figures. Start the server with `--evidence <folder>` and the Evidence tab draws all of them: overview tiles, charts with hover values and a table view, and every check with its acceptance window, grouped by evidence level and tagged as input-matched, tuned or prediction. The tab only reads files; it never runs anything.

![**The Evidence tab**, showing the Donovan et al. nucleosome study: measured and simulated association and dissociation ratios (left), the affinity ratio test (right), and residence-time curves (bottom).](figures/screen-evidence.png){#fig:evidence width=100%}

## Built-in teaching material

**How it works** and **Learn** open illustrated lessons on motion, encounters, accessibility, regulation and validation. A narrated tutorial video (18 chapters, about 12 minutes) walks through real captured states of the app. Both are separate from the simulation and never change it.
