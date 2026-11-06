# Public source release

Version **0.10.0**, prepared **2026-11-06** for [CHB-Bohrer/chromatin-studio](https://github.com/CHB-Bohrer/chromatin-studio). [Local release checks](../benchmarks/public-release-verification.json) passed. The initial GitHub upload is pending authenticated write access; the remote currently contains no files.

## Included scope

The public tool includes the Python engine, browser interface, model examples, tests, numerical/reference reports, AI development handoff, current tutorial and editable narration/capture sources. Hypothesis examples and short study summaries are included deliberately. The complete source manuscript and any full-text extraction are excluded; the tool does not need them to run.

Machine environments, credentials, caches, private continuation notes and unused backup exports are excluded. The historical `previous-run-before-cluster.json` is retained because the calibration protocol identifies it as its source model; it is a simulated model/run export, not a manuscript attachment. The linked calibration ZIP is an archived v0.8.0 proposed workload and has not been silently regenerated as v0.10.

## Software license

A software license is awaiting the maintainer's selection. Public availability alone should not be described as an open-source license grant. This section and citation metadata will be updated when the license is selected. Third-party data and dependencies retain their existing terms.

## Data and dependency attribution

The public Hansen lab reference data and derived Fbn2 curves retain the attribution and source CC-BY-4.0 records in [licenses.json](../references/hansen/licenses.json). Consult [source methods](../references/hansen/source-methods.json), [inventory](../references/hansen/inventory.json) and the original linked datasets when reusing them. A software license must not replace those data licenses.

NumPy, SciPy, polychrom and OpenMM are installed as dependencies and retain their own licenses. Their source is not vendored into the application. Optional polychrom/OpenMM versions are specified in [requirements-polychrom.txt](../requirements-polychrom.txt). Voice model weights are not included; the optional tutorial narration generator retrieves separately documented, pinned assets. Existing narration and screenshots accompany the editable tutorial sources.

Use [CITATION.cff](../CITATION.cff) for software citation. This release has no assigned DOI. Cite the original experimental/modeling studies separately where appropriate.

## Tutorial portability

The published `tutorial/render-manifest.json` and `tutorial/render-timeline.json` store narration paths relative to the `tutorial` directory. The original renderer wrote machine-specific absolute paths; the public release changes only those path strings and adds explicit path-base metadata. Audio, captures, timings, video and source checksums are preserved. [Publication metadata verification](../tutorial/publication-metadata.json) records the original/published file hashes and verifies all referenced WAVs. Earlier rendering reports remain historical evidence of the unchanged media, not a claim that their original build metadata was byte-identical after path normalization.

Future tutorial rebuilds may emit absolute build paths again. Normalize those metadata paths before committing generated records, and verify the referenced files rather than changing media hashes. Follow [tutorial/README.md](../tutorial/README.md) for narration and rendering.

## Installation and verification scope

The first public distribution supports running from a cloned or downloaded source directory. Standalone wheel packaging is unfinished. [CI](../.github/workflows/tests.yml) runs dependency-light Python CPU tests and JavaScript helper tests. Optional OpenMM/GPU tests may skip; CI success does not establish GPU execution, physical calibration or biological validity.

Current model and validation limits are recorded in [PROJECT_STATUS.md](PROJECT_STATUS.md). A 20 Mb coarse model is supported; a 20 Mb model at 200 bp per nucleosome is not. The app is a local research prototype, not a shared public simulation service.
