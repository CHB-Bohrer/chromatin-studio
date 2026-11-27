# Status

**Project:** Chromatin Studio. A local visual builder and simulator for chromatin, reacting proteins, histone marks, cohesin and transcription.
**Updated:** 2026-11-27 (project state as last reviewed 2026-11-24, application 0.10.0)

This is the short front page. The detailed record is [docs/PROJECT_STATUS.md](docs/PROJECT_STATUS.md); read it before substantive work.

## Where things stand
- Working research prototype in reduced units. No accepted empirical length/time fit yet.
- Component studies: `nucleosome_tf_kinetics` 26 of 28 checks pass, `tf_binding` (Sox2) 13 of 15; `hp1_mark_density` was redesigned after a superseded first run; `fbn2_cohesin` has not run.
- Software paper in `docs/paper/` (sections mapped to code, rebuilt with `python docs/paper/build.py`), Evidence view in the app, narrated tutorial in `tutorial/`.
- Not done: production calibration, full-model GPU validation, a 20 Mb nucleosome-resolution run, a complete published-figure reproduction.

## Last session
- 2026-11-27: added this file, `THOUGHTS.md`, `tools/thoughts.py`, `tools/next_steps.py` and the "Owner's standing preferences" section in `AGENTS.md`. No code or results changed.

## Next
- [ask] Add one HP1α hazard near 0.8 and a Sox2 hazard near 0.235 (a few runs each) to confirm the interpolations.
- [ask] Test a facilitated-dissociation rule for HP1α competition as a separate, optional mechanism.
- [ask] Run `fbn2_cohesin` overnight.
- [ask] Run the Node interface tests on a machine with Node.js.

## Waiting on the owner
- Go-ahead for the runs above (run data lives outside the repository; see `LOCAL_DATA.md` on the owner's machine).

## Where to look
- [AGENTS.md](AGENTS.md): rules and scientific distinctions. [CONTRACT.md](CONTRACT.md): configuration, API and reaction semantics.
- [docs/PROJECT_STATUS.md](docs/PROJECT_STATUS.md), [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md), [docs/COMPONENT_STUDIES.md](docs/COMPONENT_STUDIES.md).
- `chromatin_studio/`: the engine and server. `docs/paper/`: the software paper. `tutorial/`: the guided tour.
