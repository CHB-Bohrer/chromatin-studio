"""Copy the small, figure-relevant parts of validation outputs into docs/paper/data/.

Run data (raw trajectories) stays outside the repository. This script keeps only
what the paper draws: checks, figure series, summaries and provenance. It then
regenerates data/checks-table.md, the summary table included in Section 10.

    python docs/paper/analysis/collect.py --runs /path/outside/repo/val

Expected layout under --runs (any subset is fine; missing items keep their old data file):
    <study>/report.json          component studies (python -m chromatin_studio.components run ...)
    reference-results.json       python -m chromatin_studio.benchmarks --output ...
    extra.json                   docs/paper/analysis/extra_checks.py --output ...
"""

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
STUDIES = ["polymer_wlc", "protein_diffusion", "nucleosome_turnover", "telegraph_transcription",
           "cohesin_extrusion"]
KEEP = ["study", "title", "question", "scale", "analyzed_utc", "application_version", "source_matches_current",
        "changed_since_run", "complete_jobs", "expected_jobs", "check_summary", "checks", "findings", "figures",
        "summary"]

# Section each mechanism belongs to, and a plain-language name for the table.
TABLE = [
    ("benchmark", "Free particles and Rouse chain MSD", "2"),
    ("extra:ideal_chain", "Ideal chain: distances and contact scaling", "2"),
    ("polymer_wlc", "Worm-like chain bending", "2"),
    ("protein_diffusion", "Free protein diffusion (bead and IDR)", "3"),
    ("extra:binding", "Protein binding and release, two-state", "3"),
    ("nucleosome_turnover", "Nucleosome turnover, two-state", "4"),
    ("extra:marks", "Mark decay", "4"),
    ("cohesin_extrusion", "Cohesin extrusion speed, residence, barriers", "5"),
    ("extra:tethers", "Envelope tethers, two-state and mark gate", "6"),
    ("telegraph_transcription", "Promoter bursting (telegraph model)", "7"),
]


def load(path):
    return json.loads(Path(path).read_text()) if Path(path).exists() else None


def write(name, value):
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / name).write_text(json.dumps(value, indent=1))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", required=True)
    args = parser.parse_args()
    runs = Path(args.runs)
    for study in STUDIES:
        report = load(runs / study / "report.json")
        if report:
            plan = load(runs / study / "plan.json") or {}
            kept = {k: report.get(k) for k in KEEP}
            kept["parameters"] = plan.get("parameters")
            kept["environment"] = plan.get("environment")
            write(f"component-{study}.json", kept)
            print("collected", study)
    bench = load(runs / "reference-results.json")
    if bench:
        write("benchmark-reference.json", {k: bench[k] for k in ("status", "title", "scope", "created_at",
                                                                 "provenance", "parameters", "checks", "curves")})
        print("collected benchmark")
    extra = load(runs / "extra.json")
    if extra:
        old = load(DATA / "extra-checks.json") or {}
        old.update(extra)
        write("extra-checks.json", old)
        print("collected extra checks")
    write_table()


def counts(key):
    if key == "benchmark":
        d = load(DATA / "benchmark-reference.json")
        if not d:
            return None
        return sum(c["passed"] for c in d["checks"]), len(d["checks"]), d["created_at"][:10]
    if key.startswith("extra:"):
        d = load(DATA / "extra-checks.json")
        if not d or key[6:] not in d:
            return None
        checks = d[key[6:]]["checks"]
        return sum(c["passed"] for c in checks), len(checks), d.get("provenance", {}).get("created", "")
    d = load(DATA / f"component-{key}.json")
    if not d:
        return None
    s = d["check_summary"]
    return s["passed"], s["total"], (d.get("analyzed_utc") or "")[:10]


def write_table():
    rows = ["| Mechanism | Section | Checks passed | Run date |", "|---|---|---|---|"]
    for key, label, section in TABLE:
        c = counts(key)
        if c is None:
            rows.append(f"| {label} | {section} | not run | |")
        else:
            rows.append(f"| {label} | {section} | {c[0]} of {c[1]} | {c[2]} |")
    (DATA / "checks-table.md").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
