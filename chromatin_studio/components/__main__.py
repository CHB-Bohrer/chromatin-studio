"""Command line for component studies.

    python -m chromatin_studio.components list
    python -m chromatin_studio.components run STUDY --output DIR [--scale standard]
        [--replicas N] [--jobs N] [--seed N] [--set key=json ...] [--case ID ...]
    python -m chromatin_studio.components analyze DIR
    python -m chromatin_studio.components show DIR

Output directories hold run data: keep them outside the repository.
"""

import argparse
import json
from pathlib import Path
import sys

from . import catalog
from . import runner


def _print_report(report):
    summary = report["check_summary"]
    print(f"{report['study']} [{report['scale']}]  results used: {report['results_used']}"
          f" (partial {report['partial_results']}, failed jobs {len(report['failed_jobs'])})")
    if report.get("missing_jobs"):
        print(f"  missing runs ({len(report['missing_jobs'])} of {report.get('expected_jobs')}): "
              + ", ".join(report["missing_jobs"]))
    if report.get("changed_since_run"):
        print("  source changed since the run: " + ", ".join(report["changed_since_run"]))
    for check in report.get("checks", []):
        mark = "PASS" if check["passed"] else "FAIL"
        observed = "n/a" if check["observed"] is None else f"{check['observed']:.5g}"
        if check.get("low") is not None or check.get("high") is not None:
            target = f"[{check.get('low')}, {check.get('high')}]"
        else:
            target = "n/a" if check["expected"] is None else f"{check['expected']:.5g}"
        sem = f" ± {check['sem']:.2g}" if check.get("sem") else ""
        print(f"  {mark}  {check['evidence']:<22} {check['label']}: {observed}{sem} vs {target} {check['unit']}")
    for line in report.get("findings", []):
        print("  note: " + line)
    print(f"  checks passed: {summary['passed']}/{summary['total']}")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m chromatin_studio.components", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list")
    run = commands.add_parser("run")
    run.add_argument("study")
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--scale", default="standard")
    run.add_argument("--replicas", type=int)
    run.add_argument("--jobs", type=int, default=1)
    run.add_argument("--seed", type=int)
    run.add_argument("--set", action="append", default=[], metavar="KEY=JSON")
    run.add_argument("--case", action="append", default=[])
    for name in ("analyze", "show"):
        commands.add_parser(name).add_argument("output", type=Path)
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "list":
            for record in catalog():
                print(f"{record['id']:<26} {record['title']}\n{'':<26} components: {', '.join(record['components'])};"
                      f" evidence: {', '.join(record['evidence'])}; scales: {', '.join(record['scales'])}")
            return 0
        if arguments.command == "run":
            overrides = {}
            for item in arguments.set:
                key, _, value = item.partition("=")
                try:
                    overrides[key] = json.loads(value)
                except json.JSONDecodeError:
                    overrides[key] = value
            report = runner.run(arguments.study, arguments.output, scale=arguments.scale,
                                replicas=arguments.replicas, jobs=arguments.jobs, seed=arguments.seed,
                                overrides=overrides, only_cases=arguments.case or None)
        elif arguments.command == "analyze":
            report = runner.analyze(arguments.output)
        else:
            report = json.loads((arguments.output / "report.json").read_text(encoding="utf-8"))
        _print_report(report)
        return 0 if not report["failed_jobs"] else 1
    except (ValueError, FileNotFoundError) as error:
        print(f"Component study failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
