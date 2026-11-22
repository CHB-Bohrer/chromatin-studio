"""Run and analyze component studies; knows nothing about any single study.

Layout of one output directory::

    plan.json                 study, scale, parameters, cases, seeds, provenance
    raw/<case>-r<N>.json      one finished (case, replica) result
    raw/<case>-r<N>.partial.json   latest in-progress observables of a long run
    report.json               analysis of whatever raw results exist

Runs are fresh starts with recorded seeds. Finished raw files are never
overwritten: rerunning the same command resumes by skipping them. ``analyze``
can be repeated at any time and marks results that are still partial.
"""

import copy
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import sys
import time
import traceback

from .. import __version__
from . import REGISTRY, load
from .core import summarize


def _utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _write(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, allow_nan=False, indent=1), encoding="utf-8")
    os.replace(temporary, path)


def _source_hashes(study_id):
    package = Path(__file__).parent
    files = [package / "core.py", package / "runner.py", package / (REGISTRY[study_id] + ".py"),
             package.parent / "engine.py", package.parent / "transcription.py", package.parent / "schema.py"]
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files}


def parameters(study_id, scale="standard", overrides=None):
    study = load(study_id).STUDY
    if scale not in study["scales"]:
        raise ValueError(f"{study_id} has no scale '{scale}'. Available: {', '.join(sorted(study['scales']))}")
    params = copy.deepcopy(study["scales"][scale])
    for key, value in (overrides or {}).items():
        if key not in params:
            raise ValueError(f"{study_id}: unknown parameter '{key}'. Known: {', '.join(sorted(params))}")
        params[key] = value
    return params


def _cases(module, params):
    cases = module.cases(params) if hasattr(module, "cases") else [{"id": "default"}]
    identifiers = [case["id"] for case in cases]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Case ids must be unique.")
    return cases


def _job(arguments):
    study_id, params, case, replica, seed, output = arguments
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    raw_dir = Path(output) / "raw"
    stem = f"{case['id']}-r{replica}"
    final, partial = raw_dir / f"{stem}.json", raw_dir / f"{stem}.partial.json"
    if final.exists():
        return {"case": case["id"], "replica": replica, "status": "kept"}
    module = load(study_id)
    started = time.monotonic()
    envelope = {"study": study_id, "case": case, "replica": replica, "seed": seed, "started_utc": _utc()}

    def save_partial(observables):
        _write(partial, dict(envelope, complete=False, elapsed_seconds=time.monotonic() - started,
                             observables=observables))

    try:
        observables = module.simulate(params, case, seed, save_partial)
    except Exception as error:  # recorded, never hidden: analysis lists failed jobs
        _write(raw_dir / f"{stem}.failed.json", dict(envelope, error=repr(error), traceback=traceback.format_exc()))
        return {"case": case["id"], "replica": replica, "status": "failed", "error": repr(error)}
    _write(final, dict(envelope, complete=True, finished_utc=_utc(),
                       elapsed_seconds=time.monotonic() - started, observables=observables))
    if partial.exists():
        partial.unlink()
    return {"case": case["id"], "replica": replica, "status": "complete"}


def run(study_id, output, *, scale="standard", replicas=None, jobs=1, seed=None, overrides=None,
        only_cases=None):
    """Run every (case, replica) of one study, then analyze. Returns the report."""
    module = load(study_id)
    params = parameters(study_id, scale, overrides)
    if replicas is not None:
        params["replicas"] = int(replicas)
    if seed is not None:
        params["seed"] = int(seed)
    cases = _cases(module, params)
    if only_cases:
        unknown = set(only_cases) - {case["id"] for case in cases}
        if unknown:
            raise ValueError(f"Unknown case(s): {', '.join(sorted(unknown))}")
        cases = [case for case in cases if case["id"] in only_cases]
    output = Path(output)
    (output / "raw").mkdir(parents=True, exist_ok=True)
    plan_path = output / "plan.json"
    plan = {"format": 1, "study": study_id, "scale": scale, "parameters": params, "cases": cases,
            "application_version": __version__, "created_utc": _utc(),
            "source_sha256": _source_hashes(study_id),
            "environment": {"python": sys.version.split()[0], "platform": platform.platform()},
            "restart_checkpoint": False}
    if plan_path.exists():
        previous = json.loads(plan_path.read_text(encoding="utf-8"))
        if previous["study"] != study_id or previous["parameters"] != params:
            raise ValueError(f"{output} already holds a different plan; choose a new output directory.")
        if previous["source_sha256"] != plan["source_sha256"]:
            raise ValueError(f"{output} was produced by different source code; choose a new output directory.")
    else:
        _write(plan_path, plan)
    work = [(study_id, params, case, replica, int(params["seed"]) + 1000 * index + replica, str(output))
            for index, case in enumerate(_cases(module, params)) if case in cases
            for replica in range(int(params["replicas"]))]
    if jobs > 1 and len(work) > 1:
        with multiprocessing.get_context("spawn").Pool(min(int(jobs), len(work))) as pool:
            outcomes = pool.map(_job, work, chunksize=1)
    else:
        outcomes = [_job(item) for item in work]
    report = analyze(output)
    report["job_outcomes"] = outcomes
    return report


def _archive_previous_report(output, analysis_hashes):
    """Keep a report made by different analysis code instead of overwriting it.

    Reports from the same analysis source are simply superseded (they are
    reproducible from the same raw files); anything else moves to
    ``report-history/`` so it stays attributable to the code that made it.
    """
    current = output / "report.json"
    if not current.exists():
        return
    try:
        previous = json.loads(current.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        previous = {}
    if previous.get("analysis_source_sha256") == analysis_hashes:
        return
    history = output / "report-history"
    history.mkdir(exist_ok=True)
    stamp = str(previous.get("analyzed_utc", "unknown")).replace(":", "").replace("+", "Z")
    target = history / f"report-{stamp}.json"
    counter = 1
    while target.exists():
        counter += 1
        target = history / f"report-{stamp}-{counter}.json"
    os.replace(current, target)


def analyze(output):
    """(Re)build report.json from whatever raw results exist in ``output``.

    The report states which runs are complete, partial, failed or missing, the
    source that simulated them (``source_sha256``) and the source that analyzed
    them (``analysis_source_sha256``), and which files differ between the two.
    """
    output = Path(output)
    plan = json.loads((output / "plan.json").read_text(encoding="utf-8"))
    module = load(plan["study"])
    results, failed = [], []
    for path in sorted((output / "raw").glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        if path.name.endswith(".failed.json"):
            failed.append({"file": path.name, "error": record.get("error")})
        elif path.name.endswith(".partial.json"):
            if not (path.with_name(path.name.replace(".partial.json", ".json"))).exists():
                results.append(record)
        else:
            results.append(record)
    analysis = module.analyze(plan["parameters"], results)
    analysis.setdefault("figures", [])
    checks = analysis.get("checks", [])
    expected = [(case["id"], replica) for case in plan["cases"] for replica in range(int(plan["parameters"]["replicas"]))]
    finished = {(r["case"]["id"], r["replica"]) for r in results if r.get("complete")}
    started = {(r["case"]["id"], r["replica"]) for r in results}
    started |= {(case, replica) for case, replica in expected
                if any(item["file"] == f"{case}-r{replica}.failed.json" for item in failed)}
    current = _source_hashes(plan["study"])
    changed = sorted(name for name in set(current) | set(plan["source_sha256"])
                     if current.get(name) != plan["source_sha256"].get(name))
    report = {"format": 2, "study": plan["study"], "title": module.STUDY["title"],
              "components": module.STUDY["components"], "question": module.STUDY["question"],
              "scale": plan["scale"], "analyzed_utc": _utc(), "plan_created_utc": plan.get("created_utc"),
              "application_version": plan["application_version"],
              "source_sha256": plan["source_sha256"], "analysis_source_sha256": current,
              "source_matches_current": not changed, "changed_since_run": changed,
              "results_used": len(results),
              "partial_results": sum(not r.get("complete") for r in results),
              "expected_jobs": len(expected), "complete_jobs": len(finished & set(expected)),
              "missing_jobs": [f"{case}-r{replica}" for case, replica in expected if (case, replica) not in started],
              "failed_jobs": failed, "check_summary": summarize(checks), **analysis,
              "sources": module.STUDY["sources"],
              "scope": ("Component study in reduced units. Experimental comparisons depend on the unit "
                        "mapping stated in each check and do not constitute a calibration.")}
    _archive_previous_report(output, current)
    _write(output / "report.json", report)
    return report
