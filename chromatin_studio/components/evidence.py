"""Read-only access to saved component-study reports, for the interface.

Run data stays outside the repository; the server is pointed at a folder with
``--evidence DIR`` or ``CHROMATIN_STUDIO_EVIDENCE``. Only ``report.json`` files
(and their archived ``report-history/`` versions) below that folder are read.
Paths are returned relative to it, and any absolute filesystem path found
inside a report is reduced to its last component, so a capture of the
interface never shows where this machine keeps its data.
"""

import json
import os
import re
from pathlib import Path, PurePosixPath

from . import REGISTRY, catalog

REPORT_NAME = "report.json"
HISTORY = "report-history"
MAX_DEPTH = 4            # e.g. <date>/<study>/report-history/report-<stamp>.json
MAX_REPORTS = 500
MAX_BYTES = 16 * 1024 * 1024
_ABSOLUTE = re.compile(r"^(?:[A-Za-z]:[\\/]|\\\\|/(?:Users|home|mnt|c|d|Volumes|tmp|var)/)")
_SUMMARY_KEYS = ("study", "title", "components", "question", "scale", "analyzed_utc", "plan_created_utc",
                 "application_version", "source_matches_current", "changed_since_run", "results_used",
                 "partial_results", "expected_jobs", "complete_jobs", "check_summary", "format")


def configured_root(value=None):
    """The evidence folder from an explicit value or the environment, or None."""
    value = value or os.environ.get("CHROMATIN_STUDIO_EVIDENCE")
    return Path(value).expanduser().resolve() if value else None


def _hide_paths(value):
    """Copy of a JSON value with absolute filesystem paths reduced to a name."""
    if isinstance(value, dict):
        return {key: _hide_paths(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_hide_paths(item) for item in value]
    if isinstance(value, str) and _ABSOLUTE.match(value):
        return "…/" + re.split(r"[\\/]", value.rstrip("\\/"))[-1]
    return value


def _read(path):
    if path.stat().st_size > MAX_BYTES:
        raise ValueError("Report is larger than the 16 MB interface limit.")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or not isinstance(data.get("checks"), list) or "study" not in data:
        raise ValueError("Not a component-study report.")
    return data


def _candidates(root):
    for current, folders, files in os.walk(root):
        depth = len(Path(current).relative_to(root).parts)
        folders[:] = sorted(f for f in folders if not f.startswith(".") and depth < MAX_DEPTH)
        if REPORT_NAME in files and Path(current).name != HISTORY:
            yield Path(current) / REPORT_NAME


def index(root):
    """Summaries of every report below ``root``; never raises for one bad file."""
    studies = catalog()
    if root is None:
        return {"configured": False, "available": False, "reports": [], "catalog": studies,
                "hint": "Start the server with --evidence DIR or set CHROMATIN_STUDIO_EVIDENCE to the folder "
                        "that holds component-study outputs (for example the run-data folder in LOCAL_DATA.md)."}
    if not root.is_dir():
        return {"configured": True, "available": False, "reports": [], "catalog": studies,
                "hint": "The configured evidence folder does not exist on this machine."}
    reports, skipped = [], 0
    for path in _candidates(root):
        if len(reports) >= MAX_REPORTS:
            skipped += 1
            continue
        relative = path.relative_to(root).as_posix()
        try:
            data = _read(path)
        except (OSError, ValueError) as exc:
            reports.append({"path": relative, "error": str(exc)})
            continue
        record = {key: data.get(key) for key in _SUMMARY_KEYS}
        record["format"] = data.get("format", 1)
        record["figures"] = len(data.get("figures") or [])
        record["missing_jobs"] = len(data.get("missing_jobs") or [])
        record["failed_jobs"] = len(data.get("failed_jobs") or [])
        record["registered"] = data.get("study") in REGISTRY
        history = path.parent / HISTORY
        record["history"] = sorted(p.relative_to(root).as_posix() for p in history.glob("report-*.json")) \
            if history.is_dir() else []
        record["path"] = relative
        reports.append(_hide_paths(record))
    reports.sort(key=lambda r: (r.get("study") or "", r.get("analyzed_utc") or ""), reverse=False)
    return {"configured": True, "available": True, "reports": reports, "catalog": studies,
            "skipped": skipped}


def report(root, relative):
    """One full report by its path relative to ``root``."""
    if root is None:
        raise ValueError("No evidence folder is configured.")
    if not isinstance(relative, str) or not relative or "\\" in relative or "\x00" in relative:
        raise ValueError("Supply one report path relative to the evidence folder.")
    parts = PurePosixPath(relative).parts
    if PurePosixPath(relative).is_absolute() or ".." in parts or len(parts) > MAX_DEPTH + 1:
        raise ValueError("Report paths must stay inside the evidence folder.")
    name_ok = parts[-1] == REPORT_NAME or (len(parts) > 1 and parts[-2] == HISTORY
                                           and re.fullmatch(r"report-[0-9A-Za-z_.:-]+\.json", parts[-1]))
    if not name_ok:
        raise ValueError("Only report.json files and their archived versions can be opened.")
    target = (root / Path(*parts)).resolve()
    if not target.is_relative_to(root) or not target.is_file():
        raise FileNotFoundError("Report not found.")
    data = _hide_paths(_read(target))
    data["path"] = relative
    return data
