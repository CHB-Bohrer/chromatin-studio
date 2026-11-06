"""Loopback-only local UI and simulation service; no cloud service required."""

import argparse
from collections import deque
import copy
import importlib.metadata
import importlib.util
import json
import mimetypes
from pathlib import Path
import re
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, unquote, parse_qs

from . import __version__
from .schema import default_config, preset_config, PRESETS, LIMITATIONS, MAX_CHROMATIN_BEADS, MAX_PARTICLES, MAX_GENES, MAX_COHESIN, MAX_REGULATORY_REGIONS, MAX_TRANSCRIPTS, checked_config, ConfigError, validate_config

ROOT = Path(__file__).resolve().parent.parent
STATIC = ROOT / "static"
TUTORIAL = ROOT / "tutorial"
TUTORIAL_FILES = frozenset({
    "chromatin-studio-tutorial.mp4", "chromatin-studio-tutorial-poster.jpg",
    "chromatin-studio-tutorial.srt", "chromatin-studio-tutorial.vtt",
    "chromatin-studio-tutorial-chapters.json",
    "index.html", "player.js", "chromatin-studio-tutorial-2x.mp4",
    "chromatin-studio-tutorial-2x-poster.jpg", "chromatin-studio-tutorial-2x.srt",
    "chromatin-studio-tutorial-2x.vtt", "storyboard.json", "transcript.txt",
    "chromatin-studio-tutorial-2x-chapters.json",
})


def capabilities():
    available = all(importlib.util.find_spec(m) is not None for m in ("openmm", "polychrom"))
    versions = {}
    for package in ("numpy", "scipy", "openmm", "polychrom"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    polychrom_status = {"id": "polychrom", "available": available,
                       "detail": "Install optional OpenMM/polychrom dependencies to enable"}
    if available:
        from .openmm_backend import backend_status
        polychrom_status = backend_status()
    return {
        "version": __version__, "versions": versions,
        "backends": [
            {"id": "numpy", "available": True, "detail": "NumPy Brownian reference engine · local CPU"},
            polychrom_status,
        ], "limitations": LIMITATIONS,
        "units": {"length": "reduced length", "time": "reduced time", "energy": "kBT"},
        "calibration": {"available": True, "method": "fixed-length MSD time mapping",
                        "coordinates_remain_reduced": True, "experimental_reference_required": True},
        "studies": {"available": True, "engines": ["reference", "polymer"],
                    "preserves_live_simulation": True, "spatial_matched_homologs": False},
        "cluster": {"export_available": True, "scheduler": "Slurm job arrays or shell",
                    "platforms": ["CUDA", "OpenCL", "CPU", "Reference", "numpy"],
                    "remote_submission": False, "resume_supported": False,
                    "acceleration": "OpenMM integration only; forces, reactions and RNA remain on CPU."},
        "limits": {"particles": MAX_PARTICLES, "chromatin_beads": MAX_CHROMATIN_BEADS,
                   "species": 32, "genes": MAX_GENES, "cohesin": MAX_COHESIN,
                   "regulatory_regions": MAX_REGULATORY_REGIONS,
                   "transcripts": MAX_TRANSCRIPTS},
    }


class Workspace:
    def __init__(self):
        self.lock = threading.RLock()
        self.config = default_config()
        self.backend = "numpy"
        self.sim = None
        self.history = deque(maxlen=2000)
        self.created = time.time()
        self.failure = None

    def reset(self, config, backend):
        config = checked_config(config)
        if backend == "numpy":
            from .engine import Simulation
            replacement = Simulation(config)
        elif backend == "polychrom":
            if not all(importlib.util.find_spec(m) for m in ("openmm", "polychrom")):
                raise ValueError("Polychrom/OpenMM is unavailable in this Python environment. Choose the NumPy reference engine or install the optional dependencies.")
            from .openmm_backend import PolychromSimulation
            replacement = PolychromSimulation(config)
        else:
            raise ValueError("Unknown engine. Choose numpy or polychrom.")
        first = replacement.snapshot()
        # Validate output before replacing the previous working model.
        json.dumps(first, allow_nan=False)
        self.config, self.backend, self.sim = config, backend, replacement
        self.failure = None
        self.history.clear()
        self.created = time.time()
        self.record(first)
        return first

    def ensure(self):
        if self.sim is None:
            self.reset(self.config, self.backend)

    def record(self, state):
        self.history.append({"time": state["time"], "step": state["step"],
                             "metrics": copy.deepcopy(state.get("metrics", {})),
                             "genes": copy.deepcopy(state.get("genes", []))})

    def step(self, count):
        self.ensure()
        if self.failure:
            raise ValueError("Simulation stopped after numerical failure. Reduce timestep or interaction strengths and apply/reset. " + self.failure)
        try:
            self.sim.step(count)
            state = self.sim.snapshot()
            json.dumps(state, allow_nan=False)
        except Exception as exc:
            self.failure = str(exc)
            raise
        self.record(state)
        return state


WORKSPACE = Workspace()
STUDY_LOCK = threading.Lock()


class Handler(BaseHTTPRequestHandler):
    server_version = "ChromatinStudio/" + __version__

    def log_message(self, fmt, *args):
        if args and len(args) > 1 and str(args[1]).startswith(("4", "5")):
            super().log_message(fmt, *args)

    def respond(self, value, status=200, filename=None):
        raw = json.dumps(value, allow_nan=False, separators=(",", ":")).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(raw)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        if filename:
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.end_headers()
        try:
            self.wfile.write(raw)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def do_GET(self):
        path = urlparse(self.path).path
        try:
            if path.startswith("/tutorial/"):
                return self.serve_tutorial(path)
            if path == "/api/config":
                return self.respond(copy.deepcopy(WORKSPACE.config))
            if path == "/api/presets":
                return self.respond(PRESETS)
            if path.startswith("/api/presets/"):
                return self.respond(preset_config(path.split("/")[-1]))
            if path == "/api/capabilities":
                return self.respond(capabilities())
            if path == "/api/studies/catalog":
                from .studies import study_catalog
                return self.respond(study_catalog())
            if path == "/api/studies/benchmarks":
                report = ROOT / "references" / "manuscript" / "verification.json"
                if not report.is_file():
                    return self.respond({"status": "not_run", "summary": "No saved manuscript-model verification report."})
                import hashlib
                saved = json.loads(report.read_text())
                saved["current_source_matches"] = (
                    hashlib.sha256((ROOT / "chromatin_studio" / "studies.py").read_bytes()).hexdigest()
                    == saved.get("reference_result", {}).get("source_sha256"))
                return self.respond(saved)
            if path in ("/api/calibration/references", "/api/calibration/reference"):
                directory = ROOT / "references" / "hansen"
                catalog = json.loads((directory / "catalog.json").read_text())
                if path.endswith("/references"):
                    return self.respond(catalog)
                query = parse_qs(urlparse(self.path).query)
                identifiers = query.get("id", [])
                if len(identifiers) != 1 or not identifiers[0]:
                    return self.respond({"error": "Supply exactly one reference id."}, 400)
                reference = next((entry for entry in catalog["references"]
                                  if entry["id"] == identifiers[0]), None)
                if reference is None:
                    return self.respond({"error": "Calibration reference not found."}, 404)
                # Only catalog-owned relative paths are opened, never user input.
                target = (directory / reference["curve_file"]).resolve()
                if not target.is_relative_to(directory.resolve()):
                    raise ValueError("Invalid calibration reference path in catalog.")
                return self.respond({"reference": reference,
                                     "experiment": json.loads(target.read_text())})
            if path == "/api/calibration/benchmarks":
                report = ROOT / "benchmarks" / "reference-results.json"
                if not report.is_file():
                    return self.respond({"status": "not_run", "levels": [],
                                         "summary": "Run the reference benchmark suite to generate results."})
                import hashlib
                import inspect
                from .engine import Simulation
                saved = json.loads(report.read_text())
                provenance = saved.get("provenance", {})
                saved["current_source_matches"] = (
                    hashlib.sha256(inspect.getsource(Simulation._brownian_move).encode()).hexdigest()
                    == provenance.get("integrator_source_sha256")
                    and hashlib.sha256((ROOT / "chromatin_studio" / "benchmarks.py").read_bytes()).hexdigest()
                    == provenance.get("benchmark_source_sha256"))
                return self.respond(saved)
            if path in ("/api/state", "/api/export", "/api/download-config"):
                with WORKSPACE.lock:
                    WORKSPACE.ensure()
                    if path == "/api/state":
                        return self.respond(WORKSPACE.sim.snapshot())
                    if path == "/api/download-config":
                        return self.respond(WORKSPACE.config, filename="chromatin-model.json")
                    from .calibration import unit_summary
                    return self.respond({
                        "config": WORKSPACE.config, "state": WORKSPACE.sim.snapshot(),
                        "history": list(WORKSPACE.history),
                        "metadata": {"version": __version__, "backend": WORKSPACE.backend,
                                     "created_unix": WORKSPACE.created, "exported_unix": time.time(),
                                     "units": "reduced", "history_limit": 2000,
                                     "physical_mapping": unit_summary(WORKSPACE.config, WORKSPACE.config.get("calibration")),
                                     "replay": "Configuration and seed reproduce an initial run; state is an inspection snapshot, not a restart checkpoint.",
                                     "capabilities": capabilities()},
                    }, filename="chromatin-experiment.json")
            if path.startswith("/api/"):
                return self.respond({"error": "Unknown API route."}, 404)
            relative = unquote(path).lstrip("/") or "index.html"
            if relative.startswith("static/"):
                relative = relative[len("static/"):]
            target = (STATIC / relative).resolve()
            if not target.is_relative_to(STATIC) or not target.is_file():
                return self.respond({"error": "File not found."}, 404)
            raw = target.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", mimetypes.guess_type(str(target))[0] or "application/octet-stream")
            self.send_header("Content-Length", str(len(raw)))
            self.send_header("Cache-Control", "no-cache")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
            self.end_headers()
            self.wfile.write(raw)
        except KeyError:
            self.respond({"error": "Preset not found."}, 404)
        except Exception as exc:
            self.respond({"error": str(exc)}, 500)

    def serve_tutorial(self, path):
        """Stream only published tutorial assets, including MP4 byte ranges."""
        name = unquote(path[len("/tutorial/"):])
        if name not in TUTORIAL_FILES:
            return self.respond({"error": "Tutorial file not found."}, 404)
        target = (TUTORIAL / name).resolve()
        if not target.is_relative_to(TUTORIAL.resolve()) or not target.is_file():
            return self.respond({"error": "Tutorial file not found."}, 404)
        size = target.stat().st_size
        start, end, status = 0, size - 1, 200
        requested = self.headers.get("Range") if target.suffix == ".mp4" else None
        if requested is not None:
            match = re.fullmatch(r"bytes=([0-9]*)-([0-9]*)", requested.strip())
            invalid = match is None or size == 0
            if not invalid:
                first, last = match.groups()
                # Bound parsing independently of Python's integer-string limit.
                invalid = not (first or last) or max(len(first), len(last)) > 20
                if not invalid and first:
                    start = int(first)
                    end = min(int(last), size - 1) if last else size - 1
                    invalid = start >= size or end < start
                elif not invalid:
                    count = int(last)
                    invalid = count == 0
                    start = max(0, size - count)
            if invalid:
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{size}")
                self.send_header("Accept-Ranges", "bytes")
                self.send_header("Content-Length", "0")
                self.send_header("Cache-Control", "no-cache")
                self.send_header("X-Content-Type-Options", "nosniff")
                self.end_headers()
                return
            status = 206
        content_type = {".vtt": "text/vtt; charset=utf-8", ".srt": "application/x-subrip; charset=utf-8"}.get(
            target.suffix, mimetypes.guess_type(str(target))[0] or "application/octet-stream")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Cache-Control", "no-cache")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
        if target.suffix == ".mp4":
            self.send_header("Accept-Ranges", "bytes")
        if status == 206:
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.end_headers()
        try:
            with target.open("rb") as stream:
                stream.seek(start)
                remaining = end - start + 1
                while remaining:
                    chunk = stream.read(min(256 * 1024, remaining))
                    if not chunk:
                        break
                    self.wfile.write(chunk)
                    remaining -= len(chunk)
        except (BrokenPipeError, ConnectionResetError):
            # Seeking can cancel an earlier browser request partway through.
            pass

    def do_POST(self):
        origin = self.headers.get("Origin")
        if origin and urlparse(origin).netloc != self.headers.get("Host"):
            return self.respond({"error": "Cross-origin requests are not accepted."}, 403)
        if self.headers.get("Content-Type", "").split(";")[0] != "application/json":
            return self.respond({"error": "Send application/json."}, 415)
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 1024 * 1024:
                return self.respond({"error": "Request body must be between 1 byte and 1 MB."}, 413)
            data = json.loads(self.rfile.read(size))
            path = urlparse(self.path).path
            if path == "/api/validate":
                return self.respond(validate_config(data))
            if not isinstance(data, dict):
                raise ValueError("Request body must be a JSON object.")
            if path in ("/api/studies/plan", "/api/studies/run"):
                from .studies import build_study, run_study
                # Each study owns fresh simulations; never replace or step the
                # viewer's Workspace, and never hold its state lock here.
                if path.endswith("/plan"):
                    return self.respond(build_study(data.get("config"), data.get("options", {})))
                if not STUDY_LOCK.acquire(blocking=False):
                    return self.respond({"error": "A local study is already running. Wait for its result before starting another."}, 409)
                try:
                    result = run_study(data.get("config"), data.get("options", {}))
                finally:
                    STUDY_LOCK.release()
                return self.respond(result)
            if path in ("/api/calibration/fit", "/api/calibration/summary", "/api/calibration/msd"):
                from .calibration import fit_time_mapping, unit_summary, analyze_trajectory
                # Analysis never replaces or advances the live trajectory.
                if path.endswith("/msd"):
                    return self.respond(analyze_trajectory(data))
                config = checked_config(data.get("config"))
                if path.endswith("/summary"):
                    return self.respond(unit_summary(config, config.get("calibration")))
                return self.respond(fit_time_mapping(data, config=config))
            if path in ("/api/cluster/plan", "/api/cluster/export"):
                from .cluster import plan_job, build_bundle
                # Packaging a draft never initializes, replaces, or advances
                # the local run; optional GPU dependencies are not needed here.
                if path.endswith("/plan"):
                    return self.respond(plan_job(data.get("config"), data.get("job", {})))
                raw = build_bundle(data.get("config"), data.get("job", {}))
                self.send_response(200)
                self.send_header("Content-Type", "application/zip")
                self.send_header("Content-Length", str(len(raw)))
                self.send_header("Content-Disposition", 'attachment; filename="chromatin-cluster-job.zip"')
                self.send_header("Cache-Control", "no-store")
                self.send_header("X-Content-Type-Options", "nosniff")
                self.end_headers()
                try:
                    self.wfile.write(raw)
                except (BrokenPipeError, ConnectionResetError):
                    pass
                return
            with WORKSPACE.lock:
                if path == "/api/kinetics/estimate":
                    from .kinetics import estimate_binding
                    return self.respond(estimate_binding(WORKSPACE.config, data))
                if path == "/api/reset":
                    return self.respond(WORKSPACE.reset(data.get("config"), data.get("backend", "numpy")))
                if path == "/api/step":
                    count = data.get("steps", 10)
                    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 100:
                        raise ValueError("steps must be an integer between 1 and 100.")
                    return self.respond(WORKSPACE.step(count))
            return self.respond({"error": "Unknown API route."}, 404)
        except ConfigError as exc:
            self.respond({"error": "Model configuration is invalid.", "details": exc.errors}, 400)
        except (ValueError, TypeError, OverflowError) as exc:
            self.respond({"error": str(exc)}, 400)
        except Exception as exc:
            self.respond({"error": str(exc)}, 500)


def main():
    parser = argparse.ArgumentParser(description="Run Chromatin Studio locally.")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--config", type=Path, help="Optional JSON configuration to load")
    parser.add_argument("--open", action="store_true", help="Open the local app in your default browser")
    args = parser.parse_args()
    if args.config:
        WORKSPACE.config = checked_config(json.loads(args.config.read_text()))
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"Chromatin Studio {__version__}: http://127.0.0.1:{args.port}", flush=True)
    print("Local exploratory simulation. Press Ctrl+C to stop.", flush=True)
    if args.open:
        import webbrowser
        threading.Timer(0.5, lambda: webbrowser.open(f"http://127.0.0.1:{args.port}")).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
