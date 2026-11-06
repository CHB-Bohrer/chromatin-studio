"""Bounded-memory, headless simulation runs with explicit backend selection."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

import numpy as np

from . import __version__
from .cluster import canonical_json, normalize_job, sha256
from .schema import checked_config
from .calibration import unit_summary


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    raw = canonical_json(value)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="wb", dir=path.parent, prefix="." + path.name + ".", delete=False) as stream:
            temporary = stream.name
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


def verify_bundle(directory):
    """Verify bundled files if a manifest exists; return its hash or None."""
    directory = Path(directory).resolve()
    manifest_path = directory / "bundle-manifest.json"
    if not manifest_path.exists():
        return None
    raw = manifest_path.read_bytes()
    manifest = json.loads(raw)
    if manifest.get("format") != 1 or not isinstance(manifest.get("files"), dict):
        raise ValueError("Unsupported bundle manifest.")
    for name, expected in manifest["files"].items():
        path = (directory / name).resolve()
        if not path.is_relative_to(directory) or not path.is_file():
            raise ValueError(f"Bundle file is missing or outside its directory: {name}")
        content = path.read_bytes()
        if sha256(content) != expected.get("sha256") or len(content) != expected.get("bytes"):
            raise ValueError(f"Bundle integrity check failed for {name}; export a fresh bundle after model/source changes.")
    return sha256(raw)


def _environment():
    versions = {}
    for name in ("numpy", "scipy", "openmm", "polychrom"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return {"python": sys.version, "platform": platform.platform(), "packages": versions,
            "scheduler": {name: os.environ[name] for name in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK", "CUDA_VISIBLE_DEVICES", "OPENMM_CPU_THREADS", "OMP_NUM_THREADS") if name in os.environ}}


def _create_simulation(config, job):
    if job["platform"] == "numpy":
        from .engine import Simulation
        return Simulation(config)
    from .openmm_backend import PolychromSimulation
    return PolychromSimulation(config, platform=job["platform"], precision=job["precision"], device_index=job.get("device_index"))


def inspect_state(sim, state=None):
    """Fail before writing success if numerical or molecular invariants break."""
    state = sim.snapshot() if state is None else state
    for label, values in (("positions", sim.positions), ("diffusion", sim.diffusion), ("forces", sim.get_forces())):
        if not np.isfinite(np.asarray(values)).all():
            raise FloatingPointError(f"Nonfinite {label} at step {sim.step_count}.")
    if not np.isfinite(sim.potential_energy()):
        raise FloatingPointError(f"Nonfinite model energy at step {sim.step_count}.")
    if np.any(np.asarray(sim.diffusion) < 0):
        raise ValueError("Negative diffusion coefficient.")
    for molecule in state["interaction_state"]["molecules"]:
        if not 0 <= molecule["used_valence"] <= molecule["valence"]:
            raise ValueError("Protein valence accounting failed.")
    rna = state["transcription"]
    if rna["enabled"]:
        totals = rna["totals"]
        records = rna["transcripts"]
        if not (totals["live"] == len(records) <= rna["cap"] and
                totals["initiated"] == totals["live"] + totals["degraded"] and
                totals["completed"] == totals["unprocessed"] + totals["mature"] + totals["degraded"] and
                totals["live"] == totals["nascent"] + totals["unprocessed"] + totals["mature"]):
            raise ValueError("RNA mass balance or capacity accounting failed.")
        if len({record["id"] for record in records}) != len(records):
            raise ValueError("Duplicate live RNA identity.")
    # Serializability catches nonfinite RNA coordinates and summaries as well.
    json.dumps(state, allow_nan=False)
    return state


def _metric_record(state):
    return {"step": state["step"], "time": state["time"], "metrics": state["metrics"],
            "genes": state["genes"], "nucleosomes": {key: value for key, value in state["nucleosomes"].items() if key != "occupied"},
            "transcription": {"enabled": state["transcription"]["enabled"], "totals": state["transcription"]["totals"]}}


def _write_record(stream, record):
    stream.write(json.dumps(record, allow_nan=False, separators=(",", ":")) + "\n")
    stream.flush()


def _record_rna_limits(manifest, state):
    if state["transcription"]["enabled"]:
        totals = state["transcription"]["totals"]
        manifest["rna_totals"] = copy.deepcopy(totals)
        blocked = totals["blocked_initiations"]
        manifest["warnings"] = ([f"The transcript capacity blocked {blocked} initiation attempts. This numerical limit can change transcription output; increase capacity within supported limits or reduce load and compare."] if blocked else [])


def run_batch(config, job, *, replica, output, preflight=False, bundle_hash=None):
    """Execute one replica, keeping only current state in memory and refusing overwrite."""
    config = checked_config(config)
    job = normalize_job(job, config)
    if isinstance(replica, bool) or not isinstance(replica, int) or not 0 <= replica < job["replicas"]:
        raise ValueError("replica must be a zero-based index within this job's replica count.")
    seeded = copy.deepcopy(config)
    seeded["seed"] = int(config["seed"]) + replica
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    source_hashes = {path.name: sha256(path.read_bytes()) for path in sorted(Path(__file__).parent.glob("*.py"))}
    environment = _environment()
    manifest = {"format": 1, "application_version": __version__, "status": "running", "preflight_only": bool(preflight),
                "started_utc": _utc_now(), "replica": replica, "seed": seeded["seed"], "requested_job": job,
                "requested_steps": 0 if preflight else job["steps"], "completed_steps": 0,
                "model_sha256": sha256(canonical_json(config)), "seeded_model_sha256": sha256(canonical_json(seeded)),
                "job_sha256": sha256(canonical_json(job)), "bundle_manifest_sha256": bundle_hash,
                "source_sha256": source_hashes, "environment": environment,
                "environment_sha256": sha256(canonical_json(environment)), "restart_checkpoint": False,
                "coordinate_units": "reduced", "time_units": "reduced",
                "physical_mapping": unit_summary(config, config.get("calibration"))}
    started = time.monotonic()
    sim = None
    atomic_json(output / "run-manifest.json", manifest)
    atomic_json(output / "model.json", seeded)
    atomic_json(output / "job.json", job)
    atomic_json(output / "physical-units.json", manifest["physical_mapping"])
    try:
        sim = _create_simulation(seeded, job)
        manifest["backend"] = getattr(sim, "backend_metadata", {"requested_platform": "numpy", "actual_platform": "numpy", "actual_precision": "float64", "forces": "CPU", "reactions": "CPU", "rna_diffusion": "CPU"})
        atomic_json(output / "run-manifest.json", manifest)
        state = inspect_state(sim)
        _record_rna_limits(manifest, state)
        if preflight:
            atomic_json(output / "preflight-state.json", state)
        else:
            snapshots = None
            try:
                if job["snapshot_every"]:
                    snapshots = (output / "snapshots.jsonl").open("x", encoding="utf-8")
                with (output / "metrics.jsonl").open("x", encoding="utf-8") as metrics:
                    _write_record(metrics, _metric_record(state))
                    if snapshots:
                        _write_record(snapshots, state)
                    current = 0
                    while current < job["steps"]:
                        next_sample = (current // job["sample_every"] + 1) * job["sample_every"]
                        next_snapshot = ((current // job["snapshot_every"] + 1) * job["snapshot_every"] if snapshots else job["steps"])
                        target = min(next_sample, next_snapshot, job["steps"])
                        # Work chunks cap delay between safety checks and failure metadata.
                        while current < target:
                            block = min(target - current, 100)
                            sim.step(block)
                            current += block
                            if not np.isfinite(sim.positions).all():
                                raise FloatingPointError("Nonfinite coordinates during batch integration.")
                        state = inspect_state(sim)
                        _record_rna_limits(manifest, state)
                        if current % job["sample_every"] == 0 or current == job["steps"]:
                            _write_record(metrics, _metric_record(state))
                        if snapshots and (current % job["snapshot_every"] == 0 or current == job["steps"]):
                            _write_record(snapshots, state)
                        manifest["completed_steps"] = current
                        manifest["elapsed_seconds"] = time.monotonic() - started
                        atomic_json(output / "run-manifest.json", manifest)
                atomic_json(output / "final-state.json", state)
            finally:
                if snapshots:
                    snapshots.close()
        manifest.update(status="complete", completed_steps=sim.step_count, finished_utc=_utc_now(), elapsed_seconds=time.monotonic() - started)
        atomic_json(output / "run-manifest.json", manifest)
        return manifest
    except BaseException as exc:
        manifest.update(status="failed", completed_steps=sim.step_count if sim is not None else 0,
                        finished_utc=_utc_now(), elapsed_seconds=time.monotonic() - started,
                        error={"type": type(exc).__name__, "message": str(exc)})
        atomic_json(output / "run-manifest.json", manifest)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--job", required=True, type=Path)
    parser.add_argument("--replica", type=int, default=0)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--preflight", action="store_true", help="Initialize the requested backend and check the model without advancing it.")
    arguments = parser.parse_args(argv)
    try:
        bundle_hash = verify_bundle(arguments.model.parent)
        model = json.loads(arguments.model.read_text())
        job = json.loads(arguments.job.read_text())
        manifest = run_batch(model, job, replica=arguments.replica, output=arguments.output,
                             preflight=arguments.preflight, bundle_hash=bundle_hash)
        print(json.dumps({"status": manifest["status"], "replica": arguments.replica,
                          "steps": manifest["completed_steps"], "output": str(arguments.output)}))
        return 0
    except Exception as exc:
        print(f"Chromatin batch failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
