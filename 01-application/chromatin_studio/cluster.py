"""Validate and package independent, headless runs without contacting a cluster.

The job description never changes model limits and is not a restart checkpoint.
Scheduler options are deliberately a small safe subset; site-specific directives
can be edited in the exported template after inspection.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import io
import json
import math
from pathlib import Path
import re
import shlex
import zipfile

from . import __version__
from .schema import checked_config, validate_config


PLATFORMS = ("CUDA", "OpenCL", "CPU", "Reference", "numpy")
DEFAULT_JOB = {"platform": "CUDA", "precision": "mixed", "replicas": 4,
               "steps": 10000, "sample_every": 100, "snapshot_every": 1000,
               "cpus": 4, "memory_gb": 8, "walltime": "01:00:00",
               "job_name": "chromatin-studio"}
_INTEGER_LIMITS = {"replicas": (1, 1024), "steps": (1, 10**9),
                   "sample_every": (1, 10**9), "snapshot_every": (0, 10**9),
                   "cpus": (1, 256), "memory_gb": (1, 4096)}


def canonical_json(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def normalize_job(job, config):
    """Return a validated job; unknown fields fail instead of being ignored."""
    checked_config(config)
    if not isinstance(job, dict):
        raise ValueError("job must be an object.")
    unknown = set(job) - set(DEFAULT_JOB) - {"device_index", "partition", "account"}
    if unknown:
        raise ValueError("Unknown job fields: " + ", ".join(sorted(unknown)))
    result = {**DEFAULT_JOB, **job}
    if result["platform"] not in PLATFORMS:
        raise ValueError("platform must be CUDA, OpenCL, CPU, Reference, or numpy.")
    if result["precision"] not in ("single", "mixed", "double"):
        raise ValueError("precision must be single, mixed, or double.")
    for name, (low, high) in _INTEGER_LIMITS.items():
        value = result[name]
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not low <= value <= high or not math.isfinite(value) or value != int(value)):
            raise ValueError(f"{name} must be an integer between {low} and {high}.")
        result[name] = int(value)
    for name in ("job_name", "partition", "account"):
        value = result.get(name)
        if name != "job_name" and value in (None, ""):
            result.pop(name, None)
            continue
        if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", value):
            raise ValueError(f"{name} must contain 1–64 letters, digits, dots, underscores or hyphens, beginning with a letter or digit.")
    walltime = result["walltime"]
    match = re.fullmatch(r"(?:(\d{1,3})-)?(\d{1,3}):([0-5]\d):([0-5]\d)", walltime) if isinstance(walltime, str) else None
    if not match:
        raise ValueError("walltime must be HH:MM:SS or D-HH:MM:SS.")
    days, hours, minutes, seconds = (int(v or 0) for v in match.groups())
    if (match.group(1) is not None and hours > 23) or not 0 < ((days * 24 + hours) * 60 + minutes) * 60 + seconds <= 365 * 86400:
        raise ValueError("walltime must be positive and no longer than 365 days; day-format hours must be 0–23.")
    value = result.get("device_index")
    if value in (None, "") and not isinstance(value, bool):
        result.pop("device_index", None)
    else:
        if isinstance(value, bool) or not isinstance(value, (str, int)) or not re.fullmatch(r"\d{1,5}", str(value)) or int(value) > 65535:
            raise ValueError("device_index must be a single nonnegative logical device index, or omitted.")
        if result["platform"] not in ("CUDA", "OpenCL"):
            raise ValueError("device_index is only applicable to CUDA and OpenCL.")
        result["device_index"] = int(value)
    if int(config["seed"]) + result["replicas"] - 1 > 2**32 - 1:
        raise ValueError("The initial seed plus replica count exceeds the model seed range; choose a smaller seed or fewer replicas.")
    return result


def _scripts(job):
    arguments = "--platform " + shlex.quote(job["platform"]) + " --precision " + shlex.quote(job["precision"])
    if "device_index" in job:
        arguments += " --device-index " + str(job["device_index"])
    shell = f'''#!/bin/sh
set -eu
# Run from the extracted bundle; activate a prepared Python environment first.
cd -- "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPLICA=${{1:-0}}
case "$REPLICA" in ''|*[!0-9]*) echo "Replica must be a nonnegative integer." >&2; exit 2;; esac
if [ "$REPLICA" -ge {job['replicas']} ]; then echo "Replica is outside this job array." >&2; exit 2; fi
export OMP_NUM_THREADS={job['cpus']}
export OPENMM_CPU_THREADS={job['cpus']}
export OPENBLAS_NUM_THREADS={job['cpus']}
export MKL_NUM_THREADS={job['cpus']}
mkdir -p results
# A lock refuses accidental concurrent repeats and also records failed attempts.
if ! mkdir "results/replica-$REPLICA.lock"; then
  echo "Replica already attempted. Keep its results; use a fresh results directory for a new attempt." >&2
  exit 2
fi
# Verify immutable bundle files before running any scientific tests.
python -c "from chromatin_studio.batch import verify_bundle; verify_bundle('.')"
# Scientific checks run on the same allocated node and requested platform.
python -m chromatin_studio.validation {arguments} --output "results/validation-$REPLICA.json"
python -m chromatin_studio.batch --model model.json --job job.json --replica "$REPLICA" --output "results/replica-$REPLICA"
'''
    directives = ["#!/bin/sh", f"#SBATCH --job-name={job['job_name']}", f"#SBATCH --array=0-{job['replicas'] - 1}",
                  "#SBATCH --nodes=1", "#SBATCH --ntasks=1", f"#SBATCH --cpus-per-task={job['cpus']}",
                  f"#SBATCH --mem={job['memory_gb']}G", f"#SBATCH --time={job['walltime']}",
                  "#SBATCH --output=slurm-%A_%a.out"]
    if job["platform"] in ("CUDA", "OpenCL"):
        directives.append("#SBATCH --gres=gpu:1")
    for name in ("partition", "account"):
        if name in job:
            directives.append(f"#SBATCH --{name}={job[name]}")
    slurm = "\n".join(directives) + '''
set -eu
# Submit from the extracted bundle directory. Adapt GPU/queue directives to your site.
cd -- "${SLURM_SUBMIT_DIR:?Submit this script from the extracted bundle directory}"
# Do not replace CUDA_VISIBLE_DEVICES: the scheduler owns GPU visibility.
sh ./run.sh "${SLURM_ARRAY_TASK_ID:?Missing array task ID}"
'''
    return {"slurm": slurm, "shell": shell}


def plan_job(config, job):
    from .calibration import unit_summary
    config = checked_config(config)
    job = normalize_job(job, config)
    report = validate_config(config)
    warnings = list(report["warnings"])
    warnings += [
        "Each replica starts from the model and its own seed; this export does not resume the current viewer state.",
        "Forces, reactions and RNA diffusion run on CPU. CUDA or OpenCL runs only the Brownian integrator on the selected device; acceleration is not guaranteed. OpenCL device type may be CPU or GPU.",
        "Memory, walltime and GPU availability are scheduler requests, not capacity estimates. Pilot the model before a large array.",
        "Existing model particle and transcript limits still apply. A GPU request does not enable a 100,000-nucleosome model.",
        "SLURM GPU and queue directives are site-specific; inspect the template before submitting.",
    ]
    if job["platform"] in ("CUDA", "OpenCL"):
        warnings.append("One GPU is requested per independent replica. Default device selection uses the devices visible on the allocated node; requested GPU failures stop the run.")
    else:
        warnings.append("This job requests CPU execution; the precision selection only configures CUDA/OpenCL platforms.")
    metric_records = 1 + math.ceil(job["steps"] / job["sample_every"])
    full_snapshots = 0 if job["snapshot_every"] == 0 else 1 + math.ceil(job["steps"] / job["snapshot_every"])
    if full_snapshots * job["replicas"] > 10000 or metric_records * job["replicas"] > 100000:
        warnings.append(f"Large output request: {full_snapshots * job['replicas']:,} full snapshots and {metric_records * job['replicas']:,} metric records across all replicas. JSON snapshots include positions, tracks, proteins and RNA; check disk capacity and consider wider sampling intervals.")
    return {"job": job,
            "model": {"name": config["name"], "particle_count": report["particle_count"],
                      "chromatin_beads": config["chromatin"]["beads"],
                      "genomic_span_bp": report["genomic_span_bp"], "dt": config["simulation"]["dt"],
                      "steps": job["steps"], "reduced_duration": job["steps"] * config["simulation"]["dt"],
                      "physical_mapping": unit_summary(config, config.get("calibration"))},
            "replicas": [{"index": i, "seed": int(config["seed"]) + i} for i in range(job["replicas"])],
            "output": {"metric_records_per_replica": metric_records, "full_snapshots_per_replica": full_snapshots,
                       "final_states_per_replica": 1, "replicas": job["replicas"]},
            "warnings": warnings, "files": ["model.json", "job.json", "physical-units.json", "bundle-manifest.json", "run.sh", "submit.slurm", "RUN.md", "requirements.txt", "requirements-polychrom.txt", "chromatin_studio/"],
            "scripts": _scripts(job), "resume_supported": False}


def _run_instructions(job):
    dependency_file = "requirements.txt" if job["platform"] == "numpy" else "requirements-polychrom.txt"
    return f'''# Chromatin Studio cluster bundle

This is a fresh model run, not a checkpoint of the browser's current state.
Each replica uses `model.seed + replica_index`. The SLURM array requests one
independent replica per task and one GPU per task only for CUDA/OpenCL.

## Prepare the environment on your cluster

Use Python 3.10 or later and your cluster's supported GPU drivers/runtime. Follow
site instructions for environment modules and GPU allocations. From this directory:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r {dependency_file}
```

Do environment installation on a node where your site permits network access.
For offline compute nodes prepare the environment/wheels beforehand. NumPy and
SciPy minimum versions are portable requirements, not an exact environment lock;
the run manifest records actual versions. The optional file pins OpenMM and the
polychrom source commit used during development.

## Verify, then execute

Request a suitable interactive node to check GPU initialization without advancing
the model (this creates a new directory and refuses existing output):

```sh
python -m chromatin_studio.batch --model model.json --job job.json --replica 0 --preflight --output preflight-0
```

On an allocated node, `sh run.sh 0` runs scientific validation followed by replica
0. Repeat with replica indices 0 through {job['replicas'] - 1} for a generic scheduler.
Inspect `submit.slurm`, make a copy, adjust the copy's directives for your site's
GPU/partition/account rules, activate your environment, then submit from this
directory. Keep bundled files unchanged so their integrity checks remain useful:

```sh
cp submit.slurm submit-local.slurm
# Edit submit-local.slurm to match your site's requirements.
sbatch submit-local.slurm
```

GPU directives vary between clusters (for example `--gres=gpu:1` versus `--gpus=1`).
The supplied template uses `--gres=gpu:1`. The scheduler controls device visibility;
do not replace `CUDA_VISIBLE_DEVICES`. Leave `device_index` unset unless the site's
OpenMM configuration requires a logical device from the visible set. Requested
GPU failures stop the job; there is no automatic CPU fallback.

The export does not submit jobs or copy data. Transfer the ZIP or extracted folder
using your cluster's approved file-transfer tool. Keep login details out of the
model. No cluster connection or credentials are required by Chromatin Studio.

## Outputs and limits

Each `results/replica-N/` contains the seeded model, job, a run manifest, streamed
metrics JSONL, optional snapshot JSONL, and an atomically written final state.
Reports record source/configuration hashes, package versions, requested/actual
backend details and completion or failure. Bundle hashes are checked on startup.
Validation reports are stored beside replica directories. A `.lock` directory
prevents duplicate script launches. Failed attempts remain for diagnosis; run a
new attempt in a fresh results directory instead of deleting evidence blindly.
Abrupt termination can leave a manifest marked `running` and an incomplete last
JSONL line; final-state presence and manifest `complete` identify clean completion.

Metrics and sparse snapshots stream to disk; output size still grows with runtime.
Snapshots are inspection data, not restart checkpoints. Coordinates and times
remain in reduced units. A physical mapping, if attached, is preserved in
`model.json`, `physical-units.json`, and each run manifest, with its provenance
and baseline comparison. An absent, synthetic, or needs-review mapping is not
an experimentally established clock. Current particle/transcript
limits remain enforced. Forces, chemistry and RNA diffusion still run on CPU,
so GPU integration alone is not a promise of higher throughput. Pilot first.
The scientific suite detects several implementation/physics regressions; passing
it does not calibrate the biological model or establish timestep convergence for
your chosen parameters.

## Calculate chromatin MSD for calibration

Enable snapshots and choose an interval that resolves the experimental lag
window. After equilibration, compute the MSD of selected chromatin beads:

```sh
python -m chromatin_studio.calibration msd --snapshots results/replica-0/snapshots.jsonl --beads 10,20 --dimension 2 --observable relative_pair --discard-before 0 --output simulation-msd.json
```

Replace bead indices and discard time with your study's actual choices; zero
discard is only command syntax, not an equilibration recommendation. A final
partial snapshot is irregularly spaced: use a step count divisible by the
snapshot interval, or explicitly select the analyzer's final-frame drop option.
Import the resulting curve into Calibrate in the local UI. Match the observable,
dimensions, drift treatment and imaging corrections to the experimental curve.
Keep the baseline mapping fixed when comparing perturbations. Use independent
replicas, longer runs and smaller timesteps to check convergence.
'''


def build_bundle(config, job):
    """Return a portable ZIP with model, headless source, scripts, and hashes."""
    config = checked_config(config)
    plan = plan_job(config, job)
    job = plan["job"]
    root = Path(__file__).resolve().parent.parent
    files = {"model.json": canonical_json(config), "job.json": canonical_json(job),
             "physical-units.json": canonical_json(plan["model"]["physical_mapping"]),
             "run.sh": plan["scripts"]["shell"].encode(), "submit.slurm": plan["scripts"]["slurm"].encode(),
             "RUN.md": _run_instructions(job).encode()}
    for filename in ("requirements.txt", "requirements-polychrom.txt"):
        files[filename] = (root / filename).read_bytes()
    for source in sorted((root / "chromatin_studio").glob("*.py")):
        files["chromatin_studio/" + source.name] = source.read_bytes()
    if "chromatin_studio/validation.py" not in files:
        raise RuntimeError("Scientific validation source is missing; the export would not be runnable.")
    versions = {}
    for name in ("numpy", "scipy", "openmm", "polychrom"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    files["bundle-manifest.json"] = canonical_json({"format": 1, "application_version": __version__,
        "fresh_run_only": True, "export_environment_versions": versions,
        "files": {name: {"sha256": sha256(data), "bytes": len(data)} for name, data in sorted(files.items())}})
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(files.items()):
            entry = zipfile.ZipInfo(name)
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = (0o755 if name.endswith((".sh", ".slurm")) else 0o644) << 16
            archive.writestr(entry, data)
    return stream.getvalue()
