"""Measured MSD curves, explicit unit mappings, and reproducible clock fitting.

No function changes simulation dynamics. A fitted clock is conditional on a
chosen length scale, the supplied reference, observable, and fitting window;
it is not validation of unrelated mechanisms or biological parameters.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np


MAX_CURVE_POINTS = 512
OBSERVABLES = ("locus", "center_of_mass", "relative_pair")
TIME_BOUNDS = (1e-12, 1e12)
LENGTH_BOUNDS = (1e-6, 1e9)
THRESHOLDS = {"maximum_slope_difference": .15, "maximum_log_rmse": .25,
              "minimum_absolute_simulation_slope": .1}


def _number(value, name, low=0, high=1e30, *, strict_low=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number.")
    try:
        value = float(value)
    except (OverflowError, ValueError):
        raise ValueError(f"{name} must be a finite number.") from None
    if not math.isfinite(value) or value > high or (value <= low if strict_low else value < low):
        raise ValueError(f"{name} is outside its supported finite range ({low:g} to {high:g}).")
    return value


def _string(value, name, *, required=True, limit=4000):
    if not isinstance(value, str) or len(value) > limit or (required and not value.strip()):
        raise ValueError(f"{name} must be {'a nonempty' if required else 'a'} string of at most {limit} characters.")
    return value.strip()


def _conventions(value, prefix):
    dimension = value.get("dimension")
    if isinstance(dimension, bool) or dimension not in (2, 3):
        raise ValueError(f"{prefix}.dimension must be 2 or 3.")
    observable = value.get("observable")
    if observable not in OBSERVABLES:
        raise ValueError(f"{prefix}.observable must be locus, center_of_mass, or relative_pair.")
    drift = value.get("drift_correction")
    if drift not in ("none", "center_of_mass"):
        raise ValueError(f"{prefix}.drift_correction must be none or center_of_mass.")
    if observable == "relative_pair" and drift != "none":
        raise ValueError("relative_pair already cancels common translation; use drift_correction='none'.")
    return {"dimension": int(dimension), "observable": observable, "drift_correction": drift}


def config_fingerprint(config):
    """Hash model parameters, excluding mapping, display name and replica seed."""
    if not isinstance(config, dict):
        raise ValueError("config must be an object.")
    model = {key: value for key, value in config.items() if key not in ("calibration", "seed", "name")}
    try:
        data = json.dumps(model, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    except (ValueError, TypeError) as exc:
        raise ValueError("config must contain finite JSON values.") from exc
    return hashlib.sha256(data).hexdigest()


def _fingerprint(value, name, allow_none=True):
    if allow_none and value is None:
        return None
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{name} must be a lowercase SHA-256 string or null.")
    return value


def _curve(value, name):
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a curve object.")
    allowed = {"lag", "msd", "dimension", "observable", "drift_correction", "source",
               "lag_unit", "msd_unit", "counts", "analysis", "baseline_fingerprint", "warnings"}
    if set(value) - allowed:
        raise ValueError(f"Unknown {name} curve fields: {', '.join(sorted(set(value)-allowed))}.")
    convention = _conventions(value, name)
    lag, msd = value.get("lag"), value.get("msd")
    if (not isinstance(lag, list) or not isinstance(msd, list) or len(lag) != len(msd)
            or not 4 <= len(lag) <= MAX_CURVE_POINTS):
        raise ValueError(f"{name} needs 4–{MAX_CURVE_POINTS} paired lag/MSD values.")
    lag = [_number(x, name+".lag", 1e-15, 1e15) for x in lag]
    msd = [_number(x, name+".msd", 1e-30, 1e30) for x in msd]
    if any(b <= a for a, b in zip(lag, lag[1:])):
        raise ValueError(f"{name}.lag must be strictly increasing; duplicate or unordered lags are not accepted.")
    expected_units = (("reduced_time", "reduced_length_squared") if name == "simulation"
                      else ("s", "um2"))
    for key, expected in zip(("lag_unit", "msd_unit"), expected_units):
        if key in value and value[key] != expected:
            raise ValueError(f"{name}.{key} must be {expected}.")
    result = {"lag": lag, "msd": msd, **convention,
              "source": _string(value.get("source"), name+".source"),
              "lag_unit": expected_units[0], "msd_unit": expected_units[1]}
    if "baseline_fingerprint" in value:
        result["baseline_fingerprint"] = _fingerprint(value["baseline_fingerprint"], name+".baseline_fingerprint")
    if "counts" in value:
        counts = value["counts"]
        if (not isinstance(counts, list) or len(counts) != len(lag)
                or any(isinstance(x, bool) or not isinstance(x, int) or x < 1 for x in counts)):
            raise ValueError(f"{name}.counts must contain one positive integer per lag.")
        result["counts"] = counts[:]
    for key in ("analysis", "warnings"):
        if key in value:
            try:
                result[key] = json.loads(json.dumps(value[key], allow_nan=False))
            except (ValueError, TypeError) as exc:
                raise ValueError(f"{name}.{key} must contain finite JSON values.") from exc
    return result


def analyze_trajectory(payload):
    """Time-averaged MSD from regular, fixed-identity chromatin coordinates.

    Supply ``frames`` containing batch snapshots (time, chromatin), or ``time``
    and ``positions`` arrays (frame × bead × coordinate). ``beads`` selects
    indices. Locus averages selected-bead squared displacements; center_of_mass
    tracks their equal-weight center; relative_pair tracks r[bead0]-r[bead1].
    COM drift correction subtracts the *whole supplied polymer* COM per frame.
    Pair counts are overlapping time-origin counts, not independent samples.
    """
    if not isinstance(payload, dict):
        raise ValueError("Trajectory request must be an object.")
    convention = _conventions(payload, "trajectory")
    source = _string(payload.get("source"), "trajectory.source")
    if "frames" in payload:
        frames = payload["frames"]
        if not isinstance(frames, list) or not frames:
            raise ValueError("frames must be a nonempty list of snapshots.")
        try:
            times = np.asarray([f["time"] for f in frames], dtype=float)
            positions = np.asarray([f["chromatin"] for f in frames], dtype=float)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Every snapshot must have time and equally sized chromatin coordinate arrays.") from exc
    else:
        try:
            times = np.asarray(payload.get("time"), dtype=float)
            positions = np.asarray(payload.get("positions"), dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("time and positions must be regular numeric arrays.") from exc
    if (times.ndim != 1 or positions.ndim != 3 or len(positions) != len(times)
            or positions.shape[1] < 1 or positions.shape[2] not in (2, 3)
            or positions.shape[2] < convention["dimension"] or not np.isfinite(times).all()
            or not np.isfinite(positions).all() or np.any(np.diff(times) <= 0)):
        raise ValueError("Need finite, increasing times and fixed bead identities in frame × bead × 2/3 coordinates.")
    if positions.size > 12_000_000:
        raise ValueError("Trajectory exceeds the 12-million-coordinate analysis limit; select a smaller time window/file.")
    beads = payload.get("beads")
    if (not isinstance(beads, list) or not beads
            or any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < positions.shape[1] for i in beads)
            or len(set(beads)) != len(beads)):
        raise ValueError("beads must name distinct existing zero-based chromatin indices.")
    if convention["observable"] == "relative_pair" and len(beads) != 2:
        raise ValueError("relative_pair requires exactly two distinct beads.")
    discard = _number(payload.get("discard_before", 0), "discard_before", 0, 1e15)
    keep = times >= discard
    discarded = int((~keep).sum())
    times, positions = times[keep], positions[keep]
    warnings = []
    if len(times) < 9:
        raise ValueError("At least nine retained frames are needed for four lags no longer than half the track.")
    intervals = np.diff(times)
    step = float(np.median(intervals))
    regular = np.isclose(intervals, step, rtol=1e-6, atol=max(1e-12, abs(step)*1e-9))
    drop_final = payload.get("drop_final_irregular", False)
    if not isinstance(drop_final, bool):
        raise ValueError("drop_final_irregular must be true or false.")
    dropped_final = False
    if not regular.all():
        if drop_final and regular[:-1].all() and not regular[-1] and len(times) >= 10:
            times, positions = times[:-1], positions[:-1]
            dropped_final = True
            warnings.append("Explicitly dropped the final irregularly spaced frame; all preceding spacing was regular.")
        else:
            raise ValueError("MSD requires regular frame spacing. A final partial output interval may be dropped explicitly with drop_final_irregular; interior gaps are not accepted.")
    positions = positions[:, :, :convention["dimension"]].copy()
    if convention["drift_correction"] == "center_of_mass":
        positions -= positions.mean(axis=1, keepdims=True)
        warnings.append("Removed the whole supplied chromatin COM motion; this changes the observable and is not laboratory-stage drift correction.")
        if convention["observable"] == "center_of_mass" and len(beads) == positions.shape[1]:
            raise ValueError("Whole-polymer COM motion becomes identically zero after whole-polymer COM subtraction.")
    selected = positions[:, beads]
    if convention["observable"] == "center_of_mass":
        selected = selected.mean(axis=1, keepdims=True)
    elif convention["observable"] == "relative_pair":
        selected = selected[:, :1] - selected[:, 1:2]
    max_lag = payload.get("max_lag_steps", len(times)//2)
    if isinstance(max_lag, bool) or not isinstance(max_lag, int) or not 4 <= max_lag <= len(times)//2:
        raise ValueError("max_lag_steps must be an integer from 4 to half the retained frame count.")
    # Bound work/output while retaining short lags and the requested endpoint.
    offsets = (np.arange(1, max_lag+1) if max_lag <= MAX_CURVE_POINTS else
               np.unique(np.rint(np.geomspace(1, max_lag, MAX_CURVE_POINTS)).astype(int)))
    msd, counts = [], []
    for offset in offsets:
        delta = selected[offset:] - selected[:-offset]
        msd.append(float(np.mean(np.sum(delta*delta, axis=2))))
        counts.append(int(delta.shape[0]*delta.shape[1]))
    if not np.isfinite(msd).all():
        raise ValueError("MSD calculation overflowed; check coordinate units and magnitudes.")
    warnings += ["Time origins and beads can be correlated; pair counts are not independent sample sizes or confidence intervals.",
                 "No localization-noise or exposure-blur correction was applied. Match those conventions before fitting experimental data.",
                 "Choose discard_before from an equilibration analysis; retaining frames does not establish stationarity."]
    result = {"lag": (offsets*step).tolist(), "msd": msd, **convention, "source": source,
              "lag_unit": "reduced_time", "msd_unit": "reduced_length_squared", "counts": counts,
              "analysis": {"beads": beads[:], "retained_frames": len(times), "discarded_frames": discarded,
                           "discard_before": discard, "first_time": float(times[0]), "last_time": float(times[-1]),
                           "frame_interval": step, "max_lag_steps": max_lag, "dropped_final_frame": dropped_final,
                           "estimator": "overlapping time-origin mean squared displacement",
                           "projection": "first two coordinates (x,y)" if convention["dimension"] == 2 else "all three coordinates"},
              "warnings": warnings}
    if "config" in payload:
        result["baseline_fingerprint"] = config_fingerprint(payload["config"])
    return result


def _normalize_fit(payload):
    if not isinstance(payload, dict):
        raise ValueError("Fit request must be an object.")
    expected = {"simulation", "experiment", "length_nm", "length_source", "fit_min_s", "fit_max_s",
                "source_kind", "notes", "config"}
    if set(payload)-expected:
        raise ValueError("Unknown fit request fields: " + ", ".join(sorted(set(payload)-expected)))
    simulation, experiment = _curve(payload.get("simulation"), "simulation"), _curve(payload.get("experiment"), "experiment")
    for key in ("dimension", "observable", "drift_correction"):
        if simulation[key] != experiment[key]:
            raise ValueError(f"Simulation and experiment must use the same {key}; no silent dimensional or observable conversion is applied.")
    source_kind = payload.get("source_kind")
    if source_kind not in ("experimental", "synthetic"):
        raise ValueError("source_kind must explicitly be experimental or synthetic.")
    length = _number(payload.get("length_nm"), "length_nm", *LENGTH_BOUNDS)
    low = _number(payload.get("fit_min_s"), "fit_min_s", 1e-15, 1e15)
    high = _number(payload.get("fit_max_s"), "fit_max_s", 1e-15, 1e15)
    if not experiment["lag"][0] <= low < high <= experiment["lag"][-1]:
        raise ValueError("The fitting window must be increasing and fully within the experimental lag range.")
    return {"simulation": simulation, "experiment": experiment, "length_nm": length,
            "length_source": _string(payload.get("length_source"), "length_source"),
            "fit_min_s": low, "fit_max_s": high, "source_kind": source_kind,
            "notes": _string(payload.get("notes", ""), "notes", required=False)}


def _fit_log_shift(xs, ys, xe, ye, low, high):
    """Global least squares of piecewise-linear log curves, without extrapolation.

    Between interpolation knots residuals are affine in log(time scale), so
    their squared sum is quadratic. Evaluate the constrained minimum in every
    interval rather than relying on a potentially local scalar optimizer.
    """
    breaks = (xe[:, None]-xs[None, 1:-1]).ravel()
    breaks = np.unique(np.concatenate(([low, high], breaks[(breaks > low) & (breaks < high)])))
    slopes = np.diff(ys)/np.diff(xs)
    intercepts = ys[:-1]-slopes*xs[:-1]
    best_z, best_loss = low, math.inf
    if len(breaks) == 1:
        return low, float(np.mean((np.interp(xe-low, xs, ys)-ye)**2))
    for start in range(0, len(breaks)-1, 256):
        stop = min(start+256, len(breaks)-1)
        left, right = breaks[start:stop], breaks[start+1:stop+1]
        middle = (left+right)/2
        indices = np.clip(np.searchsorted(xs, xe[None, :]-middle[:, None], side="right")-1, 0, len(slopes)-1)
        a = slopes[indices]
        c = intercepts[indices]+a*xe[None, :]-ye[None, :]
        denominator = np.sum(a*a, axis=1)
        numerator = np.sum(a*c, axis=1)
        candidates = np.clip(np.divide(numerator, denominator, out=middle.copy(), where=denominator>1e-28), left, right)
        loss = np.mean((c-a*candidates[:, None])**2, axis=1)
        index = int(np.argmin(loss))
        if float(loss[index]) < best_loss:
            best_z, best_loss = float(candidates[index]), float(loss[index])
    return best_z, best_loss


def fit_time_mapping(payload, config=None):
    """Fit a single time scale at fixed length; retain diagnostics and raw input."""
    raw = _normalize_fit(payload)
    if config is None and "config" in payload:
        config = payload["config"]
    baseline = config_fingerprint(config) if config is not None else None
    sim, exp = raw["simulation"], raw["experiment"]
    mask = (np.asarray(exp["lag"]) >= raw["fit_min_s"]) & (np.asarray(exp["lag"]) <= raw["fit_max_s"])
    times, observed = np.asarray(exp["lag"])[mask], np.asarray(exp["msd"])[mask]
    if len(times) < 4:
        raise ValueError("At least four experimental lag points must fall inside the fitting window.")
    lower = max(TIME_BOUNDS[0], raw["fit_max_s"]/sim["lag"][-1])
    upper = min(TIME_BOUNDS[1], raw["fit_min_s"]/sim["lag"][0])
    if lower > upper*(1+1e-12):
        raise ValueError("No supported time scale represents the entire experimental fitting window within the simulated lag range; extend the trajectory or narrow the window.")
    lower = min(lower, upper)
    xs, ys = np.log(sim["lag"]), np.log(sim["msd"])+2*math.log(raw["length_nm"]/1000)
    xe, ye = np.log(times), np.log(observed)
    z, _ = _fit_log_shift(xs, ys, xe, ye, math.log(lower), math.log(upper))
    time_s = min(TIME_BOUNDS[1], max(TIME_BOUNDS[0], math.exp(z)))
    mapped_log = np.interp(xe-z, xs, ys)
    mapped = np.exp(mapped_log)
    log_rmse = float(np.sqrt(np.mean((mapped_log-ye)**2)))
    simulated_slope = float(np.polyfit(xe, mapped_log, 1)[0])
    experimental_slope = float(np.polyfit(xe, ye, 1)[0])
    difference = abs(simulated_slope-experimental_slope)
    boundary = min(abs(z-math.log(lower)), abs(z-math.log(upper))) <= 1e-7*max(1, abs(math.log(upper/lower)))
    identifiable = abs(simulated_slope) >= THRESHOLDS["minimum_absolute_simulation_slope"] and upper > lower*(1+1e-6)
    shape = difference <= THRESHOLDS["maximum_slope_difference"] and log_rmse <= THRESHOLDS["maximum_log_rmse"]
    curve_baseline = sim.get("baseline_fingerprint")
    provenance_matches = curve_baseline is None or baseline is None or curve_baseline == baseline
    diagnostics = {"fit_min_s": raw["fit_min_s"], "fit_max_s": raw["fit_max_s"], "n_points": len(times),
                   "log_rmse": log_rmse, "multiplicative_rmse_factor": math.exp(log_rmse),
                   "simulation_slope": simulated_slope, "experimental_slope": experimental_slope,
                   "slope_difference": difference, "shape_compatible": bool(shape),
                   "time_identifiable": bool(identifiable), "boundary_optimum": bool(boundary),
                   "time_bounds_s": [lower, upper], "thresholds": dict(THRESHOLDS),
                   "curve_baseline_matches": bool(provenance_matches),
                   "objective": "unweighted mean squared log-MSD residual; exact piecewise-linear interpolation; no extrapolation",
                   "uncertainty": "not estimated; overlapping lag estimates are correlated"}
    accepted = shape and identifiable and not boundary and provenance_matches
    status = "demo" if raw["source_kind"] == "synthetic" else "fitted" if accepted else "needs_review"
    warnings = ["A fitted clock is conditional on the chosen length, reference, observable and lag window; it does not validate unrelated mechanisms or biological predictions.",
                "MSD fitting does not calibrate the chosen length scale. bp per bead is not a physical length conversion.",
                "Bare bead diffusion differs from measured chromatin-locus diffusion. Apply the same clock to all kinetic rates and particle diffusivities.",
                "No confidence interval is supplied: use independent trajectories/reference replicates and sensitivity to window, imaging corrections and length scale."]
    if raw["source_kind"] == "synthetic":
        warnings.append("Synthetic reference: this is a demonstration and cannot establish experimental calibration.")
    if not shape:
        warnings.append("Curve slopes or log residuals exceed the declared compatibility thresholds. A clock rescaling alone does not explain this mismatch.")
    if not identifiable:
        warnings.append("The curve/window does not identify a reliable clock (near-flat slope or no time-scale search range).")
    if boundary:
        warnings.append("The best fit lies at the allowed time-scale boundary. Extend simulated lag coverage before accepting the mapping.")
    if not provenance_matches:
        warnings.append("The simulation curve was generated from a different model fingerprint; this profile needs review.")
    if baseline is None:
        warnings.append("No baseline model was supplied; this mapping is not tied to model parameters.")
    elif curve_baseline is None:
        warnings.append("The pasted simulation curve has no baseline fingerprint; its correspondence to the supplied model is user-declared.")
    profile = {"version": 1, "length_nm": raw["length_nm"], "time_s": time_s,
               "status": status, "source_kind": raw["source_kind"],
               "sources": {"simulation": sim["source"], "experiment": exp["source"], "length": raw["length_source"]},
               "notes": raw["notes"], "diagnostics": diagnostics, "baseline_fingerprint": baseline, "raw_input": raw}
    return {"profile": profile, "diagnostics": diagnostics,
            "comparison": {"lag_s": times.tolist(), "experimental_msd_um2": observed.tolist(),
                           "mapped_simulation_msd_um2": mapped.tolist()},
            "warnings": warnings, "equation": "MSD_physical(t) = (length_nm/1000)^2 * MSD_simulation(t/time_s)"}


def validate_profile(profile):
    """Strictly validate saved mappings; never promote an assumption/demo to a fit."""
    if not isinstance(profile, dict):
        raise ValueError("calibration must be a profile object.")
    required = {"version", "length_nm", "time_s", "status", "source_kind", "sources", "notes",
                "diagnostics", "baseline_fingerprint", "raw_input"}
    if set(profile) != required:
        raise ValueError("Calibration profile fields must be exactly: " + ", ".join(sorted(required)))
    if isinstance(profile["version"], bool) or profile["version"] != 1:
        raise ValueError("Unsupported calibration profile version.")
    result = copy.deepcopy(profile)
    result["length_nm"] = _number(profile["length_nm"], "calibration.length_nm", *LENGTH_BOUNDS)
    result["time_s"] = _number(profile["time_s"], "calibration.time_s", *TIME_BOUNDS)
    status, source_kind = profile["status"], profile["source_kind"]
    valid = {"fitted": ("experimental",), "needs_review": ("experimental",), "demo": ("synthetic",), "assumed": ("assumed",)}
    if status not in valid or source_kind not in valid[status]:
        raise ValueError("Calibration status is inconsistent with its explicitly declared source_kind.")
    sources = profile["sources"]
    if not isinstance(sources, dict) or set(sources) != {"simulation", "experiment", "length"}:
        raise ValueError("Calibration sources must include simulation, experiment and length.")
    result["sources"] = {key: _string(value, "calibration.sources."+key) for key, value in sources.items()}
    result["notes"] = _string(profile["notes"], "calibration.notes", required=False)
    result["baseline_fingerprint"] = _fingerprint(profile["baseline_fingerprint"], "calibration.baseline_fingerprint")
    if status == "assumed":
        if profile["raw_input"] is not None or profile["diagnostics"] is not None:
            raise ValueError("An assumed mapping has no fitted raw_input or diagnostics; use null for both.")
    else:
        raw = _normalize_fit(profile["raw_input"])
        if raw["source_kind"] != source_kind or raw["length_nm"] != result["length_nm"]:
            raise ValueError("Calibration raw inputs disagree with its source kind or length.")
        # Recompute the inexpensive one-dimensional fit to detect edited scale,
        # status or diagnostics, rather than trusting a pasted 'fitted' label.
        expected = fit_time_mapping(raw)["profile"]
        if not math.isclose(result["time_s"], expected["time_s"], rel_tol=1e-9, abs_tol=1e-15):
            raise ValueError("Calibration time_s does not reproduce its saved raw input fit.")
        diagnostics = profile["diagnostics"]
        if not isinstance(diagnostics, dict) or set(diagnostics) != set(expected["diagnostics"]):
            raise ValueError("Calibration diagnostics are missing fields or have unknown fields.")
        # Model correspondence cannot be recomputed without the baseline config;
        # verify the stored curve fingerprint directly when available.
        curve_hash = raw["simulation"].get("baseline_fingerprint")
        correspondence = curve_hash is None or result["baseline_fingerprint"] is None or curve_hash == result["baseline_fingerprint"]
        expected["diagnostics"]["curve_baseline_matches"] = correspondence
        for key, value in expected["diagnostics"].items():
            candidate = diagnostics[key]
            if isinstance(value, bool):
                if not isinstance(candidate, bool) or candidate != value:
                    raise ValueError(f"Calibration diagnostic {key} does not reproduce its inputs.")
            elif isinstance(value, int):
                if isinstance(candidate, bool) or not isinstance(candidate, int) or candidate != value:
                    raise ValueError(f"Calibration diagnostic {key} does not reproduce its inputs.")
            elif isinstance(value, float):
                if not math.isclose(_number(candidate, "diagnostics."+key, -1e30, 1e30), value, rel_tol=1e-8, abs_tol=1e-10):
                    raise ValueError(f"Calibration diagnostic {key} does not reproduce its inputs.")
            elif candidate != value:
                raise ValueError(f"Calibration diagnostic {key} does not reproduce its inputs.")
        expected_status = expected["status"]
        if source_kind == "experimental" and not correspondence:
            expected_status = "needs_review"
        if status != expected_status:
            raise ValueError("Calibration status does not agree with its reproduced fit diagnostics.")
        result["raw_input"] = raw
        result["diagnostics"] = copy.deepcopy(expected["diagnostics"])
        if result["sources"] != {"simulation": raw["simulation"]["source"], "experiment": raw["experiment"]["source"], "length": raw["length_source"]}:
            raise ValueError("Calibration sources disagree with its saved raw inputs.")
    return result


def unit_summary(config, profile=None):
    """Translate parameters without modifying rates, stiffness or dynamics.

    A changed model retains the baseline clock for controlled perturbations;
    it is flagged for interpretation rather than automatically refitted.
    """
    if profile is None:
        return {"status": "uncalibrated", "baseline_matches": False, "length_nm": None,
                "time_s": None, "warnings": ["No physical unit mapping is attached; model values remain in reduced units."]}
    p = validate_profile(profile)
    length, clock = p["length_nm"], p["time_s"]
    diffusion = (length/1000)**2/clock
    matches = p["baseline_fingerprint"] is not None and p["baseline_fingerprint"] == config_fingerprint(config)
    warnings = []
    if not matches:
        warnings.append("Current model differs from the fitted baseline or has no baseline fingerprint. The baseline clock is retained for perturbation comparisons; it has not been refitted.")
    if p["status"] != "fitted":
        warnings.append("This mapping is " + p["status"] + "; it is not an accepted experimental fit.")
    warnings.append("Converted rates are encounter hazards, not molar association constants or validated Kd values; physical mapping does not validate the mechanisms.")
    chrom = config.get("chromatin", {})
    bp = float(chrom.get("bp_per_bead", 1))
    def rates(obj):
        return {key+"_s_inv": float(value)/clock for key, value in obj.items()
                if key.endswith("_rate") and isinstance(value, (int, float)) and not isinstance(value, bool)
                and key not in ("elongation_rate", "step_rate")}
    return {"status": p["status"], "source_kind": p["source_kind"], "baseline_matches": matches,
            "length_nm": length, "time_s": clock, "warnings": warnings,
            "conversions": {"length_nm_per_unit": length, "time_s_per_unit": clock,
                            "msd_um2_per_unit": (length/1000)**2, "diffusion_um2_s_per_unit": diffusion,
                            "hazard_s_inv_per_unit": 1/clock, "spring_kbt_per_nm2_per_unit": 1/length**2},
            "simulation": {"dt_s": float(config.get("simulation", {}).get("dt", 0))*clock,
                           "confinement_radius_nm": float(config.get("simulation", {}).get("confinement_radius", 0))*length},
            "chromatin": {"bond_length_nm": float(chrom.get("bond_length", 0))*length,
                          "persistence_length_nm": float(chrom.get("persistence_length", 0))*length,
                          "bare_bead_diffusion_um2_s": float(chrom.get("diffusion", 0))*diffusion,
                          "nucleosome_rates": rates(chrom.get("nucleosomes", {}))},
            "proteins": [{"id": obj["id"], "free_molecule_diffusion_um2_s": float(obj.get("diffusion", 0))*diffusion,
                          "radius_nm": float(obj.get("radius", 0))*length, **rates(obj)} for obj in config.get("proteins", [])],
            "genes": [{"id": obj["id"], "elongation_nt_s": float(obj.get("elongation_rate", 0))*bp/clock,
                       **rates(obj)} for obj in config.get("genes", [])],
            "cohesin": {"step_rate_beads_s": float(config.get("cohesin", {}).get("step_rate", 0))/clock,
                        "step_rate_bp_s": float(config.get("cohesin", {}).get("step_rate", 0))*bp/clock,
                        **rates(config.get("cohesin", {}))},
            "transcription": {key+"_um2_s": float(config.get("transcription", {}).get(key, default))*diffusion
                              for key, default in (("rna_diffusion", .3), ("mature_rna_diffusion", .15))},
            "marks": [{"id": obj["id"], **rates(obj)} for obj in config.get("marks", [])],
            "interactions": [{"a": obj["a"], "b": obj["b"], **rates(obj)} for obj in config.get("interactions", [])]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    msd = subparsers.add_parser("msd", help="Measure chromatin MSD from a batch snapshots.jsonl file.")
    msd.add_argument("--snapshots", type=Path, required=True)
    msd.add_argument("--beads", required=True, help="Comma-separated zero-based chromatin bead indices.")
    msd.add_argument("--dimension", type=int, choices=(2, 3), default=3)
    msd.add_argument("--observable", choices=OBSERVABLES, default="locus")
    msd.add_argument("--drift-correction", choices=("none", "center_of_mass"), default="none")
    msd.add_argument("--discard-before", type=float, default=0)
    msd.add_argument("--max-lag-steps", type=int)
    msd.add_argument("--drop-final-irregular", action="store_true")
    msd.add_argument("--model", type=Path, help="Baseline config; otherwise use model.json beside snapshots when present.")
    msd.add_argument("--source", help="Provenance description; default includes the supplied trajectory path.")
    msd.add_argument("--output", type=Path, required=True)
    fit = subparsers.add_parser("fit", help="Fit a clock from a saved JSON fit request.")
    fit.add_argument("--input", type=Path, required=True)
    fit.add_argument("--model", type=Path)
    fit.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("Output exists; choose a new path to preserve previous analyses.")
    try:
        if args.command == "msd":
            if args.snapshots.stat().st_size > 256*1024*1024:
                raise ValueError("Snapshot file exceeds 256 MiB; select a smaller recorded window for this analysis.")
            payload = {"frames": [json.loads(line) for line in args.snapshots.read_text().splitlines() if line.strip()],
                       "beads": [int(x) for x in args.beads.split(",")], "dimension": args.dimension,
                       "observable": args.observable, "drift_correction": args.drift_correction,
                       "discard_before": args.discard_before, "drop_final_irregular": args.drop_final_irregular,
                       "source": args.source or "Simulation trajectory: " + str(args.snapshots)}
            if args.max_lag_steps is not None:
                payload["max_lag_steps"] = args.max_lag_steps
            model = args.model or args.snapshots.parent/"model.json"
            if model.exists():
                payload["config"] = json.loads(model.read_text())
            elif args.model is not None:
                raise ValueError("Explicit baseline model file does not exist.")
            result = analyze_trajectory(payload)
        else:
            payload = json.loads(args.input.read_text())
            result = fit_time_mapping(payload, json.loads(args.model.read_text()) if args.model else None)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, allow_nan=False)
        print("Saved " + str(args.output))
        return 0
    except (ValueError, OSError, TypeError) as exc:
        parser.exit(2, f"Calibration analysis failed: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
