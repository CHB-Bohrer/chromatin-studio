"""Shared helpers for component studies. Study modules import only from here.

A component study isolates one mechanism (bare polymer, one diffusing protein,
one binding species, cohesin alone, one gene ...) and compares the engine with
a reference. Three levels of evidence are kept apart in every record:

- ``numerical_check``: the engine reproduces its own stated rule.
- ``analytical_reference``: the engine reproduces a closed-form result.
- ``experimental_reference``: an engine output is compared with a published
  measurement. Passing means "inside the stated experimental range under the
  stated unit mapping", not that the model is calibrated or validated.

Nothing here mutates a simulation or draws from its random stream.
"""

import copy
import math

import numpy as np

from ..schema import checked_config, default_config

NUMERICAL = "numerical_check"
ANALYTICAL = "analytical_reference"
EXPERIMENTAL = "experimental_reference"
EVIDENCE_LEVELS = (NUMERICAL, ANALYTICAL, EXPERIMENTAL)

# Provisional mapping of the 2,000 x 10 kb polymer (status needs_review; see
# docs/PROJECT_STATUS.md). Studies that use it must say so in their report.
PROVISIONAL_CLOCK = {"length_nm_per_unit": 47.2, "time_s_per_unit": 0.205,
                     "time_s_per_unit_range": [0.13, 0.32], "status": "needs_review",
                     "scope": "2,000 beads x 10 kb Fbn2 baseline mechanics only"}


def bare_config(beads, *, seed=1, dt=0.002, radius=60.0, name="component study", **chromatin):
    """A validated model with chromatin only; every other mechanism is absent."""
    config = default_config()
    config.update(name=name, seed=int(seed), marks=[], proteins=[], interactions=[], bodies=[], genes=[])
    config["cohesin"] = {"enabled": False, "count": 0, "step_rate": 0, "off_rate": 0, "spring_k": 0}
    config["simulation"].update(dt=dt, confinement_radius=radius)
    config["chromatin"].update(beads=int(beads), initial_marks=[], barriers=[])
    config["chromatin"].update(chromatin)
    return config


def validated(config):
    """Schema-check a study model so studies never run an unsupported config."""
    return checked_config(copy.deepcopy(config))


def protein(identifier, count, **fields):
    """One protein species with every optional activity switched off."""
    value = dict(id=identifier, name=identifier, color="#d6a34a", count=int(count), diffusion=0.5,
                 radius=0.28, valence=1, representation="bead", idr_beads=4, sticker_fraction=1.0,
                 chromatin_affinity={"unmodified": 0.0}, on_rate=0.0, off_rate=0.0,
                 writer_mark=None, eraser_mark=None, catalytic_rate=0.0)
    value.update(fields)
    return value


def chemistry_only_step(sim, steps, interval):
    """Advance reactions with frozen coordinates (for purely 1D mechanisms).

    Uses the engine's public reaction/bookkeeping hooks, the same ones the
    optional dynamics backends call. Positions do not move, so this is valid
    only for observables that do not depend on 3D motion. ``interval`` is the
    reaction interval; it may exceed the mechanical time step because no
    mechanics is integrated. Rates are frozen within each interval, as always.
    """
    for _ in range(int(steps)):
        sim.advance_reactions(float(interval))
        sim.advance_bookkeeping(float(interval))


def mean_sem(values):
    """Mean and standard error over independent values (replicas)."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return None, None
    sem = float(np.std(values, ddof=1) / math.sqrt(values.size)) if values.size > 1 else None
    return float(np.mean(values)), sem


def check(identifier, label, observed, expected, *, evidence, tolerance=None, low=None, high=None,
          sem=None, unit="dimensionless", note=None):
    """One comparison record.

    Pass rule: inside [low, high] when an interval is given; otherwise
    |observed - expected| <= max(relative tolerance x |expected|, 3 x sem).
    The statistical allowance is stated in the record so a noisy pass is visible.
    """
    if evidence not in EVIDENCE_LEVELS:
        raise ValueError(f"Unknown evidence level: {evidence}")
    record = {"id": identifier, "label": label, "evidence": evidence, "unit": unit,
              "observed": None if observed is None else float(observed),
              "expected": None if expected is None else float(expected),
              "sem": None if sem is None else float(sem), "note": note}
    if observed is None or not math.isfinite(float(observed)):
        record.update(passed=False, rule="no finite observation")
        return record
    if low is not None or high is not None:
        lo = -math.inf if low is None else float(low)
        hi = math.inf if high is None else float(high)
        record.update(low=low, high=high, rule="inside stated interval",
                      passed=bool(lo <= float(observed) <= hi))
        return record
    allowed = abs(float(expected)) * float(tolerance if tolerance is not None else 0.05)
    statistical = 3.0 * float(sem) if sem else 0.0
    deviation = abs(float(observed) - float(expected))
    record.update(tolerance=tolerance, allowed=max(allowed, statistical), deviation=deviation,
                  relative_deviation=deviation / abs(float(expected)) if expected else None,
                  rule="|observed - expected| <= max(tolerance x |expected|, 3 x sem)",
                  passed=bool(deviation <= max(allowed, statistical)))
    return record


def source(citation, url, measured, **fields):
    """A literature source with the quantity it reports. Keep values verifiable."""
    return {"citation": citation, "url": url, "measured": measured, **fields}


def summarize(checks):
    by_level = {}
    for level in EVIDENCE_LEVELS:
        chosen = [c for c in checks if c["evidence"] == level]
        by_level[level] = {"total": len(chosen), "passed": sum(bool(c["passed"]) for c in chosen)}
    return {"total": len(checks), "passed": sum(bool(c["passed"]) for c in checks), "by_evidence": by_level}


# --- Figures: data only. The interface draws them with one generic renderer,
# so a study changes its plots by editing its own analyze() and nothing else.
# Roles fix the visual encoding: "model" = engine output, "reference" = a
# closed-form or input value, "experiment" = a published measurement.
ROLES = ("model", "reference", "experiment")
SERIES_STYLES = ("points", "connected", "line", "step")   # connected = markers joined by a line
# The interface has three hues that stay apart for every pair under color-vision
# deficiency (validated against its chart surface). Keep the hues one figure
# needs to this many; split a larger comparison into several figures.
MAX_COLORED_GROUPS = 3


def _finite_or_none(value):
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def series(label, x, y, *, role="model", low=None, high=None, style="points", group=None):
    """One plotted series. ``x`` holds numbers, or category names in a category
    figure. ``low``/``high`` are absolute error-bar ends per point; by
    convention they are +/-1 SEM over replicas (say otherwise in the note).
    Series that share a ``group`` (an engine series and the reference or
    measurement it is compared with) share a color; role sets the mark."""
    if role not in ROLES:
        raise ValueError(f"Unknown series role: {role}")
    if style not in SERIES_STYLES:
        raise ValueError(f"Unknown series style: {style}")
    xs = [value if isinstance(value, str) else _finite_or_none(value) for value in x]
    ys = [_finite_or_none(value) for value in y]
    if len(xs) != len(ys):
        raise ValueError(f"Series {label!r}: x and y differ in length.")
    record = {"label": label, "role": role, "style": style, "group": None if group is None else str(group),
              "x": xs, "y": ys}
    for key, values in (("low", low), ("high", high)):
        if values is not None:
            values = [_finite_or_none(value) for value in values]
            if len(values) != len(ys):
                raise ValueError(f"Series {label!r}: {key} differs in length.")
            record[key] = values
    return record


def band(label, *, role="experiment", x=None, y=None, category=None):
    """A shaded range. ``y=(low, high)`` spans the plot horizontally,
    ``x=(low, high)`` vertically, both give a box. In a category figure,
    ``category`` limits the range to that category's column."""
    if role not in ROLES:
        raise ValueError(f"Unknown band role: {role}")
    if x is None and y is None:
        raise ValueError(f"Band {label!r} needs x or y limits.")
    record = {"label": label, "role": role, "category": category}
    for key, limits in (("x", x), ("y", y)):
        record[key] = None if limits is None else [_finite_or_none(limits[0]), _finite_or_none(limits[1])]
    return record


def figure(identifier, title, *, x_label, y_label, series=(), bands=(), x_scale="linear", y_scale="linear",
           categories=None, note=None):
    """A chart description for the interface. ``categories`` makes the x axis categorical."""
    for scale in (x_scale, y_scale):
        if scale not in ("linear", "log"):
            raise ValueError(f"Unknown axis scale: {scale}")
    if categories is not None and x_scale != "linear":
        raise ValueError("A category axis cannot be logarithmic.")
    return {"id": identifier, "title": title, "kind": "xy" if categories is None else "category",
            "categories": None if categories is None else [str(value) for value in categories],
            "x_label": x_label, "y_label": y_label, "x_scale": x_scale, "y_scale": y_scale,
            "series": list(series), "bands": list(bands), "note": note}


def colored_groups(record):
    """Hues a figure needs in the interface: one per group that holds an engine
    series, one per ungrouped engine series. Measurements and references
    outside such a group are drawn in neutral colors."""
    return {("g:" + s["group"]) if s["group"] is not None else ("s:" + s["label"])
            for s in record["series"] if s["role"] == "model"}


def error_bar(mean, sem):
    """(low, high) of a +/-1 SEM error bar; both None when the mean is missing."""
    if mean is None:
        return None, None
    return (mean - sem, mean + sem) if sem else (mean, mean)


def survival(samples, points):
    """Empirical survival P(T > t) of dwell samples at the given times."""
    values = np.sort(np.asarray(samples, dtype=float))
    if not values.size:
        return [None] * len(points)
    return [float(1.0 - np.searchsorted(values, t, side="right") / values.size) for t in points]


class EventReader:
    """Return engine genomic events newer than the last call, in order.

    The engine keeps a bounded buffer; ``dropped`` counts events that were lost
    between two calls, which would bias kinetics. Call it often enough to keep
    ``dropped`` at zero. Reading events never touches the model state.
    """

    def __init__(self):
        self.last_id = 0
        self.dropped = 0

    def fresh(self, sim):
        events = []
        for event in reversed(sim.genomic_events):
            if event["id"] <= self.last_id:
                break
            events.append(event)
        if events:
            if events[-1]["id"] != self.last_id + 1:
                self.dropped += events[-1]["id"] - self.last_id - 1
            self.last_id = events[0]["id"]
        events.reverse()
        return events
