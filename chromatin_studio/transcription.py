"""Individual polymerases and passive RNA particles in reduced units.

This optional module adds molecular identities to the count-based model.
Nascent RNA is attached to its polymerase; released RNA is a point particle,
not an RNA chain. The same soft harmonic nuclear wall used for chromatin
confines RNA. Processing depends on the RNA position, with a piecewise-frozen
local hazard during each Brownian/reaction substep. This is an operator-split
spatial model, not an exact simulation of diffusion across a reactive boundary.
"""

import math

import numpy as np

MAX_TRANSCRIPTION_EVENTS = 10000
MAX_TRANSCRIPT_UPDATES = 1000000


def sample_promoter_interval(rng, initially_on, on_rate, off_rate,
                             initiation_rate, dt, max_events=MAX_TRANSCRIPTION_EVENTS):
    """Exact telegraph-path and initiation times for rates frozen over dt.

    Exponential switching times define the actual ON intervals. Conditional
    Poisson counts and uniform event locations sample initiation within those
    intervals. This is a two-state Markov-modulated Poisson process, rather
    than assigning a full interval's production to the final promoter state.
    See Gillespie (1977), doi:10.1021/j100540a008, and Peccoud & Ycart (1995),
    doi:10.1006/tpbi.1995.1027. Exactness is conditional on the frozen rates;
    coupling to moving chromatin remains operator split.
    """
    on = bool(initially_on)
    if any(not math.isfinite(value) or value < 0 for value in (on_rate, off_rate, initiation_rate, dt)):
        raise ValueError("Promoter rates and interval duration must be finite and nonnegative.")
    elapsed = 0.0
    occupied_time = 0.0
    switches, initiations = [], []
    while elapsed < dt:
        rate = float(off_rate if on else on_rate)
        wait = float(rng.exponential(1.0 / rate)) if rate > 0 else math.inf
        duration = min(wait, dt - elapsed)
        if on and duration > 0:
            occupied_time += duration
            mean = float(initiation_rate) * duration
            if not math.isfinite(mean) or mean > max_events:
                raise ValueError("Too many expected initiation events in one interval; reduce dt, regulatory gains, or body response multipliers.")
            count = int(rng.poisson(mean)) if initiation_rate > 0 else 0
            if len(initiations) + count + len(switches) > max_events:
                raise ValueError("Too many promoter/initiation events in one interval; reduce dt or regulatory gains.")
            if count:
                times = elapsed + duration * np.asarray(rng.random(count))
                initiations.extend(sorted(times.tolist()))
        elapsed += duration
        if wait < duration + 1e-14 and elapsed < dt:
            on = not on
            switches.append((elapsed, on))
            if len(switches) + len(initiations) > max_events:
                raise ValueError("Too many promoter/initiation events in one interval; reduce dt or switching rates.")
    return {"on": on, "on_time": occupied_time, "switches": switches,
            "initiation_times": initiations}


_TOTAL_KEYS = ("live", "nascent", "unprocessed", "mature", "initiated",
               "completed", "degraded", "blocked_initiations")


def disabled_transcription_state():
    """Disabled means individual identities/positions are not being simulated."""
    return {"enabled": False, "cap": 0,
            "promoter_timing": "within_interval_ctmc",
            "rate_coupling": "regulatory and mechanical state frozen during each reaction interval",
            "totals": dict.fromkeys(_TOTAL_KEYS, 0), "transcripts": [], "genes": []}


class ExplicitTranscription:
    """Own live transcripts and synchronize the existing per-gene summaries."""

    def __init__(self, sim):
        self.sim = sim
        self.config = sim.config.get("transcription", {})
        self.cap = int(self.config.get("max_transcripts", 2000))
        self.transcripts = {}
        self.serial = 0
        self.degraded = {g["config"]["id"]: 0 for g in sim.genes}
        self.blocked = {g["config"]["id"]: 0 for g in sim.genes}
        self.gene_by_id = {g["config"]["id"]: g for g in sim.genes}

    def _coordinate(self, bead_position, spec):
        # Each bead represents a genomic interval. The final interval is
        # represented spatially by the final gene bead, including release.
        p = float(np.clip(bead_position, int(spec["start"]), int(spec["end"])))
        left = int(math.floor(p))
        right = min(left + 1, int(spec["end"]))
        return ((1 - (p - left)) * self.sim.positions[left]
                + (p - left) * self.sim.positions[right]).copy()

    def _speed(self, spec, position):
        speed = float(spec.get("elongation_rate", 0))
        if getattr(self.sim, "nucleosomes_enabled", False):
            bead = min(int(math.floor(position)), int(spec["end"]))
            if self.sim.nucleosome_occupied[bead]:
                speed *= float(spec.get("occupied_elongation_multiplier", 1))
        return speed

    def _new_nascent(self, gene, time):
        self.serial += 1
        spec = gene["config"]
        record = {"id": f"rna-{self.serial}", "gene_id": spec["id"],
                  "polymerase_id": f"pol-{self.serial}", "state": "nascent",
                  "bead_position": float(spec["start"]), "length_nt": 0.0,
                  "position": self._coordinate(spec["start"], spec),
                  "paused": False, "initiated_time": float(time),
                  "released_time": None, "processed_time": None}
        self.transcripts[record["id"]] = record
        gene["initiated"] += 1
        return record

    def _move_polymerase(self, record, spec, active_time):
        """Traverse each bead at its own speed; return time consumed to finish."""
        remaining = float(active_time)
        end = int(spec["end"]) + 1
        while remaining > 0 and record["bead_position"] < end:
            p = record["bead_position"]
            speed = self._speed(spec, p)
            if speed <= 0:
                return active_time
            boundary = min(math.floor(p) + 1, end)
            needed = (boundary - p) / speed
            if needed <= remaining:
                record["bead_position"] = float(boundary)
                remaining -= needed
            else:
                record["bead_position"] = p + speed * remaining
                remaining = 0.0
        return active_time - remaining

    def _advance_nascent(self, record, spec, dt, start_time):
        """Exact pause/resume events and piecewise-constant elongation speeds."""
        elapsed = 0.0
        end = int(spec["end"]) + 1
        while elapsed < dt:
            rate = float(spec.get("resume_rate", 1) if record["paused"]
                         else spec.get("pause_rate", 0))
            wait = float(self.sim.rng.exponential(1 / rate)) if rate > 0 else math.inf
            segment = min(wait, dt - elapsed)
            if record["paused"]:
                elapsed += segment
            else:
                elapsed += self._move_polymerase(record, spec, segment)
                if record["bead_position"] >= end:
                    record["state"] = "unprocessed"
                    record["paused"] = False
                    record["released_time"] = start_time + elapsed
                    record["position"] = self._coordinate(end, spec)
                    record["length_nt"] = ((end - int(spec["start"]))
                                            * float(self.sim.chrom.get("bp_per_bead", 1000)))
                    record["bead_position"] = None
                    self.gene_by_id[spec["id"]]["completed"] += 1
                    self._advance_free(record, spec, dt - elapsed, start_time + elapsed)
                    return
            if wait < segment + 1e-14 and elapsed < dt:
                record["paused"] = not record["paused"]
        record["length_nt"] = ((record["bead_position"] - int(spec["start"]))
                               * float(self.sim.chrom.get("bp_per_bead", 1000)))
        record["position"] = self._coordinate(record["bead_position"], spec)

    def _processing_rate(self, record, spec):
        rate = float(spec.get("processing_rate", 0))
        # A released RNA must itself enter a speckle. Gene contact alone does
        # not change this rate. Overlapping speckles do not multiply gains.
        for body in self.sim.bodies:
            if (body.get("kind") == "speckle"
                    and np.linalg.norm(record["position"] - np.asarray(body["center"]))
                    <= float(body["radius"])):
                return rate * (1 + float(spec.get("speckle_gain", 0)))
        return rate

    def _diffusion(self, record):
        key = "mature_rna_diffusion" if record["state"] == "mature" else "rna_diffusion"
        return float(self.config.get(key, 0.15 if record["state"] == "mature" else 0.3))

    def _brownian(self, record, dt):
        diffusion = self._diffusion(record)
        if dt <= 0 or diffusion <= 0:
            return
        p = record["position"]
        radial = float(np.linalg.norm(p))
        stiffness = float(getattr(self.sim, "confinement_stiffness", 30.0))
        force = (-stiffness * (radial - self.sim.radius) / radial * p
                 if radial > self.sim.radius else np.zeros(3))
        nucleus = getattr(self.sim, "nucleus", None)
        if nucleus is not None:
            force = force + nucleus.force_energy(p[None, :])[0][0]
        record["position"] = (p + diffusion * force * dt
                              + math.sqrt(2 * diffusion * dt) * self.sim.rng.normal(size=3))
        if not np.all(np.isfinite(record["position"])):
            raise FloatingPointError("Non-finite RNA position. Reduce the time step or RNA diffusion.")

    def _advance_free(self, record, spec, dt, start_time):
        """Competing processing/decay with no repeated full-step mature decay."""
        max_diffusion = max(float(self.config.get("rna_diffusion", 0.3)),
                            float(self.config.get("mature_rna_diffusion", 0.15)))
        stiffness = float(getattr(self.sim, "confinement_stiffness", 30.0))
        if getattr(self.sim, "nucleus", None) is not None:
            stiffness += self.sim.nucleus.curvature_bound()
        substeps = max(1, int(math.ceil(dt * max_diffusion * stiffness / 0.35)))
        if substeps > 1000:
            raise ValueError("RNA diffusion requires more than 1000 integration substeps; reduce dt.")
        h = dt / substeps
        for substep in range(substeps):
            remaining = h
            elapsed = substep * h
            while remaining > 0:
                processing = self._processing_rate(record, spec) if record["state"] == "unprocessed" else 0.0
                decay = float(spec.get("rna_decay_rate", 0))
                total = processing + decay
                wait = float(self.sim.rng.exponential(1 / total)) if total > 0 else math.inf
                interval = min(remaining, wait)
                self._brownian(record, interval)
                elapsed += interval
                remaining -= interval
                if wait >= interval + 1e-14:
                    continue
                # An event exactly at the interval endpoint is harmless;
                # zero waiting times have probability zero for positive rates.
                if total > 0 and wait <= interval:
                    if processing and self.sim.rng.random() < processing / total:
                        record["state"] = "mature"
                        record["processed_time"] = start_time + elapsed
                    else:
                        self.degraded[spec["id"]] += 1
                        del self.transcripts[record["id"]]
                        return

    def _synchronize(self):
        by_gene = {key: [] for key in self.gene_by_id}
        for record in self.transcripts.values():
            by_gene[record["gene_id"]].append(record)
        for key, gene in self.gene_by_id.items():
            records = by_gene[key]
            gene["polymerases"] = [r["bead_position"] for r in records if r["state"] == "nascent"]
            gene["unprocessed"] = sum(r["state"] == "unprocessed" for r in records)
            gene["mature"] = sum(r["state"] == "mature" for r in records)

    def advance(self, dt, rate_reports):
        """Sample ON occupancy and process actual birth times chronologically.

        Molecules move for the time elapsed since their birth, including any
        release and processing in that same interval. Capacity is checked at
        each attempted birth after existing molecules have advanced to its
        time; a later degradation cannot prematurely free an earlier slot.
        Independent continuous birth times remove gene-order priority. No
        existing molecule is deleted to make space.
        """
        start = float(self.sim.time)
        births, switches, plans = [], [], []
        budget = MAX_TRANSCRIPTION_EVENTS
        previous_completed = {g["config"]["id"]: g["completed"] for g in self.sim.genes}
        for gene in self.sim.genes:
            spec = gene["config"]
            rates = rate_reports[spec["id"]]
            plan = sample_promoter_interval(self.sim.rng, gene["on"], float(rates["effective_on_rate"]),
                                            float(spec.get("off_rate", 0)),
                                            float(rates["effective_initiation_rate"]), dt, max_events=budget)
            budget -= len(plan["switches"]) + len(plan["initiation_times"])
            plans.append(plan)
            births.extend((time, spec["id"]) for time in plan["initiation_times"])
            switches.extend((time, spec["id"], on) for time, on in plan["switches"])
        estimated_updates = (len(births) + 1) * min(self.cap, len(self.transcripts) + len(births))
        if estimated_updates > MAX_TRANSCRIPT_UPDATES:
            raise ValueError("Too many transcript updates in one interval; reduce dt, initiation rates, or regulatory/body response gains.")
        for gene, plan in zip(self.sim.genes, plans):
            gene["on"] = plan["on"]
            gene["last_transcription_interval"] = {
                "start_time": start, "duration": float(dt), "on_time": plan["on_time"],
                "on_fraction": plan["on_time"] / dt if dt else 0.0,
                "switch_count": len(plan["switches"]), "initiation_attempts": len(plan["initiation_times"]),
                "accepted_initiations": 0, "blocked_initiations": 0}
        for time, gene_id, on in sorted(switches):
            self.sim._event("promoter", gene_id + (" ON" if on else " OFF"), time=start + time)

        def advance_molecules(duration, interval_start):
            if duration <= 0:
                return
            for record in list(self.transcripts.values()):
                spec = self.gene_by_id[record["gene_id"]]["config"]
                if record["state"] == "nascent":
                    self._advance_nascent(record, spec, duration, interval_start)
                else:
                    self._advance_free(record, spec, duration, interval_start)

        elapsed = 0.0
        for birth_time, gene_id in sorted(births):
            advance_molecules(birth_time - elapsed, start + elapsed)
            elapsed = birth_time
            gene = self.gene_by_id[gene_id]
            if len(self.transcripts) < self.cap:
                self._new_nascent(gene, start + birth_time)
                gene["last_transcription_interval"]["accepted_initiations"] += 1
            else:
                self.blocked[gene_id] += 1
                gene["last_transcription_interval"]["blocked_initiations"] += 1
                warning = "Transcript capacity has rejected initiation attempts; increase max_transcripts or reduce production before interpreting RNA abundance."
                if warning not in self.sim.warnings:
                    self.sim.warnings.append(warning)
        advance_molecules(dt - elapsed, start + elapsed)
        for gene in self.sim.genes:
            gene_id = gene["config"]["id"]
            completed = gene["completed"] - previous_completed[gene_id]
            if completed:
                self.sim._event("transcription", f"{gene_id}: {completed} completed", time=start + dt)
        self._synchronize()

    def snapshot(self):
        records = []
        paused_by_gene = dict.fromkeys(self.gene_by_id, 0)
        for value in self.transcripts.values():
            record = dict(value)
            spec = self.gene_by_id[record["gene_id"]]["config"]
            if record["state"] == "nascent":
                record["position"] = self._coordinate(record["bead_position"], spec)
                record["elongation_rate"] = 0.0 if record["paused"] else self._speed(spec, record["bead_position"])
                paused_by_gene[record["gene_id"]] += bool(record["paused"])
            else:
                record["elongation_rate"] = 0.0
            record["position"] = record["position"].tolist()
            records.append(record)
        genes = []
        for gene in self.sim.genes:
            spec = gene["config"]
            genes.append({"id": spec["id"], "nascent": len(gene["polymerases"]),
                          "unprocessed": gene["unprocessed"], "mature": gene["mature"],
                          "initiated": gene["initiated"], "completed": gene["completed"],
                          "degraded": self.degraded[spec["id"]], "blocked_initiations": self.blocked[spec["id"]],
                          "paused_polymerases": paused_by_gene[spec["id"]],
                          "elongation_rate": float(spec.get("elongation_rate", 0)),
                          "occupied_elongation_multiplier": float(spec.get("occupied_elongation_multiplier", 1)),
                          "pause_rate": float(spec.get("pause_rate", 0)), "resume_rate": float(spec.get("resume_rate", 1))})
        totals = {key: sum(g[key] for g in genes) for key in _TOTAL_KEYS if key != "live"}
        totals["live"] = len(records)
        return {"enabled": True, "cap": self.cap, "totals": totals, "transcripts": records, "genes": genes,
                "promoter_timing": "within_interval_ctmc",
                "rate_coupling": "regulatory and mechanical state frozen during each reaction interval"}
