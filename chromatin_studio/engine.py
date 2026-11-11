"""Small, explicit Brownian/reaction reference model in reduced units (kBT=1).

This module intentionally keeps mechanics separate from chemistry.  An optional
``dynamics(sim)`` callable may replace the Brownian move while using the same
reaction rules.  Every force is the negative gradient of ``potential_energy``.
Binding kinetics are phenomenological, not a detailed-balance calibrated model.
"""

from collections import deque
import copy
import math
import uuid

import numpy as np
from scipy.spatial import cKDTree

from .neighbors import GrowingNeighbors, neighbor_pairs


def _probability(rate, dt):
    return -np.expm1(-np.asarray(rate) * dt)


def _run_length_encode(values):
    """Encode a 1D array as full-coverage [start, end-exclusive, value] runs."""
    values = np.asarray(values)
    if not values.size:
        return []
    boundaries = np.concatenate(([0], np.flatnonzero(values[1:] != values[:-1]) + 1, [len(values)]))
    return [[int(start), int(end), values[start].item()]
            for start, end in zip(boundaries[:-1], boundaries[1:])]


class Simulation:
    """Explicit chromatin, coarse proteins, reactions, extrusion and transcription.

    Protein ``diffusion`` denotes the dilute free molecule COM diffusivity. An
    n-bead IDR therefore has bead diffusivity n*D, giving COM diffusivity D under
    independent drag. All particle indices are global, chromatin first.
    """

    def __init__(self, config, backend="numpy", dynamics=None):
        if backend != "numpy" and dynamics is None:
            raise ValueError("A non-NumPy backend requires an explicit dynamics adapter.")
        self.config = copy.deepcopy(config)
        self.rng = np.random.default_rng(int(config.get("seed", 0)))
        self.dt = float(config["simulation"]["dt"])
        if not math.isfinite(self.dt) or self.dt <= 0:
            raise ValueError("The integration time step must be positive and finite.")
        self.radius = float(config["simulation"]["confinement_radius"])
        self.confinement_stiffness = float(config["simulation"].get("confinement_stiffness", 30.0))
        if not math.isfinite(self.confinement_stiffness) or self.confinement_stiffness < 0:
            raise ValueError("Confinement stiffness must be finite and nonnegative.")
        self.n_chromatin = int(config["chromatin"]["beads"])
        self.chrom = config["chromatin"]
        self.bond_length = float(self.chrom["bond_length"])
        self.engine = backend
        self.dynamics = dynamics
        self.time = 0.0
        self.step_count = 0
        self.events = deque(maxlen=80)
        # This observation identity never consumes the model random stream.
        self.run_id = str(uuid.uuid4())
        self.genomic_events = deque(maxlen=2048)
        self.genomic_event_count = 0
        self.warnings = [
            "Reduced units: time and length have not been calibrated to experimental measurements.",
            "Independent drag omits hydrodynamics; soft repulsion permits occasional overlaps and chain crossings.",
            "Binding and catalysis are phenomenological kinetics, not a detailed-balance calibrated equilibrium model.",
            "Speckles and nucleoli are imposed stationary affinity fields; protein clusters assemble separately.",
            "IDRs are coarse domain chains, not amino-acid-resolved proteins; their input D is free molecule COM diffusion.",
        ]
        if config.get("calibration"):
            self.warnings[0] = "A physical mapping is attached (status: " + config["calibration"]["status"] + "). Coordinates and times remain reduced; inspect its baseline and fit diagnostics before biological interpretation."
        self.species = {p["id"]: p for p in config.get("proteins", [])}
        self.mark_config = {m["id"]: m for m in config.get("marks", [])}
        self.marks = {key: np.zeros(self.n_chromatin, dtype=bool) for key in self.mark_config}
        for region in self.chrom.get("initial_marks", []):
            self.marks[region["mark"]][int(region["start"]):int(region["end"]) + 1] = True
        self.bodies = [b for b in config.get("bodies", []) if b.get("enabled", True)]
        # Optional whole-nucleus envelope and tethering elements. Disabled
        # (absent) by default, so legacy models keep identical forces and RNG use.
        nucleus_config = config.get("nucleus", {})
        self.nucleus = None
        self.tether_specs = []
        if nucleus_config.get("enabled"):
            from .nucleus import NucleusEnvelope
            self.nucleus = NucleusEnvelope(nucleus_config)
            bodies_by_id = {b["id"]: b for b in self.bodies}
            for spec in nucleus_config.get("tethers", []):
                if not spec.get("enabled", True):
                    continue
                if spec.get("anchor", "envelope") == "body" and spec.get("body_id") not in bodies_by_id:
                    continue        # anchor field disabled: element inactive, never silently retargeted
                self.tether_specs.append(spec)
        self.tethers = []           # attached: {"bead", "tether", "anchor"}
        self._tether_owner = {}
        for index, spec in enumerate(self.tether_specs):
            for bead in range(int(spec["start"]), int(spec["end"]) + 1):
                self._tether_owner.setdefault(bead, index)
        self._initialize_particles()
        self.bindings = []
        self.interactions = {tuple(sorted((r["a"], r["b"]))): r for r in config.get("interactions", [])}
        self.barriers = {}
        for barrier in self.chrom.get("barriers", []):
            self.barriers.setdefault(int(barrier["bead"]), set()).add(barrier["direction"])
        self.cohesin_config = config.get("cohesin", {})
        self.cohesin = []
        if self.cohesin_config.get("enabled"):
            count = int(self.cohesin_config.get("count", 0))
            if not 0 <= count <= self.n_chromatin // 2:
                raise ValueError("The cohesin count cannot exceed half the chromatin bead count.")
            # Uniformly sample a matching of k adjacent dimers. Sampling starts
            # on the compressed lattice and expanding them avoids greedy
            # dead ends that silently underfill an otherwise valid count.
            compressed = np.sort(self.rng.choice(self.n_chromatin - count, count, replace=False))
            self.cohesin = [[int(start + offset), int(start + offset + 1)]
                            for offset, start in enumerate(compressed)]
        self.genes = []
        for gene in config.get("genes", []):
            self.genes.append({"config": gene, "on": False, "polymerases": [],
                               "unprocessed": 0, "mature": 0, "initiated": 0, "completed": 0})
        self.regulatory_regions = copy.deepcopy(config.get("regulatory_regions", []))
        self._enabled_regions = [region for region in self.regulatory_regions if region.get("enabled", True)]
        self._regions_by_gene = {}
        self._region_affinity_multipliers = {}
        coverage = {}
        for region in self._enabled_regions:
            self._regions_by_gene.setdefault(region["target_gene"], []).append(region)
            species = region.get("tf_species")
            if not species:
                continue
            if species not in coverage:
                coverage[species] = np.zeros(self.n_chromatin, dtype=bool)
                self._region_affinity_multipliers[species] = np.ones(self.n_chromatin)
            interval = slice(int(region["start"]), int(region["end"]) + 1)
            values = self._region_affinity_multipliers[species][interval]
            covered = coverage[species][interval]
            multiplier = float(region.get("tf_affinity_multiplier", 1))
            values[:] = np.where(covered, np.maximum(values, multiplier), multiplier)
            coverage[species][interval] = True
        self.nucleosome_config = copy.deepcopy(self.chrom.get("nucleosomes", {}))
        self.nucleosomes_enabled = bool(self.nucleosome_config.get("enabled", False))
        self.nucleosome_occupied = np.ones(self.n_chromatin, dtype=bool)
        self.nucleosome_counters = {"endpoint_evictions": 0, "endpoint_reassemblies": 0,
                                    "mark_reset_intervals": 0, "slides": 0}
        if self.nucleosomes_enabled:
            occupancy = float(self.nucleosome_config.get("initial_occupancy", 1))
            if occupancy == 0:
                self.nucleosome_occupied[:] = False
            elif occupancy < 1:
                self.nucleosome_occupied[:] = self.rng.random(self.n_chromatin) < occupancy
            # Explicit intervals override the sampled initial state, in order.
            for region in self.nucleosome_config.get("initial_regions", []):
                self.nucleosome_occupied[int(region["start"]):int(region["end"]) + 1] = bool(region["occupied"])
            self._clear_vacant_marks()
        self.last_substeps = 1
        self.explicit_transcription = None
        if config.get("transcription", {}).get("enabled", False):
            from .transcription import ExplicitTranscription
            self.explicit_transcription = ExplicitTranscription(self)

    def _random_direction(self):
        direction = self.rng.normal(size=3)
        return direction / np.linalg.norm(direction)

    def _random_position(self, fraction=0.8):
        return self._random_direction() * self.radius * fraction * self.rng.random() ** (1.0 / 3.0)

    def _initialize_particles(self):
        positions = [np.zeros(3)]
        radii = [self.bond_length / 2] * self.n_chromatin
        diffusion = [float(self.chrom.get("diffusion", 0.1))] * self.n_chromatin
        max_radius = max([self.bond_length / 2] + [float(p.get("radius", 0.3)) for p in self.species.values()])
        growing = GrowingNeighbors(cell_size=1.3 * max_radius, max_radius=max_radius)
        growing.add(0, positions[0], self.bond_length / 2)
        tangent = self._random_direction()
        persistence = float(self.chrom.get("persistence_length", 0.0))
        keep = math.exp(-self.bond_length / persistence) if persistence > 0 else 0.0
        backtracks = 0
        while len(positions) < self.n_chromatin:
            candidate = None
            for attempt in range(120):
                direction = keep * tangent + (1 - keep + 0.2) * self._random_direction()
                direction /= np.linalg.norm(direction)
                point = positions[-1] + self.bond_length * direction
                if np.linalg.norm(point) > self.radius * 0.78 or (
                        self.nucleus is not None and self.nucleus.rho(point)[0] > 0.97):
                    direction = -positions[-1] / max(np.linalg.norm(positions[-1]), 1e-12)
                    direction += 0.5 * self._random_direction()
                    direction /= np.linalg.norm(direction)
                    point = positions[-1] + self.bond_length * direction
                if growing.is_clear(point, self.bond_length / 2):
                    candidate = point
                    tangent = direction
                    break
            if candidate is None:
                # Escape a locally trapped growing tip without introducing a
                # severe overlap or changing the configured contour length.
                backtracks += 1
                if backtracks > max(100, self.n_chromatin // 2):
                    raise ValueError("Could not grow chromatin without severe overlaps. Increase confinement radius or reduce bead count/bond length.")
                trim = max(1, len(positions) - 20)
                for index in range(trim, len(positions)):
                    growing.remove(index)
                del positions[trim:]
                tangent = self._random_direction()
                continue
            growing.add(len(positions), candidate, self.bond_length / 2)
            positions.append(candidate)
        self.protein_particles = []
        self.molecules = []
        self.protein_chains = []
        self.idr_bonds = []
        for species_id, protein in self.species.items():
            size = int(protein.get("idr_beads", 4)) if protein.get("representation") == "idr" else 1
            radius = float(protein.get("radius", 0.3))
            for _ in range(int(protein.get("count", 0))):
                molecule = len(self.molecules)
                indices = []
                sites = int(round(size * float(protein.get("sticker_fraction", 1)))) if size > 1 else 1
                sites = max(0, min(size, sites))
                stickers = set(self.rng.choice(size, sites, replace=False).tolist()) if sites else set()
                last = None
                for bead in range(size):
                    point = None
                    for _attempt in range(150):
                        proposal = self._random_position() if last is None else last + 2 * radius * self._random_direction()
                        if np.linalg.norm(proposal) > self.radius * 0.92:
                            continue
                        if self.nucleus is not None and self.nucleus.rho(proposal)[0] > 0.97:
                            continue
                        if growing.is_clear(proposal, radius):
                            point = proposal
                            break
                    if point is None:
                        # Fail visibly instead of distorting a chain or silently clipping.
                        raise ValueError("Could not initialize particles without severe overlaps. Increase confinement radius or reduce counts/radii.")
                    index = len(positions)
                    indices.append(index)
                    growing.add(index, point, radius)
                    positions.append(point)
                    radii.append(radius)
                    diffusion.append(float(protein.get("diffusion", 0.5)) * size)
                    self.protein_particles.append({"index": index, "species": species_id, "molecule": molecule,
                                                   "sticker": bead in stickers, "capacity": (1 if size > 1 else int(protein.get("valence", 1)))})
                    if last is not None:
                        self.idr_bonds.append((index - 1, index, 20.0, 2 * radius))
                    last = point
                self.molecules.append({"species": species_id, "indices": indices,
                                       "valence": int(protein.get("valence", 1))})
                if size > 1:
                    self.protein_chains.append(indices)
        self.positions = np.asarray(positions, dtype=float)
        self.diffusion = np.asarray(diffusion, dtype=float)
        self.radii = np.asarray(radii, dtype=float)
        self._particle_info = {p["index"]: p for p in self.protein_particles}
        excluded = {(i, i + 1) for i in range(self.n_chromatin - 1)}
        excluded.update((a, b) for a, b, _, _ in self.idr_bonds)
        self._excluded_pair_codes = np.array(sorted(a * len(self.positions) + b for a, b in excluded), dtype=np.int64)
        self._repulsion_cutoff = 2 * float(np.max(self.radii))
        self.initialization_backtracks = backtracks
        self._body_affinities = []
        for body in self.bodies:
            affinity = np.zeros(len(self.positions))
            affinity[:self.n_chromatin] = float(body.get("chromatin_affinity", 0))
            for particle in self.protein_particles:
                # Divide the molecular affinity across its IDR domain beads.
                n = len(self.molecules[particle["molecule"]]["indices"])
                affinity[particle["index"]] = float(body.get("protein_affinity", {}).get(particle["species"], 0)) / n
            self._body_affinities.append(affinity * float(body.get("strength", 1)))

    def local_persistence(self):
        values = np.full(self.n_chromatin, float(self.chrom.get("persistence_length", 0)))
        for gene in self.genes:
            if gene["polymerases"]:
                spec = gene["config"]
                values[int(spec["start"]):int(spec["end"]) + 1] *= float(spec.get("persistence_multiplier", 1))
        if self.nucleosomes_enabled:
            values *= np.where(self.nucleosome_occupied,
                               float(self.nucleosome_config.get("occupied_persistence_multiplier", 1)),
                               float(self.nucleosome_config.get("vacant_persistence_multiplier", 0.5)))
        return values

    def _repulsion_pairs(self):
        """Actual nonbonded soft-core contacts, with no stored N-by-N arrays."""
        pairs = neighbor_pairs(self.positions, self._repulsion_cutoff)
        if not pairs.size:
            return pairs, np.empty(0)
        codes = pairs[:, 0] * len(self.positions) + pairs[:, 1]
        if self._excluded_pair_codes.size:
            slots = np.searchsorted(self._excluded_pair_codes, codes)
            in_range = slots < len(self._excluded_pair_codes)
            excluded = np.zeros(len(pairs), dtype=bool)
            excluded[in_range] = self._excluded_pair_codes[slots[in_range]] == codes[in_range]
            pairs = pairs[~excluded]
        sigma2 = (self.radii[pairs[:, 0]] + self.radii[pairs[:, 1]]) ** 2
        delta = self.positions[pairs[:, 0]] - self.positions[pairs[:, 1]]
        inside = np.einsum("ij,ij->i", delta, delta) < sigma2
        return pairs[inside], sigma2[inside]

    def _mechanics(self):
        x = self.positions
        f = np.zeros_like(x)
        energy = 0.0

        def springs(a, b, stiffness, length):
            nonlocal energy
            a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
            if a.size == 0:
                return
            delta = x[b] - x[a]
            distance = np.linalg.norm(delta, axis=1)
            if np.any(distance < 1e-10):
                raise FloatingPointError("A bonded pair coincided. Reduce dt or relax the initial density.")
            extension = distance - length
            force = (stiffness * extension / distance)[:, None] * delta
            np.add.at(f, a, force)
            np.add.at(f, b, -force)
            energy += float(np.sum(0.5 * stiffness * extension ** 2))

        backbone = np.arange(self.n_chromatin - 1)
        springs(backbone, backbone + 1, float(self.chrom["bond_k"]), self.bond_length)
        if self.idr_bonds:
            a, b, k, length = map(np.asarray, zip(*self.idr_bonds))
            springs(a, b, k, length)
        if self.bindings:
            springs([b["a"] for b in self.bindings], [b["b"] for b in self.bindings],
                    np.array([b["k"] for b in self.bindings]), np.array([b["rest"] for b in self.bindings]))
        if self.cohesin:
            springs([p[0] for p in self.cohesin], [p[1] for p in self.cohesin],
                    float(self.cohesin_config.get("spring_k", 0)), self.bond_length)
        if self.n_chromatin >= 3:
            a = x[1:self.n_chromatin - 1] - x[:self.n_chromatin - 2]
            b = x[2:self.n_chromatin] - x[1:self.n_chromatin - 1]
            la, lb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
            if np.any(la < 1e-10) or np.any(lb < 1e-10):
                raise FloatingPointError("A backbone segment collapsed. Reduce dt.")
            u, v = a / la[:, None], b / lb[:, None]
            cosine = np.einsum("ij,ij->i", u, v)
            stiffness = self.local_persistence()[1:-1] / self.bond_length
            ga = -stiffness[:, None] / la[:, None] * (v - cosine[:, None] * u)
            gb = -stiffness[:, None] / lb[:, None] * (u - cosine[:, None] * v)
            f[:self.n_chromatin - 2] += ga
            f[1:self.n_chromatin - 1] += -ga + gb
            f[2:self.n_chromatin] -= gb
            energy += float(np.sum(stiffness * (1 - cosine)))
        # Smooth, bounded soft repulsion: U=eps/2*(1-r^2/sigma^2)^2.
        eps = float(self.chrom.get("repulsion", 0))
        if eps:
            pairs, sigma2 = self._repulsion_pairs()
            delta = x[pairs[:, 0]] - x[pairs[:, 1]]
            r2 = np.einsum("ij,ij->i", delta, delta)
            q = 1 - r2 / sigma2
            force = (2 * eps * q / sigma2)[:, None] * delta
            np.add.at(f, pairs[:, 0], force)
            np.add.at(f, pairs[:, 1], -force)
            energy += float(0.5 * eps * np.sum(q ** 2))
        if self.chrom.get("state_interactions"):
            contacts = self._chromatin_contact_state()
            f[:self.n_chromatin] += contacts["forces"]
            energy += contacts["energy"]
        # Confinement is a soft harmonic wall, not an undocumented coordinate clamp.
        radial = np.linalg.norm(x, axis=1)
        outside = radial > self.radius
        if np.any(outside):
            penetration = radial[outside] - self.radius
            f[outside] -= (self.confinement_stiffness * penetration / radial[outside])[:, None] * x[outside]
            energy += float(0.5 * self.confinement_stiffness * np.sum(penetration ** 2))
        if self.nucleus is not None:
            wall, wall_energy = self.nucleus.force_energy(x)
            f += wall
            energy += wall_energy
            if self.tethers:
                beads = np.array([t["bead"] for t in self.tethers], dtype=int)
                anchors = np.array([t["anchor"] for t in self.tethers], dtype=float)
                k = np.array([float(self.tether_specs[t["tether"]].get("spring_k", 2.0)) for t in self.tethers])
                rest = np.array([float(self.tether_specs[t["tether"]].get("rest_length", 0.5)) for t in self.tethers])
                delta = x[beads] - anchors
                distance = np.linalg.norm(delta, axis=1)
                extension = distance - rest
                safe = np.maximum(distance, 1e-12)
                force = -(k * extension / safe)[:, None] * delta
                force[distance < 1e-12] = 0.0
                np.add.at(f, beads, force)
                energy += float(np.sum(0.5 * k * extension ** 2))
        # Fixed Gaussian affinity wells; signed amplitudes also support exclusion.
        for body, affinity in zip(self.bodies, self._body_affinities):
            delta = x - np.asarray(body["center"], dtype=float)
            r2 = np.einsum("ij,ij->i", delta, delta)
            scale2 = float(body["radius"]) ** 2
            potential = -affinity * np.exp(-0.5 * r2 / scale2)
            energy += float(np.sum(potential))
            f += potential[:, None] * delta / scale2
        if not np.all(np.isfinite(f)) or not math.isfinite(energy):
            raise FloatingPointError("Non-finite mechanics. Reduce dt, stiffness, or particle density.")
        return f, energy

    def forces(self):
        return self._mechanics()[0]

    get_forces = forces

    def potential_energy(self):
        return self._mechanics()[1]

    def _chromatin_contact_state(self):
        from .chromatin_contacts import compute_chromatin_contacts
        return compute_chromatin_contacts(self.positions[:self.n_chromatin], self.marks,
                                          self.chrom.get("state_interactions", []))

    def integration_substeps(self):
        """Conservative stability subdivision; chemistry still advances each dt.

        This is a numerical safeguard, not an accuracy guarantee. Convergence
        with smaller dt remains necessary, especially for binding kinetics.
        """
        dmax = float(np.max(self.diffusion)) if self.diffusion.size else 0.0
        if dmax == 0:
            return 1
        stiffness = np.zeros(len(self.positions))
        stiffness[:self.n_chromatin] = 2 * float(self.chrom["bond_k"])
        for a, b, k, _ in self.idr_bonds:
            stiffness[a] += k
            stiffness[b] += k
        for bond in self.bindings:
            stiffness[bond["a"]] += bond["k"]
            stiffness[bond["b"]] += bond["k"]
        for left, right in self.cohesin:
            stiffness[left] += float(self.cohesin_config.get("spring_k", 0))
            stiffness[right] += float(self.cohesin_config.get("spring_k", 0))
        stiffness[:self.n_chromatin] += 8 * self.local_persistence() / self.bond_length ** 3
        stiffness += 2 * float(self.chrom.get("repulsion", 0)) / (2 * self.radii) ** 2
        if self.chrom.get("state_interactions"):
            stiffness[:self.n_chromatin] += self._chromatin_contact_state()["stiffness"]
        # The spectral magnitude of the Gaussian-well Hessian is bounded by
        # |amplitude|/radius^2. Both attraction and exclusion can be stiff.
        for body, amplitude in zip(self.bodies, self._body_affinities):
            stiffness += np.abs(amplitude) / float(body["radius"]) ** 2
        # The radial wall's maximum curvature is its spring constant.
        stiffness += self.confinement_stiffness
        if self.nucleus is not None:
            stiffness += self.nucleus.curvature_bound()
            for tether in self.tethers:
                stiffness[tether["bead"]] += float(self.tether_specs[tether["tether"]].get("spring_k", 2.0))
        limit = float(np.max(self.diffusion * stiffness))
        subdivisions = max(1, int(math.ceil(self.dt * limit / 0.35)))
        if subdivisions > 1000:
            raise ValueError("This parameter combination needs more than 1000 integration substeps per step. Reduce dt, diffusion, bond stiffness or body affinity/strength, or increase the body radius.")
        return subdivisions

    def _brownian_move(self):
        subdivisions = self.integration_substeps()
        self.last_substeps = subdivisions
        dt = self.dt / subdivisions
        noise = np.sqrt(2 * self.diffusion * dt)[:, None]
        for _ in range(subdivisions):
            force = self.get_forces()
            proposal = self.positions + self.diffusion[:, None] * force * dt + noise * self.rng.normal(size=self.positions.shape)
            self.validate_positions(proposal)
            self.positions[:] = proposal

    def validate_positions(self, positions=None):
        positions = self.positions if positions is None else positions
        nucleus = getattr(self, "nucleus", None)
        if (not np.all(np.isfinite(positions)) or
                (getattr(self, "confinement_stiffness", 30.0) > 0
                 and np.max(np.linalg.norm(positions, axis=1)) > self.radius * 50) or
                (nucleus is not None and nucleus.stiffness > 0 and np.max(nucleus.rho(positions)) > 50)):
            raise FloatingPointError("Integration became unstable; coordinates were not clipped. Reduce dt, stiffness, or diffusion.")

    def _event(self, kind, detail, time=None):
        self.events.append({"time": float(self.time if time is None else time), "type": kind, "detail": detail})

    def _genomic_event(self, kind, bead, dt, *, time=None, timing="interval_endpoint", **fields):
        """Observe accepted changes, with bounded retention and no random draw.

        IDs order observations, not a globally interleaved reaction timeline:
        reaction families are operator split. Turnover endpoints do not resolve
        individual microscopic eviction/reassembly times within an interval.
        """
        self.genomic_event_count += 1
        self.genomic_events.append({"id": self.genomic_event_count, "type": kind,
            "time": float(self.time + dt if time is None else time),
            "step": self.step_count + 1, "bead": int(bead), "timing": timing,
            "interval_start": float(self.time), "interval_end": float(self.time + dt), **fields})

    def genomic_event_state(self):
        return {"run_id": self.run_id, "capacity": self.genomic_events.maxlen,
                "total_emitted": self.genomic_event_count, "last_id": self.genomic_event_count,
                "first_retained_id": self.genomic_events[0]["id"] if self.genomic_events else None,
                "dropped": self.genomic_event_count - len(self.genomic_events),
                "events": copy.deepcopy(list(self.genomic_events))}

    def _chromatin_off_rates(self, species_id, base_rate=None):
        protein = self.species[species_id]
        rates = np.full(self.n_chromatin, float(protein.get("off_rate", 0) if base_rate is None else base_rate))
        if self.nucleosomes_enabled:
            factors = protein.get("nucleosome_off_rate_multiplier", {})
            rates *= np.where(self.nucleosome_occupied, float(factors.get("occupied", 1)),
                              float(factors.get("vacant", 1)))
        return rates

    def _bond_off_rate(self, bond):
        particle = self._particle_info.get(bond["b"])
        fallback = self.species[particle["species"]].get("off_rate", 0) if bond["kind"] == "chromatin" else 0
        rate = float(bond.get("off_rate", fallback))
        if bond["kind"] == "chromatin" and self.nucleosomes_enabled:
            state = "occupied" if self.nucleosome_occupied[bond["a"]] else "vacant"
            rate *= float(self.species[particle["species"]].get("nucleosome_off_rate_multiplier", {}).get(state, 1))
        return rate

    def _chromatin_affinity(self, species_id):
        mapping = self.species[species_id].get("chromatin_affinity", {})
        result = np.full(self.n_chromatin, float(mapping.get("unmodified", 0)))
        mapped = np.zeros(self.n_chromatin, dtype=bool)
        for mark, flags in self.marks.items():
            if mark in mapping:
                value = float(mapping[mark])
                result[flags & ~mapped] = value
                result[flags & mapped] = np.maximum(result[flags & mapped], value)
                mapped |= flags
        if species_id in self._region_affinity_multipliers:
            result *= self._region_affinity_multipliers[species_id]
        if self.nucleosomes_enabled:
            factors = self.species[species_id].get("nucleosome_affinity", {})
            result *= np.where(self.nucleosome_occupied, float(factors.get("occupied", 1)),
                               float(factors.get("vacant", 1)))
        return result

    def _clear_vacant_marks(self):
        if self.nucleosomes_enabled:
            for flags in self.marks.values():
                flags[~self.nucleosome_occupied] = False

    def _remodeling_targets(self):
        """Unique molecule/anchor pairs from direct or configured recruitment."""
        context = self._regulatory_binding_context()
        targets = {(molecule, bead) for molecule, sites in context["direct_sites"].items() for bead in sites}
        for region in self._enabled_regions:
            _, recruited = self._region_recruitment(region, context)
            targets.update((molecule, bead) for molecule, sites in recruited.items() for bead in sites)
        return sorted(targets)

    def _nucleosome_rates(self):
        """Per-site remodeling hazards from the reservoir and actual anchors.

        A protein contributes once per molecule/bead, including a recruited
        coactivator's TF anchor. A paused polymerase remains engaged and can
        contribute the gene's optional remodeling hazard at its present bead.
        """
        eviction = np.full(self.n_chromatin, float(self.nucleosome_config.get("eviction_rate", 0)))
        reassembly = np.full(self.n_chromatin, float(self.nucleosome_config.get("reassembly_rate", 0)))
        if any(float(p.get("nucleosome_eviction_rate", 0)) or
               float(p.get("nucleosome_reassembly_rate", 0)) for p in self.species.values()):
            for molecule, bead in self._remodeling_targets():
                protein = self.species[self.molecules[molecule]["species"]]
                eviction[bead] += float(protein.get("nucleosome_eviction_rate", 0))
                reassembly[bead] += float(protein.get("nucleosome_reassembly_rate", 0))
        for gene in self.genes:
            rate = float(gene["config"].get("nucleosome_eviction_rate", 0))
            if rate and gene["polymerases"]:
                beads = np.clip(np.floor(gene["polymerases"]).astype(int), 0, self.n_chromatin - 1)
                np.add.at(eviction, beads, rate)
        return eviction, reassembly

    def _advance_nucleosomes(self, dt):
        """Sample exact two-state endpoints and retain eviction memory.

        For an initially occupied site, P(final vacant) is the standard
        two-state transition probability, while P(any eviction)=1-exp(-E*dt).
        Sampling these nested events with one uniform also clears histone
        marks when eviction and reassembly occur within the same interval.
        Counters describe endpoint changes or mark-reset intervals, not the
        unobserved number of microscopic events inside an interval.
        """
        if not self.nucleosomes_enabled:
            return
        self._clear_vacant_marks()
        eviction, reassembly = self._nucleosome_rates()
        occupied = self.nucleosome_occupied.copy()
        if not np.any(np.where(occupied, eviction, reassembly) > 0):
            return
        total = eviction + reassembly
        positive = total > 0
        transition = np.zeros(self.n_chromatin)
        transition[positive] = _probability(total[positive], dt) / total[positive]
        final_vacant = eviction * transition
        final_occupied = reassembly * transition
        uniforms = self.rng.random(self.n_chromatin)
        self.nucleosome_occupied[:] = np.where(occupied, uniforms >= final_vacant,
                                               uniforms < final_occupied)
        # Analytically final_vacant <= any_eviction; max protects the nesting
        # from tiny floating point discrepancies at extreme rate ratios.
        any_eviction = np.maximum(final_vacant, _probability(eviction, dt))
        reset = occupied & (uniforms < any_eviction)
        for flags in self.marks.values():
            flags[reset | ~self.nucleosome_occupied] = False
        counts = {"endpoint_evictions": int(np.count_nonzero(occupied & ~self.nucleosome_occupied)),
                  "endpoint_reassemblies": int(np.count_nonzero(~occupied & self.nucleosome_occupied)),
                  "mark_reset_intervals": int(np.count_nonzero(reset))}
        for key, value in counts.items():
            self.nucleosome_counters[key] += value
        if any(counts.values()):
            self._event("nucleosome", f"{counts['endpoint_evictions']} occupied→vacant; "
                        f"{counts['endpoint_reassemblies']} vacant→occupied; "
                        f"{counts['mark_reset_intervals']} occupied-site mark-reset intervals")
        for bead in np.flatnonzero(occupied != self.nucleosome_occupied):
            self._genomic_event("nucleosome_eviction" if occupied[bead] else "nucleosome_reassembly",
                bead, dt, reason="turnover_endpoint", occupied=bool(self.nucleosome_occupied[bead]),
                eviction_rate=float(eviction[bead]), reassembly_rate=float(reassembly[bead]))
        for bead in np.flatnonzero(reset & self.nucleosome_occupied):
            self._genomic_event("nucleosome_reset", bead, dt,
                reason="evicted_and_reassembled_within_interval", occupied=True,
                eviction_rate=float(eviction[bead]), reassembly_rate=float(reassembly[bead]))

    def _advance_nucleosome_sliding(self, dt):
        """Coarse remodeling hops, exact conditional on frozen protein anchors.

        Each unique bound/recruited molecule+anchor supplies its configured
        total sliding hazard. Valid hops share that rate uniformly: an occupied
        anchor or immediate neighbor can move one lattice site into a vacancy.
        Histone flags move with the nucleosome. This conserves nucleosome number
        and represents neither ATP turnover nor sub-bead wrapped-DNA geometry.
        Sliding is operator split from eviction, binding and other chemistry.
        """
        if not self.nucleosomes_enabled or not any(float(p.get("nucleosome_sliding_rate", 0)) > 0
                                                   for p in self.species.values()):
            return
        targets = [(molecule, anchor, float(self.species[self.molecules[molecule]["species"]]
                                            .get("nucleosome_sliding_rate", 0)))
                   for molecule, anchor in self._remodeling_targets()]
        targets = [target for target in targets if target[2] > 0]
        elapsed, count = 0.0, 0
        while targets:
            candidates, hazards = [], []
            for molecule, anchor, rate in targets:
                moves = [(source, destination)
                         for source in range(max(0, anchor - 1), min(self.n_chromatin, anchor + 2))
                         if self.nucleosome_occupied[source]
                         for destination in (source - 1, source + 1)
                         if 0 <= destination < self.n_chromatin and not self.nucleosome_occupied[destination]]
                if not moves:
                    continue
                candidates.extend((molecule, anchor, source, destination, rate) for source, destination in moves)
                hazards.extend([rate / len(moves)] * len(moves))
            if not candidates:
                return
            total = float(sum(hazards))
            elapsed += float(self.rng.exponential(1.0 / total))
            if elapsed >= dt:
                return
            index = min(int(np.searchsorted(np.cumsum(hazards), float(self.rng.random()) * total, side="right")),
                        len(candidates) - 1)
            molecule, anchor, source, destination, rate = candidates[index]
            if count >= 10000:
                raise ValueError("Nucleosome sliding exceeded 10000 events in one interval; reduce dt or sliding rates.")
            self.nucleosome_occupied[source] = False
            self.nucleosome_occupied[destination] = True
            for flags in self.marks.values():
                flags[destination], flags[source] = flags[source], False
            self.nucleosome_counters["slides"] += 1
            count += 1
            self._genomic_event("nucleosome_slide", source, dt, time=self.time + elapsed,
                timing="sampled_slide", reason="anchored_remodeler", to_bead=int(destination),
                molecule=int(molecule), species=self.molecules[molecule]["species"],
                anchor_bead=int(anchor), rate=rate)

    def nucleosome_state(self):
        """Read occupancy without treating disabled occupancy as a measurement."""
        return {"enabled": self.nucleosomes_enabled,
                "occupied": self.nucleosome_occupied.astype(int).tolist() if self.nucleosomes_enabled else [],
                "occupied_fraction": float(np.mean(self.nucleosome_occupied)) if self.nucleosomes_enabled else None,
                **self.nucleosome_counters}

    def _advance_binding(self, dt, chromatin_tree=None):
        retained = []
        for bond in self.bindings:
            off_rate = self._bond_off_rate(bond)
            uniform = self.rng.random()
            if uniform < _probability(off_rate, dt):
                self._event("unbind", bond["kind"])
                if bond["kind"] == "chromatin":
                    particle = self._particle_info[bond["b"]]
                    self._genomic_event("chromatin_unbind", bond["a"], dt,
                        time=self.time - math.log1p(-uniform) / off_rate,
                        timing="sampled_dissociation", reason="off_rate",
                        particle=int(particle["index"]), molecule=int(particle["molecule"]),
                        species=particle["species"], off_rate=off_rate)
            else:
                retained.append(bond)
        self.bindings = retained
        used_molecule = np.zeros(len(self.molecules), dtype=int)
        used_site = np.zeros(len(self.positions), dtype=int)
        existing = set()
        molecule_chromatin = set()
        for bond in self.bindings:
            existing.add(tuple(sorted((bond["a"], bond["b"]))))
            for index in (bond["a"], bond["b"]):
                if index in self._particle_info:
                    used_molecule[self._particle_info[index]["molecule"]] += 1
                    used_site[index] += 1
            if bond["kind"] == "chromatin":
                molecule_chromatin.add((self._particle_info[bond["b"]]["molecule"], bond["a"]))

        def available(particle):
            return (particle["sticker"] and used_site[particle["index"]] < particle["capacity"]
                    and used_molecule[particle["molecule"]] < self.molecules[particle["molecule"]]["valence"])

        candidates = []
        active = [p for p in self.protein_particles if available(p)]
        if not active:
            return
        chromatin_tree = cKDTree(self.positions[:self.n_chromatin]) if chromatin_tree is None else chromatin_tree
        binding_rates = {key: float(value.get("on_rate", 0)) * self._chromatin_affinity(key)
                         for key, value in self.species.items()}
        active_indices = np.array([p["index"] for p in active], dtype=int)
        captures = self.radii[active_indices] + self.bond_length
        neighborhoods = chromatin_tree.query_ball_point(self.positions[active_indices], captures, return_sorted=True)
        for particle, capture, neighborhood in zip(active, captures, neighborhoods):
            index = particle["index"]
            rates = binding_rates[particle["species"]]
            neighbors = np.asarray(neighborhood, dtype=int)
            if not neighbors.size:
                continue
            delta = self.positions[neighbors] - self.positions[index]
            distance2 = np.einsum("ij,ij->i", delta, delta)
            neighbors = neighbors[(distance2 < capture * capture) & (rates[neighbors] > 0)]
            for chromatin in neighbors:
                if (particle["molecule"], int(chromatin)) in molecule_chromatin:
                    continue
                waiting = self.rng.exponential(1.0 / rates[chromatin])
                if waiting < dt:
                    candidates.append((waiting, int(chromatin), index, "chromatin", None))
        cutoff = max([0] + [float(r["radius"]) for r in self.interactions.values() if float(r.get("on_rate", 0)) > 0])
        encounter_pairs = neighbor_pairs(self.positions[active_indices], cutoff)
        for left, right in encounter_pairs:
            a, b = active[left], active[right]
            if a["molecule"] == b["molecule"]:
                continue
            pair = (a["index"], b["index"])
            if pair in existing:
                continue
            interaction = self.interactions.get(tuple(sorted((a["species"], b["species"]))))
            if not interaction or float(interaction.get("on_rate", 0)) <= 0:
                continue
            delta = self.positions[a["index"]] - self.positions[b["index"]]
            if float(np.dot(delta, delta)) >= float(interaction["radius"]) ** 2:
                continue
            waiting = self.rng.exponential(1.0 / float(interaction["on_rate"]))
            if waiting < dt:
                candidates.append((waiting, a["index"], b["index"], "protein", interaction))
        # Exponential clocks order competing encounters and prevent overbooking sites.
        for waiting, a, b, kind, interaction in sorted(candidates, key=lambda item: item[0]):
            pa, pb = self._particle_info.get(a), self._particle_info[b]
            if not available(pb) or (pa is not None and not available(pa)):
                continue
            if kind == "chromatin" and (pb["molecule"], a) in molecule_chromatin:
                continue
            if tuple(sorted((a, b))) in existing:
                continue
            if kind == "chromatin":
                off_rate = float(self.species[pb["species"]].get("off_rate", 0))
                stiffness = 12.0
                molecule_chromatin.add((pb["molecule"], a))
            else:
                off_rate = float(interaction.get("off_rate", 0))
                stiffness = float(interaction.get("strength", 1))
            self.bindings.append({"a": a, "b": b, "kind": kind, "off_rate": off_rate,
                                  "k": stiffness, "rest": float(self.radii[a] + self.radii[b])})
            existing.add(tuple(sorted((a, b))))
            for particle in (pa, pb):
                if particle is not None:
                    used_site[particle["index"]] += 1
                    used_molecule[particle["molecule"]] += 1
            self._event("bind", kind)
            if kind == "chromatin":
                self._genomic_event("chromatin_bind", a, dt, time=self.time + waiting,
                    timing="sampled_encounter", reason="capture_encounter", particle=int(b),
                    molecule=int(pb["molecule"]), species=pb["species"],
                    on_rate=float(binding_rates[pb["species"]][a]),
                    off_rate=self._bond_off_rate(self.bindings[-1]))

    def _advance_marks(self, dt, chromatin_tree=None):
        if not self.marks:
            return
        self._clear_vacant_marks()
        write = {key: np.zeros(self.n_chromatin) for key in self.marks}
        erase = {key: np.full(self.n_chromatin, float(self.mark_config[key].get("decay_rate", 0))) for key in self.marks}
        for bond in self.bindings:
            if bond["kind"] != "chromatin":
                continue
            protein = self.species[self._particle_info[bond["b"]]["species"]]
            rate = float(protein.get("catalytic_rate", 0))
            if protein.get("writer_mark") in write:
                write[protein["writer_mark"]][bond["a"]] += rate
            if protein.get("eraser_mark") in erase:
                erase[protein["eraser_mark"]][bond["a"]] += rate
        # Recruitment extends catalysis only to the actual TF anchor of an
        # existing one-hop protein bond. One coactivator can reach several
        # anchors, but duplicate IDR links/regions do not multiply its rate at
        # the same bead, nor does its own direct attachment add a second rate.
        if self._enabled_regions:
            context = self._regulatory_binding_context()
            recruited_targets = set()
            for region in self._enabled_regions:
                _, targets = self._region_recruitment(region, context)
                for molecule, beads in targets.items():
                    own_sites = context["direct_sites"].get(molecule, set())
                    recruited_targets.update((molecule, bead) for bead in beads if bead not in own_sites)
            for molecule, bead in sorted(recruited_targets):
                protein = self.species[self.molecules[molecule]["species"]]
                rate = float(protein.get("catalytic_rate", 0))
                if protein.get("writer_mark") in write:
                    write[protein["writer_mark"]][bead] += rate
                if protein.get("eraser_mark") in erase:
                    erase[protein["eraser_mark"]][bead] += rate
        contact_radius = max([0] + [float(m["contact_radius"]) for m in self.mark_config.values() if float(m.get("spread_rate", 0)) > 0])
        contact_pairs = neighbor_pairs(self.positions[:self.n_chromatin], contact_radius, tree=chromatin_tree)
        delta = self.positions[contact_pairs[:, 0]] - self.positions[contact_pairs[:, 1]]
        contact_distance2 = np.einsum("ij,ij->i", delta, delta)
        for key, state in self.marks.items():
            spec = self.mark_config[key]
            spread = float(spec.get("spread_rate", 0))
            if spread:
                pairs = contact_pairs[contact_distance2 < float(spec["contact_radius"]) ** 2]
                # Each undirected contact contributes its donor state to the
                # opposite target. This exactly matches a dense adjacency dot
                # product, with no N-by-N adjacency or displacement arrays.
                donors = np.bincount(pairs[:, 0], weights=state[pairs[:, 1]], minlength=self.n_chromatin)
                donors += np.bincount(pairs[:, 1], weights=state[pairs[:, 0]], minlength=self.n_chromatin)
                write[key] += spread * donors
            if self.nucleosomes_enabled:
                write[key][~self.nucleosome_occupied] = 0
            total = write[key] + erase[key]
            positive = total > 0
            flip = np.zeros(self.n_chromatin)
            hazard = np.where(state, erase[key], write[key])
            flip[positive] = hazard[positive] / total[positive] * _probability(total[positive], dt)
            changed = self.rng.random(self.n_chromatin) < flip
            if np.any(changed):
                state[changed] = ~state[changed]
                self._event("mark", f"{key}: {int(changed.sum())} bead changes")

    def _load_cohesin(self, ignore=None):
        occupied = {bead for pair in self.cohesin if pair is not ignore for bead in pair}
        candidates = [i for i in range(self.n_chromatin - 1) if i not in occupied and i + 1 not in occupied]
        if not candidates:
            return None
        left = int(self.rng.choice(candidates))
        return [left, left + 1]

    def _advance_cohesin(self, dt):
        if not self.cohesin_config.get("enabled"):
            return
        rate = float(self.cohesin_config.get("step_rate", 0))
        off_probability = _probability(float(self.cohesin_config.get("off_rate", 0)), dt)
        retained = []
        for pair in self.cohesin:
            if self.rng.random() < off_probability:
                self._event("cohesin", "unloaded")
            else:
                retained.append(pair)
        self.cohesin = retained
        # Unloaded extruders form an implicit free pool and load immediately
        # when a pair of adjacent unoccupied sites is available. No kinetic
        # loading rate is claimed because the configuration has none.
        while len(self.cohesin) < int(self.cohesin_config.get("count", 0)):
            replacement = self._load_cohesin()
            if replacement is None:
                break
            self.cohesin.append(replacement)
            self._event("cohesin", "loaded on adjacent sites")
        for index in self.rng.permutation(len(self.cohesin)):
            pair = self.cohesin[index]
            moves = [(side, int(self.rng.poisson(rate * dt))) for side in (0, 1)]
            self.rng.shuffle(moves)
            for side, count in moves:
                direction = "left" if side == 0 else "right"
                increment = -1 if side == 0 else 1
                for _ in range(count):
                    current = pair[side]
                    barrier = self.barriers.get(current, set())
                    if "both" in barrier or direction in barrier:
                        break
                    target = current + increment
                    if not 0 <= target < self.n_chromatin:
                        break
                    occupied = {bead for other in self.cohesin if other is not pair for bead in other}
                    if target in occupied:
                        break
                    pair[side] = target

    def _tether_hazard(self, spec, bead, kind):
        """Current attach or detach hazard for one bead (reduced time).

        Attachment = attach_rate x occupancy factor x mark gate; detachment =
        detach_rate x occupancy factor. Occupancy factors use the current
        occupied/vacant state when nucleosomes are enabled, else 1.
        """
        rate = float(spec.get("attach_rate" if kind == "attach" else "detach_rate", 0))
        if self.nucleosomes_enabled:
            factors = spec.get("occupancy_attach_multiplier" if kind == "attach" else "occupancy_detach_multiplier", {})
            rate *= float(factors.get("occupied" if self.nucleosome_occupied[bead] else "vacant", 1))
        if kind == "attach" and spec.get("required_mark") and not self.marks[spec["required_mark"]][bead]:
            return 0.0
        return rate

    def _tether_geometry(self, spec, points):
        """(|distance to anchor surface|, surface points) for candidate beads."""
        from .nucleus import body_surface_distance, body_surface_point
        if spec.get("anchor", "envelope") == "body":
            body = next(b for b in self.bodies if b["id"] == spec["body_id"])
            return (np.abs(body_surface_distance(body, points)),
                    np.array([body_surface_point(body, p) for p in points]))
        return np.abs(self.nucleus.signed_distance(points)), self.nucleus.surface_point(points)

    def _advance_tethers(self, dt):
        if self.nucleus is None or not self.tether_specs:
            return
        kept = []
        for tether in self.tethers:
            spec = self.tether_specs[tether["tether"]]
            if self.rng.random() < _probability(self._tether_hazard(spec, tether["bead"], "detach"), dt):
                self._genomic_event("tether_detach", tether["bead"], dt, tether=spec["id"],
                                    anchor_kind=spec.get("anchor", "envelope"))
            else:
                kept.append(tether)
        self.tethers = kept
        attached = {t["bead"] for t in kept}
        candidates = [bead for bead in sorted(self._tether_owner) if bead not in attached]
        if not candidates:
            return
        # One uniform per candidate keeps the draw count independent of geometry.
        draws = self.rng.random(len(candidates))
        by_spec = {}
        for bead, u in zip(candidates, draws):
            by_spec.setdefault(self._tether_owner[bead], []).append((bead, u))
        for index, items in sorted(by_spec.items()):
            spec = self.tether_specs[index]
            beads = [bead for bead, _ in items]
            distance, surface = self._tether_geometry(spec, self.positions[beads])
            capture = float(spec.get("capture_distance", 1.5))
            for (bead, u), d, anchor in zip(items, distance, surface):
                if d > capture:
                    continue
                if u < _probability(self._tether_hazard(spec, bead, "attach"), dt):
                    self.tethers.append({"bead": int(bead), "tether": index, "anchor": [float(v) for v in anchor]})
                    self._genomic_event("tether_attach", bead, dt, tether=spec["id"],
                                        anchor_kind=spec.get("anchor", "envelope"))

    def nucleus_state(self):
        """Read the envelope and tether state without RNG or mutation."""
        if self.nucleus is None:
            return {"enabled": False, "tethers": [], "attachments": []}
        chromatin = self.positions[:self.n_chromatin]
        distance = self.nucleus.signed_distance(chromatin)
        records = []
        for index, spec in enumerate(self.tether_specs):
            beads = np.arange(int(spec["start"]), int(spec["end"]) + 1)
            attached = sum(1 for t in self.tethers if t["tether"] == index)
            geometry, _ = self._tether_geometry(spec, chromatin[beads])
            records.append({"id": spec["id"], "name": spec.get("name", spec["id"]), "start": int(spec["start"]),
                            "end": int(spec["end"]), "anchor": spec.get("anchor", "envelope"),
                            "body_id": spec.get("body_id"), "color": spec.get("color"),
                            "attached": attached, "beads": int(len(beads)),
                            "within_capture": int(np.count_nonzero(geometry <= float(spec.get("capture_distance", 1.5)))),
                            "mean_anchor_surface_distance": float(np.mean(geometry))})
        return {"enabled": True, **self.nucleus.state(self.positions),
                "chromatin_envelope_distance": {"min": float(np.min(-distance)), "mean": float(np.mean(-distance)),
                                                "within_1": int(np.count_nonzero(np.abs(distance) <= 1.0))},
                "tethers": records,
                "attachments": [{"bead": t["bead"], "tether": self.tether_specs[t["tether"]]["id"],
                                 "anchor": list(t["anchor"])} for t in self.tethers]}

    def _legacy_gene_on_rate(self, spec):
        promoter = int(spec["start"])
        required, repressive = spec.get("required_mark"), spec.get("repressive_mark")
        if required and not self.marks[required][promoter]:
            return 0.0
        if repressive and self.marks[repressive][promoter]:
            return 0.0
        activators = set()
        for bond in self.bindings:
            if bond["kind"] == "chromatin" and int(spec["start"]) <= bond["a"] <= int(spec["end"]):
                protein = self._particle_info[bond["b"]]
                if protein["species"] == spec.get("activator_species"):
                    activators.add(protein["molecule"])
        return float(spec.get("on_rate", 0)) * (1 + len(activators))

    def _regulatory_binding_context(self):
        """Index actual bonds once for a chemistry step or read-only report."""
        direct_sites, by_species, partners = {}, {}, {}
        for bond in self.bindings:
            if bond["kind"] == "chromatin":
                particle = self._particle_info[bond["b"]]
                molecule = particle["molecule"]
                sites = direct_sites.setdefault(molecule, set())
                sites.add(int(bond["a"]))
                by_species.setdefault(particle["species"], {})[molecule] = sites
            elif bond["kind"] == "protein":
                left = self._particle_info[bond["a"]]["molecule"]
                right = self._particle_info[bond["b"]]["molecule"]
                partners.setdefault(left, set()).add(right)
                partners.setdefault(right, set()).add(left)
        return {"direct_sites": direct_sites, "by_species": by_species, "partners": partners}

    def _region_recruitment(self, region, context):
        """Return distinct bound TFs and one-hop coactivator anchor sets."""
        start, end = int(region["start"]), int(region["end"])
        tf_sites, recruited = {}, {}
        for molecule, sites in context["by_species"].get(region.get("tf_species"), {}).items():
            anchors = {bead for bead in sites if start <= bead <= end}
            if not anchors:
                continue
            tf_sites[molecule] = anchors
            coactivator = region.get("coactivator_species")
            if not coactivator:
                continue
            for partner in context["partners"].get(molecule, set()):
                if self.molecules[partner]["species"] == coactivator:
                    recruited.setdefault(partner, set()).update(anchors)
        return tf_sites, recruited

    def _regulatory_region_records(self, *, include_marks=False):
        """Measure regions from present positions and bonds, never cached state."""
        context = self._regulatory_binding_context()
        genes = {gene["config"]["id"]: gene["config"] for gene in self.genes}
        promoters = {region["target_gene"]: region for region in self._enabled_regions
                     if region["kind"] == "promoter"}
        promoter_trees = {}
        records = []
        for region in self.regulatory_regions:
            start, end = int(region["start"]), int(region["end"])
            spec = genes[region["target_gene"]]
            promoter = promoters.get(region["target_gene"])
            target_start = int(promoter["start"]) if promoter else int(spec["start"])
            target_end = int(promoter["end"]) if promoter else target_start
            if region["kind"] == "promoter":
                distance, contact = 0.0, True
            else:
                key = (target_start, target_end)
                if key not in promoter_trees:
                    promoter_trees[key] = cKDTree(self.positions[target_start:target_end + 1])
                distances, _ = promoter_trees[key].query(self.positions[start:end + 1], k=1)
                distance = float(np.min(distances))
                contact = distance <= float(region.get("contact_radius", 0))
            tf_sites, recruited = self._region_recruitment(region, context)
            tf_count, recruited_count = len(tf_sites), len(recruited)
            signal = tf_count + float(region.get("coactivator_gain", 0)) * recruited_count
            enabled = bool(region.get("enabled", True))
            active = enabled and contact and signal > 0
            record = {"id": region["id"], "name": region.get("name", region["id"]),
                      "kind": region["kind"], "start": start, "end": end,
                      "target_gene": region["target_gene"], "enabled": enabled,
                      "tf_species": region.get("tf_species"),
                      "coactivator_species": region.get("coactivator_species"),
                      "tf_molecules": sorted(tf_sites), "tf_count": tf_count,
                      "recruited_molecules": sorted(recruited), "recruited_count": recruited_count,
                      "tf_bound_beads": sorted({bead for sites in tf_sites.values() for bead in sites}),
                      "recruited_anchors": [{"molecule": molecule, "beads": sorted(recruited[molecule])}
                                            for molecule in sorted(recruited)],
                      "contact": bool(contact), "distance": distance,
                      "target_promoter": {"start": target_start, "end": target_end,
                                          "region_id": promoter["id"] if promoter else None},
                      "signal": signal, "active": bool(active),
                      "activation_contribution": float(region.get("activation_gain", 0)) * signal if active else 0.0,
                      "initiation_contribution": float(region.get("initiation_gain", 0)) * signal if active else 0.0}
            if include_marks:
                record["mark_fractions"] = {mark: float(np.mean(flags[start:end + 1]))
                                            for mark, flags in self.marks.items()}
            records.append(record)
        return records

    def _gene_regulatory_rates(self, spec, region_records):
        required, repressive = spec.get("required_mark"), spec.get("repressive_mark")
        promoter = int(spec["start"])
        required_present = bool(self.marks[required][promoter]) if required else None
        repressive_present = bool(self.marks[repressive][promoter]) if repressive else None
        gate = (required_present is not False) and (repressive_present is not True)
        base_on, base_initiation = float(spec.get("on_rate", 0)), float(spec.get("initiation_rate", 0))
        regional = spec["id"] in self._regions_by_gene
        contributions = [{key: region[key] for key in ("contact", "signal", "activation_contribution", "initiation_contribution")}
                         | {"region_id": region["id"]}
                         for region in region_records if region["target_gene"] == spec["id"] and region["enabled"]]
        legacy_activators = set()
        if regional:
            activation = 1 + sum(region["activation_contribution"] for region in contributions)
            initiation = 1 + sum(region["initiation_contribution"] for region in contributions)
            on_rate = base_on * activation if gate else 0.0
        else:
            for bond in self.bindings:
                if bond["kind"] == "chromatin" and int(spec["start"]) <= bond["a"] <= int(spec["end"]):
                    particle = self._particle_info[bond["b"]]
                    if particle["species"] == spec.get("activator_species"):
                        legacy_activators.add(particle["molecule"])
            activation, initiation = 1 + len(legacy_activators), 1.0
            on_rate = self._legacy_gene_on_rate(spec)
        # Optional environmental responses are explicit hypotheses independent
        # of body forces and RNA processing. Log interpolation makes factor 1
        # neutral and treats promoting/repressive multipliers symmetrically.
        body_responses = []
        activation_log, initiation_log = 0.0, 0.0
        for body in self.bodies:
            response = spec.get("body_responses", {}).get(body.get("id"))
            if response is None:
                continue
            distance = float(np.linalg.norm(self.positions[promoter] - np.asarray(body["center"])))
            weight = math.exp(-0.5 * (distance / float(body["radius"])) ** 2)
            am = float(response.get("activation_multiplier", 1))
            im = float(response.get("initiation_multiplier", 1))
            if not (math.isfinite(am) and math.isfinite(im) and am > 0 and im > 0):
                raise ValueError("Body response multipliers must be positive and finite.")
            af, inf = weight * math.log(am), weight * math.log(im)
            activation_log += af
            initiation_log += inf
            body_responses.append({"body_id": body["id"], "kind": body.get("kind"),
                                   "distance": distance, "weight": weight,
                                   "activation_multiplier": am, "initiation_multiplier": im,
                                   "activation_factor": math.exp(af), "initiation_factor": math.exp(inf)})
        if max(activation_log, initiation_log) > 600:
            raise ValueError("Combined nuclear-body response is too large; reduce response multipliers or overlapping bodies.")
        environment_activation = math.exp(activation_log)
        environment_initiation = math.exp(initiation_log)
        on_rate *= environment_activation
        initiation_rate = base_initiation * initiation * environment_initiation
        if not math.isfinite(on_rate + initiation_rate):
            raise ValueError("Transcription rates overflowed; reduce regulatory gains or body response multipliers.")
        explicit = self.explicit_transcription is not None
        causal_terms = {
            "on_rate": {"base": base_on, "regulatory_multiplier": activation,
                        "environment_multiplier": environment_activation,
                        "mark_gate": bool(gate), "effective": on_rate},
            "initiation_rate": {"base": base_initiation, "regulatory_multiplier": initiation,
                                "environment_multiplier": environment_initiation, "effective": initiation_rate},
            "off_rate": {"effective": float(spec.get("off_rate", 0))},
            "elongation": {"base_beads_per_time": float(spec.get("elongation_rate", 0)),
                           "occupied_multiplier": float(spec.get("occupied_elongation_multiplier", 1)),
                           "pause_rate": float(spec.get("pause_rate", 0)) if explicit else 0.0,
                           "resume_rate": float(spec.get("resume_rate", 1)) if explicit else 0.0,
                           "occupancy_applies": bool(explicit and self.nucleosomes_enabled)},
            "processing": {"base": float(spec.get("processing_rate", 0)),
                           "speckle_gain": float(spec.get("speckle_gain", 0)),
                           "contact_mode": "released_rna_position" if explicit else "gene_contact"},
            "nucleosomes": {"enabled": self.nucleosomes_enabled,
                            "promoter_occupied": bool(self.nucleosome_occupied[promoter]) if self.nucleosomes_enabled else None,
                            "gene_occupied_fraction": float(np.mean(self.nucleosome_occupied[int(spec["start"]):int(spec["end"]) + 1])) if self.nucleosomes_enabled else None}}
        return {"id": spec["id"], "mode": "regions" if regional else "legacy",
                "base_on_rate": base_on, "effective_on_rate": on_rate,
                "base_initiation_rate": base_initiation,
                "effective_initiation_rate": initiation_rate,
                "mark_gate_passed": bool(gate), "required_mark_present": required_present,
                "repressive_mark_present": repressive_present,
                "activation_multiplier": activation, "initiation_multiplier": initiation,
                "environment_activation_multiplier": environment_activation,
                "environment_initiation_multiplier": environment_initiation,
                "body_responses": body_responses, "causal_terms": causal_terms,
                "legacy_activator_count": len(legacy_activators), "contributions": contributions}

    def regulation_state(self):
        """Read current regulatory mechanisms and rates without consuming RNG.

        Gain parameters are explicit modeling hypotheses. Contact only removes
        regional rate gains when lost; it does not forcibly turn a promoter
        OFF. Existing promoter switching and chromatin marks retain memory.
        """
        regions = self._regulatory_region_records(include_marks=True) if self.regulatory_regions else []
        return {"regions": regions,
                "genes": [self._gene_regulatory_rates(gene["config"], regions) for gene in self.genes]}

    def _gene_on_rate(self, spec):
        if spec["id"] not in self._regions_by_gene and not spec.get("body_responses"):
            return self._legacy_gene_on_rate(spec)
        return self._gene_regulatory_rates(spec, self._regulatory_region_records())["effective_on_rate"]

    def _gene_speckle_contact(self, spec):
        locus = self.positions[int(spec["start"]):int(spec["end"]) + 1]
        for body in self.bodies:
            if body.get("kind") == "speckle":
                if np.any(np.linalg.norm(locus - np.asarray(body["center"]), axis=1) <= float(body["radius"]) + 0.5 * self.bond_length):
                    return True
        return False

    def _advance_genes(self, dt):
        from .transcription import MAX_TRANSCRIPTION_EVENTS, sample_promoter_interval
        regions = self._regulatory_region_records() if self._enabled_regions else []
        rates = {gene["config"]["id"]: self._gene_regulatory_rates(gene["config"], regions)
                 for gene in self.genes}
        if self.explicit_transcription is not None:
            self.explicit_transcription.advance(dt, rates)
            return
        # Sample first, so a pathological event count cannot partially update
        # the gene population. Rates are fixed during this reaction interval.
        plans, budget = [], MAX_TRANSCRIPTION_EVENTS
        for gene in self.genes:
            spec = gene["config"]
            report = rates[spec["id"]]
            plan = sample_promoter_interval(self.rng, gene["on"], report["effective_on_rate"],
                                            float(spec.get("off_rate", 0)),
                                            report["effective_initiation_rate"], dt, max_events=budget)
            budget -= len(plan["switches"]) + len(plan["initiation_times"])
            plans.append(plan)
        for gene, plan in zip(self.genes, plans):
            spec = gene["config"]
            gene["on"] = plan["on"]
            for offset, on in plan["switches"]:
                self._event("promoter", spec["id"] + (" ON" if on else " OFF"), time=self.time + offset)
            births = plan["initiation_times"]
            gene["last_transcription_interval"] = {
                "start_time": float(self.time), "duration": float(dt), "on_time": plan["on_time"],
                "on_fraction": plan["on_time"] / dt if dt else 0.0,
                "switch_count": len(plan["switches"]), "initiation_attempts": len(births),
                "accepted_initiations": len(births), "blocked_initiations": 0}
            speed = float(spec.get("elongation_rate", 0))
            end = int(spec["end"]) + 1
            releases, polymerases = [], []
            for position, birth_time in ([(p, 0.0) for p in gene["polymerases"]]
                                         + [(float(spec["start"]), time) for time in births]):
                available = dt - birth_time
                completion_time = (end - position) / speed if speed > 0 else math.inf
                if completion_time <= available:
                    releases.append(birth_time + completion_time)
                else:
                    polymerases.append(position + speed * available)
            gene["polymerases"] = polymerases
            gene["initiated"] += len(births)
            gene["completed"] += len(releases)
            if releases:
                self._event("transcription", f"{spec['id']}: {len(releases)} completed", time=self.time + dt)
            processing = float(spec.get("processing_rate", 0))
            if self._gene_speckle_contact(spec):
                processing *= 1 + float(spec.get("speckle_gain", 0))
            decay = float(spec.get("rna_decay_rate", 0))
            # With equal decay hazards before/after maturation, independent
            # survival and processing draws give the exact constant-rate
            # count transition. Newly matured RNA does not decay twice.
            def survivors(count, duration):
                if not count or duration <= 0:
                    return count, 0
                alive = int(self.rng.binomial(count, math.exp(-decay * duration))) if decay else count
                unprocessed = int(self.rng.binomial(alive, math.exp(-processing * duration))) if processing and alive else alive
                return unprocessed, alive - unprocessed
            if gene["mature"] and decay:
                gene["mature"] = int(self.rng.binomial(gene["mature"], math.exp(-decay * dt)))
            gene["unprocessed"], matured = survivors(gene["unprocessed"], dt)
            gene["mature"] += matured
            for release_time in releases:
                unprocessed, matured = survivors(1, dt - release_time)
                gene["unprocessed"] += unprocessed
                gene["mature"] += matured

    def advance_reactions(self, dt):
        # Binding and marks see the same post-dynamics coordinates. Reuse this
        # tree within the reaction interval, never across a Brownian move.
        need_tree = bool(self.protein_particles) or any(float(m.get("spread_rate", 0)) > 0 for m in self.mark_config.values())
        chromatin_tree = cKDTree(self.positions[:self.n_chromatin]) if need_tree else None
        self._advance_binding(dt, chromatin_tree)
        self._advance_nucleosomes(dt)
        self._advance_nucleosome_sliding(dt)
        self._advance_marks(dt, chromatin_tree)
        self._advance_cohesin(dt)
        self._advance_tethers(dt)
        self._advance_genes(dt)

    def advance_bookkeeping(self, dt):
        self.time += float(dt)
        self.step_count += 1

    def step(self, steps):
        if isinstance(steps, bool) or not isinstance(steps, (int, np.integer)) or steps < 0:
            raise ValueError("steps must be a nonnegative integer")
        for _ in range(int(steps)):
            if self.dynamics is None:
                self._brownian_move()
            else:
                self.dynamics(self)
                self.validate_positions()
            self.advance_reactions(self.dt)
            self.advance_bookkeeping(self.dt)

    def _largest_cluster(self):
        parent = list(range(len(self.molecules)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for bond in self.bindings:
            if bond["kind"] == "protein":
                a = self._particle_info[bond["a"]]["molecule"]
                b = self._particle_info[bond["b"]]["molecule"]
                parent[find(a)] = find(b)
        if not parent:
            return 0
        counts = np.bincount([find(i) for i in range(len(parent))])
        return int(np.max(counts))

    def genome_tracks(self):
        """Return authoritative, read-only one-dimensional chromatin tracks.

        ``affinity_runs`` and ``occupancy_runs`` map every configured species ID
        to runs ``[start_bead, end_bead_exclusive, value]``. Runs cover [0, N),
        including zero-valued intervals. Bead i represents the genomic interval
        [i*bp_per_bead, (i+1)*bp_per_bead); span_bp is N*bp_per_bead.

        Affinity values are the exact dimensionless encounter-rate multipliers
        used by binding chemistry, including the max-present-mapped-mark rule
        and unmodified fallback. They are neither occupancies nor dissociation
        constants. Occupancy counts distinct directly bound protein molecules
        of that species at each chromatin bead. Multiple IDR domain links from
        one molecule to the same bead count once. Protein-protein links do not
        imply direct chromatin occupancy. No state, event log, or RNG is changed.
        """
        occupancy = {species_id: np.zeros(self.n_chromatin, dtype=np.int64)
                     for species_id in self.species}
        binding_sites = []
        seen = set()
        for bond in self.bindings:
            if bond["kind"] != "chromatin":
                continue
            particle = self._particle_info[bond["b"]]
            binding_sites.append({"bead": int(bond["a"]), "particle": int(bond["b"]),
                                  "molecule": int(particle["molecule"]), "species": particle["species"],
                                  "off_rate": self._bond_off_rate(bond)})
            key = (bond["a"], particle["molecule"])
            if key not in seen:
                occupancy[particle["species"]][bond["a"]] += 1
                seen.add(key)
        bp_per_bead = self.chrom.get("bp_per_bead", 1)
        return {"bp_per_bead": bp_per_bead,
                "span_bp": self.n_chromatin * bp_per_bead,
                "nucleosome_occupancy_runs": (_run_length_encode(self.nucleosome_occupied.astype(int))
                                             if self.nucleosomes_enabled else []),
                "persistence_runs": _run_length_encode(self.local_persistence()),
                "affinity_runs": {species_id: _run_length_encode(self._chromatin_affinity(species_id))
                                  for species_id in self.species},
                "on_rate_runs": {species_id: _run_length_encode(float(protein.get("on_rate", 0)) *
                                                                self._chromatin_affinity(species_id))
                                 for species_id, protein in self.species.items()},
                "off_rate_runs": {species_id: _run_length_encode(self._chromatin_off_rates(species_id))
                                  for species_id in self.species},
                "binding_sites": binding_sites,
                "tether_runs": _run_length_encode(self._tether_attached_array()),
                "tether_region_runs": _run_length_encode(self._tether_region_array()),
                "occupancy_runs": {species_id: _run_length_encode(values)
                                   for species_id, values in occupancy.items()}}

    def _tether_attached_array(self):
        values = np.zeros(self.n_chromatin, dtype=np.int64)
        for tether in self.tethers:
            values[tether["bead"]] = 1
        return values

    def _tether_region_array(self):
        """0 = no active tethering element; i+1 = owned by active element i."""
        values = np.zeros(self.n_chromatin, dtype=np.int64)
        for bead, index in self._tether_owner.items():
            values[bead] = index + 1
        return values

    def interaction_state(self):
        """Summarize the current molecular binding graph without changing state.

        A molecule is ``direct`` when it has its own chromatin bond,
        ``indirect`` when it only reaches chromatin through protein bonds, and
        ``free`` when its protein-connected component has no chromatin bond.
        Protein pair contexts use *direct* binding only: ``free_free`` means
        neither partner has its own chromatin bond, even if indirectly tethered.
        Link counts count domain bonds; molecule-pair counts deduplicate all
        domain bonds between the same two molecules. Permanent IDR backbone
        bonds are not intermolecular links. Shared chromatin beads do not merge
        protein components. Component IDs are minimum member molecule IDs in
        this frame, not persistent assembly lineages.
        """
        count = len(self.molecules)
        parents = list(range(count))
        ranks = [0] * count
        chromatin_sites = [set() for _ in range(count)]
        used_valence = [0] * count
        molecular_links = {}

        def find(index):
            while parents[index] != index:
                parents[index] = parents[parents[index]]
                index = parents[index]
            return index

        def union(left, right):
            left, right = find(left), find(right)
            if left == right:
                return
            if ranks[left] < ranks[right]:
                left, right = right, left
            parents[right] = left
            if ranks[left] == ranks[right]:
                ranks[left] += 1

        for bond in self.bindings:
            for particle in (bond["a"], bond["b"]):
                if particle in self._particle_info:
                    used_valence[self._particle_info[particle]["molecule"]] += 1
            if bond["kind"] == "chromatin":
                molecule = self._particle_info[bond["b"]]["molecule"]
                chromatin_sites[molecule].add(int(bond["a"]))
            elif bond["kind"] == "protein":
                left = self._particle_info[bond["a"]]["molecule"]
                right = self._particle_info[bond["b"]]["molecule"]
                pair = tuple(sorted((left, right)))
                molecular_links[pair] = molecular_links.get(pair, 0) + 1
                union(left, right)

        by_root = {}
        for molecule in range(count):
            by_root.setdefault(find(molecule), []).append(molecule)
        components = []
        component_by_molecule = {}
        for members in sorted(by_root.values(), key=lambda group: group[0]):
            species_counts = {}
            sites = set()
            direct = 0
            for molecule in members:
                species = self.molecules[molecule]["species"]
                species_counts[species] = species_counts.get(species, 0) + 1
                sites.update(chromatin_sites[molecule])
                direct += bool(chromatin_sites[molecule])
            component = {"id": members[0], "molecule_ids": members, "size": len(members),
                         "species_counts": species_counts, "direct_bound": direct,
                         "indirect_bound": len(members) - direct if direct else 0,
                         "chromatin_bound": bool(direct), "protein_links": 0,
                         "unique_molecule_pairs": 0, "chromatin_beads": sorted(sites)}
            components.append(component)
            for molecule in members:
                component_by_molecule[molecule] = component

        species_counts = {species: {"id": species, "total": 0, "direct_bound": 0,
                                   "indirect_bound": 0, "free": 0} for species in self.species}
        molecules = []
        for index, molecule in enumerate(self.molecules):
            component = component_by_molecule[index]
            status = ("direct" if chromatin_sites[index] else
                      "indirect" if component["chromatin_bound"] else "free")
            species_counts[molecule["species"]]["total"] += 1
            species_counts[molecule["species"]][status + "_bound" if status != "free" else status] += 1
            molecules.append({"id": index, "species": molecule["species"],
                              "particle_ids": list(molecule["indices"]),
                              "component_id": component["id"], "binding_status": status,
                              "chromatin_beads": sorted(chromatin_sites[index]),
                              "used_valence": used_valence[index], "valence": molecule["valence"]})

        def empty_counts():
            return {"total": 0, "free_free": 0, "bound_free": 0, "bound_bound": 0}

        def pair_record(pair):
            return {"a": pair[0], "b": pair[1], "configured": pair in self.interactions,
                    "links": empty_counts(), "molecule_pairs": empty_counts()}

        pairs = {pair: pair_record(pair) for pair in self.interactions}
        links = []
        for (left, right), link_count in sorted(molecular_links.items()):
            direct = int(bool(chromatin_sites[left])) + int(bool(chromatin_sites[right]))
            context = ("free_free", "bound_free", "bound_bound")[direct]
            species_pair = tuple(sorted((self.molecules[left]["species"], self.molecules[right]["species"])))
            if species_pair not in pairs:
                pairs[species_pair] = pair_record(species_pair)
            record = pairs[species_pair]
            record["links"][context] += link_count
            record["links"]["total"] += link_count
            record["molecule_pairs"][context] += 1
            record["molecule_pairs"]["total"] += 1
            component = component_by_molecule[left]
            component["protein_links"] += link_count
            component["unique_molecule_pairs"] += 1
            links.append({"a": left, "b": right, "links": link_count, "context": context})

        totals = {"molecules": count,
                  "direct_bound": sum(item["direct_bound"] for item in species_counts.values()),
                  "indirect_bound": sum(item["indirect_bound"] for item in species_counts.values()),
                  "free": sum(item["free"] for item in species_counts.values()),
                  "protein_links": sum(molecular_links.values()), "molecule_pairs": len(molecular_links),
                  "components": len(components), "assemblies": sum(item["size"] > 1 for item in components)}
        return {"definitions": {
                    "binding_status": "direct: own chromatin bond; indirect: protein-bond path to a directly bound molecule; free: no such path",
                    "pair_context": "free_free / bound_free / bound_bound count zero / one / two partners with their own chromatin bond; indirect tethering is not direct binding",
                    "component_id": "minimum member molecule ID in this frame; not an assembly lineage"},
                "species": list(species_counts.values()), "molecules": molecules,
                "pairs": [pairs[key] for key in sorted(pairs)], "molecule_links": links,
                "components": components, "totals": totals}

    def confinement_state(self):
        """Read soft-wall occupancy without moving particles or consuming RNG.

        Coordinates are never projected back inside the reference sphere. A
        zero stiffness disables its force. Released RNA is counted separately;
        attached nascent RNA is excluded to avoid double counting chromatin.
        """
        def measure(positions):
            positions = np.asarray(positions, dtype=float).reshape((-1, 3))
            radial = np.linalg.norm(positions, axis=1)
            penetration = np.maximum(radial - self.radius, 0)
            outside = int(np.count_nonzero(penetration))
            return {"total": len(radial), "outside_count": outside,
                    "outside_fraction": outside / len(radial) if len(radial) else 0.0,
                    "max_radius": float(np.max(radial)) if len(radial) else 0.0,
                    "max_penetration": float(np.max(penetration)) if len(radial) else 0.0,
                    "energy": float(0.5 * self.confinement_stiffness * np.sum(penetration ** 2))}
        particles = measure(self.positions)
        rna_positions = ([record["position"] for record in self.explicit_transcription.transcripts.values()
                          if record["state"] != "nascent"] if self.explicit_transcription else [])
        rna = measure(rna_positions)
        return {"kind": "soft_harmonic_sphere", "enabled": self.confinement_stiffness > 0,
                "radius": self.radius, "stiffness": self.confinement_stiffness,
                "energy": particles["energy"] + rna["energy"], "particles": particles, "rna": rna}

    def snapshot(self):
        chromatin = self.positions[:self.n_chromatin]
        bound = {self._particle_info[b["b"]]["molecule"] for b in self.bindings if b["kind"] == "chromatin"}
        marked = np.any(np.asarray(list(self.marks.values())), axis=0) if self.marks else np.zeros(self.n_chromatin, bool)
        rg = float(np.sqrt(np.mean(np.sum((chromatin - chromatin.mean(axis=0)) ** 2, axis=1))))
        bonds = [{"a": self.positions[b["a"]].tolist(), "b": self.positions[b["b"]].tolist(), "kind": b["kind"]}
                 for b in self.bindings]
        bonds.extend({"a": self.positions[a].tolist(), "b": self.positions[b].tolist(), "kind": "cohesin"} for a, b in self.cohesin)
        bonds.extend({"a": self.positions[t["bead"]].tolist(), "b": list(t["anchor"]), "kind": "tether"} for t in self.tethers)
        regulation = self.regulation_state()
        gene_regulation = {gene["id"]: gene for gene in regulation["genes"]}
        chromatin_contacts = (self._chromatin_contact_state()["report"] if self.chrom.get("state_interactions")
                              else {"rules": [], "total_contacts": 0, "energy": 0.0})
        if self.explicit_transcription is not None:
            transcription = self.explicit_transcription.snapshot()
        else:
            from .transcription import disabled_transcription_state
            transcription = disabled_transcription_state()
        return {"time": float(self.time), "step": self.step_count, "engine": self.engine,
                "chromatin": chromatin.tolist(),
                "proteins": [{"id": p["index"], "species": p["species"], "molecule": p["molecule"],
                              "position": self.positions[p["index"]].tolist(), "sticker": p["sticker"]} for p in self.protein_particles],
                "protein_chains": copy.deepcopy(self.protein_chains), "bonds": bonds,
                "marks": {key: flags.astype(int).tolist() for key, flags in self.marks.items()},
                "genome_tracks": self.genome_tracks(),
                "interaction_state": self.interaction_state(),
                "regulation": regulation, "chromatin_contacts": chromatin_contacts,
                "nucleosomes": self.nucleosome_state(),
                "genomic_events": self.genomic_event_state(),
                "transcription": transcription, "confinement": self.confinement_state(),
                "nucleus": self.nucleus_state(),
                "cohesin": copy.deepcopy(self.cohesin),
                "genes": [{"id": g["config"]["id"], "name": g["config"].get("name", g["config"]["id"]),
                           "on": bool(g["on"]), "nascent": len(g["polymerases"]), "unprocessed": g["unprocessed"],
                           "mature": g["mature"], "initiated": g["initiated"], "completed": g["completed"],
                           "polymerase_positions": list(g["polymerases"]),
                           "last_transcription_interval": copy.deepcopy(g.get("last_transcription_interval")),
                           "regulation": gene_regulation[g["config"]["id"]],
                           "speckle_contact": self._gene_speckle_contact(g["config"])} for g in self.genes],
                "metrics": {"bound_fraction": len(bound) / len(self.molecules) if self.molecules else 0,
                            "largest_cluster": self._largest_cluster(), "marked_fraction": float(np.mean(marked)),
                            "rg": rg, "mean_persistence": float(np.mean(self.local_persistence())),
                            "protein_count": len(self.molecules), "particle_count": len(self.positions),
                            "cohesin_bound": len(self.cohesin),
                            "cohesin_free": max(0, int(self.cohesin_config.get("count", 0)) - len(self.cohesin)) if self.cohesin_config.get("enabled") else 0,
                            "integration_substeps": self.last_substeps},
                "events": list(self.events), "warnings": list(self.warnings)}
