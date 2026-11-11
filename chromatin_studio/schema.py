"""Versioned model configuration and strict validation for local simulations.

Chromatin regions use zero-based, inclusive bead indices. At resolution b bp,
bead i covers genomic coordinates [i*b, (i+1)*b); a region [start, end] covers
[start*b, (end+1)*b). Genomic coordinates label the model, not its time units.
"""

import copy
import math
import re


MAX_CHROMATIN_BEADS = 4000
MAX_PARTICLES = 6000
MAX_GENES = 256
MAX_COHESIN = 256
MAX_REGULATORY_REGIONS = 128
MAX_TRANSCRIPTS = 5000


LIMITATIONS = [
    "Dynamics are stored in reduced units. Physical interpretation requires a documented length scale and an MSD-derived clock; attaching a mapping does not validate the biological model.",
    "Speckles and nucleoli are imposed stationary affinity fields. Protein assemblies can form independently; nucleolar subphases and rRNA processing are not modeled.",
    "IDRs are explicit domain-bead sticker chains, not residue-level or sequence-parametrized CALVADOS proteins.",
    "Independent particle drag omits hydrodynamic coupling; soft excluded volume can allow chain crossings.",
    "Histone marks are independent binary flags per bead; a nucleosome-unit bead can carry these flags but does not resolve histone tails, individual residues, or wrapped DNA.",
    "Regulatory gains and state-selective chromatin attraction are configurable hypotheses; protein preferences are encounter multipliers, not measured dissociation constants.",
    "Transcription-dependent persistence and processing enhancement are configurable hypotheses, not universal biological laws.",
    "Assembly and gene parameters are illustrative. Seeing a cluster is not evidence of calibrated liquid-liquid phase separation.",
]


def _protein(identifier, name, color, count, **overrides):
    value = dict(id=identifier, name=name, color=color, count=count, diffusion=0.5,
                 radius=0.28, valence=3, representation="bead", idr_beads=4,
                 sticker_fraction=0.75, chromatin_affinity={"unmodified": 0.15},
                 on_rate=1.2, off_rate=0.4, writer_mark=None, eraser_mark=None,
                 catalytic_rate=0.8)
    value.update(overrides)
    return value


def default_config():
    return {
        "version": "0.1", "name": "Chromatin × nuclear bodies", "seed": 42,
        "simulation": {"dt": 0.002, "steps_per_frame": 10, "confinement_radius": 11.0},
        "chromatin": {"beads": 96, "bp_per_bead": 1000, "bond_length": 0.85,
                      "bond_k": 45.0, "diffusion": 0.08, "persistence_length": 2.0,
                      "repulsion": 8.0,
                      "initial_marks": [{"start": 8, "end": 25, "mark": "H3K9me3"},
                                        {"start": 40, "end": 58, "mark": "H3K27ac"}],
                      "barriers": [{"bead": 30, "direction": "right"},
                                   {"bead": 69, "direction": "left"}]},
        "marks": [
            {"id": "H3K9me3", "name": "H3K9me3", "color": "#8064b0", "decay_rate": 0.015, "spread_rate": 0.04, "contact_radius": 1.25},
            {"id": "H3K27ac", "name": "H3K27ac", "color": "#dfaa49", "decay_rate": 0.02, "spread_rate": 0.01, "contact_radius": 1.25},
        ],
        "proteins": [
            _protein("hp1", "HP1-like reader", "#8d71bc", 12,
                     chromatin_affinity={"unmodified": 0.04, "H3K9me3": 4.0, "H3K27ac": 0.04}),
            _protein("writer", "H3K9 writer", "#b196cf", 6, writer_mark="H3K9me3",
                     chromatin_affinity={"unmodified": 0.1, "H3K9me3": 3.0}),
            _protein("srrm2", "SRRM2-like scaffold", "#399c97", 8, representation="idr", idr_beads=4,
                     valence=5, diffusion=0.35, chromatin_affinity={"unmodified": 0.08, "H3K27ac": 0.5}),
            _protein("npm1", "NPM1-like scaffold", "#c67c68", 7, representation="idr", idr_beads=4,
                     valence=5, diffusion=0.35, chromatin_affinity={"unmodified": 0.05}),
            _protein("activator", "Activator", "#d6a34a", 10, valence=2, diffusion=0.8,
                     chromatin_affinity={"unmodified": 0.3, "H3K27ac": 3.0, "H3K9me3": 0.02}),
            _protein("p300", "Acetylation writer", "#bfa357", 4, writer_mark="H3K27ac",
                     chromatin_affinity={"unmodified": 0.25, "H3K27ac": 2.0}),
        ],
        "interactions": [
            {"a": "hp1", "b": "hp1", "strength": 4.0, "on_rate": 1.2, "off_rate": 0.35, "radius": 1.1},
            {"a": "srrm2", "b": "srrm2", "strength": 4.0, "on_rate": 1.8, "off_rate": 0.3, "radius": 1.1},
            {"a": "npm1", "b": "npm1", "strength": 4.0, "on_rate": 1.8, "off_rate": 0.3, "radius": 1.1},
            {"a": "srrm2", "b": "activator", "strength": 2.0, "on_rate": 1.0, "off_rate": 0.5, "radius": 1.1},
        ],
        "bodies": [
            {"id": "speckle", "name": "Nuclear speckle", "kind": "speckle", "color": "#4caaa3",
             "center": [-3.0, 2.0, 0.0], "radius": 2.6, "strength": 1.0,
             "chromatin_affinity": 0.3, "protein_affinity": {"srrm2": 4.0, "activator": 2.0, "p300": 1.0}, "enabled": True},
            {"id": "nucleolus", "name": "Nucleolus", "kind": "nucleolus", "color": "#c17b67",
             "center": [3.4, -2.1, 0.0], "radius": 2.7, "strength": 1.0,
             "chromatin_affinity": -0.4, "protein_affinity": {"npm1": 4.0, "hp1": 0.5}, "enabled": True},
        ],
        "cohesin": {"enabled": True, "count": 3, "step_rate": 1.5, "off_rate": 0.08, "spring_k": 2.5},
        "genes": [{"id": "gene_a", "name": "Gene A", "start": 40, "end": 58,
                   "on_rate": 0.3, "off_rate": 0.2, "initiation_rate": 1.0,
                   "elongation_rate": 3.0, "processing_rate": 0.5, "rna_decay_rate": 0.03,
                   "activator_species": "activator", "required_mark": "H3K27ac",
                   "repressive_mark": "H3K9me3", "speckle_gain": 2.0,
                   "persistence_multiplier": 0.65}],
    }


PRESETS = [
    {"id": "nuclear", "name": "Nuclear neighborhood", "description": "Six species, two nuclear landmarks, writers, IDRs and cohesin."},
    {"id": "assembly", "name": "Multivalent assembly", "description": "Explicit flexible sticker chains with higher-order binding; no imposed bodies."},
    {"id": "epigenetic", "name": "Reader–writer memory", "description": "Mark-dependent binding and spreading on a chromatin polymer."},
    {"id": "extrusion", "name": "Loop extrusion", "description": "Cohesin movement with directional barriers and a transcribing gene."},
    {"id": "twenty_mb", "name": "20 Mb chromatin", "description": "2,000 beads at 10 kb each, six protein species, six gene regions and 24 cohesins; illustrative reduced-unit dynamics."},
    {"id": "regulatory_20mb", "name": "20 Mb regulatory landscape", "description": "Enhancers, promoters, TF–coactivator recruitment and mark-selective attraction at 10 kb per bead; sub-bead elements are unresolved."},
    {"id": "nucleosome_locus", "name": "80 kb nucleosome-unit locus", "description": "400 beads at a chosen 200 bp repeat, with enhancers, promoters and coactivators; coarse nucleosome units, without wrapped-DNA geometry."},
    {"id": "nucleosome_transcription", "name": "Nucleosome access & RNA · 80 kb", "description": "Dynamic occupancy, TF accessibility, local stiffness, individual polymerases and diffusing RNA at 200 bp per bead."},
    {"id": "nucleosome_access", "name": "Live nucleosome access · 25.6 kb", "description": "128 nucleosome units, open regions, distinct TF access and residence preferences, and a recruited sliding/evicting remodeler. Illustrative kinetics."},
    {"id": "transcription_20mb", "name": "20 Mb transcription landscape", "description": "Individual polymerases and RNA with effective occupancy/accessibility states in 10 kb bins; individual nucleosomes remain unresolved."},
    {"id": "manuscript_locus", "name": "Manuscript hypotheses · 80 kb", "description": "Coarse nucleosome units, TF/coactivators, acetylation and Polycomb hypotheses, and optional nuclear-environment responses. Illustrative parameters, not a fit to the manuscript."},
    {"id": "manuscript_20mb", "name": "Manuscript hypotheses · 20 Mb", "description": "The same candidate mechanisms across a 20 Mb coarse chain; occupancy represents 10 kb bins, not individual nucleosomes. Use matched ablations."},
]


def preset_config(identifier):
    c = default_config()
    if identifier == "nuclear":
        return c
    if identifier == "assembly":
        c["name"] = "Multivalent IDR assembly"
        c["chromatin"]["beads"] = 64
        c["chromatin"]["barriers"] = []
        c["chromatin"]["initial_marks"] = []
        c["genes"] = []
        c["cohesin"]["enabled"] = False
        c["bodies"] = []
        c["proteins"] = [p for p in c["proteins"] if p["id"] in ("srrm2", "npm1")]
        for p in c["proteins"]:
            p["count"] = 16
            p["idr_beads"] = 5
            p["valence"] = 6
        c["interactions"] = [p for p in c["interactions"] if p["a"] in ("srrm2", "npm1") and p["b"] in ("srrm2", "npm1")]
        c["simulation"]["confinement_radius"] = 7.0
    elif identifier == "epigenetic":
        c["name"] = "Reader–writer chromatin"
        c["proteins"] = [p for p in c["proteins"] if p["id"] in ("hp1", "writer")]
        c["interactions"] = [p for p in c["interactions"] if p["a"] == "hp1"]
        c["bodies"] = []
        c["genes"] = []
        c["cohesin"]["enabled"] = False
    elif identifier == "extrusion":
        c["name"] = "Cohesin and gene mechanics"
        c["proteins"] = []
        c["interactions"] = []
        c["bodies"] = []
        c["cohesin"]["count"] = 5
        c["genes"][0]["activator_species"] = None
    elif identifier == "twenty_mb":
        c["name"] = "20 Mb chromatin"
        c["simulation"]["confinement_radius"] = 20.0
        c["simulation"]["steps_per_frame"] = 10
        c["chromatin"].update(
            beads=2000, bp_per_bead=10000,
            initial_marks=[
                {"start": 80, "end": 249, "mark": "H3K9me3"},
                {"start": 720, "end": 899, "mark": "H3K9me3"},
                {"start": 1450, "end": 1619, "mark": "H3K9me3"},
                {"start": 320, "end": 359, "mark": "H3K27ac"},
                {"start": 500, "end": 539, "mark": "H3K27ac"},
                {"start": 1050, "end": 1099, "mark": "H3K27ac"},
                {"start": 1250, "end": 1289, "mark": "H3K27ac"},
                {"start": 1750, "end": 1789, "mark": "H3K27ac"},
            ],
            barriers=[
                {"bead": 300, "direction": "right"}, {"bead": 650, "direction": "left"},
                {"bead": 1000, "direction": "right"}, {"bead": 1350, "direction": "left"},
                {"bead": 1700, "direction": "right"}, {"bead": 1900, "direction": "left"},
            ],
        )
        counts = {"hp1": 60, "writer": 24, "srrm2": 20, "npm1": 20,
                  "activator": 50, "p300": 20}
        for protein in c["proteins"]:
            protein["count"] = counts[protein["id"]]
        c["cohesin"]["count"] = 24
        c["bodies"][0].update(center=[-6.0, 3.0, 0.0], radius=4.0)
        c["bodies"][1].update(center=[6.0, -3.0, 0.0], radius=4.0)
        gene_template = c["genes"][0]
        c["genes"] = []
        for letter, start, end in [("a", 330, 341), ("b", 508, 515),
                                   ("c", 1060, 1074), ("d", 1260, 1269),
                                   ("e", 1760, 1779), ("f", 760, 774)]:
            gene = copy.deepcopy(gene_template)
            gene.update(id="gene_" + letter, name="Gene " + letter.upper(),
                        start=start, end=end)
            if letter == "f":
                gene.update(name="Gene F · repressed region", required_mark=None)
            c["genes"].append(gene)
    elif identifier in ("regulatory_20mb", "nucleosome_locus"):
        if identifier == "regulatory_20mb":
            c = preset_config("twenty_mb")
            c["name"] = "20 Mb regulatory landscape"
            c["chromatin"]["representation"] = "coarse"
        else:
            c["name"] = "80 kb nucleosome-unit locus"
            c["simulation"]["confinement_radius"] = 13.0
            c["chromatin"].update(
                beads=400, bp_per_bead=200, representation="nucleosome_unit",
                initial_marks=[{"start": 20, "end": 69, "mark": "H3K9me3"},
                               {"start": 110, "end": 159, "mark": "H3K27ac"},
                               {"start": 260, "end": 299, "mark": "H3K27ac"}],
                barriers=[{"bead": 80, "direction": "right"},
                          {"bead": 210, "direction": "left"},
                          {"bead": 240, "direction": "right"},
                          {"bead": 350, "direction": "left"}])
            template = c["genes"][0]
            c["genes"] = [dict(copy.deepcopy(template), id="gene_a", name="Gene A", start=140, end=199),
                          dict(copy.deepcopy(template), id="gene_b", name="Gene B", start=280, end=329)]
            c["cohesin"]["count"] = 6
        c["chromatin"]["state_interactions"] = [
            {"a": "H3K9me3", "b": "H3K9me3", "strength": 0.5, "radius": 1.5, "min_separation": 3}]
        c["interactions"].append({"a": "activator", "b": "p300", "strength": 2.0,
                                   "on_rate": 1.2, "off_rate": 0.4, "radius": 1.1})
        c["regulatory_regions"] = []
        for gene in c["genes"]:
            for kind, start, end in (("enhancer", gene["start"] - 20, gene["start"] - 18),
                                     ("promoter", gene["start"], gene["start"])):
                c["regulatory_regions"].append({
                    "id": gene["id"] + "_" + kind, "name": gene["name"] + " " + kind,
                    "color": "#648dc2" if kind == "enhancer" else "#d6a34a",
                    "kind": kind, "start": start, "end": end, "target_gene": gene["id"],
                    "tf_species": "activator", "coactivator_species": "p300",
                    "tf_affinity_multiplier": 3.0, "activation_gain": 0.5,
                    "initiation_gain": 0.25, "coactivator_gain": 1.0,
                    "contact_radius": 1.5, "enabled": True})
    elif identifier in ("nucleosome_transcription", "transcription_20mb"):
        c = preset_config("nucleosome_locus" if identifier == "nucleosome_transcription" else "regulatory_20mb")
        c["name"] = "Nucleosome access & RNA · 80 kb" if identifier == "nucleosome_transcription" else "20 Mb transcription landscape"
        c["chromatin"]["nucleosomes"] = {
            "enabled": True, "initial_occupancy": 1.0,
            "initial_regions": [{"start": r["start"], "end": r["end"], "occupied": False}
                                for r in c["regulatory_regions"]],
            "eviction_rate": 0.01, "reassembly_rate": 0.05,
            "occupied_persistence_multiplier": 1.0, "vacant_persistence_multiplier": 0.5}
        c["transcription"] = {"enabled": True, "max_transcripts": 2000,
                              "rna_diffusion": 0.3, "mature_rna_diffusion": 0.15}
        for protein in c["proteins"]:
            protein["nucleosome_affinity"] = {"occupied": 1.0, "vacant": 1.0}
            protein["nucleosome_eviction_rate"] = 0.0
            protein["nucleosome_reassembly_rate"] = 0.0
            if protein["id"] == "activator":
                protein["nucleosome_affinity"] = {"occupied": 0.1, "vacant": 1.0}
            if protein["id"] in ("hp1", "writer"):
                protein["nucleosome_affinity"]["vacant"] = 0.0
        # A separate illustrative remodeler avoids attributing eviction to p300 itself.
        c["proteins"].append(_protein("remodeler", "TF-recruited remodeler", "#5e9dbe", 8,
            valence=2, chromatin_affinity={"unmodified": 0.05},
            nucleosome_affinity={"occupied": 1.0, "vacant": 0.2},
            nucleosome_eviction_rate=1.0, nucleosome_reassembly_rate=0.0))
        c["interactions"].append({"a": "activator", "b": "remodeler", "strength": 2,
                                  "on_rate": 1.2, "off_rate": .4, "radius": 1.1})
        # Enhancers recruit the remodeler; promoters retain the acetylating
        # coactivator. Each element still describes one physical interval.
        for r in c["regulatory_regions"]:
            if r["kind"] == "enhancer":
                r["coactivator_species"] = "remodeler"
        for gene in c["genes"]:
            gene.update(on_rate=1.5, off_rate=.2, initiation_rate=1.5,
                        required_mark=None, elongation_rate=20.0,
                        pause_rate=.1, resume_rate=1.0, occupied_elongation_multiplier=.5,
                        nucleosome_eviction_rate=.3)
    elif identifier == "nucleosome_access":
        c = default_config()
        c["name"] = "Live nucleosome access · 25.6 kb"
        c["simulation"].update(confinement_radius=6.5, steps_per_frame=25)
        c["chromatin"].update(beads=128, bp_per_bead=200, representation="nucleosome_unit",
            bond_length=.65, initial_marks=[{"start": 0, "end": 25, "mark": "H3K9me3"},
                                           {"start": 72, "end": 103, "mark": "H3K27ac"}],
            barriers=[{"bead": 42, "direction": "right"}, {"bead": 110, "direction": "left"}],
            nucleosomes={"enabled": True, "initial_occupancy": 1,
                "initial_regions": [{"start": 10, "end": 13, "occupied": False},
                                    {"start": 36, "end": 39, "occupied": False},
                                    {"start": 72, "end": 74, "occupied": False}],
                "eviction_rate": .015, "reassembly_rate": .08,
                "occupied_persistence_multiplier": 1, "vacant_persistence_multiplier": .5})
        c["proteins"] = [
            _protein("activator", "Open-site TF", "#d6a34a", 40, valence=2, on_rate=3, off_rate=1,
                     chromatin_affinity={"unmodified": 1, "H3K27ac": 2, "H3K9me3": .2},
                     nucleosome_affinity={"occupied": .05, "vacant": 1},
                     nucleosome_off_rate_multiplier={"occupied": 2, "vacant": .5}),
            _protein("pioneer", "Pioneer-like TF", "#5daf9a", 16, valence=1, on_rate=2, off_rate=.8,
                     chromatin_affinity={"unmodified": .8, "H3K27ac": 1},
                     nucleosome_affinity={"occupied": .8, "vacant": 1}),
            _protein("reader", "Nucleosome reader", "#9b81c7", 12, valence=2, on_rate=2, off_rate=.6,
                     chromatin_affinity={"unmodified": .2, "H3K9me3": 3},
                     nucleosome_affinity={"occupied": 1, "vacant": 0},
                     nucleosome_off_rate_multiplier={"occupied": 1, "vacant": 4}),
            _protein("remodeler", "Sliding / evicting remodeler", "#649dcc", 8, valence=2, on_rate=2, off_rate=.5,
                     chromatin_affinity={"unmodified": .3},
                     nucleosome_affinity={"occupied": 1, "vacant": .5},
                     nucleosome_eviction_rate=.3, nucleosome_reassembly_rate=.05,
                     nucleosome_sliding_rate=1.5)]
        c["interactions"] = [{"a": "activator", "b": "remodeler", "strength": 2,
                              "on_rate": 2, "off_rate": .5, "radius": 1.1}]
        c["bodies"] = []
        c["cohesin"].update(count=2)
        c["genes"] = [dict(c["genes"][0], start=72, end=103, required_mark=None,
                           elongation_rate=10, occupied_elongation_multiplier=.5,
                           nucleosome_eviction_rate=.1)]
        c["regulatory_regions"] = [
            {"id": kind, "name": "Gene A " + kind, "kind": kind, "start": start, "end": end,
             "color": "#649dcc" if kind == "enhancer" else "#d6a34a", "target_gene": "gene_a",
             "tf_species": "activator", "coactivator_species": "remodeler",
             "tf_affinity_multiplier": 3, "activation_gain": .5, "initiation_gain": .25,
             "coactivator_gain": 1, "contact_radius": 1.5, "enabled": True}
            for kind, start, end in (("enhancer", 36, 39), ("promoter", 72, 74))]
        c["transcription"] = {"enabled": True, "max_transcripts": 2000,
                              "rna_diffusion": .3, "mature_rna_diffusion": .15}
    elif identifier in ("manuscript_locus", "manuscript_20mb"):
        c = preset_config("nucleosome_transcription" if identifier == "manuscript_locus" else "transcription_20mb")
        c["name"] = "Manuscript hypotheses · " + ("80 kb nucleosome units" if identifier == "manuscript_locus" else "20 Mb coarse landscape")
        c["simulation"]["confinement_stiffness"] = 30.0
        # These named modules are controllable hypotheses, not a conversion of
        # SIF associations into measured rate constants or causal evidence.
        for identifier, name, color in [
                ("H3K4me1", "H3K4me1 · enhancer-associated", "#7fb486"),
                ("H3K36me3", "H3K36me3 · transcription-associated", "#4caaa3"),
                ("H3K27me3", "H3K27me3 · Polycomb-associated", "#9b82c3"),
                ("H2AK119ub1", "H2AK119ub1 · Polycomb-associated", "#bd86ac")]:
            c["marks"].append(dict(id=identifier, name=name, color=color,
                                   decay_rate=.02, spread_rate=.02, contact_radius=1.25))
        for region in c["regulatory_regions"]:
            if region["kind"] == "enhancer":
                c["chromatin"]["initial_marks"].append(dict(start=region["start"], end=region["end"], mark="H3K4me1"))
        n = c["chromatin"]["beads"]
        # A separate terminal domain provides a controllable repressive module.
        for mark in ("H3K27me3", "H2AK119ub1"):
            c["chromatin"]["initial_marks"].append(dict(start=int(.91*n), end=int(.97*n)-1, mark=mark))
        for gene in c["genes"]:
            c["chromatin"]["initial_marks"].append(dict(start=gene["start"], end=gene["end"], mark="H3K36me3"))
            gene["body_responses"] = {
                "speckle": {"activation_multiplier": 1.5, "initiation_multiplier": 1.0},
                "nucleolus": {"activation_multiplier": .7, "initiation_multiplier": 1.0}}
        c["proteins"].extend([
            _protein("hdac", "HDAC-like eraser", "#d18d6b", 4, valence=2,
                     eraser_mark="H3K27ac", catalytic_rate=.3,
                     chromatin_affinity={"unmodified": .1, "H3K27ac": 1.5}),
            _protein("prc2", "PRC2-like writer", "#9b82c3", 4, valence=2,
                     writer_mark="H3K27me3", catalytic_rate=.3,
                     chromatin_affinity={"unmodified": .05, "H3K27me3": 2}),
            _protein("prc1", "PRC1-like writer", "#bd86ac", 4, valence=3,
                     writer_mark="H2AK119ub1", catalytic_rate=.3,
                     chromatin_affinity={"unmodified": .05, "H3K27me3": 2, "H2AK119ub1": 1})])
        c["chromatin"]["state_interactions"].append(
            dict(a="H3K27me3", b="H3K27me3", strength=.3, radius=1.5, min_separation=3))
    else:
        raise KeyError(identifier)
    return c


class ConfigError(ValueError):
    def __init__(self, errors):
        self.errors = errors
        super().__init__("; ".join(errors))


def validate_config(config):
    """Return a complete report; malformed input never reaches the integrator."""
    errors = []
    warnings = list(LIMITATIONS)
    if not isinstance(config, dict):
        return {"valid": False, "errors": ["Configuration must be a JSON object."],
                "warnings": [], "particle_count": 0, "genomic_span_bp": None}

    if "calibration" in config:
        from .calibration import validate_profile, config_fingerprint
        try:
            profile = validate_profile(config["calibration"])
            if profile["status"] != "fitted":
                warnings.append("Physical mapping status: " + profile["status"] + ". This is not an accepted experimental clock.")
            if profile.get("baseline_fingerprint") != config_fingerprint(config):
                warnings.append("This model differs from the calibration baseline. The original clock is retained for comparison; review the changes before interpreting physical values.")
        except (ValueError, TypeError, OverflowError) as exc:
            errors.append("calibration: " + str(exc))
    else:
        warnings.append("No physical length/time mapping is attached. Times and diffusion coefficients cannot yet be read as biological seconds or micrometers squared per second.")

    def obj(value, path):
        if not isinstance(value, dict):
            errors.append(path + " must be an object.")
            return {}
        return value

    def arr(value, path, limit):
        if not isinstance(value, list):
            errors.append(path + " must be a list.")
            return []
        if len(value) > limit:
            errors.append(f"{path} exceeds the prototype limit of {limit} entries.")
        return value[:limit]

    def number(value, path, low=0.0, high=100.0, integer=False):
        try:
            finite = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
        except OverflowError:
            finite = False
        if not finite:
            errors.append(path + " must be a finite number.")
            return low
        if not low <= value <= high or (integer and value != int(value)):
            errors.append(f"{path} must be {'an integer ' if integer else ''}between {low} and {high}.")
        return min(high, max(low, value))

    def text(value, path, length=120):
        if not isinstance(value, str) or not value.strip() or len(value) > length:
            errors.append(f"{path} must be nonempty text of at most {length} characters.")

    def identity(items, path):
        ids = set()
        for k, item in enumerate(items):
            item = obj(item, f"{path}[{k}]")
            key = item.get("id")
            if not isinstance(key, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,47}", key):
                errors.append(f"{path}[{k}].id must be a short identifier beginning with a letter.")
                continue
            if key in ids:
                errors.append(f"Duplicate {path} identifier: {key}.")
            ids.add(key)
            text(item.get("name"), f"{path}[{k}].name")
            if path != "genes" and not re.fullmatch(r"#[0-9a-fA-F]{6}", str(item.get("color", ""))):
                errors.append(f"{path}[{k}].color must be a six-digit hexadecimal color.")
        return ids

    if config.get("version") != "0.1":
        errors.append("Unsupported version; expected 0.1.")
    text(config.get("name"), "name")
    number(config.get("seed"), "seed", 0, 2**32 - 1, True)
    sim = obj(config.get("simulation"), "simulation")
    number(sim.get("dt"), "simulation.dt", 0.00001, 0.01)
    number(sim.get("steps_per_frame"), "simulation.steps_per_frame", 1, 100, True)
    boundary = number(sim.get("confinement_radius"), "simulation.confinement_radius", 3, 100)
    wall_stiffness = number(sim.get("confinement_stiffness", 30), "simulation.confinement_stiffness", 0, 200)
    if wall_stiffness == 0:
        warnings.append("Confinement force is disabled. The radius remains an initialization/reference scale; particles can diffuse away from it.")
    chrom = obj(config.get("chromatin"), "chromatin")
    if chrom.get("representation", "coarse") not in ("coarse", "nucleosome_unit"):
        errors.append("chromatin.representation must be coarse or nucleosome_unit.")
    n = number(chrom.get("beads"), "chromatin.beads", 8, MAX_CHROMATIN_BEADS, True)
    bp_per_bead = number(chrom.get("bp_per_bead"), "chromatin.bp_per_bead", 1, 1000000, True)
    raw_dimensions = (chrom.get("beads"), chrom.get("bp_per_bead"))
    genomic_span = None
    try:
        if all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
               and v >= 1 and v == int(v) for v in raw_dimensions):
            genomic_span = int(raw_dimensions[0]) * int(raw_dimensions[1])
    except (ValueError, OverflowError):
        pass
    if bp_per_bead > 1000:
        warnings.append(
            f"Coarse genomic resolution: each bead represents {bp_per_bead:,.0f} bp. "
            "Individual nucleosomes, promoters, and other sub-bead features are unresolved. "
            "Changing bp per bead does not calibrate or rescale diffusion, mechanics, or reaction rates."
        )
    if chrom.get("representation") == "nucleosome_unit":
        warnings.append("Nucleosome-unit representation labels one isotropic bead per chosen repeat. Wrapped DNA, core orientation, linker DNA and histone tails are not resolved; mechanics and rates require calibration.")
        if not 147 <= bp_per_bead <= 300:
            warnings.append("The chosen bp per nucleosome unit lies outside the illustrative 147–300 bp repeat range; check your intended mapping.")
    for field, low, high in [("bond_length", 0.2, 2), ("bond_k", 1, 200),
                             ("diffusion", 0, 2), ("persistence_length", 0, 20), ("repulsion", 0, 100)]:
        number(chrom.get(field), "chromatin." + field, low, high)
    marks = arr(config.get("marks"), "marks", 16)
    species = arr(config.get("proteins"), "proteins", 32)
    genes = arr(config.get("genes"), "genes", MAX_GENES)
    bodies = arr(config.get("bodies"), "bodies", 12)
    mark_ids, species_ids = identity(marks, "marks"), identity(species, "proteins")
    if "unmodified" in mark_ids:
        errors.append("marks: unmodified is reserved for the absence of mark flags.")
    gene_ids = identity(genes, "genes")
    body_ids = identity(bodies, "bodies")

    def reference(value, choices, path, optional=False):
        if optional and value in (None, ""):
            return
        if not isinstance(value, str) or value not in choices:
            errors.append(f"{path} references an unknown identifier: {value!r}.")

    def affinities(value, choices, path, signed=False):
        for key, val in obj(value, path).items():
            reference(key, choices, path)
            number(val, path + "." + key, -20 if signed else 0, 20)

    def region(value, path):
        start = number(value.get("start"), path + ".start", 0, max(n-1, 0), True)
        end = number(value.get("end"), path + ".end", 0, max(n-1, 0), True)
        if start > end:
            errors.append(path + ": start must not exceed end.")

    if "nucleosomes" in chrom:
        nuc = obj(chrom["nucleosomes"], "chromatin.nucleosomes")
        if not isinstance(nuc.get("enabled", False), bool):
            errors.append("chromatin.nucleosomes.enabled must be boolean.")
        number(nuc.get("initial_occupancy", 1), "chromatin.nucleosomes.initial_occupancy", 0, 1)
        for field in ("eviction_rate", "reassembly_rate"):
            number(nuc.get(field, 0), "chromatin.nucleosomes." + field, 0, 20)
        for field, default in (("occupied_persistence_multiplier", 1), ("vacant_persistence_multiplier", .5)):
            number(nuc.get(field, default), "chromatin.nucleosomes." + field, 0, 5)
        for k, r in enumerate(arr(nuc.get("initial_regions", []), "chromatin.nucleosomes.initial_regions", 512)):
            path = f"chromatin.nucleosomes.initial_regions[{k}]"
            r = obj(r, path)
            region(r, path)
            if not isinstance(r.get("occupied"), bool):
                errors.append(path + ".occupied must be boolean.")
        if nuc.get("enabled"):
            warnings.append("Nucleosome eviction clears histone marks; reassembly starts unmarked. Occupancy changes affect new binding hazards, while existing protein tethers retain their off-rates.")
            if bp_per_bead > 300:
                warnings.append("Occupancy is an effective binary state of a coarse genomic bin at this resolution, not the occupancy of one physical nucleosome.")
    if "transcription" in config:
        tr = obj(config["transcription"], "transcription")
        if not isinstance(tr.get("enabled", False), bool):
            errors.append("transcription.enabled must be boolean.")
        number(tr.get("max_transcripts", 2000), "transcription.max_transcripts", 1, MAX_TRANSCRIPTS, True)
        for field, default in (("rna_diffusion", .3), ("mature_rna_diffusion", .15)):
            number(tr.get(field, default), "transcription." + field, 0, 5)
        if tr.get("enabled"):
            warnings.append("Individual RNA molecules are passive Brownian points. Nascent transcripts follow polymerases; released molecules diffuse without RNA–protein binding, steric interactions or RNA-chain forces. The live-transcript cap blocks new initiation and must remain inactive for quantitative interpretation.")

    for k, m in enumerate(marks):
        m = obj(m, f"marks[{k}]")
        for field in ("decay_rate", "spread_rate"):
            number(m.get(field), f"marks[{k}].{field}", 0, 10)
        number(m.get("contact_radius"), f"marks[{k}].contact_radius", 0.1, 5)
    for k, r in enumerate(arr(chrom.get("initial_marks"), "chromatin.initial_marks", 512)):
        r = obj(r, f"initial_marks[{k}]")
        region(r, f"initial_marks[{k}]")
        reference(r.get("mark"), mark_ids, f"initial_marks[{k}].mark")
    for k, b in enumerate(arr(chrom.get("barriers"), "chromatin.barriers", 512)):
        b = obj(b, f"barriers[{k}]")
        number(b.get("bead"), f"barriers[{k}].bead", 0, n-1, True)
        if b.get("direction") not in ("left", "right", "both"):
            errors.append(f"barriers[{k}].direction must be left, right or both.")
    state_pairs = set()
    for k, rule in enumerate(arr(chrom.get("state_interactions", []), "chromatin.state_interactions", 153)):
        path = f"chromatin.state_interactions[{k}]"
        rule = obj(rule, path)
        for side in ("a", "b"):
            reference(rule.get(side), mark_ids | {"unmodified"}, path + "." + side)
        pair = tuple(sorted([str(rule.get("a")), str(rule.get("b"))]))
        if pair in state_pairs:
            errors.append(path + " duplicates an existing symmetric chromatin-state pair.")
        state_pairs.add(pair)
        number(rule.get("strength"), path + ".strength", 0, 20)
        number(rule.get("radius"), path + ".radius", .2, 5)
        number(rule.get("min_separation", 3), path + ".min_separation", 1, MAX_CHROMATIN_BEADS, True)
    particle_count = n
    for k, p in enumerate(species):
        p = obj(p, f"proteins[{k}]")
        path = f"proteins[{k}]"
        count = number(p.get("count"), path + ".count", 0, 500, True)
        beads = number(p.get("idr_beads"), path + ".idr_beads", 2, 16, True)
        number(p.get("valence"), path + ".valence", 0, 16, True)
        number(p.get("diffusion"), path + ".diffusion", 0, 5)
        number(p.get("radius"), path + ".radius", 0.1, 1)
        number(p.get("sticker_fraction"), path + ".sticker_fraction", 0, 1)
        for field in ("on_rate", "off_rate", "catalytic_rate"):
            number(p.get(field), path + "." + field, 0, 20)
        if p.get("representation") not in ("bead", "idr"):
            errors.append(path + ".representation must be bead or idr.")
        particle_count += count * (beads if p.get("representation") == "idr" else 1)
        affinities(p.get("chromatin_affinity"), mark_ids | {"unmodified"}, path + ".chromatin_affinity")
        for field in ("nucleosome_affinity", "nucleosome_off_rate_multiplier"):
            if field in p:
                affinity = obj(p[field], path + "." + field)
                for key, value in affinity.items():
                    if key not in ("occupied", "vacant"):
                        errors.append(path + "." + field + " has an unknown occupancy state.")
                    # Release factors span the measured nucleosome effects on TF
                    # dissociation (up to ~1000-fold faster; Luo et al. 2014).
                    number(value, path + "." + field + "." + str(key), 0,
                           1000 if field == "nucleosome_off_rate_multiplier" else 20)
        for field in ("nucleosome_eviction_rate", "nucleosome_reassembly_rate", "nucleosome_sliding_rate"):
            number(p.get(field, 0), path + "." + field, 0, 20)
        for field in ("writer_mark", "eraser_mark"):
            reference(p.get(field), mark_ids, path + "." + field, True)
    if particle_count > MAX_PARTICLES:
        errors.append(f"Model has {int(particle_count)} particles; this local prototype supports at most {MAX_PARTICLES}. Reduce counts or IDR chain lengths.")
    pair_keys = set()
    for k, r in enumerate(arr(config.get("interactions"), "interactions", 528)):
        r = obj(r, f"interactions[{k}]")
        path = f"interactions[{k}]"
        for side in ("a", "b"):
            reference(r.get(side), species_ids, path + "." + side)
        pair = tuple(sorted([str(r.get("a")), str(r.get("b"))]))
        if pair in pair_keys:
            errors.append(path + " duplicates an existing symmetric species pair.")
        pair_keys.add(pair)
        for field, low, high in [("strength", 0, 20), ("on_rate", 0, 20), ("off_rate", 0, 20), ("radius", 0.2, 5)]:
            number(r.get(field), path + "." + field, low, high)
    for k, b in enumerate(bodies):
        b = obj(b, f"bodies[{k}]")
        path = f"bodies[{k}]"
        center = arr(b.get("center"), path + ".center", 3)
        if len(center) != 3:
            errors.append(path + ".center must contain three coordinates.")
        for coord in center:
            number(coord, path + ".center", -100, 100)
        radius = number(b.get("radius"), path + ".radius", 0.2, 20)
        if len(center) == 3 and all(isinstance(v, (int, float)) and -100 <= v <= 100 for v in center):
            if math.sqrt(sum(v*v for v in center)) + radius > boundary:
                warnings.append(b.get("name", path) + " extends beyond the confinement boundary.")
        number(b.get("strength"), path + ".strength", 0, 20)
        number(b.get("chromatin_affinity"), path + ".chromatin_affinity", -20, 20)
        affinities(b.get("protein_affinity"), species_ids, path + ".protein_affinity", True)
        if b.get("kind") not in ("speckle", "nucleolus"):
            errors.append(path + ".kind must be speckle or nucleolus.")
        if not isinstance(b.get("enabled"), bool):
            errors.append(path + ".enabled must be boolean.")
    coh = obj(config.get("cohesin"), "cohesin")
    if not isinstance(coh.get("enabled"), bool):
        errors.append("cohesin.enabled must be boolean.")
    number(coh.get("count"), "cohesin.count", 0, min(MAX_COHESIN, int(n)//2), True)
    for field in ("step_rate", "off_rate", "spring_k"):
        number(coh.get(field), "cohesin." + field, 0, 20)
    for k, g in enumerate(genes):
        g = obj(g, f"genes[{k}]")
        path = f"genes[{k}]"
        region(g, path)
        for field in ("on_rate", "off_rate", "initiation_rate", "processing_rate", "rna_decay_rate"):
            number(g.get(field), path + "." + field, 0, 20)
        number(g.get("elongation_rate"), path + ".elongation_rate", 0.01, 100)
        number(g.get("speckle_gain"), path + ".speckle_gain", 0, 20)
        number(g.get("persistence_multiplier"), path + ".persistence_multiplier", 0.1, 5)
        for field, default in (("pause_rate", 0), ("resume_rate", 1), ("nucleosome_eviction_rate", 0)):
            number(g.get(field, default), path + "." + field, 0, 20)
        number(g.get("occupied_elongation_multiplier", 1), path + ".occupied_elongation_multiplier", 0, 5)
        reference(g.get("activator_species"), species_ids, path + ".activator_species", True)
        reference(g.get("required_mark"), mark_ids, path + ".required_mark", True)
        reference(g.get("repressive_mark"), mark_ids, path + ".repressive_mark", True)
        for body_id, response in obj(g.get("body_responses", {}), path + ".body_responses").items():
            response_path = path + ".body_responses." + str(body_id)
            reference(body_id, body_ids, response_path)
            response = obj(response, response_path)
            if set(response) - {"activation_multiplier", "initiation_multiplier"}:
                errors.append(response_path + " accepts only activation_multiplier and initiation_multiplier.")
            for field in ("activation_multiplier", "initiation_multiplier"):
                number(response.get(field, 1), response_path + "." + field, .01, 20)
        if g.get("body_responses"):
            warnings.append(path + ": nuclear-body responses are prescribed proximity-dependent hypotheses, not measured SIF-to-rate mappings.")
    if "nucleus" in config:
        nucleus = obj(config["nucleus"], "nucleus")
        if not isinstance(nucleus.get("enabled", False), bool):
            errors.append("nucleus.enabled must be boolean.")
        radii = arr(nucleus.get("radii"), "nucleus.radii", 3)
        center = arr(nucleus.get("center", [0, 0, 0]), "nucleus.center", 3)
        if len(radii) != 3 or len(center) != 3:
            errors.append("nucleus.radii and nucleus.center need three values.")
        semi = [number(v, f"nucleus.radii[{k}]", 1, 2000) for k, v in enumerate(radii)]
        mid = [number(v, f"nucleus.center[{k}]", -2000, 2000) for k, v in enumerate(center)]
        number(nucleus.get("envelope_stiffness", 30), "nucleus.envelope_stiffness", 0, 200)
        if nucleus.get("enabled") and len(semi) == 3 and len(mid) == 3:
            scaled = math.sqrt(sum((m / a) ** 2 for m, a in zip(mid, semi)))
            if scaled > 0.9:
                errors.append("nucleus: the model origin (chromatin seed point) must lie inside the envelope, at most 0.9 of the way to the wall.")
            warnings.append("Nuclear envelope is a soft ellipsoidal wall; tethers are phenomenological anchor springs. No lamina proteins, pores or membrane mechanics are represented.")
        tethers = arr(nucleus.get("tethers", []), "nucleus.tethers", 128)
        identity(tethers, "nucleus.tethers")
        for k, t in enumerate(tethers):
            path = f"nucleus.tethers[{k}]"
            t = obj(t, path)
            region(t, path)
            if not isinstance(t.get("enabled", True), bool):
                errors.append(path + ".enabled must be boolean.")
            if t.get("anchor", "envelope") not in ("envelope", "body"):
                errors.append(path + ".anchor must be envelope or body.")
            if t.get("anchor") == "body":
                reference(t.get("body_id"), body_ids, path + ".body_id")
            for field, low, high, default in (("attach_rate", 0, 20, None), ("detach_rate", 0, 20, None),
                                              ("capture_distance", 0.1, 10, 1.5), ("spring_k", 0, 50, 2.0),
                                              ("rest_length", 0, 5, 0.5)):
                number(t.get(field, default), path + "." + field, low, high)
            reference(t.get("required_mark"), mark_ids, path + ".required_mark", True)
            for field in ("occupancy_attach_multiplier", "occupancy_detach_multiplier"):
                if field in t:
                    factors = obj(t[field], path + "." + field)
                    for key, value in factors.items():
                        if key not in ("occupied", "vacant"):
                            errors.append(path + "." + field + " has an unknown occupancy state.")
                        number(value, path + "." + field + "." + str(key), 0, 1000)
    regions = arr(config.get("regulatory_regions", []), "regulatory_regions", MAX_REGULATORY_REGIONS)
    identity(regions, "regulatory_regions")
    gene_lookup = {g.get("id"): g for g in genes if isinstance(g, dict) and isinstance(g.get("id"), str)}
    promoters = set()
    for k, r in enumerate(regions):
        path = f"regulatory_regions[{k}]"
        r = obj(r, path)
        region(r, path)
        if r.get("kind") not in ("enhancer", "promoter"):
            errors.append(path + ".kind must be enhancer or promoter.")
        reference(r.get("target_gene"), gene_ids, path + ".target_gene")
        for field in ("tf_species", "coactivator_species"):
            reference(r.get(field), species_ids, path + "." + field, True)
        for field in ("tf_affinity_multiplier", "activation_gain", "initiation_gain", "coactivator_gain"):
            number(r.get(field), path + "." + field, 0, 20)
        number(r.get("contact_radius"), path + ".contact_radius", .2, 5)
        if not isinstance(r.get("enabled"), bool):
            errors.append(path + ".enabled must be boolean.")
        target = r.get("target_gene")
        if r.get("kind") == "promoter" and r.get("enabled") and isinstance(target, str):
            if target in promoters:
                errors.append(path + ": only one enabled promoter region is allowed per target gene.")
            promoters.add(target)
            promoter_start = gene_lookup.get(target, {}).get("start")
            values = (r.get("start"), r.get("end"), promoter_start)
            if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values):
                if not values[0] <= values[2] <= values[1]:
                    errors.append(path + ": promoter region must contain the target gene start bead.")
        if r.get("enabled") and r.get("coactivator_species") and r.get("tf_species"):
            pair = tuple(sorted([str(r["tf_species"]), str(r["coactivator_species"])]))
            if pair not in pair_keys:
                warnings.append(path + ": recruitment requires a TF–coactivator protein interaction rule, which is currently absent.")
    return {"valid": not errors, "errors": errors, "warnings": warnings,
            "particle_count": int(particle_count), "genomic_span_bp": genomic_span}


def checked_config(config):
    report = validate_config(config)
    if not report["valid"]:
        raise ConfigError(report["errors"])
    return copy.deepcopy(config)
