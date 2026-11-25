"""Component studies: one mechanism at a time, compared with a reference.

Each study is one self-contained module in this package. Adding, changing or
removing a study touches only its own file and its line in ``REGISTRY``.
See docs/COMPONENT_STUDIES.md for the module contract.
"""

import importlib

# study id -> module name inside this package. Order is the suggested run order:
# cheap single-mechanism checks first, coupled experimental comparisons last.
REGISTRY = {
    "polymer_wlc": "polymer_wlc",
    "protein_diffusion": "protein_diffusion",
    "nucleosome_turnover": "nucleosome_turnover",
    "telegraph_transcription": "telegraph_transcription",
    "cohesin_extrusion": "cohesin_extrusion",
    "tf_binding": "tf_binding",
    "nucleosome_tf_kinetics": "nucleosome_tf_kinetics",
    "hp1_mark_density": "hp1_mark_density",
    "fbn2_cohesin": "fbn2_cohesin",
}

REQUIRED_STUDY_KEYS = ("id", "title", "components", "question", "evidence", "sources", "scales")


def load(study_id):
    """Import one study module and verify the module contract."""
    if study_id not in REGISTRY:
        raise ValueError(f"Unknown component study: {study_id}. Known: {', '.join(REGISTRY)}")
    module = importlib.import_module(f"{__name__}.{REGISTRY[study_id]}")
    study = getattr(module, "STUDY", None)
    if not isinstance(study, dict) or any(key not in study for key in REQUIRED_STUDY_KEYS):
        raise ValueError(f"{study_id}: STUDY must define {', '.join(REQUIRED_STUDY_KEYS)}.")
    if study["id"] != study_id:
        raise ValueError(f"{study_id}: STUDY['id'] does not match its registry key.")
    for name in ("simulate", "analyze"):
        if not callable(getattr(module, name, None)):
            raise ValueError(f"{study_id}: module must define {name}().")
    if "smoke" not in study["scales"] or "standard" not in study["scales"]:
        raise ValueError(f"{study_id}: scales must include 'smoke' and 'standard'.")
    return module


def catalog():
    """Read-only description of every registered study."""
    records = []
    for study_id in REGISTRY:
        study = load(study_id).STUDY
        records.append({key: study[key] for key in ("id", "title", "components", "question", "evidence")}
                       | {"scales": sorted(study["scales"]), "sources": study["sources"]})
    return records
