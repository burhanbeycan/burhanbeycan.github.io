"""Record which simulated population and experimental specimen a result refers to."""
from copy import deepcopy

FIELDS = {
    "population": "Chain population (DP, distribution, sequence)",
    "termini": "End groups and attachment sites",
    "protonation": "Ionization state and counterions",
    "force_field": "Force field, charge method and parameter source",
    "conditions": "Geometry, solvent and thermodynamic conditions",
    "replicas": "Number of independently initialized systems",
    "target": "Observable and experimental specimen / source",
}
METHODS = ("calculation", "planned_md", "executed_md", "trajectory_analysis", "published_md")


def validate_study(value):
    if value is None:
        return None
    if not isinstance(value, dict) or value.get("method") not in METHODS:
        raise ValueError("Select a valid study method.")
    if set(value) - {"method", *FIELDS}:
        raise ValueError("Unknown study record fields.")
    for key, v in value.items():
        if key == "replicas":
            if v is not None and (type(v) is not int or not 1 <= v <= 100000):
                raise ValueError("Independent systems must be a positive integer or omitted.")
        elif not isinstance(v, str) or not v.strip() or len(v) > 1000:
            raise ValueError("Study fields require nonempty text up to 1,000 characters.")
    return deepcopy(value)


def study_audit(study):
    if not study or study["method"]=="calculation":
        return []
    missing=[f for f in FIELDS if not study.get(f)]
    notes=["Record "+FIELDS[f].lower()+"." for f in missing]
    if study.get("replicas")==1:
        notes.append("One initial system was declared. Between-initialization variability is unavailable.")
    if study["method"]=="planned_md":
        notes.append("This linked model is planned. Simulation evidence requires a completed trajectory.")
    return notes


def study_snapshot(material):
    context=validate_study(material.get("study_context"))
    return {"specimen": material["sample_id"], "declaration": context, "audit":study_audit(context),
            "basis":"User-entered context for this specimen. Calculation inputs and the producing laboratory identify the actual result."}
