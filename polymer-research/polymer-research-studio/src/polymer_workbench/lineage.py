"""Provenance for physical specimens and evidence, independent of model run history."""
from copy import deepcopy

FORMS = ("batch", "cast_film", "electrospun_mat", "solution", "other")
ORIGINS = ("own_measurement", "primary_article", "secondary_citation", "model_assumption", "illustrative")


def required_text(record, key, limit=1000):
    value = record.get(key)
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ValueError(f"{key} must contain 1–{limit} characters.")
    return value


def check_lineage(material):
    if "lineage" in material:
        node = material["lineage"]
        if not isinstance(node, dict) or node.get("form") not in FORMS:
            raise ValueError("Choose a supported specimen form.")
        required_text(node, "process")
        ancestors = node.get("ancestors")
        if not isinstance(ancestors, list) or len(ancestors) > 50:
            raise ValueError("A specimen lineage supports up to 50 ancestor identifiers.")
        if any(not isinstance(v, str) or not v.strip() or len(v) > 200 for v in ancestors):
            raise ValueError("Ancestor identifiers must be nonempty text.")
        if len(set([*ancestors, material["sample_id"]])) != len(ancestors) + 1:
            raise ValueError("A specimen cannot be its own ancestor.")
    entries = material.get("evidence", [])
    if not isinstance(entries, list) or len(entries) > 200:
        raise ValueError("A specimen supports at most 200 evidence records.")
    for e in entries:
        if not isinstance(e, dict) or e.get("origin") not in ORIGINS:
            raise ValueError("Select an evidence origin.")
        for key in ("quantity", "value", "unit", "conditions", "locator"):
            required_text(e, key)
        if e["origin"] in ("primary_article", "secondary_citation"):
            required_text(e, "reported_in")
        if e["origin"] == "secondary_citation":
            required_text(e, "attributed_to")  # 'unresolved' is explicit and remains secondary.


def derive_material(parent, sample_id, name, form, process):
    from .exchange import default_material, validate_material
    parent = validate_material(parent)
    child = default_material()
    child.update(sample_id=sample_id, name=name)
    child["chemistry"] = deepcopy(parent["chemistry"])
    child["conditions"] = {"temperature_K": None, "state": "unspecified", "solvent": "unspecified",
                           "notes": "Record preparation and measurement conditions for this specimen."}
    child["lineage"] = {"form": form, "process": process,
                        "ancestors": [*parent.get("lineage", {}).get("ancestors", []), parent["sample_id"]]}
    child["evidence"] = []
    return validate_material(child)


def check_registry(samples, active):
    from .exchange import validate_material
    if not isinstance(samples, list) or not 1 <= len(samples) <= 100:
        raise ValueError("A workspace supports 1–100 registered specimens.")
    by_id = {}
    for entry in samples:
        item = validate_material(entry)
        if item["sample_id"] in by_id:
            raise ValueError("Specimen identifiers must be unique within the workspace.")
        by_id[item["sample_id"]] = item
    if by_id.get(active["sample_id"]) != active:
        raise ValueError("The active specimen must match its registered material record.")
    for item in by_id.values():
        ancestors = item.get("lineage", {}).get("ancestors", [])
        for index, ancestor_id in enumerate(ancestors):
            ancestor = by_id.get(ancestor_id)
            # A single imported specimen may refer to an external ancestor.
            if ancestor and ancestor.get("lineage", {}).get("ancestors", []) != ancestors[:index]:
                raise ValueError("Conflicting specimen ancestry in this workspace.")
