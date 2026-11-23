"""Versioned, unit-explicit material and calculation exchange.

Checksums detect accidental edits; they are not digital signatures or proof of
experimental validity. All imported measurements retain their stated origin.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
from typing import Any
import uuid

from .lineage import check_lineage, check_registry
from .study import validate_study, study_snapshot

MATERIAL_SCHEMA = "polymer-material/1.0"
RESULT_SCHEMA = "polymer-result/1.0"
SESSION_SCHEMA = "polymer-session/1.0"
MAX_JSON_BYTES = 4_000_000


def ensure_json(value: Any) -> Any:
    """Reject non-finite numbers, excessive payloads and non-JSON objects."""
    def check(v: Any, depth: int = 0) -> None:
        if depth > 35:
            raise ValueError("The exchange document is nested too deeply.")
        if isinstance(v, float) and not math.isfinite(v):
            raise ValueError("NaN and infinite values are not valid scientific inputs.")
        if isinstance(v, dict):
            if any(not isinstance(k, str) for k in v):
                raise ValueError("JSON object keys must be strings.")
            for item in v.values():
                check(item, depth + 1)
        elif isinstance(v, (list, tuple)):
            for item in v:
                check(item, depth + 1)
        elif v is not None and not isinstance(v, (str, bool, int, float)):
            raise ValueError("Only ordinary JSON values can be exchanged.")
    check(value)
    if len(json.dumps(value, allow_nan=False).encode("utf-8")) > MAX_JSON_BYTES:
        raise ValueError("Exchange document exceeds the 4 MB size limit.")
    return value


def default_material() -> dict:
    return {
        "schema": MATERIAL_SCHEMA,
        "sample_id": "demonstration-01",
        "name": "Model polymer · demonstration",
        "chemistry": {
            "repeat_unit": "unspecified",
            "architecture": "linear",
            "tacticity": "unspecified",
            "sequence": "unspecified",
            "end_groups": "unspecified",
            "repeat_smiles": "",
            "stereochemistry_source": "unspecified",
            "branch_topology": "unspecified",
            "functional_groups": "unspecified",
            "synthesis_route": "unspecified",
            "conversion_basis": "unspecified",
            "notes": "Illustrative inputs. No material-specific calibration is implied.",
        },
        "conditions": {
            "temperature_K": 298.15,
            "state": "unspecified",
            "solvent": "unspecified",
            "notes": "Each calculation records its own model conditions.",
        },
        "parameters": [],
        "history": [],
    }


def validate_material(material: dict | None) -> dict:
    m = deepcopy(default_material() if material is None else material)
    ensure_json(m)
    if not isinstance(m, dict) or m.get("schema") != MATERIAL_SCHEMA:
        raise ValueError("Expected a polymer-material/1.0 record.")
    for key in ("sample_id", "name"):
        if not isinstance(m.get(key), str) or not 1 <= len(m[key]) <= 200:
            raise ValueError(f"Material {key} must be text with 1–200 characters.")
    for key in ("chemistry", "conditions"):
        if not isinstance(m.get(key), dict):
            raise ValueError(f"Material {key} must be an object.")
    for key, value in m["chemistry"].items():
        if not isinstance(value, str) or len(value) > 4000:
            raise ValueError("Chemistry metadata must be text, at most 4000 characters per field.")
    for key in ("state", "solvent", "notes"):
        if key in m["conditions"] and not isinstance(m["conditions"][key], str):
            raise ValueError("State, solvent and condition notes must be text.")
    temp = m["conditions"].get("temperature_K")
    if temp is not None and (isinstance(temp, bool) or not isinstance(temp, (float, int)) or temp <= 0):
        raise ValueError("Material temperature must be a positive value in kelvin.")
    for key in ("parameters", "history"):
        if not isinstance(m.get(key, []), list):
            raise ValueError(f"Material {key} must be a list.")
        m.setdefault(key, [])
    if any(not isinstance(p, dict) for p in m["parameters"]):
        raise ValueError("Material parameters must be objects with explicit quantity and unit metadata.")
    for h in m["history"]:
        if not isinstance(h, dict) or any(not isinstance(h.get(k), str) for k in ("run_id", "application", "version")):
            raise ValueError("Calculation history requires run_id, application and version strings.")
    validate_study(m.get("study_context"))
    check_lineage(m)
    distribution = m.get("molar_mass_distribution")
    if distribution is not None:
        if not isinstance(distribution, dict):
            raise ValueError("Molar-mass distribution must be an object.")
        masses = distribution.get("molar_mass_g_mol", [])
        weights = distribution.get("number_fraction", [])
        if not isinstance(masses, list) or not isinstance(weights, list) or len(masses) != len(weights) or not 1 <= len(masses) <= 10000:
            raise ValueError("Molar masses and number fractions need matching, nonempty lists.")
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or x <= 0 for x in masses):
            raise ValueError("Molar masses must be positive, in g/mol.")
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or x < 0 for x in weights) or not math.isclose(sum(weights), 1.0, rel_tol=1e-7, abs_tol=1e-9):
            raise ValueError("Number fractions must be nonnegative and sum to one.")
        if distribution.get("basis", "number_fraction") != "number_fraction":
            raise ValueError("Exchanged distributions use explicitly normalized number fractions.")
    return m


def digest(document: dict) -> str:
    clean = {k: v for k, v in document.items() if k != "checksum_sha256"}
    return hashlib.sha256(json.dumps(clean, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def validate_output(output: dict) -> None:
    """Check the public result contract before rendering or accepting a file."""
    if not isinstance(output, dict):
        raise ValueError("Calculation results must be an object.")
    for key in ("metrics", "plots", "tables", "notes", "warnings", "handoffs"):
        if not isinstance(output.get(key, []), list):
            raise ValueError(f"Calculation returned an invalid {key} field.")
    for key in ("notes", "warnings"):
        if any(not isinstance(value, str) for value in output.get(key, [])):
            raise ValueError(f"Calculation {key} must contain text.")
    views = output.get("molecular_views", [])
    if not isinstance(views, list) or len(views) > 8:
        raise ValueError("A calculation may contain at most eight molecular views.")
    for view in views:
        if not isinstance(view, dict) or not isinstance(view.get("title"), str) or not isinstance(view.get("caption", ""), str):
            raise ValueError("Molecular views need a title and text caption.")
        atoms, bonds = view.get("atoms"), view.get("bonds")
        if not isinstance(atoms, list) or not 1 <= len(atoms) <= 2500 or not isinstance(bonds, list) or len(bonds) > 4000:
            raise ValueError("A molecular view needs 1–2500 atoms and at most 4000 bonds.")
        ids = set()
        for atom in atoms:
            if not isinstance(atom, dict) or isinstance(atom.get("index"), bool) or not isinstance(atom.get("index"), int) or atom["index"] in ids:
                raise ValueError("Molecular atom indices must be unique integers.")
            ids.add(atom["index"])
            if not isinstance(atom.get("element"), str) or not 1 <= len(atom["element"]) <= 3:
                raise ValueError("Each molecular atom needs a chemical element label.")
            if any(isinstance(atom.get(k), bool) or not isinstance(atom.get(k), (int, float)) or not math.isfinite(atom[k]) or abs(atom[k]) > 1e6 for k in ("x", "y")):
                raise ValueError("Molecular drawing coordinates must be bounded finite numbers.")
            if not isinstance(atom.get("label", ""), str) or len(atom.get("label", "")) > 40:
                raise ValueError("Molecular atom labels must be text of at most 40 characters.")
            if "unit_index" in atom and atom["unit_index"] is not None and (isinstance(atom["unit_index"], bool) or not isinstance(atom["unit_index"], int)):
                raise ValueError("Repeat-unit memberships must be integer indices or null.")
        edges = set()
        for bond in bonds:
            if not isinstance(bond, dict) or any(isinstance(bond.get(k), bool) or not isinstance(bond.get(k), int) or bond[k] not in ids for k in ("begin", "end")) or bond["begin"] == bond["end"]:
                raise ValueError("Every molecular bond must connect two existing distinct atom indices.")
            edge = tuple(sorted((bond["begin"], bond["end"])))
            if edge in edges:
                raise ValueError("Duplicate molecular graph bonds are invalid.")
            edges.add(edge)
            if isinstance(bond.get("order"), bool) or bond.get("order") not in (1, 1.5, 2, 3):
                raise ValueError("Molecular bond orders must be 1, 1.5, 2 or 3.")
            if not isinstance(bond.get("direction", "NONE"), str) or not isinstance(bond.get("stereo", ""), str):
                raise ValueError("Molecular bond stereochemistry must be text.")
    for metric in output.get("metrics", []):
        if not isinstance(metric, dict) or not isinstance(metric.get("label"), str):
            raise ValueError("Each metric needs a text label.")
        value = metric.get("value")
        if value is not None and (isinstance(value, bool) or not isinstance(value, (str, int, float))):
            raise ValueError("A metric must be a number, text, or null.")
    for plot in output.get("plots", []):
        if not isinstance(plot, dict) or any(not isinstance(plot.get(k), str) for k in ("title", "x_label", "y_label")):
            raise ValueError("Each plot needs a title and labelled axes.")
        if any(plot.get(k, "linear") not in ("linear", "log") for k in ("x_scale", "y_scale")):
            raise ValueError("Plot scales must be linear or log.")
        if not isinstance(plot.get("series"), list):
            raise ValueError("Each plot needs a series list.")
        for series in plot["series"]:
            if not isinstance(series, dict) or not isinstance(series.get("name"), str):
                raise ValueError("Each curve needs a name.")
            x, y = series.get("x"), series.get("y")
            if not isinstance(x, list) or not isinstance(y, list) or len(x) != len(y) or len(x) > 1200:
                raise ValueError("Curve x/y arrays must match and contain at most 1200 points.")
            if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in x + y):
                raise ValueError("Curve coordinates must be finite numbers.")
    for table in output.get("tables", []):
        if not isinstance(table, dict) or not isinstance(table.get("title"), str):
            raise ValueError("Each table needs a title.")
        columns, rows = table.get("columns"), table.get("rows")
        if not isinstance(columns, list) or any(not isinstance(c, str) for c in columns) or not isinstance(rows, list):
            raise ValueError("Tables need text columns and a row list.")
        if any(not isinstance(row, list) or len(row) != len(columns) for row in rows):
            raise ValueError("Table rows must match their columns.")
    for link in output.get("handoffs", []):
        if not isinstance(link, dict) or not isinstance(link.get("target"), str) or not isinstance(link.get("inputs"), dict):
            raise ValueError("A handoff needs a target application and explicit input map.")
        if not isinstance(link.get("notes", []), list) or any(not isinstance(n, str) for n in link.get("notes", [])):
            raise ValueError("Handoff notes must be text.")


def make_result(app: dict, inputs: dict, output: dict, material: dict | None) -> dict:
    ensure_json(output)
    validate_output(output)
    mat = validate_material(material)
    for key in ("metrics", "plots", "tables", "notes", "warnings", "handoffs"):
        if not isinstance(output.get(key, []), list):
            raise ValueError(f"Calculation returned an invalid {key} field.")
        output.setdefault(key, [])
    updates = output.get("material_updates", {})
    if not isinstance(updates, dict):
        raise ValueError("Material updates must be an object.")
    for key, value in updates.items():
        if key in ("schema", "sample_id", "name", "history", "lineage", "evidence", "study_context"):
            raise ValueError("A model may not overwrite the material identity.")
        if key == "parameters":
            if not isinstance(value, list):
                raise ValueError("Derived parameters must be a list.")
            mat["parameters"].extend(deepcopy(value))
        elif key in ("chemistry", "conditions"):
            mat[key].update(deepcopy(value))
        else:
            mat[key] = deepcopy(value)
    run_id = str(uuid.uuid4())
    mat["history"].append({"run_id": run_id, "application": app["id"], "version": app["version"], "kind": "model_calculation"})
    mat["history"] = mat["history"][-100:]
    validate_material(mat)
    result = {
        "schema": RESULT_SCHEMA,
        "run_id": run_id,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "application": {k: app[k] for k in ("id", "name", "version")},
        "material": mat,
        "inputs": deepcopy(inputs),
        "results": output,
        "references": deepcopy(app.get("references", [])),
        "study": study_snapshot(mat),
        "interpretation": "Model calculation. Numerical verification does not establish experimental validity for a named polymer.",
    }
    result["checksum_sha256"] = digest(result)
    return ensure_json(result)


def validate_exchange(document: dict) -> dict:
    ensure_json(document)
    if not isinstance(document, dict):
        raise ValueError("Choose a material, calculation, or workbench session JSON document.")
    schema = document.get("schema")
    if schema == MATERIAL_SCHEMA:
        validate_material(document)
    elif schema == RESULT_SCHEMA:
        if document.get("checksum_sha256") != digest(document):
            raise ValueError("Calculation checksum does not match. Re-export the original result; edited values must be entered as new inputs.")
        if not isinstance(document.get("material"), dict):
            raise ValueError("Calculation document has no material record.")
        validate_material(document["material"])
        if not isinstance(document.get("inputs"), dict) or not isinstance(document.get("results"), dict) or not isinstance(document.get("application"), dict):
            raise ValueError("Incomplete calculation document.")
        if any(not isinstance(document["application"].get(k), str) or not document["application"][k] for k in ("id", "name", "version")):
            raise ValueError("Calculation application requires id, name and version strings.")
        if any(not isinstance(document.get(k), str) for k in ("run_id", "created_utc")):
            raise ValueError("Calculation run id and creation time must be text.")
        if "study" in document and document["study"] != study_snapshot(document["material"]):
            raise ValueError("Study snapshot must agree with the specimen and declared model context.")
        validate_output(document["results"])
        references = document.get("references", [])
        if not isinstance(references, list) or any(not isinstance(r, dict) or not isinstance(r.get("title"), str) or not isinstance(r.get("url"), str) for r in references):
            raise ValueError("Calculation references require titles and URLs.")
    elif schema == SESSION_SCHEMA:
        if not isinstance(document.get("material"), dict):
            raise ValueError("Session document has no material record.")
        validate_material(document["material"])
        if "samples" in document:
            check_registry(document["samples"], document["material"])
        runs = document.get("runs")
        if not isinstance(runs, list) or len(runs) > 50:
            raise ValueError("A session may contain at most 50 calculation runs.")
        for run in runs:
            if not isinstance(run, dict) or run.get("schema") != RESULT_SCHEMA:
                raise ValueError("Session runs must be complete calculation records.")
            validate_exchange(run)
    else:
        raise ValueError("Unsupported exchange schema. This version reads polymer-material, polymer-result and polymer-session version 1.0.")
    return deepcopy(document)


def handoff(document: dict, target: str, target_fields: list[dict], destination: dict | None = None) -> dict:
    source = validate_exchange(document)
    if source.get("schema") != RESULT_SCHEMA:
        raise ValueError("A parameter handoff requires a completed calculation.")
    if destination is not None:
        destination = validate_material(destination)
        if destination["sample_id"] != source["material"]["sample_id"]:
            return {"material": destination, "inputs": {}, "notes": [
                "The selected specimen has a different identifier. Enter coefficients and geometry measured for this specimen."]}
        if destination.get("lineage") != source["material"].get("lineage"):
            raise ValueError("Source and destination disagree on specimen lineage.")
    links = source["results"].get("handoffs", [])
    link = next((x for x in links if isinstance(x, dict) and x.get("target") == target), None)
    if link is None:
        return {"material": destination if destination is not None else source["material"], "inputs": {}, "notes": ["Material identity and provenance transferred. No physical coefficient was inferred or mapped."]}
    allowed = {f["key"] for f in target_fields}
    inputs = link.get("inputs", {})
    if not isinstance(inputs, dict) or set(inputs) - allowed:
        raise ValueError("This handoff contains inputs unsupported by the target application.")
    return {"material": destination if destination is not None else source["material"], "inputs": deepcopy(inputs), "notes": link.get("notes", []), "source_run_id": source["run_id"]}

