"""Academic GUI adapter; chemistry is implemented in the model module."""
from __future__ import annotations

from . import experiments, microstructure
from .model import TEMPLATES, build_oligomer, molecular_view, sequence_statistics

REFERENCES = [
    {"id": "rdkit_graph", "title": "RDKit: Getting Started in Python — molecule graphs, stereochemistry and descriptors", "url": "https://www.rdkit.org/docs/GettingStartedInPython.html"},
    {"id": "rdkit_zip", "title": "RDKit rdmolops: molzip and attachment-point labels", "url": "https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html"},
    {"id": "rdkit_mass", "title": "RDKit rdMolDescriptors: molecular formula and exact mass", "url": "https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html"},
    {"id": "rdkit_stereo", "title": "The RDKit Book — stereochemistry representation", "url": "https://www.rdkit.org/docs/RDKit_Book.html"},
    {"id": "petox_chemistry", "title": "Chiu, Thill and Fairchok: Poly(2-ethyl-2-oxazoline), ACS Advances in Chemistry 213, Chapter 23", "url": "https://doi.org/10.1021/ba-1986-0213.ch023", "year": 1986},
]

APP = {
    "id": "chemistry", "name": "Polymer Chemistry Lab", "version": "0.4.0", "accent": "#936326",
    "tagline": "Connect repeat units, resolve end groups, and inspect molecular sequence.",
    "description": "Build a finite linear oligomer from explicit polymerized repeat connectivity; inspect atom/bond structure, formula, end-group-aware molar mass and A/B sequence dyads.",
    "scope": "A graph constructor and chemistry ledger. H/H graph caps are stated assumptions; no reaction mechanism, conformational energy, tacticity distribution or material-property accuracy is inferred.",
    "fields": [
        {"key": "template_a", "label": "Repeat A · mapped SMILES", "type": "textarea", "default": TEMPLATES["petox"]["smiles"], "unit": "polymerized connectivity", "group": "Repeat chemistry", "help": "Exactly two neutral single-bond dummy ports [*:1] and [*:2]. This is the polymerized repeat, not the cyclic oxazoline monomer. Real atoms may not have map numbers."},
        {"key": "template_b", "label": "Repeat B · mapped SMILES", "type": "textarea", "default": TEMPLATES["pmetox"]["smiles"], "unit": "polymerized connectivity", "group": "Repeat chemistry", "help": "Used only if the explicit sequence contains B. Every unit is oriented port 1 → port 2; no reaction compatibility is inferred."},
        {"key": "sequence", "label": "Finite chain sequence", "type": "textarea", "default": "AAAAAA", "unit": "A/B units", "group": "Sequence and end chemistry", "help": "1–64 explicit A/B symbols (at most 512 non-port atoms). AAABBB defines one block oligomer; ABABAB defines one alternating oligomer. It is not an ensemble of stochastic copolymers."},
    ],
    "presets": [
        {"name": "PEtOx · N-propionyl backbone", "description": "Six polymerized N-propionylethyleneimine units. H/H caps are graph assumptions, not actual CROP initiator/terminator groups.", "inputs": {"template_a": TEMPLATES["petox"]["smiles"], "sequence": "AAAAAA"}},
        {"name": "PEtOx / PMeOx · explicit block", "description": "Three EtOx-derived then three MeOx-derived repeats; sequence is explicit, without a statistical sequence model.", "inputs": {"template_a": TEMPLATES["petox"]["smiles"], "template_b": TEMPLATES["pmetox"]["smiles"], "sequence": "AAABBB"}},
        {"name": "PEtOx / PMeOx · alternating", "description": "Same unit counts as the block example, different AB/BA dyads; molar masses agree but connectivity differs.", "inputs": {"template_a": TEMPLATES["petox"]["smiles"], "template_b": TEMPLATES["pmetox"]["smiles"], "sequence": "ABABAB"}},
        *[{"name": item["name"] + " · six repeats", "description": "Connectivity and elemental composition template. Stereochemistry remains unspecified unless encoded explicitly in the SMILES; no atactic distribution is generated.", "inputs": {"template_a": item["smiles"], "sequence": "AAAAAA"}} for key, item in TEMPLATES.items() if key in ("pe", "peo", "pva", "pmma")],
    ],
    "models": [
        {"name": "Oriented repeat graph assembly", "equation": "tail(i)–head(i+1); cap(head(1))=H; cap(tail(n))=H", "assumptions": ["Already-polymerized repeat connectivity", "Two neutral single-bond dummy attachment points per repeat", "One explicit finite sequence, with each repeat oriented identically", "H/H model caps; no initiator or termination reaction inferred"], "references": ["rdkit_graph", "rdkit_zip"]},
        {"name": "Repeat and finite-chain molar masses", "equation": "M(chain)=Σ_i M(repeat_i)+2M(H); homopolymer DP=[M(chain)−2M(H)]/M(repeat)", "assumptions": ["Average molar masses use RDKit's isotopically averaged atomic-weight convention except explicitly labelled isotopes", "Exact mass uses explicit isotopes or RDKit's most abundant isotope masses", "A finite single molecule is not a sample Mn or Mw measurement", "End-group correction is 2H only for the H/H model in this module"], "references": ["rdkit_mass"]},
        {"name": "Sequence dyads", "equation": "f_A=N_A/n; f_XY=N_XY/(n−1)", "assumptions": ["Directed noncyclic chain; n−1 adjacent pairs", "No random-copolymer kinetics or reactivity ratios inferred", "Explicit input stereochemistry is retained; no tacticity assignment is generated"], "references": ["rdkit_stereo"]},
    ],
    "references": REFERENCES,
}

# Preserve v0.2 input replay and default graph behaviour. Presets reset the
# mode explicitly when a researcher switches back from a comparison.
for field in APP["fields"]:
    field["show_if"] = {"key": "chemistry_mode", "value": "graph"}
for preset in APP["presets"]:
    preset["inputs"]["chemistry_mode"] = "graph"
APP["fields"] = [experiments.MODE_FIELD, *APP["fields"], *experiments.FIELDS]
APP["fields"][0]["options"].extend([
    {"value":"pva_benchmark","label":"PVA · independent NMR grade holdout"},
    {"value":"petox_hydrolysis","label":"WBPU blend partner · explicit PEtOx/EI ion states"},
])
APP["fields"].extend(microstructure.FIELDS)
for f in APP["fields"]:
    if f["key"] in ("ensemble_dp","ensemble_count","ensemble_seed"):
        f["show_if"] = {"key":"chemistry_mode", "value":[*experiments.PVA_MODES, "pva_benchmark"]}
APP["presets"] = [*microstructure.PRESETS, *experiments.PRESETS, *APP["presets"]]
APP["models"].extend([*experiments.MODELS, *microstructure.MODELS])
APP["references"].extend([*experiments.REFERENCES, *microstructure.REFERENCES])
APP["description"] = "Build finite oligomers, examine PEtOx/EI hydrolysis and charge, and test PVA sequence models against experimental NMR grades."
APP["scope"] = "Graph identity, stated positive-ion mass targets and stationary sequence models. Chain-specific assignments require independent structural evidence."


def calculate(inputs, material=None):
    if not isinstance(inputs, dict):
        raise ValueError("Calculation inputs must be an object.")
    keys = {field["key"] for field in APP["fields"]}
    if set(inputs) - keys:
        raise ValueError("Unknown chemistry input(s): " + ", ".join(sorted(set(inputs) - keys)))
    values = {field["key"]: field["default"] for field in APP["fields"]}
    values.update(inputs)
    mode = values["chemistry_mode"]
    if mode == "pva_benchmark":
        return microstructure.pva_benchmark(values)
    if mode == "petox_hydrolysis":
        return microstructure.hydrolysis(values)
    if mode == "petox":
        return experiments.compare_ends(values)
    if mode in ("pva_nmr", "pva_random"):
        return experiments.model_pva(values)
    if mode != "graph":
        raise ValueError("Unknown chemistry analysis mode.")
    oligomer = build_oligomer(values["template_a"], values["sequence"], values["template_b"])
    summary = oligomer.summary
    statistics = sequence_statistics(oligomer.sequence)
    warnings = []
    if summary["unassigned_tetrahedral_centres"]:
        warnings.append("Some possible tetrahedral stereocentres are unspecified. An unspecified structure does not represent a calculated atactic, isotactic or syndiotactic population.")
    if summary["formal_charge"]:
        warnings.append("This is a charged molecular graph without separate counterions. No solution speciation, electroneutral formulation or ion association model is implied.")
    notes = [
        "The displayed structure is a finite linear model with explicitly stated H/H caps. For PEtOx, the N-propionyl repeat is distinct from the cyclic 2-ethyl-2-oxazoline monomer; the N-terminal model cap forms a secondary amide.",
        "Atom coordinates are a 2D drawing in arbitrary units. They do not describe an energy-minimized chain, chain dimensions, solvent conformation or force-field topology.",
        "The mass of one constructed molecule is not an experimental number-average molar mass (Mn), mass-average molar mass (Mw), or molar-mass distribution.",
        "Formula and mass encode connectivity and isotope/charge input. They do not establish synthesis feasibility or property-prediction accuracy. RDKit version: " + summary["rdkit_version"],
        "Every repeat is oriented from port 1 to port 2. A/B sequence statistics describe only this finite directed chain; no stochastic ensemble is generated.",
    ]
    handoffs = []
    if len(set(oligomer.sequence)) == 1:
        symbol = oligomer.sequence[0]
        repeat_mass = oligomer.repeats[symbol]["average_mass_g_mol"]
        handoffs.append({"target": "chain", "label": "Transfer repeat and end-group masses", "inputs": {"repeat_mass_g_mol": repeat_mass, "end_group_mass_g_mol": summary["end_group_mass_g_mol"]}, "notes": ["Transferred the average polymerized repeat mass and explicit 2H cap correction for this homopolymer model.", "A constructed finite molecule is not a measured MWD. Supply an independent sample MWD, contour increment and chain stiffness in Chain Lab. Use a different end correction for real initiator/terminator groups."]})
    else:
        notes.append("No single-repeat mass is transferred for this mixed-repeat chain. A copolymer sample requires composition/sequence-resolved mass-to-contour mapping.")
    chem = {
        "repeat_unit": "; ".join(key + ": " + info["smiles"] for key, info in oligomer.repeats.items()),
        "architecture": "linear", "sequence": oligomer.sequence,
        "end_groups": "H/H computational graph caps (2.016 g/mol with current RDKit atomic weights); not inferred synthesis end groups",
        "tacticity": "not inferred; stereochemistry only as explicitly encoded in the molecular graph",
        "finite_oligomer_smiles": summary["canonical_isomeric_smiles"],
        "end_group_evidence": "Illustrative graph caps.",
    }
    metrics = [
        {"key": "dp", "label": "Explicit repeat count", "value": len(oligomer.sequence), "unit": "units", "note": "Single finite chain, not DPn of a sample."},
        {"key": "formula", "label": "H-capped chain formula", "value": summary["formula"], "unit": "", "note": "Includes the two H model caps."},
        {"key": "mass", "label": "Average molecular molar mass", "value": summary["average_molar_mass_g_mol"], "unit": "g/mol", "note": "One finite graph; not experimental Mn."},
        {"key": "exact_mass", "label": "Exact molecular mass", "value": summary["exact_mass_Da"], "unit": "Da", "note": "RDKit isotopic mass convention."},
        {"key": "cap_mass", "label": "Explicit end correction", "value": summary["end_group_mass_g_mol"], "unit": "g/mol", "note": "Two protium caps, average atomic-weight convention."},
        {"key": "heavy_atoms", "label": "Heavy atoms", "value": summary["heavy_atoms"], "unit": "atoms", "note": "Atom/bond graph constructed and sanitized."},
    ]
    return {"metrics": metrics, "plots": [], "tables": [
        {"title": "Polymerized repeat contributions · no terminal H caps", "columns": ["Symbol", "Formula", "Average mass (g/mol)", "Exact mass (Da)", "Count", "Unit mole fraction"], "rows": [[key, info["formula"], info["average_mass_g_mol"], info["exact_mass_Da"], statistics["counts"][key], statistics["fractions"][key]] for key, info in oligomer.repeats.items()]},
        {"title": "Explicit directed sequence dyads", "columns": ["Dyad", "Count", "Fraction of n−1 bonds"], "rows": [[key, value, statistics["dyad_fractions"][key]] for key, value in statistics["dyad_counts"].items()]},
        {"title": "Structure and stereochemistry ledger", "columns": ["Quantity", "Value"], "rows": [["Canonical isomeric SMILES", summary["canonical_isomeric_smiles"]], ["Explicit sequence", oligomer.sequence], ["Assigned tetrahedral centres", summary["assigned_tetrahedral_centres"]], ["Unassigned possible tetrahedral centres", summary["unassigned_tetrahedral_centres"]], ["Formal charge", summary["formal_charge"]], ["RDKit version", summary["rdkit_version"]]]},
    ], "molecular_views": [molecular_view(oligomer)], "notes": notes, "warnings": warnings, "handoffs": handoffs,
        "structure": {**summary, "repeat_information": oligomer.repeats}, "material_updates": {"chemistry": chem}}
