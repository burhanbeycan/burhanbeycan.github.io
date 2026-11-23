from __future__ import annotations
import math
import numpy as np
from .model import finite, nernst_einstein, onsager
from .md_views import DOI, FIELDS, OPTIONS, PRESETS, md_result
from .records import FORM_GROUPS, PESKO, RECORDS, gui_inputs, peo_results

APP = {
    "id": "electrolyte", "name": "Polymer Electrolyte Lab", "version": "0.4.0",
    "tagline": "Ion motion, collective transport & reference frames",
    "description": "Compare the ideal Nernst–Einstein estimate with a user-supplied Onsager matrix. Preserve the difference between tracer diffusion and collective charge transport.",
    "accent": "#6453a0",
    "scope": "Isothermal linear response of a homogeneous binary 1:1 electrolyte in a neutral polymer-fixed frame. User-supplied transport parameters; no inferred material constants, morphology or battery lifetime.",
    "fields": [
        {"key":"model","label":"Transport description","type":"select","default":"nernst_einstein","group":"Model","help":"The collective model requires independently supplied coefficients in the same polymer-fixed reference frame.","options":[{"value":"nernst_einstein","label":"Nernst–Einstein · ideal baseline"},{"value":"onsager","label":"Onsager · collective coefficients"}]},
        {"key":"temperature_K","label":"Temperature","type":"number","default":298.15,"min":1,"max":5000,"unit":"K","group":"Tracer baseline","help":"Temperature of the supplied self-diffusion coefficients. No temperature extrapolation of D is performed."},
        {"key":"concentration_mol_L","label":"Salt concentration","type":"number","default":0.5,"min":1e-9,"max":100,"unit":"mol/L","group":"Tracer baseline","help":"For a 1:1 salt, c+ = c− = c_salt. Converted to mol/m³ internally."},
        {"key":"d_plus_m2_s","label":"Cation self-diffusivity","type":"number","default":1e-12,"min":1e-30,"max":1e-4,"unit":"m²/s","group":"Tracer baseline","help":"Tracer/self coefficient at the specified state; not a salt chemical diffusivity."},
        {"key":"d_minus_m2_s","label":"Anion self-diffusivity","type":"number","default":2e-12,"min":1e-30,"max":1e-4,"unit":"m²/s","group":"Tracer baseline","help":"No named-polymer default: the initial value is an illustrative input."},
        {"key":"Lpp","label":"Collective L++","type":"number","default":2e-13,"min":0,"max":1,"unit":"mol²/(J·m·s)","group":"Collective transport","show_if":{"key":"model","value":"onsager"},"help":"Includes same-species distinct correlations; do not substitute cD/RT without stating the ideal approximation."},
        {"key":"Lmm","label":"Collective L−−","type":"number","default":8e-13,"min":0,"max":1,"unit":"mol²/(J·m·s)","group":"Collective transport","show_if":{"key":"model","value":"onsager"},"help":"Must use the same normalization and reference frame as L++."},
        {"key":"Lpm","label":"Cross coefficient L+−","type":"number","default":3e-13,"min":-1,"max":1,"unit":"mol²/(J·m·s)","group":"Collective transport","show_if":{"key":"model","value":"onsager"},"help":"Symmetric Onsager matrix; |L+−| cannot exceed √(L++L−−)."},
    ],
    "presets": [
        {"name":"Ideal binary electrolyte","description":"Synthetic demonstration; coefficients are not assigned to PEO or another material.","inputs":{"model":"nernst_einstein","temperature_K":298.15,"concentration_mol_L":0.5,"d_plus_m2_s":1e-12,"d_minus_m2_s":2e-12}},
        {"name":"Negative t+ · admissible matrix","description":"Synthetic positive-definite Onsager matrix illustrating that a collective transference number can be negative.","inputs":{"model":"onsager","Lpp":2e-13,"Lmm":8e-13,"Lpm":3e-13}},
        {"name":"Symmetric collective response","description":"Synthetic symmetric diagonal coefficients: t+ = 1/2 when conductivity is nonzero.","inputs":{"model":"onsager","Lpp":4e-13,"Lmm":4e-13,"Lpm":1e-13}},
    ],
    "models": [
        {"name":"Nernst–Einstein estimate","equation":"σ_NE = F² c (D+ + D−)/(R T); t+_ideal = D+/(D+ + D−)","assumptions":["Binary monovalent 1:1 salt; electroneutral bulk.","Independent ions: collective correlations are omitted.","Input D and c refer to one specified state."],"references":["ion_pairing","codata"]},
        {"name":"Collective Onsager response","equation":"J_i = −Σ_j L_ij ∇μ̃_j; σ = F²(L++ + L−− − 2L+−); t+ = (L++ − L+−)/(L++ + L−− − 2L+−)","assumptions":["Neutral polymer-fixed reference frame, identical for all supplied L values.","Symmetric positive-semidefinite L; isothermal linear response.","t+ is the current fraction under an electrical force with no chemical-potential gradients; it is not automatically a Bruce–Vincent measurement."],"references":["frame_gap"]},
    ],
    "references": [
        {"id":"ion_pairing","title":"France-Lanord & Grossman (2019). Correlations from Ion Pairing and the Nernst–Einstein Equation.","url":"https://doi.org/10.1103/PhysRevLett.122.136001","year":2019},
        {"id":"frame_gap","title":"Transference Number in Polymer Electrolytes: Mind the Reference-Frame Gap (2022).","url":"https://doi.org/10.1021/jacs.2c02389","year":2022},
        {"id":"codata","title":"NIST CODATA 2022: exact defining constants and derived R and F.","url":"https://physics.nist.gov/cuu/pdf/wall_2022.pdf","year":2022},
    ],
}


APP["fields"][0]["options"].append({"value": "peo_matched", "label": "PEO/LiTFSI · matched tracer and EIS measurements"})
for field in APP["fields"]:
    if field["group"] == "Tracer baseline":
        field["show_if"] = {"key": "model", "value": ["nernst_einstein", "onsager"]}
for key, label, help_text in FORM_GROUPS:
    APP["fields"].append({"key": f"peo_{key}", "label": label, "type": "textarea",
                          "default": gui_inputs(PESKO)[f"peo_{key}"], "group": "PEO/LiTFSI specimen",
                          "show_if": {"key": "model", "value": "peo_matched"}, "help": help_text})
APP["presets"].extend({"name": name, "description": desc, "inputs": gui_inputs(record)}
                      for name, desc, record in RECORDS)
APP["models"].append({"name": "Dry PEO/LiTFSI matched-state comparison",
                      "equation": "r = Li/EO; c (mol/L) = 1000 ρ (g/mL) r/(M_EO + r M_LiTFSI) on the long-chain repeat basis; σ_EIS = l/(A Rb).",
                      "assumptions": ["Dry binary formulation, anhydrous salt, a single reported composition and thermal state.",
                                      "An independent mixture volume, density or concentration at the same temperature is required for the volumetric NE estimate."],
                      "references": ["pesko2017", "zheng2018", "ion_pairing", "frame_gap"]})
APP["references"].extend([
    {"id": "pesko2017", "title": "Pesko et al. Negative Transference Numbers in Poly(ethylene oxide)-Based Electrolytes.",
     "url": "https://doi.org/10.1149/2.0581711jes", "year": 2017},
    {"id": "zheng2018", "title": "Zheng et al. Optimizing Ion Transport in Polyether-Based Electrolytes for Lithium Batteries.",
     "url": "https://doi.org/10.1021/acs.macromol.7b02706", "year": 2018},
])


APP["fields"][0]["options"].extend({"value":value, "label":label} for value,label in OPTIONS)
APP["fields"].extend(FIELDS)
APP["presets"].extend(PRESETS)
APP["scope"] += " The MD modes compare published first-shell coordination, atom-pair persistence, and formal charge displacement in dry PEO/LiTFSI."
APP["references"].append({"id":"force_scaling","title":"Charge Scaling in Classical Force Fields for Lithium Ions in Polymers (2024).","year":2024,"url":"https://doi.org/"+DOI})
APP["models"].append({"name":"Oxygen-resolved Li shells and collective transport",
    "equation":"CN = contacts/(frames × Li); σ = e² slope⟨|Σ zᵢ Δrᵢ|²⟩/(6VkBT).",
    "assumptions":["Count distinct oxygen atoms within the declared first-shell cutoff.","Use unwrapped positions and all ions of the neutral box for the collective term."],
    "references":["force_scaling","ion_pairing","frame_gap"]})


def calculate(inputs, material=None):
    if not isinstance(inputs, dict):
        raise ValueError("Inputs must be a dictionary.")
    if material is not None and not isinstance(material, dict):
        raise ValueError("Material metadata must be a dictionary.")
    p = {f["key"]: f["default"] for f in APP["fields"]}
    unknown = set(inputs) - set(p)
    if unknown:
        raise ValueError("Unknown input fields: " + ", ".join(sorted(str(key) for key in unknown)))
    p.update(inputs)
    if not isinstance(p["model"], str) or p["model"] not in {option["value"] for option in APP["fields"][0]["options"]}:
        raise ValueError("Choose a listed electrolyte model.")
    # The backend enforces the declared GUI ranges even for hidden fields.
    for field in APP["fields"]:
        if field["type"] == "number":
            p[field["key"]] = finite(p[field["key"]], field["label"],
                                     lower=field.get("min"), upper=field.get("max"))
    if p["model"] in {value for value,_ in OPTIONS}:
        return md_result(p)
    if p["model"] == "peo_matched":
        return peo_results(p)
    T = p["temperature_K"]
    c_L = p["concentration_mol_L"]
    ideal = nernst_einstein(c_L * 1000, p["d_plus_m2_s"], p["d_minus_m2_s"], T)
    metrics = [
        {"key":"sigma_ne","label":"Nernst–Einstein estimate","value":ideal["conductivity_S_m"],"unit":"S/m","note":"Correlation-free baseline"},
        {"key":"sigma_ne_display","label":"Nernst–Einstein estimate","value":ideal["conductivity_S_m"]*10,"unit":"mS/cm","note":"1 S/m = 10 mS/cm"},
        {"key":"t_ideal","label":"Ideal cation current fraction","value":ideal["cation_current_fraction_ideal"],"unit":"1","note":"Derived only under the independent-ion assumption"},
    ]
    warnings = ["The initial coefficients are synthetic demonstration inputs, not experimentally validated polymer parameters."]
    notes = ["The chemical material record is preserved. Molar mass, tacticity, sequence, end groups and crystallinity do not generate transport coefficients in this model.", "Ion self-diffusion is never handed to Transport Lab as a Fickian penetrant or salt chemical diffusion coefficient.", "The neutral polymer defines the velocity reference. A solvent-fixed frame is equivalent only if the solvent has the same mean velocity; no reference-frame conversion is performed."]
    if p["model"] == "onsager":
        collective = onsager(p["Lpp"], p["Lmm"], p["Lpm"])
        metrics.extend([
            {"key":"sigma_collective","label":"Collective conductivity","value":collective["conductivity_S_m"],"unit":"S/m","note":"From the supplied Onsager matrix"},
            {"key":"t_collective","label":"Collective t+","value":collective["t_plus_polymer_fixed"],"unit":"1","note":"Polymer-fixed frame; undefined if charge conductivity vanishes"},
        ])
        if collective["t_plus_polymer_fixed"] is None:
            warnings.append("The charge-conduction mode is zero or numerically singular. t+ is undefined and is not reported as zero.")
        pp, mm = float(p["Lpp"]), float(p["Lmm"])
        bound = math.sqrt(pp)*math.sqrt(mm)
        plots = []
        if bound > 0:
            rho = np.linspace(-0.95, 0.95, 191)
            sigmas, ts = [], []
            for r in rho:
                value = onsager(pp, mm, float(r)*bound)
                sigmas.append(value["conductivity_S_m"])
                ts.append(value["t_plus_polymer_fixed"])
            plots.append({"title":"Cross-correlation sensitivity · fixed diagonal coefficients","x_label":"L+− / √(L++ L−−)","y_label":"Conductivity (S/m)","x_scale":"linear","y_scale":"linear","series":[{"name":"Admissible symmetric matrices","x":rho.tolist(),"y":sigmas}]})
            if all(v is not None for v in ts):
                plots.append({"title":"Reference-frame-specific current fraction","x_label":"L+− / √(L++ L−−)","y_label":"t+ (dimensionless)","x_scale":"linear","y_scale":"linear","series":[{"name":"Polymer-fixed t+","x":rho.tolist(),"y":ts}]})
            notes.append("The correlation sweep is a controlled mathematical sensitivity study, not a prediction of how salt concentration changes the matrix.")
        else:
            notes.append("The cross-correlation ratio is undefined because a diagonal coefficient is zero. Positive semidefiniteness then requires L+−=0; no correlation sweep is plotted.")
        rows = [["L++",float(p["Lpp"]),"mol²/(J·m·s)","User supplied"],["L−−",float(p["Lmm"]),"mol²/(J·m·s)","User supplied"],["L+−",float(p["Lpm"]),"mol²/(J·m·s)","User supplied"],["Minimum eigenvalue",collective["minimum_eigenvalue"],"mol²/(J·m·s)","Numerical matrix check"]]
    else:
        concentrations = np.linspace(max(1e-9, c_L * 0.1), min(100, c_L * 2), 150)
        plots = [{"title":"Concentration sensitivity · D held fixed","x_label":"Salt concentration (mol/L)","y_label":"Ideal conductivity (S/m)","x_scale":"linear","y_scale":"linear","series":[{"name":"Nernst–Einstein baseline","x":concentrations.tolist(),"y":[ideal["conductivity_S_m"]*float(x)/c_L for x in concentrations]}]}]
        rows = [["L++ ideal",ideal["Lpp_ideal"],"mol²/(J·m·s)","cD+/(RT)"],["L−− ideal",ideal["Lmm_ideal"],"mol²/(J·m·s)","cD−/(RT)"],["L+− ideal",0,"mol²/(J·m·s)","Independent-ion assumption"]]
        notes.append("The concentration curve holds D+ and D− fixed. Real concentration dependence requires independent coefficients at each state.")
    return {"metrics":metrics,"plots":plots,"tables":[{"title":"Transport coefficient definitions","columns":["Quantity","Value","Unit","Origin"],"rows":rows}],"notes":notes,"warnings":warnings,"material_updates":{"parameters":[{"key":"electrolyte_sigma_ne","value":ideal["conductivity_S_m"],"unit":"S/m","quantity":"ideal_ionic_conductivity","origin":"model_calculation","model_scope":"Independent monovalent 1:1 ions at user-specified c, D, T"}]},"handoffs":[]}

