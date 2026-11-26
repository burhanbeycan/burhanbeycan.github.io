"""GUI contracts for end-group comparisons and NMR-constrained PVA ensembles."""
from .end_groups import END_PAIRS, petox_candidates
from .model import TEMPLATES, build_oligomer, molecular_view
from .sequences import pva_ensemble

NMR_SOURCE = "Budhlall et al. (2000), Macromol. Symp. 155, 63–84, eqs. 3–5, ¹³C NMR, Poval 217."
NMR_URL = "https://bpb-us-w2.wpmucdn.com/sites.uml.edu/dist/4/136/files/2018/09/Budhlall_2000_Macromol.Symposia_PVOHviaNMRandHPLC-22fn4h7.pdf"
PVAC_REPEAT = "[*:1]CC(OC(=O)C)[*:2]"
MODES = [("graph", "Explicit A/B graph"), ("petox", "PEtOx end groups & ion masses"),
         ("pva_random", "PVA · hydrolysis fraction only"), ("pva_nmr", "PVA · NMR dyad constraints")]


def _number(key, label, value, limits, modes, unit, help_text, step="any"):
    return {"key": key, "label": label, "type": "number", "default": value,
            "min": limits[0], "max": limits[1], "step": step, "unit": unit,
            "show_if": {"key": "chemistry_mode", "value": modes},
            "group": "Model inputs", "help": help_text}


MODE_FIELD = {"key": "chemistry_mode", "label": "Chemistry analysis", "type": "select",
              "default": "graph", "group": "Analysis",
              "options": [{"value": key, "label": label} for key, label in MODES]}
PVA_MODES = ["pva_random", "pva_nmr"]
FIELDS = [
    _number("petox_dp", "PEtOx repeat count", 6, (1, 64), "petox", "units",
            "All three candidates have the same N-propionylethyleneimine backbone.", 1),
    {"key": "selected_ends", "label": "Assign a candidate to the model record", "type": "select",
     "group": "End-group assignment", "default": "none", "show_if": {"key": "chemistry_mode", "value": "petox"},
     "help": "Compare only preserves the current sample identity. Selecting a model enables an explicit mass handoff.",
     "options": [{"value": "none", "label": "Compare only"}] +
                [{"value": k, "label": vals[0]} for k, vals in END_PAIRS.items()]},
    {"key": "end_evidence", "label": "Sample-specific end-group evidence (optional)", "type": "textarea",
     "default": "", "group": "End-group assignment", "show_if": {"key": "chemistry_mode", "value": "petox"},
     "help": "Record a spectrum identifier and assigned end-retaining product ions if available. Free text records an interpretation."},
    _number("oh_fraction", "Vinyl alcohol units after hydrolysis", 0.881, (0,1), "pva_random", "mole fraction",
            "Bernoulli reference. Hydrolysis degree alone leaves sequence correlations undetermined."),
    *[_number("nmr_"+k, label, value, (0,1), "pva_nmr", "bond fraction",
              "Use normalized, quantitative dyad areas. Mixed includes both orientations.")
      for k, label, value in [("oo", "OH–OH", .837), ("mixed", "OH–OAc (both directions)", .088),
                              ("aa", "OAc–OAc", .075)]],
    {"key": "dyad_source", "label": "NMR spectrum / data source", "type": "textarea",
     "default": "", "group": "Model inputs", "show_if": {"key": "chemistry_mode", "value": "pva_nmr"},
     "help": "Give the grade/lot, nucleus and spectrum, or paper and table. Published default fractions are identified in the result."},
    _number("ensemble_dp", "Length of each simulated sequence", 30, (2,64), PVA_MODES, "units",
            "A finite illustration. The reported NMR dyads characterize the sample population.", 1),
    _number("ensemble_count", "Number of simulated chains", 400, (1,4096), PVA_MODES, "chains",
            "The analytical target and finite sampled population are both reported.", 1),
    _number("ensemble_seed", "Sequence random seed", 23, (0,2**32-1), PVA_MODES, "integer",
            "Reuse this value to reproduce the same sample with the same software environment.", 1),
]
PRESETS = [
    {"name": "PEtOx · same DP, three end groups", "description": "Compare neutral graphs and explicit positive-ion m/z targets.",
     "inputs": {"chemistry_mode": "petox", "petox_dp": 6, "selected_ends": "none", "end_evidence": ""}},
    {"name": "PVA · Poval 217 (¹³C NMR dyads)", "description": "Reproduce the published long-chain OH fraction, η and mean run lengths from three measured dyad fractions.",
     "inputs": {"chemistry_mode": "pva_nmr", "nmr_oo": .837, "nmr_mixed": .088, "nmr_aa": .075,
                "dyad_source": "", "ensemble_dp": 30, "ensemble_count": 400, "ensemble_seed": 23}},
    {"name": "PVA · equal-hydrolysis Bernoulli reference", "description": "Compare with the Poval 217 preset at the same OH mole fraction.",
     "inputs": {"chemistry_mode": "pva_random", "oh_fraction": .881, "ensemble_dp": 30, "ensemble_count": 400, "ensemble_seed": 23}},
]
REFERENCES = [
    {"id":"petox_ends", "title":"Morgan et al. (2022), end-resolved PEtOx / EI oligomer tandem MS", "url":"https://doi.org/10.1039/D2PY00437B"},
    {"id":"pva_dyads", "title":"Budhlall et al. (2000), measured PVA dyads and sequence statistics", "url":NMR_URL},
    {"id":"pva_context", "title":"Budhlall et al. (2003), PVA architecture and aqueous conformation", "url":"https://doi.org/10.1021/ma030027d"},
    {"id":"adducts", "title":"UC Davis Fiehn Lab, MS Adduct Calculator (conventional positive-ion shifts)", "url":"https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator"},
]
MODELS = [
    {"name":"PEtOx end-group assignment","equation":"M = n M(C5H9NO) + M(α + ω); m/z = (m_neutral + Δm_ion)/|z|",
     "assumptions":["α-Me bonds to backbone N; ω-OH or ω-N3 bonds to backbone C", "Linear unhydrolysed neutral chains; theoretical positive adducts"],
     "references":["petox_ends","adducts","rdkit_mass"]},
    {"name":"Partially hydrolysed PVAc · stationary binary ensemble","equation":"h = f(OH,OH) + f(mixed)/2; η = f(mixed) / [2h(1−h)]",
     "assumptions":["Quantitative dyad fractions with negligible end correction", "Stationary first-order Markov sequences; no stereochemical sequence model"],
     "references":["pva_dyads","pva_context"]},
]


def _text(value, label):
    if not isinstance(value, str) or len(value) > 2000:
        raise ValueError(f"{label} must be text up to 2000 characters.")
    return value.strip()


def compare_ends(values):
    dp = values["petox_dp"]
    candidates = petox_candidates(dp)
    selected = values["selected_ends"]
    evidence = _text(values["end_evidence"], "End-group evidence")
    if not isinstance(selected, str) or (selected != "none" and selected not in END_PAIRS):
        raise ValueError("Choose a supported end-group assignment.")
    result = {"metrics":[{"key":"dp", "label":"Common PEtOx repeat count", "value":dp,"unit":"units"}],
              "molecular_views":[c["view"] for c in candidates], "plots":[], "tables":[
                  {"title":"Neutral molecules and theoretical positive-ion targets",
                   "columns":["Capping model","Formula","Average M (g/mol)","Neutral exact (Da)","[M+H]⁺ (m/z)","[M+Na]⁺ (m/z)","[M+2H]²⁺ (m/z)"],
                   "rows":[[c[k] for k in ("name","formula","avg_g_mol","neutral_Da","mz_H","mz_Na","mz_2H")] for c in candidates]},
                  {"title":"Covalent terminal assignments", "columns":["Model","α, attached to N","ω, attached to C","End correction (g/mol)","Canonical SMILES"],
                   "rows":[[c[k] for k in ("name","alpha","omega","end_mass_g_mol","smiles")] for c in candidates]}],
              "notes":["Methyl tosylate initiation places methyl on the α backbone N. The hydroxyl and azide candidates correspond to the termination assignments reported by Morgan et al. (2022). H/H is an illustrative graph reference.",
                       "The table gives prospective m/z values using a stated ion composition and charge. Precursor mass agreement alone cannot establish a terminal assignment. End-retaining MS/MS product ions, with independent assignment and sample identity, supply structural evidence.",
                       "Graph comparison uses unhydrolysed PEtOx. Peaks from partially hydrolysed EtOx/EI, fragments, isotopologues, salts and other adducts require their own compositions."],
              "warnings":[], "handoffs":[], "material_updates":{},
              "end_group_models":[{k:v for k,v in c.items() if k!="view"} for c in candidates]}
    if evidence:
        result["tables"].append({"title":"User-supplied evidence", "columns":["Evidence record","Evaluation"],
                                 "rows":[[evidence,"Recorded as supplied; assignment requires spectrum review."]]})
    if selected != "none":
        c = next(c for c in candidates if c["key"] == selected)
        repeat_mass = build_oligomer(TEMPLATES["petox"]["smiles"],"A").repeats["A"]["average_mass_g_mol"]
        result["handoffs"].append({"target":"chain","label":"Transfer selected end model",
                                   "inputs":{"repeat_mass_g_mol":repeat_mass,"end_group_mass_g_mol":c["end_mass_g_mol"]},
                                   "notes":["Selected PEtOx model: " + c["name"] + ". The user supplies the independent MWD."]})
        result["material_updates"] = {"chemistry":{"architecture":"linear","repeat_unit":TEMPLATES["petox"]["smiles"],
                                            "sequence":"A"*int(dp), "end_groups":"Selected model: "+c["name"],
                                            "finite_oligomer_smiles":c["smiles"],
                                            "end_group_evidence":evidence or "No sample-specific evidence supplied."}}
    return result


def model_pva(values):
    mode = values["chemistry_mode"]
    keys = ("nmr_oo","nmr_mixed","nmr_aa")
    target = [values[k] for k in keys]
    r = pva_ensemble(mode, values["oh_fraction"], target, values["ensemble_dp"], values["ensemble_count"], values["ensemble_seed"])
    source = "Bernoulli model at user-specified OH mole fraction."
    if mode == "pva_nmr":
        provided = _text(values["dyad_source"], "Dyad source")
        if not provided and target != [.837, .088, .075]:
            raise ValueError("Give the source of the supplied dyads, or label them as a hypothetical scenario.")
        source = provided or NMR_SOURCE
    seq = r["examples"][0]
    example = build_oligomer(TEMPLATES["pva"]["smiles"], seq, PVAC_REPEAT)
    view = molecular_view(example)
    view.update(title="Simulated PVA/PVAc chain 1", caption="A is a vinyl alcohol unit; B is a residual vinyl acetate unit. The first seeded draw is shown with illustrative H/H caps and unspecified tacticity.")
    return {"metrics":[{"key":"oh","label":"Target OH repeat fraction","value":r["oh_fraction"],"unit":"mol/mol"},
                       {"key":"eta","label":"Blockiness index η","value":r["eta"],"unit":"",
                        "note":"1 for independent units; undefined at pure composition."},
                       {"key":"sample_oh","label":"Realized OH repeat fraction","value":r["sample_oh_fraction"],"unit":"mol/mol"}],
            "plots":[{"title":"Variation of hydrolysis fraction across simulated chains",
                      "x_label":"OH units per finite chain","y_label":"Simulated chain count",
                      "series":[{"name":"Equal-DP model population",
                                 "x":list(range(int(values["ensemble_dp"])+1)),
                                 "y":[dict(r["histogram"]).get(i,0) for i in range(int(values["ensemble_dp"])+1)]}]}],
            "tables":[{"title":"Dyad evidence and simulated population","columns":["Pair","Target bond fraction","Simulated bond fraction"],
                       "rows":[[k,a,b] for k,a,b in zip(("OH–OH","OH–OAc (both directions)","OAc–OAc"),r["target"],r["observed_dyads"])]},
                      {"title":"Long-chain run statistics","columns":["Unit","P(switch | unit)","Mean run length (units)"],
                       "rows":[["OH",r["switch_oh_to_oac"],r["mean_run_oh"]],["OAc",r["switch_oac_to_oh"],r["mean_run_oac"]]]},
                      {"title":"First seeded model sequences (A = OH, B = OAc)","columns":["Draw","Sequence"],
                       "rows":[[idx+1,s] for idx,s in enumerate(r["examples"])]},
                      {"title":"Input provenance","columns":["Basis","Record"],"rows":[[mode,source]]}],
            "molecular_views":[view],"sequence_ensemble":r,"handoffs":[],"material_updates":{},
            "warnings":["Hydrolysis degree and dyad averages do not uniquely reconstruct a sample chain, tacticity, branching or higher-order correlations. The single drawing is a simulated example."],
            "notes":["The NMR mode uses normalized, quantitatively assigned OH–OH, total mixed and OAc–OAc dyads. Both mixed directions have equal expectation in the stationary model.",
                     "Run lengths are bulk geometric means; a zero switch probability gives an unbounded run for a present state. An absent state has no defined run length. Both cases display a blank in the table.",
                     "Simulated chains have finite length and independent stationary starting units. The empirical dyad frequencies fluctuate around the analytical target.",
                     "No simulated sequence or example mass is assigned as an experimentally identified sample chain."]}
