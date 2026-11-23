"""Scientific GUI adapter for polymer reaction balances."""
import numpy as np
from .model import finite, living_chain_growth, step_growth, terminal_copolymer


def number(key, label, default, unit, lower, upper, group, help_text, model):
    return {"key": key, "label": label, "default": default, "type": "number", "unit": unit,
            "min": lower, "max": upper, "group": group, "help": help_text,
            "show_if": {"key": "model", "value": model}}


REFERENCES = [
    {"id": "carothers", "title": "Meaurio (2024). Generalization of the Carothers equation for linear step growth polymers.", "year": 2024, "url": "https://doi.org/10.1038/s41598-024-78476-2"},
    {"id": "iupac_chain", "title": "IUPAC Gold Book: chain polymerization", "url": "https://old.goldbook.iupac.org/html/C/C00958.html"},
    {"id": "iupac_copolymer", "title": "IUPAC Gold Book: copolymerisation composition equation", "url": "https://goldbook.iupac.org/terms/view/15391"},
    {"id": "iupac_ratios", "title": "IUPAC recommended experimental methods and data evaluation procedures for determination of radical copolymerization reactivity ratios (2024).", "year": 2024, "url": "https://doi.org/10.1039/D4PY00270A"},
]

APP = {
    "id": "synthesis", "name": "Polymer Synthesis Lab", "version": "0.2.0", "accent": "#b66b33",
    "tagline": "From reactive equivalents to chain count, end groups and composition.",
    "description": "Track functional-group conversion and mass conservation; distinguish theoretical feed bookkeeping from measured molar mass and kinetic prediction.",
    "scope": "Ideal acyclic bifunctional step growth, ideal chain-count bookkeeping and fixed-feed terminal copolymerization. Conversion, initiation efficiency and reactivity ratios are supplied by the user.",
    "fields": [
        {"key": "model", "label": "Polymer chemistry model", "type": "select", "default": "step_growth", "group": "Model", "unit": "", "help": "Select a bounded model with explicit reaction and population assumptions.", "options": [
            {"value": "step_growth", "label": "AA + BB · functional-group and mass balances"},
            {"value": "living", "label": "Ideal chain growth · initiated-chain bookkeeping"},
            {"value": "copolymer", "label": "Terminal copolymer · instantaneous composition"}]},
        number("n_aa_mol", "AA monomer amount", .01, "mol", 1e-12, 1e6, "Bifunctional feed", "Moles of molecules with exactly two A groups. Functionality is not a fit parameter.", "step_growth"),
        number("n_bb_mol", "BB monomer amount", .01, "mol", 1e-12, 1e6, "Bifunctional feed", "Moles of molecules with exactly two B groups. A reacts only with B.", "step_growth"),
        number("p_groups", "Limiting-group conversion", .99, "fraction", 0, 1, "Reaction extent", "Fraction of the initially limiting A or B groups consumed. This is supplied, never predicted from time or temperature.", "step_growth"),
        number("mass_aa", "AA monomer molar mass", 146.14, "g/mol", 1e-6, 1e7, "Mass balance", "Input molecule molar mass, including reactive groups before reaction.", "step_growth"),
        number("mass_bb", "BB monomer molar mass", 116.21, "g/mol", 1e-6, 1e7, "Mass balance", "Input molecule molar mass. Unreacted monomer remains in the theoretical population.", "step_growth"),
        number("mass_loss", "Byproduct mass per A–B bond", 18.015, "g/mol of bonds", 0, 1e7, "Mass balance", "Use 0 for polyaddition; use the chemically correct loss per bond for condensation. Each bond forms one byproduct equivalent in this bookkeeping.", "step_growth"),
        number("n_monomer", "Initial monomer amount", .01, "mol", 1e-12, 1e6, "Chain-growth feed", "Initial moles of molecules contributing one specified residue on incorporation.", "living"),
        number("n_initiator", "Initial initiator amount", .0001, "mol", 1e-12, 1e6, "Chain-growth feed", "Initiated chains = initiator × efficiency × chains per initiator. Initiator chemistry is not inferred.", "living"),
        number("p_monomer", "Monomer conversion", .9, "fraction", 0, 1, "Chain-growth feed", "User-supplied or independently measured conversion. Temperature and time cannot generate this value here.", "living"),
        number("initiation_efficiency", "Initiation efficiency", 1, "fraction", 1e-9, 1, "Chain count", "Fraction of the specified initiator capacity producing chains; must be independently justified.", "living"),
        number("chains_per_initiator", "Chains per initiator molecule", 1, "chains/molecule", 1, 100, "Chain count", "Positive integer. Counts separate chains, not growing arms that remain joined to one multifunctional core.", "living"),
        number("repeat_mass", "Incorporated residue molar mass", 100, "g/mol", 1e-6, 1e7, "Chain mass", "Use residue mass after incorporation; for noncondensative ring opening it can equal monomer mass. Initial value is illustrative.", "living"),
        number("end_mass", "Combined end-group/core mass per chain", 30, "g/mol", 0, 1e7, "Chain mass", "Mass outside DP × residue mass, including initiator-derived and termination-derived atoms without double counting. Supply a chemically derived value.", "living"),
        number("f_A", "Instantaneous free-monomer fraction A", .5, "mole fraction", 1e-9, 1-1e-9, "Copolymer feed", "One fixed feed state. No cumulative composition or batch composition drift is calculated.", "copolymer"),
        number("r_A", "Reactivity ratio rA = kAA/kAB", .5, "1", 0, 1e6, "Terminal reactivity", "Supplied for the relevant mechanism and conditions; no named-monomer coefficients are assumed.", "copolymer"),
        number("r_B", "Reactivity ratio rB = kBB/kBA", .5, "1", 0, 1e6, "Terminal reactivity", "Stationary terminal model; penultimate effects and depropagation are omitted.", "copolymer"),
    ],
    "presets": [
        {"name": "Ideal nylon-6,6 mass balance", "description": "Adipic acid + hexamethylenediamine approximate input masses; p=0.99 is stipulated, not an experimental outcome.", "inputs": {"model": "step_growth", "n_aa_mol": .01, "n_bb_mol": .01, "p_groups": .99, "mass_aa": 146.14, "mass_bb": 116.21, "mass_loss": 18.015}},
        {"name": "Stoichiometric imbalance · finite ceiling", "description": "Synthetic bifunctional polyaddition. r=0.98 limits DPn to 99 even at full limiting-group conversion.", "inputs": {"model": "step_growth", "n_aa_mol": .0098, "n_bb_mol": .01, "p_groups": 1, "mass_aa": 1000, "mass_bb": 250, "mass_loss": 0}},
        {"name": "Ideal initiated-chain balance", "description": "Synthetic residue/end-group masses; no kinetic model and no SEC calibration.", "inputs": {"model": "living", "n_monomer": .01, "n_initiator": .0001, "p_monomer": .9, "initiation_efficiency": 1, "chains_per_initiator": 1, "repeat_mass": 100, "end_mass": 30}},
        {"name": "Bernoulli terminal limit", "description": "rA=rB=1 gives F_A=f_A and independent stationary additions at this fixed feed.", "inputs": {"model": "copolymer", "f_A": .3, "r_A": 1, "r_B": 1}},
        {"name": "Alternating terminal limit", "description": "rA=rB=0 prohibits AA and BB propagation; stationary interior-feed composition is 1:1.", "inputs": {"model": "copolymer", "f_A": .3, "r_A": 0, "r_B": 0}},
    ],
    "models": [
        {"name": "Carothers and molecule balance", "equation": "r=min(nAA,nBB)/max(nAA,nBB); bonds=2min(nAA,nBB)p; n=nAA+nBB−bonds; Xn=(1+r)/(1+r−2rp)", "assumptions": ["Exactly difunctional AA and BB; A reacts only with B", "Every bond is intermolecular; no cyclization, branching or side reactions", "Xn counts original monomer molecules, not AA+BB repeat pairs", "Whole retained population includes unreacted monomer"], "references": ["carothers"]},
        {"name": "Ideal initiated-chain mass balance", "equation": "nchains=η s nI; DPn=p nM/nchains; Mn=DPn Mresidue+Mends", "assumptions": ["Supplied conversion and initiation efficiency", "Chain count remains fixed after initiation; no transfer, termination or coupling", "All initiated chains counted, including zero-addition chains", "End/core mass is supplied with a consistent residue convention"], "references": ["iupac_chain"]},
        {"name": "Mayo–Lewis terminal model", "equation": "F_A=(rA f_A²+f_A f_B)/(rA f_A²+2f_A f_B+rB f_B²); PAA=rA f_A/(rA f_A+f_B)", "assumptions": ["Stationary terminal reactivity, irreversible growth and fixed free-monomer feed", "No penultimate effects or composition drift", "Dyads and run lengths describe the stationary infinite-chain limit, not a finite sampled oligomer"], "references": ["iupac_copolymer", "iupac_ratios"]},
    ],
    "references": REFERENCES,
}


def metric(key, label, value, unit="1", note=""):
    return {"key": key, "label": label, "value": value, "unit": unit, "note": note}


def plot(title, xlabel, ylabel, name, x, y):
    return {"title": title, "x_label": xlabel, "y_label": ylabel, "series": [{"name": name, "x": list(map(float, x)), "y": list(map(float, y))}]}


def calculate(inputs, material=None):
    if not isinstance(inputs, dict):
        raise ValueError("Calculation inputs must be an object.")
    p = {f["key"]: f["default"] for f in APP["fields"]}
    unknown = set(inputs) - set(p)
    if unknown:
        raise ValueError("Unknown inputs: " + ", ".join(sorted(map(str, unknown))))
    p.update(inputs)
    if p["model"] not in ("step_growth", "living", "copolymer"):
        raise ValueError("Unknown synthesis model.")
    for f in APP["fields"]:
        if f["type"] == "number":
            p[f["key"]] = finite(p[f["key"]], f["label"], f["min"], f["max"])
    notes = ["All results are conditional on supplied conversion and coefficients. Reaction time, solvent and temperature do not predict conversion here.", "A theoretical number-average mass is not a calibrated SEC/GPC measurement. Isolation or removal of residual monomer changes the measured population."]
    warnings = []
    if p["model"] == "step_growth":
        v = step_growth(p["n_aa_mol"], p["n_bb_mol"], p["p_groups"], p["mass_aa"], p["mass_bb"], p["mass_loss"])
        metrics = [metric("DPn", "Original-monomer units per molecule", v["DPn_original_monomer_units"]), metric("Mn", "Whole-population Mn", v["Mn_population_g_mol"], "g/mol"), metric("r", "Limiting/excess group ratio", v["r"]), metric("chain_mol", "Retained molecule amount", v["remaining_molecules_mol"], "mol")]
        rows = [[k, v[k], unit] for k, unit in [("bonds_mol", "mol"), ("remaining_A_mol", "mol"), ("remaining_B_mol", "mol"), ("conversion_A", "fraction"), ("conversion_B", "fraction"), ("initial_mass_g", "g"), ("byproduct_mass_g", "g"), ("retained_population_mass_g", "g")]]
        grid = np.linspace(0, .999, 200)
        r = v["r"]
        plots = [plot("Conversion and stoichiometric ceiling", "Limiting-group conversion", "Original-monomer Xn", "Current AA/BB feed", grid, (1+r)/(1+r-2*r*grid))]
        notes += ["Each condensation bond removes the specified byproduct mass. For polyaddition, set byproduct mass to zero.", "The Carothers Xn counts each AA or BB molecule as one unit. Do not multiply it directly by an AA+BB paired-repeat mass."]
        if v["balanced_Flory_length_dispersity"] is not None:
            metrics.append(metric("flory_length_D", "Balanced Flory length dispersity", v["balanced_Flory_length_dispersity"], "1", "Reference Xw/Xn=1+p; not generally the molar-mass dispersity when residue masses differ or end masses matter."))
        if v["DPn_original_monomer_units"] is None:
            warnings.append("Exact balance with p=1 is a singular ideal infinite-chain limit. Finite Xn and Mn are undefined; no finite acyclic population is claimed.")
        tables = [{"title": "Functional-group and mass balances", "columns": ["Quantity", "Value", "Unit"], "rows": rows}]
    elif p["model"] == "living":
        args = (p["n_monomer"], p["n_initiator"], p["p_monomer"], p["initiation_efficiency"], p["chains_per_initiator"], p["repeat_mass"], p["end_mass"])
        v = living_chain_growth(*args)
        metrics = [metric("DPn", "Residues per initiated chain", v["DPn"]), metric("Mn", "Conditional theoretical Mn", v["Mn_g_mol"], "g/mol"), metric("chain_mol", "Initiated chain amount", v["initiated_chains_mol"], "mol"), metric("end_fraction", "End/core mass fraction", v["end_mass_fraction"], "1")]
        grid = np.linspace(0, 1, 200)
        plots = [plot("Conversion sensitivity · fixed chain count", "Supplied conversion", "Theoretical Mn (g/mol)", "Fixed initiation efficiency", grid, grid*p["n_monomer"]/v["initiated_chains_mol"]*p["repeat_mass"]+p["end_mass"])]
        tables = [{"title": "Initiated-chain population balance", "columns": ["Quantity", "Value", "Unit"], "rows": [[k, v[k], "mol" if k.endswith("mol") else "g"] for k in ["initiated_chains_mol", "reacted_monomer_mol", "unreacted_monomer_mol", "chain_population_mass_g"]]}]
        notes += ["The initial residue and end-group masses are illustrative. Use a chemically derived finite-chain formula to avoid double counting terminal atoms.", "Chains per initiator means distinct molecules. A multiarm star attached to one core remains one molecule and requires a different molecular bookkeeping convention."]
        if v["DPn"] < 10:
            warnings.append("The average chain contains fewer than ten incorporated residues. End groups and the treatment of zero-addition chains can dominate; a long-chain polymer approximation needs scrutiny.")
    else:
        v = terminal_copolymer(p["f_A"], p["r_A"], p["r_B"])
        metrics = [metric("F_A", "Instantaneous copolymer fraction A", v["F_A"]), metric("F_B", "Instantaneous copolymer fraction B", v["F_B"]), metric("run_A", "Stationary mean A run length", v["mean_A_run"], "residues"), metric("run_B", "Stationary mean B run length", v["mean_B_run"], "residues")]
        grid = np.linspace(.001, .999, 200)
        plots = [plot("Instantaneous composition response", "Free-monomer mole fraction A", "Instantaneous incorporated fraction A", "Terminal model", grid, [terminal_copolymer(x, p["r_A"], p["r_B"])["F_A"] for x in grid])]
        tables = [{"title": "Stationary ordered dyad fractions", "columns": ["Dyad", "Fraction"], "rows": [[name, v["dyad_"+name]] for name in ["AA", "AB", "BA", "BB"]]}, {"title": "Chain-end transition probabilities", "columns": ["Transition", "Probability"], "rows": [[name, v["P_"+name]] for name in ["AA", "AB", "BA", "BB"]]}]
        notes += ["Dyads are oriented stationary adjacent-residue fractions; AB+BA is the combined unlike-dyad fraction.", "No finite sequence, batch conversion profile or cumulative copolymer composition is inferred. Reactivity ratios must be valid for the actual reaction mechanism and state."]
    return {"metrics": metrics, "plots": plots, "tables": tables, "notes": notes, "warnings": warnings, "handoffs": []}
