"""Workbench views for the WBPU formulation, characterization and comparison tools."""

from .wbpu import (
    DEFAULT_RECIPE, DMA_DEMO, EXTRACTION_DEMO, PAPER_DMA, PAPER_DOI, PAPER_SERIES, PAPER_SOURCE, dma_sweeps,
    extraction, recipe_balance, swelling_folds,
)


def _choice(key, title, options, default, group, help_text, mode=None):
    field = {"key": key, "label": title, "type": "select", "options": [
        {"value": value, "label": label} for value, label in options], "default": default,
        "group": group, "help": help_text}
    if mode:
        field["show_if"] = {"key": "network_mode", "value": mode}
    return field


def _text(key, title, default, mode, group, help_text):
    return {"key": key, "label": title, "default": default, "type": "textarea",
            "show_if": {"key": "network_mode", "value": mode}, "group": group, "help": help_text}


def extend_app(app):
    for field in app["fields"]:
        if "show_if" not in field:
            field["show_if"] = {"key": "network_mode", "value": "equilibrium"}
    # Existing conditional fields depend on another equilibrium field; the runtime
    # resolves their parent visibility. Keep their local condition intact.
    for preset in app["presets"]:
        preset["inputs"]["network_mode"] = "equilibrium"
    app["fields"][:0] = [
        _choice("network_mode", "Analysis", [
            ("equilibrium", "Neutral network equilibrium"), ("recipe", "WBPU functional-group ledger"),
            ("evidence", "Film extraction and DMA"), ("holdout", "Withheld-formulation comparison")],
            "equilibrium", "Analysis", "Use a common swelling protocol to compare experimental formulations."),
        _text("recipe_csv", "Formulation components (CSV)", DEFAULT_RECIPE, "recipe", "Stoichiometry",
              "Choose molar_mass (g/mol plus functionality), OH_number (mg KOH/g), NCO_wt_pct, or equivalent_weight (g/eq). For reacted water only, choose net_urea and leave value and functionality blank."),
        _text("recipe_source", "Recipe source and assay basis", "Illustrative PEG/PCL/BDO feed using nominal dry-component molar masses.",
              "recipe", "Stoichiometry", "Record the lot, solids basis, measured OH number, NCO assay and amine functionality as available."),
        _text("extraction_csv", "Extraction mass record (CSV)", EXTRACTION_DEMO, "evidence", "Extraction",
              "Use neat unfilled films, dry to constant mass and record solvent, duration and temperature. The blotted sample may weigh less than its initial mass after dissolution."),
        _text("dma_csv", "Optional tensile DMA sweep (CSV)", DMA_DEMO, "evidence", "Mechanical evidence",
              "Paste E′ at each temperature and frequency. Leave this field empty for extraction alone."),
        _text("evidence_context", "Sample history, controls and analytical evidence", "Synthetic mass and DMA example; solvent quality and thermal transitions have not been established.",
              "evidence", "Mechanical evidence", "Record drying, matrix solubility control, DSC/WAXS transitions, thermal cycling, FTIR or post-cure chemistry. Use identical specimen IDs to associate observations."),
        _text("series_csv", "Matched formulation series (CSV)", PAPER_SERIES, "holdout", "Predictive comparison",
              "Hard segment is recipe wt%. Target is net mass increase per 100 g initial dry film. Use a single protocol including identical hard segment definition. Formulation IDs group replicates."),
        _text("series_source", "Source, table or dataset identifier", PAPER_SOURCE,
              "holdout", "Predictive comparison", "Identify original measurements, their uncertainty and any extraction or digitization."),
        _text("dma_protocol", "Optional independent modulus protocol", PAPER_DMA, "holdout", "Predictive comparison",
              "For an augmented model enter pre-immersion E′ in MPa in every row and identify the frequency, temperature, strain regime and conditioning. Use separate dry-film measurements."),
    ]
    app["presets"].extend([
        {"name": "PEG 20k / PCL 1k WBPU feed", "description": "A dry-component example with equal PEG/PCL masses and a BDO chain-extension stage.",
         "inputs": {"network_mode": "recipe", "recipe_csv": DEFAULT_RECIPE}},
        {"name": "Extracted film and DMA sweep", "description": "Synthetic masses on three weighing bases with a separate tensile DMA trace.",
         "inputs": {"network_mode": "evidence", "extraction_csv": EXTRACTION_DEMO, "dma_csv": DMA_DEMO}},
        {"name": "Jatropha WBPU / toluene benchmark", "description": "Three experimental WBPU formulations from Saalah et al. (2021), tested in turn against the other two.",
         "inputs": {"network_mode": "holdout", "series_csv": PAPER_SERIES, "dma_protocol": PAPER_DMA,
                    "series_source": PAPER_SOURCE}},
    ])
    app["references"].extend([
        {"id": "polymeg", "title": "POLYMEG formulation calculations — hydroxyl and isocyanate equivalent weights",
         "url": "https://www.dollmar.com/wp-content/uploads/2022/09/Brochure-Polymeg.pdf"},
        {"id": "saalah2021", "title": "Saalah et al. — Waterborne polyurethane dispersion derived from jatropha oil",
         "url": "https://doi.org/" + PAPER_DOI, "year": 2021},
        {"id": "segmented2019", "title": "Tunable Structure and Properties of Segmented Thermoplastic Polyurethanes as a Function of Flexible Segment",
         "url": "https://doi.org/10.3390/polym11121910", "year": 2019},
    ])
    app["models"].extend([
        {"name": "Reactive feed ledger", "equation": "n_group = m f / Mn; index = 100 nNCO / (nOH + namine + 2 nreacted water)",
         "assumptions": ["One ordinary primary or secondary amine site consumes one NCO to form a urea. Bulk dispersion water sits outside the ledger."], "references": ["polymeg"]},
        {"name": "Extraction mass bases", "equation": "gel = md/m0; sol = 1 − gel; net swell = (mw−m0)/m0; gel uptake = (mw−md)/md",
         "assumptions": ["m0 is initially dry film, mw is blotted wet residue and md is the same residue redried. Evidence depends on extraction and thermal history."],
         "references": ["saalah2021", "segmented2019"]},
        {"name": "Matched swelling series", "equation": "S = a + b w_hard [ + c ln(E′ / 1 MPa) ]",
         "assumptions": ["Equal-weight formulation means, least-squares fit inside each fold; metrics use percentage points."],
         "references": ["saalah2021"]},
    ])


def _metric(key, label, value, unit, note):
    return {"key": key, "label": label, "value": value, "unit": unit, "note": note}


def _table(title, columns, rows):
    return {"title": title, "columns": columns, "rows": rows}


def calculate_wbpu(v):
    result = {"metrics": [], "tables": [], "plots": [], "material_updates": {}, "handoffs": [], "warnings": [], "notes": []}
    if v["network_mode"] == "recipe":
        data = recipe_balance(v["recipe_csv"])
        totals = data["totals"]
        metrics = result["metrics"]
        for key in ("OH", "NCO", "amine"):
            metrics.append(_metric(key, f"{key} groups", totals[key], "mmol", "Across all included feed stages."))
        if totals["OH"]:
            metrics.append(_metric("nco_oh", "Feed NCO/OH", totals["NCO"] / totals["OH"], "mol/mol", "Only hydroxyl groups enter the denominator."))
        metrics.extend([
            _metric("nco_index", "Overall reactive index", data["index"], "", "100 × NCO equivalents / included NCO demand."),
            _metric("excess", "Net NCO equivalents", data["balance"], "mmol", "Positive implies nominal excess in the feed.")])
        result["tables"] = [
            _table("Component equivalents", ["Component", "Stage", "Group", "Mass (g)", "Basis", "Reactive equivalents (mmol)"], data["entries"]),
            _table("Equivalents added by stage", ["Stage", "NCO (mmol)", "OH (mmol)", "Amine (mmol)", "Water demand (mmol NCO)"],
                   [[s, g["NCO"], g["OH"], g["amine"], g["reacted_water"]] for s, g in data["stages"].items()])]
        result["notes"] = [v["recipe_source"], "The ledger counts feed groups. Stage totals record additions. Residual NCO or conversion requires a measurement.",
                           "For commercial polyols, lot-specific hydroxyl number or equivalent weight accounts for deviations from ideal functionality. Use dry active mass with a matched assay basis.",
                           "Primary and secondary amine sites contribute one ordinary urea-forming reaction each. Reacted water, when independently quantified, consumes two NCO equivalents in the net urea pathway. Biuret, allophanate, isocyanurate and other side reactions are outside this balance.",
                           "A feed of diols and diisocyanates can yield linear segmented PU. PEG/PCL crystallites, ionic clusters and hard domains can support transient junctions. Independent structural evidence is required to specify elastically active strands."]
    elif v["network_mode"] == "evidence":
        rows = extraction(v["extraction_csv"])
        series = dma_sweeps(v["dma_csv"])
        result["tables"] = [_table("Operational extraction on neat, unfilled films",
            ["Sample", "Solvent", "°C", "h", "Insoluble (gel) %", "Extractable (sol) %", "Net swell / initial dry %", "Solvent / retained dry gel %"], rows)]
        if series:
            result["plots"] = [{"title": "Frequency-resolved dry-film tensile DMA", "x_label": "Temperature (°C)",
                                "y_label": "Storage modulus E′ (MPa)", "x_scale": "linear", "y_scale": "log", "series": series}]
        result["notes"] = [v["evidence_context"], "A retained dry residue establishes operational insolubility under the reported solvent, time and temperature. Check matrix solubility and completion of extraction; insoluble fillers require a separate mass balance.",
                           "A finite-frequency E′ plateau can include hard-domain association, crystallites and entanglements. Compare thermal cycling, DSC or WAXS and solvent controls when assigning permanent junctions. Conversion of the DMA trace to a covalent strand density is unsupported.",
                           "A zero retained mass gives an undefined gel-normalized uptake. Initial-mass net swelling combines uptake and extractable mass loss."]
    else:
        if not v["series_source"].strip():
            raise ValueError("Identify the experimental data source, including a specimen record or citation.")
        if v["series_csv"].strip() != PAPER_SERIES.strip() and v["series_source"].strip() == PAPER_SOURCE:
            raise ValueError("Update the dataset citation or transformation note when modifying the published series.")
        data = swelling_folds(v["series_csv"], v["dma_protocol"])
        labels = data["models"]
        for index, (label, mae) in enumerate(zip(labels, data["mae"])):
            result["metrics"].append(_metric(f"mae_{index}", f"{label} MAE", mae, "pp", "Equal weight per withheld formulation."))
        folds = data["folds"]
        result["tables"] = [_table("One formulation withheld per fold",
            ["Formulation", "Observation rows", "Hard segment wt%", "Observed swell %", *[s + " prediction %" for s in labels], "Predictor range"],
            [[row["sample"], row["replicates_or_rows"], row["hard"], row["observed"], *row["predictions"],
              "outside marginal ranges" if row["outside_training_range"] else "inside marginal ranges"] for row in folds])]
        sorted_folds = sorted(folds, key=lambda r: r["hard"])
        series = [{"name": "Measured", "x": [r["hard"] for r in sorted_folds], "y": [r["observed"] for r in sorted_folds]}]
        for i, label in enumerate(labels):
            series.append({"name": label, "x": series[0]["x"], "y": [r["predictions"][i] for r in sorted_folds]})
        result["plots"] = [{"title": "Measured and withheld-formulation mass swelling", "x_label": "Recipe hard segment (wt%)",
                            "y_label": "Net swell (% initial dry mass)", "x_scale": "linear", "y_scale": "linear", "series": series}]
        family, solvent, temp, hours, protocol = data["conditions"]
        result["notes"] = [v["series_source"], f"{family}; {solvent}, {temp:g} °C, {hours:g} h; {protocol}.",
                           "Replicate rows are averaged within formulations before fitting each fold. Endpoints of a formulation series can require extrapolation. All folds share the supplied specimen family and exposure protocol.",
                           "The regression is a within-series empirical comparison. Other recipe components can covary with hard-segment fraction. Generalization to another laboratory, formulation family or solvent requires independent data."]
        if data["augmented"] and len(folds) == 3:
            result["notes"].append("Three formulations permit separate one-variable regressions. A joint recipe + DMA fit would have only two training rows for three coefficients and is omitted.")
        if data["augmented"]:
            result["notes"].append("Independent DMA protocol: " + v["dma_protocol"])
        if v["series_csv"].strip() == PAPER_SERIES.strip():
            result["notes"].append("The reported targets are three formulation means from a single jatropha study. Tables 3 and 4 use JPU 188; Table 1 prints JPU 181. The benchmark joins only Tables 3 and 4. The study includes TEA, DMPA, HEMA and IPDI in its hard-segment mass definition. The source reports triplicate extractions; published swell uncertainties in table order are ±7.6, ±17.9 and ±14.9 percentage points. Toluene gel fractions are 51.2%, 55.4% and 74.2%.")
    return result
