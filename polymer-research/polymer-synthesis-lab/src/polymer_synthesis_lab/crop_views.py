"""Input forms and result tables for recipe predictions and experimental checks."""
from .crop import ENDS, Recipe, audit, load_cohort

BATCH = load_cohort("batch_azide")
FLOW = load_cohort("flow_feed")
REFERENCES = [
    {"id": "crop_rates", "title": "Vergaelen et al. Arrhenius coefficients in MeCN, Table 2 (2023).",
     "url": "https://doi.org/10.1021/acs.macromol.2c01930"},
    {"id": "crop_origin", "title": "Vergaelen et al. Original kinetic experiments (2015).",
     "url": "https://doi.org/10.1021/acsmacrolett.5b00392"},
    {"id": "crop_branch", "title": "Arráez et al. Chain transfer and SEC constraints on branching (2019).",
     "url": "https://doi.org/10.1021/acs.macromol.9b00544"},
    {"id": "crop_batch", "title": "PEtOx–N3 batch experiment, SI Fig. S2D (2019).",
     "url": "https://doi.org/10.1039/C9CC00407F"},
    {"id": "crop_flow", "title": "EtOx continuous-flow feed series, SI Table S1 (2015).",
     "url": "https://doi.org/10.1039/C5CC04319K"},
]


def select(key, label, default, options, model, help_text):
    return {"key": key, "label": label, "default": default, "type": "select", "group": "PEtOx study",
            "unit": "", "help": help_text, "options": [{"value": k, "label": v} for k,v in options],
            "show_if": {"key": "model", "value": model}}


def add_models(app, number):
    predict, benchmark = "petox_recipe", "petox_benchmark"
    app["fields"][0]["options"].extend([
        {"value": predict, "label": "PEtOx · CROP from feed, temperature and time"},
        {"value": benchmark, "label": "PEtOx · measured Mn and group holdouts"}])
    app["fields"].extend([
        number("crop_m", "Starting EtOx concentration", 3, "mol/L", 1e-6, 12, "EtOx feed", "Total initial reaction volume basis.", predict),
        number("crop_i", "Starting MeOTs concentration", .03, "mol/L", 1e-8, 12, "EtOx feed", "Total initial reaction volume basis.", predict),
        number("crop_T", "Isothermal reaction temperature", 100, "°C", 0, 200, "Kinetics", "Use actual material temperature during hold.", predict),
        number("crop_minutes", "Isothermal reaction time", 30, "min", 0, 100000, "Kinetics", "Exclude ramp and cooling, or assess their contribution separately.", predict),
        select("crop_ends", "Quenched terminal groups", "Me/OH", [(v,v) for v in ENDS], predict, "Combined terminal residue masses assume quantitative capping."),
        select("crop_split", "Withhold complete groups", "paper", [
            ("paper","Publication"), ("series","Reaction series"), ("concentration","Initial monomer concentration"),
            ("feed","Monomer/initiator feed ratio"), ("concentration_and_feed","Both new concentration and new feed")], benchmark,
            "In the joint split, only training rows sharing neither the held-out concentration nor the held-out feed are eligible."),
        {"key": "crop_csv", "label": "Experimental PEtOx observations (CSV)", "type": "textarea",
         "default": BATCH, "group": "Experimental dataset", "unit": "", "help": "One reactor, assay and terminal chemistry per cohort. Use the original experiment as paper ID. X_observed may be empty.",
         "show_if": {"key": "model", "value": benchmark}},
    ])
    app["presets"].extend([
        {"name": "PEtOx · predict from an isothermal recipe", "description": "MeOTs in MeCN, 3 M EtOx, 100:1 feed. Illustrative unknown sample.",
         "inputs": {"model": predict, "crop_m": 3, "crop_i": .03, "crop_T": 100, "crop_minutes": 30, "crop_ends": "Me/OH"}},
        {"name": "PEtOx · batch SEC time course", "description": "Four crude azide-terminated PEtOx samples, 82 °C, PMMA-calibrated DMAc SEC, SI Fig. S2D. Nominal concentrations.",
         "inputs": {"model": benchmark, "crop_csv": BATCH, "crop_split": "paper"}},
        {"name": "PEtOx · four withheld flow feeds", "description": "DP targets 15–120 at 4 M, 160 °C and 5 min residence; one PMMA-SEC assay, SI Table S1.",
         "inputs": {"model": benchmark, "crop_csv": FLOW, "crop_split": "feed"}},
    ])
    app["models"].append({
        "name": "EtOx propagation and a conditional Mn model",
        "equation": "kp=A exp(−Ea/RT); X=1−exp(−kp[I]₀t); DPn=([M]₀/[I]₀)X; Mn=Mends+99.133 DPn",
        "assumptions": ["EtOx + MeOTs in MeCN, isothermal hold, rapid quantitative initiation",
                        "Constant active site population and volume, no transfer, branching or termination during propagation",
                        "One nominal terminal chemistry. SEC calibration and sample preparation can shift apparent Mn."],
        "references": [r["id"] for r in REFERENCES]})
    app["references"].extend(REFERENCES)


def recipe_view(p, metric, plot):
    r = Recipe(p["crop_m"], p["crop_i"], p["crop_T"], p["crop_minutes"], p["crop_ends"])
    x, mn, full = r.predict()
    times = [k * r.minutes / 99 for k in range(100)] if r.minutes else [0]
    curves = [Recipe(r.M, r.I, r.T, t, r.ends).predict()[1] for t in times]
    warnings = []
    if abs(r.M-3) > 1e-6 or abs(r.feed-100) > .001:
        warnings.append("Transfer away from the published kinetic reference feed (3 M, 100:1).")
    if not 60 <= r.T <= 140:
        warnings.append("Temperature outside the 60–140 °C range tabulated in the 2023 paper's Table 1.")
    if x * r.feed < 10:
        warnings.append("DPn <10: low-DP populations, actual end groups and isolation strongly affect Mn.")
    return {
        "metrics": [metric("X", "Predicted conversion", x, "fraction"), metric("Mn", "Kinetic theoretical Mn", mn, "g/mol"),
                    metric("Mn_full", "Full-conversion reference Mn", full, "g/mol"), metric("DPn", "Predicted mean residue count", x*r.feed)],
        "plots": [plot("Isothermal propagation", "Time (min)", "Theoretical Mn (g/mol)", "Published kinetics", times, curves)],
        "tables": [{"title":"Kinetic and structural inputs","columns":["Quantity","Value","Unit"],
                    "rows":[["Monomer/initiator",r.feed,"mol/mol"],["Apparent kp",r.kp,"L/(mol s)"],["kp[I]₀",r.kp*r.I,"1/s"],
                            ["Average EtOx residue",99.133,"g/mol"],[r.ends,ENDS[r.ends],"g/mol of chains"]]}],
        "notes": ["Rates from Vergaelen et al., Table 2. MeCN kinetics in that paper reuse the 2015 study.",
                  "Theoretical Mn assumes no chain transfer and quantitative initiation. SEC-apparent Mn also depends on calibration and purification.",
                  "At t=0 the mass term denotes the formal capped zero-addition limit, not an isolated polymer."],
        "warnings": warnings, "handoffs": []}


def benchmark_view(p, metric):
    raw = p["crop_csv"]
    rows, folds, results, summary = audit(raw, p["crop_split"])
    methods = {"kinetic": "Fixed published kinetics", "full": "Full incorporation", "scaled": "Training-only assay scale × kinetic Mn"}
    notes = [
        "Each held-out prediction uses feed, temperature, time and terminal groups. Measured conversion appears only for comparison.",
        "The third model learns one multiplicative Mn assay correction from training-series median ratios. It does not identify initiation efficiency or transfer.",
        "Errors are sample-weighted across tested rows of this cohort. P90 is the empirical linearly interpolated 90th percentile. The project target is median absolute relative error below 10%.",
        "Each publication ID must identify the original experiments; deduplicate reprinted data. A split by feed alone may share a laboratory and concentration with training."]
    warnings = []
    if len({r.paper for r in rows}) < 2:
        warnings.append("Single experimental paper. Cross-paper transfer of a trained model cannot be established here.")
    if len({round(r.recipe.M, 6) for r in rows}) < 2:
        warnings.append("Single initial monomer concentration. Cross-concentration Mn performance is unresolved.")
    if any(f[1] == 0 for f in folds):
        warnings.append("One or more folds have no eligible training samples. Assay-scaled predictions are omitted there; fixed baselines remain computable.")
    if any(not 60 <= r.recipe.T <= 140 for r in rows):
        warnings.append("At least one temperature is outside the 60–140 °C interval in the kinetic source's Table 1.")
    if any(r.paper == "10.1039/c9cc00407f" for r in rows):
        notes.append("SI Fig. S2D batch cohort: nominal I0=(24.6 mg/186.23 g mol⁻¹)/(1.3+1.9 mL), M0=113 I0. Assumes additive volumes and negligible MeOTs volume; exact solution concentration was not reported.")
    if any(r.paper == "10.1039/c5cc04319k" for r in rows):
        notes.append("SI Table S1 flow cohort: 5 min residence at 20 bar. Me/OH is a nominal end-group assumption after KOH/MeOH quench. Terminal identity was not established per sample. A single isothermal residence time ignores axial dispersion.")
    plots = []
    # Show data on a physical axis when the cohort varies along one coordinate.
    selected = None
    if len({(r.recipe.M,r.recipe.I,r.recipe.T) for r in rows}) == 1:
        selected = ("Isothermal time (min)", lambda r: r.recipe.minutes)
    elif len({(r.recipe.M,r.recipe.T,r.recipe.minutes) for r in rows}) == 1:
        selected = ("Initial EtOx / MeOTs (mol/mol)", lambda r: r.recipe.feed)
    if selected and len(rows) > 1:
        label, axis = selected
        lookup = {(r.paper,r.series,r.sample): axis(r) for r in rows}
        points = sorted(results, key=lambda v: lookup[(v["paper"],v["series"],v["sample"])])
        curves = []
        for key, name in [("Mn","Reported Mn"),("kinetic","Fixed kinetic prediction"),("scaled","Training-only assay correction")]:
            available = [r for r in points if r[key] is not None]
            if available:
                curves.append({"name":name, "x":[lookup[(r["paper"],r["series"],r["sample"])] for r in available],
                               "y":[r[key] for r in available]})
        plots.append({"title":"Cohort observations and recipe predictions", "x_label":label,
                      "y_label":"Mn (g/mol)", "series":curves})
    n, med, tail = summary["kinetic"]
    return {
        "metrics": [metric("samples","Measured Mn records",n,"records"), metric("median_ape","Fixed kinetics · median absolute relative error",med,"%"),
                    metric("p90_ape","Fixed kinetics · P90 absolute relative error",tail,"%"),metric("folds","Held-out groups",len(folds),"groups")],
        "tables": [
            {"title":"Measured cohort","columns":["Original paper","Reactor","Mn assay","Termini","Concentration or feed basis"],"rows":[[paper,rows[0].reactor,rows[0].assay,rows[0].recipe.ends,"Initial molar concentrations"] for paper in sorted({r.paper for r in rows})]},
            {"title":"Model comparison","columns":["Model","Scored rows","Median absolute error (%)","P90 absolute error (%)"],"rows":[[v,*summary[k]] for k,v in methods.items()]},
            {"title":"Holdout composition","columns":["Excluded group","Training rows","Test rows","Training assay factor"],"rows":folds},
            {"title":"Individual predictions","columns":["Paper","Series","Sample","Predicted X","Observed X","Measured Mn","Kinetic Mn","Full-incorporation Mn","Assay-scaled Mn","Kinetic error (%)"],
             "rows":[[r[k] for k in ("paper","series","sample","X","X_observed","Mn","kinetic","full","scaled","ape_kinetic")] for r in results]}],
        "plots": plots, "notes": notes, "warnings": warnings, "handoffs": []}
