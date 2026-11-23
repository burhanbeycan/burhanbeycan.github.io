"""Workbench controls and evidence views for water transport in PVA."""
from importlib.resources import files

import numpy as np

from . import pva
from .model import finite_number, integer_number

FORWARD = "pva_forward"
HOLDOUT = "pva_holdout"
PAPER = "pva_published"
DATA = files("polymer_transport_lab").joinpath("data")
EQ_DEMO = DATA.joinpath("pva_eq_demo.csv").read_text(encoding="utf-8")
KIN_DEMO = DATA.joinpath("pva_kinetics_demo.csv").read_text(encoding="utf-8")


def add_controls(app, select, number):
    both = {"key": "mode", "value": [FORWARD, HOLDOUT]}
    forward = {"key": "mode", "value": FORWARD}
    holdout = {"key": "mode", "value": HOLDOUT}
    legacy = {"key": "mode", "value": ["release", "uptake"]}
    controls = app["fields"]
    controls[0]["options"].extend([{"value": key, "label": label} for key, label in [
        (FORWARD, "PVA film · independent isotherm + humidity step"),
        (HOLDOUT, "PVA film · compare models at a held-out thickness"),
        (PAPER, "PVA · inspect published equilibrium and apparent D"),
    ]])
    for field in controls[1:]:
        if "show_if" not in field:
            field["show_if"] = legacy

    def csv_field(key, label, default, help_text, visible):
        return dict(key=key, label=label, type="textarea", default=default, unit="CSV",
                    help=help_text, group="PVA water", show_if=visible)
    controls.extend([
        csv_field("pva_isotherm", "Independent equilibrium record", EQ_DEMO,
                  "Water activity 0–1 and g water per g dry polymer on one branch. material_id denotes a grade and processing history; include a traceable record_id.", both),
        select("pva_origin", "Origin of input records", "synthetic_demonstration",
               [("synthetic_demonstration","Synthetic demonstration"),("user_supplied_unverified","User-supplied records")],
               "Check the origin when replacing either example. The project records its input CSVs with each run.", "PVA water", show_if=both),
        number("pva_thickness_um", "Full, fixed film thickness", 90,"µm","Use measured film thickness. A backed film has one exposed face.", "PVA film",0.01,10000,show_if=forward),
        select("pva_faces","Exposed faces","one",[("one","One exposed face; impermeable backing"),("two","Two faces at identical humidity")],"Both faces see the same step when exposed.","PVA film",show_if=forward),
        number("pva_start","Initial water activity",0,"","Initial film is equilibrated at this activity.","Humidity step",0,1,show_if=forward),
        number("pva_end","New water activity",0.7,"","Surface instantaneously equilibrates on the specified isotherm branch.","Humidity step",0,1,show_if=forward),
        number("pva_dref","D at midpoint of input loading range",1.5e-12,"m²/s","The constant-D baseline and variable-D scenario use the same mid-loading D.","PVA film",1e-19,1e-7,show_if=forward),
        number("pva_contrast","ln[D(high loading)/D(low loading)]",2,"","Contrast spans the full input isotherm range. Zero gives constant diffusivity.","PVA film",-6,6,show_if=forward),
        number("pva_duration","Elapsed time",3600,"s","Time since humidity step.","Humidity step",0.001,1e9,show_if=forward),
        csv_field("pva_curves","Independent time-course records",KIN_DEMO,
                  "One curve_id per physical film humidity step, separate from equilibrium record IDs. Supply elapsed seconds and g/g dry polymer, full nominal thickness and 1 or 2 faces.",holdout),
        number("pva_holdout_um","Nominal held-out full thickness",90,"µm","All curves at this exact nominal thickness are excluded from both fits. Use the same nominal grouping across replicates.","Validation",0.01,10000,show_if=holdout),
    ])
    app["presets"].extend([
        {"name":"Published PVA / HDI · 25 °C","description":"El Shafee and Naguib, Tables 2 and 3. Published apparent diffusivities and fitted equilibrium isotherm.", "inputs":{"mode":PAPER}},
        {"name":"PVA / water · one exposed face","description":"Synthetic fixed-film humidity step with independent equilibrium input.", "inputs":{"mode":FORWARD}},
        {"name":"PVA / water · withhold 90 µm","description":"Synthetic 50 µm training set and 90 µm holdout, two humidity steps each. Fixture uses a finer grid.", "inputs":{"mode":HOLDOUT}},
    ])
    app["references"].extend([
        {"id":"jeck","title":"Jeck et al. (2011), crosslinked PVA water sorption kinetics and independent equilibrium","url":"https://doi.org/10.1016/j.cep.2010.09.004"},
        {"id":"shafee","title":"El Shafee & Naguib (2003), HDI-crosslinked PVA, Tables 2 and 3","url":"https://doi.org/10.1016/S0032-3861(02)00865-0"},
    ])
    app["models"].append({"name":"Humidity step in a fixed PVA film", "equation":"∂w/∂t = ∂x[D(w) ∂xw]; D(w) = Dref exp[B (w − wmid)/(wmax − wmin)]",
                          "assumptions":["w in g water per g dry polymer; uniform dry mass density and fixed geometry",
                                         "Measured branch supplies surface equilibrium by interpolation; instantaneous humidity step",
                                         "Reflected domain has length L for one face and L/2 for two exposed faces"],
                          "references":["jeck","shafee","bdf"]})


def run(values, metric, plot):
    mode = values["mode"]
    result = dict(metrics=[], tables=[], plots=[], warnings=[], notes=[], material_updates={"parameters":[]}, handoffs=[])
    if mode == PAPER:
        rows = pva.literature_rows()
        result["tables"].append(dict(title="Nominal 5 mol% HDI per vinyl alcohol repeat · 25 °C", columns=[
            "Vapour activity", "GAB-calculated g/g dry PVA network", "Reported apparent D × 10⁷ (cm²/s)", "Apparent D (m²/s)"], rows=rows))
        result["plots"].append(plot("Reported early-time apparent D", "Water activity", "Apparent D (m²/s)",
                                    [dict(name="El Shafee & Naguib · Table 3", x=[r[0] for r in rows], y=[r[3] for r in rows])]))
        result["metrics"] = [metric("D_ratio", "Apparent D ratio at aw 0.74 versus 0.38", rows[-1][3]/rows[0][3], "×")]
        result["notes"] = ["Table 3 reports early-time apparent coefficients for initially dry films from reduced sorption curves. Values depend on the humidity step.",
                           "Water loadings shown are evaluations of the article’s fitted GAB relation (Table 2, hm=0.062 g/g, f=0.83, c=1.70). They are derived equilibrium values.",
                           "This source provides no thickness-labelled raw time series in these tables. Its three apparent coefficients cannot be substituted for the local D(w)."]
        check = pva.paper_humidity_check()
        result["tables"].append(dict(title="One reserved activity at 25 °C · apparent D (m²/s)",
                                     columns=["Model", "Train aw", "Held-out aw", "Observed", "Prediction", "Absolute relative error (%)"],
                                     rows=[[label,"0.38, 0.74",check["aw"],check["observed"],check[key],
                                            100*abs(check[key]/check["observed"]-1)] for label,key in
                                           [("Constant D","constant"),("Exponential local D(w)","predicted")]]))
        result["notes"].append("For this one-activity interpolation, fit the reported apparent coefficients at aw 0.38 and 0.74, "
                               "predict aw 0.61 and compare to its reported coefficient. The exponential model matches each endpoint "
                               "through its predicted early-time gravimetric slope for a dry film. The source used an initial finite fitting window; "
                               "the numerical comparator uses its t→0 limit. This within-study, single-point check has no thickness transfer information.")
        return result

    iso = pva.read_isotherm(values["pva_isotherm"])
    origin = values["pva_origin"]
    if origin not in {"synthetic_demonstration", "user_supplied_unverified"}:
        raise ValueError("Choose a valid PVA data origin.")
    if origin != "synthetic_demonstration" and (values["pva_isotherm"].strip() == EQ_DEMO.strip() or
         (mode == HOLDOUT and values["pva_curves"].strip() == KIN_DEMO.strip())):
        raise ValueError("The supplied example contains synthetic records; replace all demo CSVs or retain the synthetic origin.")
    result["notes"] = [
        f"Input origin: {origin}. Condition: {iso.condition[0]}, {iso.condition[1]:g} K, {iso.condition[2]}. Equilibrium record(s): {', '.join(sorted(iso.record_ids))}.",
        "w means g water per g dry polymer. Surface water loading is interpolated within the independent equilibrium branch; the initial bulk loading is its initial-activity value.",
        "The model assumes a uniform dry mass per volume, negligible thickness change and negligible gas-side resistance. Use fixed-film sorption or desorption. Drying a freshly cast solution requires evolving geometry and interfacial mass transfer.",
    ]
    result["plots"].append(plot("Equilibrium branch", "Water activity", "g water / g dry polymer",
                                [dict(name="Input isotherm", x=iso.aw.tolist(), y=iso.water.tolist())]))
    if mode == FORWARD:
        face_choice = values["pva_faces"]
        if face_choice not in {"one","two"}:
            raise ValueError("Select one or two exposed faces.")
        times = np.r_[0, np.geomspace(1e-6, 1, 80) * finite_number(values["pva_duration"],"Elapsed time",0.001,1e9)]
        curve = pva.Curve("forward",values["pva_thickness_um"],1 if face_choice=="one" else 2,
                          values["pva_start"],values["pva_end"],times,np.array([]))
        constant = pva.loading_curve(iso, curve, values["pva_dref"], 0)
        variable = pva.loading_curve(iso, curve, values["pva_dref"], values["pva_contrast"], 64)
        wi, wf = iso.step(curve.start,curve.end)
        result["plots"].insert(0, plot("Fixed-film humidity step", "Elapsed time (s)", "g water / g dry polymer",
                                      [dict(name=name,x=times.tolist(),y=line.tolist()) for name,line in
                                       [("Constant D",constant),("Exponential D(w)",variable)]]))
        result["metrics"] = [metric("initial_w","Initial equilibrium water", wi, "g/g dry polymer"),
                             metric("final_eq_w","Boundary equilibrium water",wf,"g/g dry polymer"),
                             metric("last_w","Mean water at selected time",float(variable[-1]),"g/g dry polymer")]
        return result

    curves = pva.read_curves(values["pva_curves"], iso)
    models, train, test = pva.compare_models(iso, curves, values["pva_holdout_um"])
    result["notes"].append(f"Fit curves: {', '.join(c.key for c in train)}. Held-out curves: {', '.join(c.key for c in test)}. "
                           "Both models fit only training kinetic observations with equal total weight per curve. The independent isotherm is fixed across the split.")
    result["notes"].append("Fit bounds: Dref = 10⁻¹⁹ to 10⁻⁷ m²/s; −6 ≤ ln(Dmax-loading/Dmin-loading) ≤ 6. Calibration uses 32 finite-volume cells for the variable-D model.")
    names = ["Constant D", "Exponential D(w)"]
    for name, model in zip(names, models):
        if model["weak_sensitivity"] or model["near_bound"]:
            result["warnings"].append(f"{name}: low local sensitivity or fit near a search bound; parameters may be poorly resolved.")
    for c in test:
        result["plots"].insert(0,plot(f"Held-out film · {c.key} · {c.thickness:g} µm","Elapsed time (s)","g water / g dry polymer",
                                     [dict(name="Input observations", x=c.times.tolist(), y=c.observed.tolist())]+
                                     [dict(name=name,x=c.times.tolist(),y=model["prediction"][c.key]) for name,model in zip(names,models)]))
    result["tables"].append(dict(title="Equal-curve RMSE, g water per g dry polymer",columns=["Model","Dref (m²/s)","Log D contrast","Fit RMSE","Held-out RMSE"],
                                 rows=[[name,m["dref"],m["contrast"],m["train_rmse"],m["test_rmse"]] for name,m in zip(names,models)]))
    result["metrics"]=[metric("pva_holdout_const","Held-out RMSE · constant D",models[0]["test_rmse"],"g/g dry polymer"),
                       metric("pva_holdout_var","Held-out RMSE · exponential D(w)",models[1]["test_rmse"],"g/g dry polymer")]
    return result
