"""Workbench adapter; the mathematical model does not depend on this module."""

from __future__ import annotations

import numpy as np

from . import __version__
from .model import (
    MODEL_ID, characteristic_times, finite_number, finite_volume_slab,
    fit_diffusivity, fraction_transferred, integer_number, parse_fraction_data,
    steady_permeation,
)


REFERENCES = [
    {"id": "slab", "title": "PolyKin official documentation: plane-sheet fractional uptake; Crank (1975), p.48", "url": "https://hugomvale.github.io/polykin/reference/hmt/diffusion/uptake_constc_sheet/"},
    {"id": "bdf", "title": "SciPy solve_ivp: stiff initial-value integration and BDF", "url": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html"},
    {"id": "fit", "title": "SciPy least_squares: nonlinear least squares with bounds", "url": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html"},
    {"id": "permeation", "title": "Sandia National Laboratories: Hydrogen Compatibility of Materials — Polymers, equations 1–3", "url": "https://www.sandia.gov/app/uploads/sites/158/2021/12/8100TechRef_polymers.pdf"},
]


def _number(key, label, default, unit, help_text, group, minimum, maximum, **extra):
    return {"key": key, "label": label, "type": "number", "default": default,
            "unit": unit, "help": help_text, "group": group,
            "min": minimum, "max": maximum, **extra}


def _select(key, label, default, options, help_text, group, **extra):
    return {"key": key, "label": label, "type": "select", "default": default,
            "unit": "", "help": help_text, "group": group,
            "options": [{"value": value, "label": label} for value, label in options], **extra}


_FIT_VISIBLE = {"key": "fit_mode", "value": "fit_diffusivity"}
_PERMEATION_VISIBLE = {"key": "permeation_mode", "value": "henry"}
_SYNTHETIC_TIMES = np.array([0, 30, 60, 120, 240, 480, 900, 1800, 3600], dtype=float)
_SYNTHETIC_FRACTIONS = fraction_transferred(_SYNTHETIC_TIMES, 2e-12, 100e-6)
_SYNTHETIC_CSV = "time_s,fraction\n" + "\n".join(
    f"{t:g},{f:.12g}" for t, f in zip(_SYNTHETIC_TIMES, _SYNTHETIC_FRACTIONS)
)

APP = {
    "id": "transport", "name": "Polymer Transport Lab", "version": __version__,
    "tagline": "Resolve diffusion through a polymer film, with an independent mass-balance check.",
    "description": "Own analytical and finite-volume implementations for uniform, constant-D slab release and uptake. Fit D at known full thickness; inspect sensitivity and optional Henry-law permeation. Default scenarios are synthetic demonstrations.",
    "accent": "#0891b2",
    "scope": "One-dimensional fixed full slab; constant positive effective diffusivity; two faces at imposed concentration; no swelling kinetics, binding, degradation, finite bath or chemical inference of D.",
    "fields": [
        _select("mode", "Experiment", "release", [("release", "Release to two perfect sinks"), ("uptake", "Uptake from two constant-concentration faces")], "The normalized transferred fraction has the same mathematical form in both modes; concentration profiles differ.", "Model"),
        _number("thickness_um", "Full film thickness", 100, "µm", "Full distance between the two exposed faces, not half-thickness. For a network handoff use fully equilibrated wet thickness.", "Model", 0.001, 1e6),
        _number("diffusivity_m2_s", "Effective diffusivity D", 1e-12, "m²/s", "Explicit user-supplied value or demonstration input. No value is inferred from polymer identity, molecular weight or mesh size.", "Model", 1e-24, 1e-4),
        _number("temperature_K", "Temperature label", 298.15, "K", "Metadata for the supplied D. No Arrhenius correction is applied.", "Conditions", 1, 2000),
        _number("duration_s", "Observation duration", 3600, "s", "The displayed observation window; characteristic times may lie outside it.", "Numerics", 1e-9, 1e12),
        _number("n_points", "Output time points", 240, "", "Uniform output times including t=0; the BDF integrator uses its own adaptive steps.", "Numerics", 20, 1000, step=1),
        _number("n_cells", "Finite-volume cells", 80, "", "Independent numerical verification grid; finer cells resolve the initial surface layer better.", "Numerics", 8, 160, step=1),
        _select("fit_mode", "Fit an effective D", "off", [("off", "Forward calculation"), ("fit_diffusivity", "Fit user time/fraction data")], "Only D is fitted; full thickness, time origin and normalization are fixed.", "Fit"),
        {"key": "fit_data", "label": "Time / transferred-fraction data", "type": "textarea", "default": "", "unit": "s, fraction", "help": "CSV header time_s,fraction is optional. At least five rows; 0≤fraction≤1. Supply released fraction or normalized uptake, never fraction remaining.", "group": "Fit", "show_if": _FIT_VISIBLE},
        _select("fit_data_origin", "Data provenance", "user_supplied_unverified", [("user_supplied_unverified", "User supplied; unverified"), ("synthetic_demonstration", "Synthetic demonstration")], "This label does not certify an experiment or a material.", "Fit", show_if=_FIT_VISIBLE),
        _number("fit_lower_d", "Lower fit bound", 1e-16, "m²/s", "Positive lower bound on D, not a material database value.", "Fit", 1e-24, 1e-4, show_if=_FIT_VISIBLE),
        _number("fit_upper_d", "Upper fit bound", 1e-8, "m²/s", "Must exceed the lower bound. The objective scan reveals dependence on the selected bounds.", "Fit", 1e-24, 1e-4, show_if=_FIT_VISIBLE),
        _select("permeation_mode", "Separate steady permeation estimate", "off", [("off", "Disabled"), ("henry", "Henry solution-diffusion P = D S")], "A distinct steady gas experiment; not the two-face transient release boundary condition.", "Permeation"),
        _number("solubility_mol_m3_Pa", "Henry solubility S = c/p", 1e-5, "mol m⁻³ Pa⁻¹", "User supplied on an explicit molar, pressure-based basis. Requires c=S p and pressure≈fugacity; no concentration-independent claim for a real polymer.", "Permeation", 1e-20, 1e6, show_if=_PERMEATION_VISIBLE),
        _number("delta_pressure_Pa", "Upstream minus downstream pressure", 1e5, "Pa", "Signed pressure difference for the separate steady-state calculation.", "Permeation", -1e12, 1e12, show_if=_PERMEATION_VISIBLE),
    ],
    "presets": [
        {"name": "Synthetic film release", "description": "Synthetic 100 µm film with chosen D=1e-12 m²/s; no named polymer or measured calibration.", "inputs": {"mode": "release", "thickness_um": 100, "diffusivity_m2_s": 1e-12, "duration_s": 3600, "fit_mode": "off", "permeation_mode": "off"}},
        {"name": "Synthetic uptake", "description": "Same ideal slab under normalized uptake boundary conditions; demonstrates release/uptake complementarity.", "inputs": {"mode": "uptake", "thickness_um": 100, "diffusivity_m2_s": 1e-12, "duration_s": 3600, "fit_mode": "off", "permeation_mode": "off"}},
        {"name": "Synthetic D recovery", "description": "Noiseless data generated from this analytical model at D=2e-12 m²/s; a numerical recovery example, not external validation.", "inputs": {"mode": "release", "thickness_um": 100, "diffusivity_m2_s": 1e-12, "duration_s": 3600, "fit_mode": "fit_diffusivity", "fit_data": _SYNTHETIC_CSV, "fit_data_origin": "synthetic_demonstration", "fit_lower_d": 1e-15, "fit_upper_d": 1e-9, "permeation_mode": "off"}},
        {"name": "Synthetic gas permeation", "description": "Chosen D and Henry S demonstrate consistent SI permeability units; this is not a real membrane dataset.", "inputs": {"mode": "uptake", "thickness_um": 100, "diffusivity_m2_s": 1e-12, "duration_s": 3600, "fit_mode": "off", "permeation_mode": "henry", "solubility_mol_m3_Pa": 1e-5, "delta_pressure_Pa": 1e5}},
    ],
    "models": [
        {"name": "Two-face full-slab Fickian transport", "equation": "Fo = D t/L²; f = 1 - (8/π²) Σ[n≥0] exp[-(2n+1)²π²Fo]/(2n+1)²", "assumptions": ["L is full thickness; one-dimensional fixed geometry", "Uniform initial loading or concentration deficit", "Constant scalar D; both faces at imposed concentration", "No external resistance, binding, reaction, erosion, finite bath or dynamic swelling"], "references": ["slab"]},
        {"name": "Conservative finite-volume check", "equation": "dr_i/dFo = N²(r_(i-1)-2r_i+r_(i+1)); boundary diagonal -3N²; dq/dFo = 2N(r_0+r_(N-1))", "assumptions": ["Cell-centred equal volumes; boundary at half-cell distance", "Generic stiff BDF integration; cumulative boundary transfer is a separate state", "At Fo>8 the numerically negligible terminal profile is held"], "references": ["bdf"]},
        {"name": "Bounded effective-D fit", "equation": "min over log10(D) Σ_i [f(t_i;D,L)-f_i]²", "assumptions": ["Known fixed L, correct time zero and normalization", "Only D/L² is identifiable if thickness is also unknown", "Equal weights; any local standard error is conditional on independent equal-variance errors"], "references": ["fit", "slab"]},
        {"name": "Optional steady Henry permeation", "equation": "c = S p; P = D S; J_ss = P Δp/L", "assumptions": ["A separate steady experiment with different concentrations at the two faces", "Constant D and Henry S; ideal gas pressure approximates fugacity", "No defects, convection, sorption nonlinearity or surface resistance"], "references": ["permeation"]},
    ],
    "references": REFERENCES,
}

DEFAULTS = {field["key"]: field["default"] for field in APP["fields"]}


def _choice(inputs, key, allowed):
    value = inputs[key]
    if not isinstance(value, str) or value not in allowed:
        raise ValueError(f"{key} must be one of: {', '.join(allowed)}.")
    return value


def _metric(key, label, value, unit="", note=""):
    return {"key": key, "label": label, "value": value, "unit": unit, "note": note}


def _plot(title, x_label, y_label, series):
    return {"title": title, "x_label": x_label, "y_label": y_label,
            "x_scale": "linear", "y_scale": "linear", "series": series}


def calculate(inputs: dict, material: dict | None = None) -> dict:
    if not isinstance(inputs, dict):
        raise ValueError("Inputs must be a dictionary.")
    unknown = set(inputs) - {field["key"] for field in APP["fields"]}
    if unknown:
        raise ValueError("Unknown calculation input(s): " + ", ".join(sorted(map(str, unknown))))
    if material is not None and not isinstance(material, dict):
        raise ValueError("Material metadata must be a dictionary.")
    values = {**DEFAULTS, **inputs}
    # Validate every supplied numeric field, including currently hidden fields.
    numbers = {}
    for field in APP["fields"]:
        if field["type"] == "number":
            numbers[field["key"]] = finite_number(values[field["key"]], field["label"], field["min"], field["max"])
    mode = _choice(values, "mode", ["release", "uptake"])
    fit_mode = _choice(values, "fit_mode", ["off", "fit_diffusivity"])
    data_origin = _choice(values, "fit_data_origin", ["user_supplied_unverified", "synthetic_demonstration"])
    permeation_mode = _choice(values, "permeation_mode", ["off", "henry"])
    points = integer_number(values["n_points"], "Output time points", 20, 1000)
    cells = integer_number(values["n_cells"], "Finite-volume cells", 8, 160)
    if not isinstance(values["fit_data"], str):
        raise ValueError("Fit data must be text.")
    if numbers["fit_upper_d"] <= numbers["fit_lower_d"]:
        raise ValueError("Upper D fit bound must exceed lower bound.")
    length = numbers["thickness_um"] * 1e-6
    d = numbers["diffusivity_m2_s"]
    duration = numbers["duration_s"]
    notes = [
        "Full thickness L is the distance between the two exposed faces. Both faces are perfect sinks (release) or held at constant concentration (uptake).",
        "Default parameter values and presets are synthetic demonstrations, not measurements or a polymer property database.",
        "D is constant and user supplied. Temperature is recorded without a temperature correction; chemistry, molar-mass distribution and mesh size do not determine D here.",
        "The analytical expression and an independently discretized finite-volume/BDF solution are compared. Agreement verifies implementation under the stated model; it does not validate a real material.",
        "For a network handoff, thickness must be the fully equilibrated wet thickness. Transport then uses fixed geometry; swelling kinetics are not coupled.",
    ]
    if material:
        notes.append("Imported chemistry, condition and population metadata are retained by the workbench but are not used to infer transport coefficients. An imported record is not experimental verification.")
    warnings = []
    metrics = []
    tables = []
    plots = []
    d_origin = "user_input_or_demonstration"
    if fit_mode == "fit_diffusivity":
        fit_t, fit_f = parse_fraction_data(values["fit_data"])
        fitted = fit_diffusivity(fit_t, fit_f, length, numbers["fit_lower_d"], numbers["fit_upper_d"])
        d = fitted.diffusivity_m2_s
        d_origin = f"fitted_to_{data_origin}"
        warnings.extend(fitted.warnings)
        notes.extend([
            "The forward solution uses fitted D. Thickness, time origin and normalization were fixed; fractional curves alone determine D/L², not D and L separately.",
            "Fit RMSE is an equal-weight residual statistic. The objective scan is not a confidence interval. Local log10(D) standard error, when shown, assumes independent equal-variance errors and ignores thickness, boundary-condition and model uncertainty.",
            "A relative thickness error propagates approximately twice as strongly to fitted D (D ∝ L²).",
        ])
        if data_origin == "synthetic_demonstration":
            notes.append("Fit data are explicitly labelled synthetic. Recovering a generating parameter is numerical verification only.")
        if fitted.log10_d_standard_error is None:
            notes.append("Local log10(D) standard error is undefined or uninformative for this fit; no confidence claim is made.")
        metrics.extend([
            _metric("fit_rmse", "Fit RMSE", fitted.rmse, "fraction", "Equal-weight residuals; not predictive validation."),
            _metric("fit_log10_d_se", "Local log10(D) standard error", fitted.log10_d_standard_error, "decades", "Conditional on fixed L and independent equal-variance errors; null when unresolved."),
            _metric("fit_sensitivity", "Local D sensitivity", fitted.sensitivity_rms, "fraction / decade", "RMS derivative of predicted fractions with respect to log10(D)."),
        ])
        plots.extend([
            _plot("Fit observations and predictions", "Time (s)", "Transferred fraction", [
                {"name": f"Data: {data_origin}", "x": fit_t.tolist(), "y": fit_f.tolist()},
                {"name": "Fitted slab prediction", "x": fit_t.tolist(), "y": fitted.predicted.tolist()},
            ]),
            _plot("Bounded D objective scan", "log10(D / (m²/s))", "Sum of squared fraction residuals", [
                {"name": "Objective; not a confidence interval", "x": fitted.scan_log10_d.tolist(), "y": fitted.scan_sse.tolist()},
            ]),
        ])
        tables.append({"title": "Fit residuals (prediction − observation)", "columns": ["time_s", "observed_fraction", "predicted_fraction", "residual"], "rows": [[float(t), float(y), float(p), float(r)] for t, y, p, r in zip(fit_t, fit_f, fitted.predicted, fitted.residuals)]})
    times = np.linspace(0, duration, points)
    analytic = fraction_transferred(times, d, length)
    numerical = finite_volume_slab(times, d, length, cells)
    scales = characteristic_times(d, length)
    difference = float(np.max(np.abs(analytic - numerical.fraction_transferred)))
    fo_end = d * duration / length**2
    metrics[0:0] = [
        _metric("diffusivity_m2_s", "D used", d, "m²/s", d_origin),
        _metric("final_fraction", "Final transferred fraction", float(analytic[-1]), "fraction", "Released mass / initial mass, or uptake / equilibrium uptake."),
        _metric("first_mode_time_s", "First-mode relaxation time", scales["first_mode_time_s"], "s", "L²/(π²D), using full thickness L; not the time to complete release."),
        _metric("t50_s", "Time to 50% transfer", scales["t50_s"], "s"),
        _metric("t90_s", "Time to 90% transfer", scales["t90_s"], "s"),
        _metric("fourier_number_final", "Final Fourier number", float(fo_end), "D t/L²"),
        _metric("analytic_numeric_max_error", "Analytical–numerical maximum difference", difference, "fraction", f"At requested output times; {cells} finite-volume cells."),
        _metric("mass_balance_max_error", "Maximum mass-balance residual", numerical.mass_balance_max_error, "fraction", "Raw mean remaining concentration + independently integrated two-face flux − 1."),
    ]
    plots.insert(0, _plot("Two-face release" if mode == "release" else "Two-face uptake", "Time (s)", "Transferred fraction", [
        {"name": "Analytical slab", "x": times.tolist(), "y": analytic.tolist()},
        {"name": f"Finite volume ({cells} cells)", "x": times.tolist(), "y": numerical.fraction_transferred.tolist()},
        {"name": "Integrated boundary transfer", "x": times.tolist(), "y": numerical.integrated_boundary_transfer.tolist()},
    ]))
    profile_series = []
    indices = sorted(set([0, max(1, points // 20), points // 4, points // 2, points - 1]))
    x_um = np.r_[0.0, numerical.cell_centres_m * 1e6, length * 1e6]
    for idx in indices:
        remaining = numerical.remaining_profiles[idx]
        concentration = remaining if mode == "release" else 1 - remaining
        boundary = (1.0 if mode == "release" else 0.0) if idx == 0 else (0.0 if mode == "release" else 1.0)
        profile_series.append({"name": f"t = {times[idx]:.5g} s", "x": x_um.tolist(), "y": np.r_[boundary, concentration, boundary].tolist()})
    plots.append(_plot("Finite-volume concentration profiles", "Position across full thickness (µm)", "c/c₀" if mode == "release" else "c/c∞", profile_series))
    notes.append("The t=0 profile shows the uniform initial condition. The imposed face concentrations jump at t=0+; an arbitrarily thin initial boundary layer cannot be resolved by a finite grid.")
    if difference > 0.005:
        warnings.append("Finite-volume discrepancy exceeds 0.005 fraction at the sampled times. Refine the grid and inspect the early boundary layer; the analytical solution remains the reference within this model.")
    if numerical.time_tail_held_at_equilibrium:
        notes.append("For Fo>8 the finite-volume terminal profile is held to bound computation. The exact continuum remaining fraction is below exp(−8π²), approximately 5.2e−35; solver roundoff is larger than this tail.")
    if duration / (points - 1) > scales["t50_s"]:
        warnings.append("The output time spacing exceeds t50. Shorten the observation window or increase output points to display the initial kinetics.")
    if analytic[-1] < 0.05:
        warnings.append("Less than 5% transfer occurs in this window under the model. Extend duration to inspect the later kinetics.")
    tables.append({"title": "Characteristic times (full thickness)", "columns": ["Quantity", "Time (s)", "Meaning"], "rows": [
        ["L²/D", scales["diffusion_scale_s"], "Diffusion scaling time"],
        ["L²/(π²D)", scales["first_mode_time_s"], "First-mode e-folding time"],
        ["t50", scales["t50_s"], "50% transferred"],
        ["t90", scales["t90_s"], "90% transferred"],
        ["t99", scales["t99_s"], "99% transferred"],
    ]})
    parameters = [
        {"quantity": "effective_solute_diffusivity", "value": d, "unit": "m2/s", "origin": d_origin, "model_scope": MODEL_ID, "temperature_K": numbers["temperature_K"]},
        {"quantity": "full_transport_slab_thickness", "value": length, "unit": "m", "origin": "user_input_or_equilibrium_network_handoff", "model_scope": MODEL_ID},
    ]
    if permeation_mode == "henry":
        steady = steady_permeation(d, length, numbers["solubility_mol_m3_Pa"], numbers["delta_pressure_Pa"])
        metrics.extend([
            _metric("permeability", "Steady permeability P = D S", steady["permeability_mol_m_s_Pa"], "mol m⁻¹ s⁻¹ Pa⁻¹", "Explicit molar Henry basis; separate steady gas experiment."),
            _metric("steady_flux", "Steady molar flux", steady["steady_flux_mol_m2_s"], "mol m⁻² s⁻¹", "J = P Δp/L; positive from upstream to downstream."),
        ])
        notes.append("P=D S assumes constant D and Henry solubility c=S p, local interfacial equilibrium, and pressure≈fugacity. This separate steady calculation does not change the transient two-face boundary conditions. No nonlinear sorption is fitted.")
        parameters.extend([
            {"quantity": "henry_solubility_molar_pressure_basis", "value": numbers["solubility_mol_m3_Pa"], "unit": "mol/m3/Pa", "origin": "user_input_or_demonstration", "model_scope": "constant-D-Henry-steady-permeation/1.0"},
            {"quantity": "permeability_molar_pressure_basis", "value": steady["permeability_mol_m_s_Pa"], "unit": "mol/m/s/Pa", "origin": "derived_from_explicit_D_and_S", "model_scope": "constant-D-Henry-steady-permeation/1.0"},
        ])
    return {"metrics": metrics, "plots": plots, "tables": tables, "notes": notes,
            "warnings": warnings, "material_updates": {"parameters": parameters}, "handoffs": []}
