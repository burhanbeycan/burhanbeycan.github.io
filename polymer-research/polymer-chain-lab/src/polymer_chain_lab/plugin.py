"""Workbench adapter for the independently usable scientific kernels."""

from __future__ import annotations

import math

import numpy as np

from .md_workflow import md_output

from .viscosity import CSV_FIELDS, benchmark, literature_series, mixture, read_observations
from .model import (degrees_of_polymerization, finite_number, gaussian_dimensions, integer,
                    make_distribution, parse_distribution, population_form_factor,
                    sample_ris, wormlike_dimensions)

POPULATION_MODELS = {"key": "model", "value": ["gaussian", "wlc"]}
MD_MODEL = {"key": "model", "value": "petox_md"}
FIT_MODEL = {"key": "model", "value": "intrinsic_viscosity_fit"}
MIX_MODEL = {"key": "model", "value": "intrinsic_viscosity_mwd"}
MWD_MODELS = {"key": "model", "value": ["gaussian", "wlc", "intrinsic_viscosity_mwd"]}
RIS_MODEL = {"key": "model", "value": "ris"}


def _number(key, label, default, unit, minimum, maximum, help_text, group, show_if=None, step=None):
    field = {"key": key, "label": label, "type": "number", "default": default,
             "unit": unit, "min": minimum, "max": maximum, "help": help_text, "group": group}
    if show_if is not None:
        field["show_if"] = show_if
    if step is not None:
        field["step"] = step
    return field


REFERENCES = [
    {"id": "galata2024", "title":"Galata and Kröger, aqueous linear/cyclic PEtOx simulation", "url":"https://doi.org/10.1021/acs.macromol.4c00174"},
    {"id": "gubarev2018", "title": "Gubarev et al., PEtOx hydrodynamics in PBS", "url": "https://doi.org/10.1039/C8PY00255J", "year": 2018},
    {"id": "gubarev_si", "title": "Experimental SI, Tables S1 and S2 (SEC–MALS and intrinsic viscosity)", "url": "https://www.rsc.org/suppdata/c8/py/c8py00255j/c8py00255j1.pdf", "year": 2018},
    {"id": "iupac_mw_viscosity", "title": "IUPAC Gold Book: viscosity-average molar mass", "url": "https://goldbook.iupac.org/terms/view/12220"},
    {"id": "iupac_dp", "title": "IUPAC Gold Book: degree of polymerization", "url": "https://goldbook.iupac.org/terms/view/D01569"},
    {"id": "iupac_mn", "title": "IUPAC Gold Book: number-average molar mass", "url": "https://goldbook.iupac.org/terms/view/12216", "year": 2025},
    {"id": "iupac_dispersity", "title": "IUPAC Gold Book: molar-mass dispersity", "url": "https://goldbook.iupac.org/terms/view/12224", "year": 2025},
    {"id": "sasview_debye", "title": "SasView: monodisperse Gaussian coil and Debye form factor", "url": "https://www.sasview.org/docs/user/models/mono_gauss_coil.html"},
    {"id": "wlc_dimensions", "title": "AFRC documentation: worm-like chain and Benoit–Doty gyration radius", "url": "https://afrc.readthedocs.io/en/stable/polymer_models/worm_like_chain_obrien.html"},
    {"id": "iupac_ris", "title": "IUPAC Gold Book: rotational isomeric state", "url": "https://goldbook.iupac.org/terms/view/12187", "year": 2025},
]

DEFAULT_CSV = "molar_mass_g_mol,amount\n20000,0.5\n80000,0.5"
TAIL_CSV = "molar_mass_g_mol,amount\n10000,0.09375\n50000,0.8875\n250000,0.01875"

APP = {
    "id": "chain", "name": "Polymer Chain Lab", "version": "0.4.0", "accent": "#287f80",
    "tagline": "Resolve populations, chain stiffness, and local conformational statistics.",
    "description": "Analyze user-supplied molar-mass distributions; compare Gaussian and worm-like chain dimensions; sample an explicit finite-chain rotational-isomeric-state model; evaluate PEtOx dilute-solution viscometry by held-out molar-mass class.",
    "scope": "Equilibrium single-chain statistics and population weighting. PEtOx viscometry benchmarks use measured literature data. Conformational presets and the component-law MWD example use illustrative parameters.",
    "fields": [
        {"key": "model", "label": "Physical model", "type": "select", "default": "gaussian", "unit": "", "help": "Population models use MWD; RIS describes independent chains with an explicit local torsion Hamiltonian.", "group": "Model", "options": [
            {"value":"petox_md","label":"PEtOx · end-group MD trajectories"},
            {"value": "gaussian", "label": "Gaussian chain + MWD + scattering"},
            {"value": "wlc", "label": "Worm-like chain + MWD dimensions"},
            {"value": "ris", "label": "Finite-chain RIS ensemble"},
            {"value": "intrinsic_viscosity_fit", "label": "PEtOx · measured intrinsic-viscosity benchmark"},
            {"value": "intrinsic_viscosity_mwd", "label": "Intrinsic viscosity from a complete MWD"}]},
        _number("pilot_dp","Degree of polymerization",4,"repeats",1,10000,"The bundled finite-cluster trajectories cover DP 2 and 4.","End-group MD",MD_MODEL,1),
        {"key":"trajectory_csv","label":"Analyzed trajectory CSV (optional)","type":"textarea","default":"","group":"End-group MD","show_if":MD_MODEL,
         "help":"Leave empty for the executed MMFF94 cluster pilot, or paste the output from polymer-petox-md gro. One sample series and matched conditions per comparison."},
        _number("temperature_K", "Temperature", 298.15, "K", 1, 5000, "RIS uses temperature in Boltzmann weights. Gaussian/WLC stiffness is user-supplied at this temperature; no temperature law is inferred.", "Conditions", {"key": "model", "value": ["gaussian", "wlc", "ris"]}),
        {"key": "mwd_source", "label": "Population source", "type": "select", "default": "csv", "unit": "", "group": "Population", "help": "Use the table below or the number-fraction MWD in an imported material record.", "show_if": MWD_MODELS, "options": [{"value": "csv", "label": "Entered CSV"}, {"value": "material", "label": "Imported material MWD"}]},
        {"key": "distribution_basis", "label": "CSV amount basis", "type": "select", "default": "number", "unit": "", "group": "Population", "help": "Discrete fractions/counts are normalized. Log-density means density per d(log10 M), integrated over the supplied range by trapezoidal quadrature. This field does not alter imported material number fractions.", "show_if": MWD_MODELS, "options": [
            {"value": "number", "label": "Discrete number amounts"},
            {"value": "mass", "label": "Discrete mass amounts"},
            {"value": "number_log_density", "label": "Number density per log10 M"},
            {"value": "mass_log_density", "label": "Mass density per log10 M"}]},
        {"key": "distribution_csv", "label": "Molar mass, amount", "type": "textarea", "default": DEFAULT_CSV, "unit": "g/mol, relative amount", "group": "Population", "help": "Two comma-separated columns; optional molar_mass_g_mol,amount header. Up to 512 rows. Masses must be positive; amounts nonnegative. Raw SEC retention times are not molar masses.", "show_if": MWD_MODELS},
        _number("viscosity_temp_C", "Viscometry temperature", 37, "°C", -273.14, 600, "The published series covers 15, 25, 37, 45 and 55 °C. Each fit uses one measured temperature.", "Dilute-solution conditions", FIT_MODEL),
        {"key": "calibration_csv", "label": "Custom PEtOx measurement CSV (optional)", "type": "textarea", "default": "", "unit": "", "group": "Dilute-solution conditions", "show_if": FIT_MODEL, "help": "Leave blank to use Gubarev et al. Experimental linear PEtOx and SEC–MALS Mw only. Keep the same chemistry, architecture and medium within a fit. Header: " + ",".join(CSV_FIELDS)},
        _number("target_mw_g_mol", "Target sample Mw", 50000, "g/mol", 1, 1e10, "Interpolate for a sample of comparable dispersity and end chemistry measured under the same conditions.", "Dilute-solution conditions", FIT_MODEL),
        _number("component_k", "Component-law K", 0.04, "(mL/g)/(g/mol)^a", 1e-12, 1e5, "Enter an independently justified law for individual mass components. The default is illustrative.", "Component viscosity law", MIX_MODEL),
        _number("component_a", "Component-law exponent a", 0.65, "", 0.001, 2, "All populations must share the same linear chemistry, solvent and temperature.", "Component viscosity law", MIX_MODEL),
        _number("component_low_g_mol", "Lower calibrated component mass", 10000, "g/mol", 1, 1e12, "Smallest mass for which the component law applies.", "Component viscosity law", MIX_MODEL),
        _number("component_high_g_mol", "Upper calibrated component mass", 150000, "g/mol", 1, 1e12, "All populated mass bins must lie within the calibration interval.", "Component viscosity law", MIX_MODEL),
        {"key": "component_medium", "label": "Polymer, solvent and temperature", "type": "text", "default": "Illustrative linear homopolymer in solvent S at 25 °C", "unit": "", "group": "Component viscosity law", "show_if": MIX_MODEL, "help": "Identify the actual polymer, solvent composition and temperature used for the constants."},
        {"key": "component_source", "label": "Source of component law", "type": "text", "default": "Synthetic example", "unit": "", "group": "Component viscosity law", "show_if": MIX_MODEL, "help": "Cite a calibration record and mass range for measured constants."},
        _number("tail_threshold_g_mol", "High-mass tail threshold", 100000, "g/mol", 1e-6, 1e12, "Report fractions strictly above this molar mass. Useful when Mn and dispersity hide different tails.", "Population", POPULATION_MODELS),
        _number("repeat_mass_g_mol", "Repeat-unit molar mass", 100, "g/mol", 0.001, 1e6, "Common repeat contribution M0. Effective DP=(M−Mend)/M0; do not substitute the molar mass of a capped monomer. The default is synthetic.", "Chain mapping", POPULATION_MODELS),
        _number("end_group_mass_g_mol", "Net end-group mass correction", 0, "g/mol per chain", 0, 1e6, "Mend in M=DP×M0+Mend: terminal atoms not already counted in the repeat contribution. Zero explicitly neglects this correction. Supply the value appropriate to your chain chemistry.", "Chain mapping", POPULATION_MODELS),
        _number("contour_increment_nm", "Contour increment per repeat", 0.25, "nm/repeat", 1e-6, 100, "L=[(M−Mend)/M0]×increment. This mapping must come from a justified chain geometry; terminal geometric contributions are not included.", "Chain mapping", POPULATION_MODELS),
        _number("kuhn_length_nm", "Kuhn length", 1.0, "nm", 1e-8, 1e8, "Long-chain Gaussian parameter b; valid for many Kuhn segments under ideal-chain conditions.", "Chain mapping", {"key": "model", "value": "gaussian"}),
        _number("persistence_length_nm", "Persistence length", 5.0, "nm", 1e-8, 1e8, "WLC tangent-correlation length lp. The long-chain Kuhn length is 2lp; do not confuse these quantities.", "Chain mapping", {"key": "model", "value": "wlc"}),
        _number("q_min_nm_inverse", "Minimum scattering vector", 0.001, "nm⁻¹", 1e-8, 1e4, "Positive lower edge of the logarithmic scattering grid.", "Scattering", {"key": "model", "value": "gaussian"}),
        _number("q_max_nm_inverse", "Maximum scattering vector", 0.2, "nm⁻¹", 1e-8, 1e4, "Continuum Gaussian scattering loses microscopic validity when qb approaches unity.", "Scattering", {"key": "model", "value": "gaussian"}),
        _number("n_points", "Scattering grid points", 180, "", 16, 512, "Bounded logarithmic grid; applies only to Gaussian scattering.", "Scattering", {"key": "model", "value": "gaussian"}, 1),
        _number("n_bonds", "Skeletal bonds", 80, "", 3, 500, "An open chain has n+1 equal-mass sites and n−2 torsions. No molar mass is assigned to these sites.", "RIS geometry", RIS_MODEL, 1),
        _number("bond_length_nm", "Fixed bond length", 0.154, "nm", 1e-6, 100, "Synthetic geometry; no polymer identity is inferred.", "RIS geometry", RIS_MODEL),
        _number("bond_angle_deg", "Internal bond angle", 112, "degrees", 1, 179, "Angle between consecutive bonds at their shared vertex; fixed for every site.", "RIS geometry", RIS_MODEL),
        _number("gauche_energy_kj_mol", "Gauche state energy above trans", 1.5, "kJ/mol", -50, 50, "E(trans)=0; E(gauche+)=E(gauche−)=this value. Values are synthetic, not assigned to a named polymer.", "RIS Hamiltonian", RIS_MODEL),
        _number("opposite_gauche_penalty_kj_mol", "Adjacent opposite-gauche penalty", 4, "kJ/mol", -50, 100, "Extra energy for consecutive gauche+/gauche− or gauche−/gauche+; all other pair energies are zero.", "RIS Hamiltonian", RIS_MODEL),
        _number("n_chains", "Independent chain samples", 200, "", 2, 2000, "Exact finite-chain torsion sampling. Total sites across all samples must not exceed 200,000.", "RIS sampling", RIS_MODEL, 1),
        _number("seed", "Random seed", 1729, "", 0, 2**32 - 1, "Reproducible with the same code and NumPy environment; no Markov-chain burn-in is required.", "RIS sampling", RIS_MODEL, 1),
    ],
    "presets": [
        {"name":"Executed MD · PEtOx termini with 12 waters","description":"Short atomistic MMFF94 microsolvation trajectories for two chemical ends; inspect three independent packing seeds per end.","inputs":{"model":"petox_md","pilot_dp":4}},
        {"name": "Measured PEtOx · PBS at 37 °C", "description": "Seven SEC–MALS mass classes with independent dilute-solution viscometry. Inspect withheld-class errors and the separate SEC conditions.", "inputs": {"model": "intrinsic_viscosity_fit", "viscosity_temp_C": 37, "calibration_csv": "", "target_mw_g_mol": 50000}},
        {"name": "Synthetic MWD · viscosity-average molar mass", "description": "Apply a specified single-component law to a two-bin polymer mass distribution.", "inputs": {"model": "intrinsic_viscosity_mwd", "mwd_source": "csv", "distribution_basis": "number", "distribution_csv": DEFAULT_CSV}},
        {"name": "Synthetic A · two mass populations", "description": "Mn=50,000 g/mol, Mw=68,000 g/mol, Đ=1.36; compare with Synthetic B to expose high-tail differences.", "inputs": {"model": "gaussian", "mwd_source": "csv", "distribution_basis": "number", "distribution_csv": DEFAULT_CSV, "repeat_mass_g_mol": 100, "contour_increment_nm": 0.25, "kuhn_length_nm": 1, "tail_threshold_g_mol": 100000}},
        {"name": "Synthetic B · same Mn and Đ, different tail", "description": "Matches Synthetic A in Mn and Mw but changes Mz and the scattering-weighted radius. This is an exact discrete demonstration.", "inputs": {"model": "gaussian", "mwd_source": "csv", "distribution_basis": "number", "distribution_csv": TAIL_CSV, "repeat_mass_g_mol": 100, "contour_increment_nm": 0.25, "kuhn_length_nm": 1, "tail_threshold_g_mol": 100000}},
        {"name": "Synthetic WLC · rod-to-coil span", "description": "A discrete log-spaced population explores finite contour length relative to persistence length; no WLC scattering curve is substituted.", "inputs": {"model": "wlc", "mwd_source": "csv", "distribution_basis": "number", "distribution_csv": "molar_mass_g_mol,amount\n1000,1\n10000,1\n100000,1", "repeat_mass_g_mol": 100, "contour_increment_nm": 0.25, "persistence_length_nm": 25}},
        {"name": "Synthetic RIS · local steric coupling", "description": "Finite-chain trans/gauche model with an explicit adjacent opposite-gauche penalty. Parameters are illustrative and have no assigned chemistry.", "inputs": {"model": "ris", "temperature_K": 298.15, "n_bonds": 80, "bond_length_nm": 0.154, "bond_angle_deg": 112, "gauche_energy_kj_mol": 1.5, "opposite_gauche_penalty_kj_mol": 4, "n_chains": 200, "seed": 1729}},
    ],
    "models": [
        {"name":"Finite-cluster MD and aqueous trajectory analysis","equation":"Rg² = Σmi |ri−rCOM|²/Σmi; n = Σ I(rsite,water ≤ cutoff)/Nsites",
         "assumptions":["Matched chain length and α end", "Bundled simulations are isolated clusters; PBC only for user-supplied orthorhombic liquid trajectories", "Replica is the statistical unit"],"references":["galata2024"]},
        {"name": "Intrinsic viscosity versus sample Mw", "equation": "ln[η] = ln K + a ln Mw", "assumptions": ["Matched linear PEtOx, SEC–MALS protocol, viscometry medium and temperature", "Similar within-series dispersities", "Testing with complete nominal mass classes held out"], "references": ["gubarev2018", "gubarev_si"]},
        {"name": "Dilute-solution MWD viscosity average", "equation": "[η]mix = Σwi K Mi^a; Mv = (Σwi Mi^a)^(1/a)", "assumptions": ["Noninteracting chains at infinite dilution", "Mass fractions; uniform chemistry, architecture and conditions", "Component-level constants supplied independently"], "references": ["iupac_mw_viscosity"]},
        {"name": "Chemically explicit mass-to-contour mapping", "equation": "DPi=(Mi−Mend)/M0; Li=DPi a", "assumptions": ["One common repeat contribution and terminal correction across the population", "Open linear chains; fractional effective DP for continuum bins", "Terminal geometry is not included in contour length", "No sequence-resolved copolymer mapping is inferred"], "references": ["iupac_dp"]},
        {"name": "Molar-mass population", "equation": "Mn=ΣnM; Mw=ΣnM²/ΣnM; Mz=ΣnM³/ΣnM²; Đ=Mw/Mn", "assumptions": ["Normalized number fractions", "Each row represents a discrete bin or a documented quadrature node", "Scattering weights nM² require equal contrast per mass"], "references": ["iupac_mn", "iupac_dispersity"]},
        {"name": "Gaussian dimensions and Debye scattering", "equation": "<Re²>=bL; <Rg²>=bL/6; P(q)=2[exp(−x)−1+x]/x², x=q²<Rg²>", "assumptions": ["Ideal linear chains with many Kuhn segments", "Unperturbed dimensions; no excluded-volume expansion", "Population scattering assumes independent chains and uniform contrast per mass"], "references": ["sasview_debye"]},
        {"name": "Worm-like dimensions", "equation": "<Re²>=2lp²[L/lp−1+exp(−L/lp)]; <Rg²>=Llp/3−lp²+2lp³/L−2lp⁴[1−exp(−L/lp)]/L²", "assumptions": ["Continuous inextensible open chain", "Uniform local persistence length", "No nonlocal excluded volume", "No WLC form factor in this release"], "references": ["wlc_dimensions"]},
        {"name": "Finite-chain RIS", "equation": "Pr(s1,…,sm)=exp{−[ΣEi(si)+ΣJ(si,si+1)]/(RT)}/Z", "assumptions": ["Three torsions: trans 180°, gauche+ 60°, gauche− −60°", "Uniform bond length, bond angle and torsion energies", "Exact finite-chain transfer messages; independent samples", "Equal site masses; no nonlocal interactions"], "references": ["iupac_ris"]},
    ],
    "references": REFERENCES,
}


def _validated_inputs(inputs):
    if not isinstance(inputs, dict):
        raise ValueError("Calculation inputs must be an object.")
    unknown = set(inputs) - {field["key"] for field in APP["fields"]}
    if unknown:
        raise ValueError("Unknown calculation input(s): " + ", ".join(sorted(map(str, unknown))))
    values = {field["key"]: field["default"] for field in APP["fields"]}
    values.update(inputs)
    integers = {"pilot_dp","n_points", "n_bonds", "n_chains", "seed"}
    for field in APP["fields"]:
        key, value = field["key"], values[field["key"]]
        if field["type"] == "number":
            validator = integer if key in integers else finite_number
            values[key] = validator(value, field["label"], field["min"], field["max"])
        elif field["type"] == "select":
            if value not in [option["value"] for option in field["options"]]:
                raise ValueError(f"Unknown option for {field['label']}.")
        elif not isinstance(value, str):
            raise ValueError(f"{field['label']} must be text.")
    return values


def _metric(key, label, value, unit="", note=""):
    return {"key": key, "label": label, "value": None if value is None else float(value), "unit": unit, "note": note}


def _parameter(quantity, value, unit, model_scope):
    return {"quantity": quantity, "value": float(value), "unit": unit, "origin": "user_input_or_demonstration", "model_scope": model_scope}


def _plot(title, x_label, y_label, series, x_scale="linear", y_scale="linear"):
    return {"title": title, "x_label": x_label, "y_label": y_label, "series": series, "x_scale": x_scale, "y_scale": y_scale}


def _series(name, x, y):
    return {"name": name, "x": np.asarray(x).tolist(), "y": np.asarray(y).tolist()}


def _material_notes(material):
    notes = ["Chemistry names, sequence, tacticity, architecture, end groups, solvent, and imported parameter metadata do not infer model constants. Supply parameters justified for the selected physical model."]
    if material:
        notes.append("An imported material record is context, not evidence of experimental validation. Only its explicit MWD is used when Population source is Imported material MWD; all model parameters and temperature come from the inputs shown here.")
    return notes


def _distribution(values, material):
    if values["mwd_source"] == "material":
        record = material.get("molar_mass_distribution")
        if not isinstance(record, dict) or record.get("basis") != "number_fraction":
            raise ValueError("The imported material needs an MWD with basis='number_fraction'.")
        distribution = make_distribution(record.get("molar_mass_g_mol", []), record.get("number_fraction", []), "number")
    else:
        distribution = parse_distribution(values["distribution_csv"], values["distribution_basis"])
    return distribution


def _population(values, material):
    distribution = _distribution(values, material)
    mass, number = distribution.molar_mass_g_mol, distribution.number_fraction
    mass_fraction, scattering = distribution.mass_fraction, distribution.scattering_fraction
    moments = distribution.moments()
    dp = degrees_of_polymerization(mass, values["repeat_mass_g_mol"], values["end_group_mass_g_mol"])
    contour = dp * values["contour_increment_nm"]
    gaussian = values["model"] == "gaussian"
    if gaussian:
        re2, rg2 = gaussian_dimensions(contour, values["kuhn_length_nm"])
        segment_ratio = contour / values["kuhn_length_nm"]
    else:
        re2, rg2 = wormlike_dimensions(contour, values["persistence_length_nm"])
        segment_ratio = contour / values["persistence_length_nm"]
    metrics = [_metric(key, label, moments[key], unit) for key, label, unit in [
        ("Mn", "Number-average molar mass Mn", "g/mol"), ("Mw", "Mass-average molar mass Mw", "g/mol"),
        ("Mz", "z-average molar mass Mz", "g/mol"), ("dispersity", "Molar-mass dispersity Đ", "")]]
    metrics.extend([
        _metric("dp_number", "Number-average effective DP", float(number @ dp), "repeats/chain", "Σn(M−Mend)/M0 for the common supplied chemical mapping"),
        _metric("rg_number_nm", "Number-weighted RMS Rg", math.sqrt(float(number @ rg2)), "nm", "sqrt(Σn<Rg²>), not the mean of Rg"),
        _metric("rg_mass_nm", "Mass-weighted RMS Rg", math.sqrt(float(mass_fraction @ rg2)), "nm", "sqrt(Σw<Rg²>)"),
        _metric("rg_scattering_nm", "Zero-q intensity-weighted RMS Rg", math.sqrt(float(scattering @ rg2)), "nm", "sqrt(Σs<Rg²>); s∝nM² at uniform contrast per mass"),
        _metric("re_number_nm", "Number-weighted RMS end-to-end distance", math.sqrt(float(number @ re2)), "nm"),
    ])
    tail = mass > values["tail_threshold_g_mol"]
    for key, label, weights in [("tail_number", "Number fraction above threshold", number), ("tail_mass", "Mass fraction above threshold", mass_fraction), ("tail_scattering", "Zero-q intensity fraction above threshold", scattering)]:
        metrics.append(_metric(key, label, float(weights[tail].sum()), "fraction", f"M > {values['tail_threshold_g_mol']:g} g/mol"))
    plots = [
        _plot("Population weights by molar mass", "Molar mass (g/mol)", "Component/bin fraction", [
            _series("Number fraction", mass, number), _series("Mass fraction", mass, mass_fraction), _series("Zero-q intensity fraction", mass, scattering)], "log"),
        _plot("Single-component dimensions", "Molar mass (g/mol)", "RMS distance (nm)", [
            _series("Radius of gyration", mass, np.sqrt(rg2)), _series("End-to-end distance", mass, np.sqrt(re2))], "log", "log"),
    ]
    notes = _material_notes(material) + [
        "The contour mapping DP=(M−Mend)/M0 and L=DP×contour increment uses a common repeat contribution and net terminal mass correction. It is not a composition-resolved copolymer mapping and does not add terminal geometric length.",
        "Number, mass, and scattering-weighted RMS sizes average squared dimensions before taking a square root. They are different observables.",
        "The nM² scattering weights require identical contrast per unit mass, dilute independent chains, and a common composition. No absolute intensity or interparticle structure factor is calculated.",
        "Population plots show discrete component/bin fractions. Connected points are visual guides, not a recovered continuous SEC distribution.",
        "Temperature does not automatically change Gaussian/WLC stiffness: supply a Kuhn or persistence length applicable at the entered temperature.",
    ]
    warnings = []
    if distribution.input_basis.endswith("log_density"):
        notes.append("Log-density input is integrated by trapezoidal node weights in log10(M), then converted to number fractions. Only the entered range is included; refine the grid and check missing tails for convergence.")
    if np.min(dp) < 10:
        warnings.append("Some components contain fewer than ten effective repeats; end-group masses and discrete geometry can materially affect the continuum mapping.")
    if values["end_group_mass_g_mol"] == 0:
        notes.append("The net terminal mass correction is set to zero. This is an explicit approximation unless the selected repeat convention already includes the terminal contribution.")
    parameters = [_parameter("repeat_unit_molar_mass", values["repeat_mass_g_mol"], "g/mol", "effective_homopolymer_contour_mapping"),
                  _parameter("net_end_group_molar_mass_correction", values["end_group_mass_g_mol"], "g/mol", "effective_homopolymer_contour_mapping"),
                  _parameter("contour_increment_per_repeat", values["contour_increment_nm"], "nm", "effective_homopolymer_contour_mapping")]
    if gaussian:
        if values["q_max_nm_inverse"] <= values["q_min_nm_inverse"]:
            raise ValueError("Maximum scattering vector must exceed the minimum.")
        q = np.geomspace(values["q_min_nm_inverse"], values["q_max_nm_inverse"], values["n_points"])
        intensity = population_form_factor(q, rg2, distribution)
        plots.append(_plot("Gaussian population form factor", "q (nm⁻¹)", "I(q) / I(0)", [_series("M²-weighted Debye mixture", q, intensity)], "log", "log"))
        if np.any(segment_ratio < 10):
            warnings.append("At least one component has fewer than ten Kuhn segments. Gaussian dimensions and Debye scattering are long-chain approximations; inspect the WLC model where appropriate.")
        if values["q_max_nm_inverse"] * values["kuhn_length_nm"] >= 0.3:
            warnings.append("The scattering range reaches qb ≥ 0.3. The Debye curve is a continuum prediction; local geometry becomes important as qb approaches unity. This 0.3 flag is a conservative diagnostic, not a universal boundary.")
        parameters.append(_parameter("kuhn_length", values["kuhn_length_nm"], "nm", "ideal_gaussian_chain"))
    else:
        notes.append("WLC dimensions retain finite contour length and the rigid-rod limit. A WLC scattering model is not implemented; no Gaussian Debye form factor is assigned to this mode.")
        parameters.append(_parameter("persistence_length", values["persistence_length_nm"], "nm", "continuous_worm_like_chain"))
    rows = [[float(mass[i]), float(number[i]), float(mass_fraction[i]), float(scattering[i]), float(dp[i]), float(contour[i]), float(segment_ratio[i]), float(np.sqrt(rg2[i])), float(np.sqrt(re2[i]))] for i in range(mass.size)]
    return {"metrics": metrics, "plots": plots, "tables": [{"title": "Component-resolved population and chain dimensions", "columns": ["M (g/mol)", "Number fraction", "Mass fraction", "Zero-q intensity fraction", "Effective DP", "L (nm)", "L/b" if gaussian else "L/lp", "RMS Rg (nm)", "RMS Re (nm)"], "rows": rows}], "notes": notes, "warnings": warnings,
            "material_updates": {"molar_mass_distribution": {"molar_mass_g_mol": mass.tolist(), "number_fraction": number.tolist(), "basis": "number_fraction", "origin": "user_input_or_demonstration"}, "parameters": parameters}, "handoffs": []}


def _ris(values, material):
    keys = ("n_bonds", "bond_length_nm", "bond_angle_deg", "temperature_K", "gauche_energy_kj_mol", "opposite_gauche_penalty_kj_mol", "n_chains", "seed")
    result = sample_ris(**{key: values[key] for key in keys})
    metrics = [
        _metric("rg_rms_nm", "Ensemble RMS Rg", math.sqrt(result["mean_rg2_nm2"]), "nm", "sqrt(mean Rg²) over independent sampled conformations"),
        _metric("re_rms_nm", "Ensemble RMS end-to-end distance", math.sqrt(result["mean_re2_nm2"]), "nm"),
        _metric("mean_rg2_nm2", "Mean-square gyration radius", result["mean_rg2_nm2"], "nm²"),
        _metric("sem_rg2_nm2", "Monte Carlo SE of mean Rg²", result["sem_rg2_nm2"], "nm²", "Sampling uncertainty only; excludes uncertainty in the Hamiltonian and chemical parameters"),
        _metric("sem_re2_nm2", "Monte Carlo SE of mean Re²", result["sem_re2_nm2"], "nm²", "Independent-chain sample standard deviation divided by sqrt(number of chains)"),
        _metric("characteristic_ratio", "Finite-chain characteristic ratio Cn", result["finite_characteristic_ratio"], "", "<Re²>/(n l²); this is not an extrapolated C∞"),
        _metric("trans_fraction", "Exact mean trans fraction", result["exact_state_fraction"][0], "fraction", "Finite-chain transfer-matrix marginal, including end effects"),
        _metric("opposite_gauche_fraction", "Exact opposite-gauche pair fraction", result["exact_opposite_gauche_fraction"], "fraction", "Undefined for a chain with only one torsion"),
    ]
    rg = np.sqrt(result["rg2_nm2"])
    if float(np.ptp(rg)) <= 1e-12 * max(float(rg.mean()), 1e-12):
        # NumPy's default equal-value histogram adds an arbitrary +/-0.5 range,
        # which can create unphysical negative radius bins for short chains.
        centers = np.array([float(rg.mean())])
        histogram = np.array([values["n_chains"]])
    else:
        histogram, edges = np.histogram(rg, bins=min(32, max(5, int(math.sqrt(values["n_chains"])))))
        centers = (edges[:-1] + edges[1:]) / 2
    plots = [_plot("Sampled conformational spread", "Rg of a sampled conformation (nm)", "Sample fraction per bin", [_series("RIS ensemble", centers, histogram / values["n_chains"])])]
    projections = []
    for i in range(min(3, values["n_chains"])):
        coords = result["coordinates_nm"][i]
        projections.append(_series(f"Sample {i + 1}", coords[:, 0], coords[:, 1]))
    plots.append(_plot("Equal-mass skeletal sites · first three XY projections", "x (nm)", "y (nm)", projections))
    notes = _material_notes(material) + [
        "RIS defaults are synthetic and are not a validated parameter set for any named polymer. Local torsion energies and pair penalties need independent chemical calibration.",
        "The torsion distribution is sampled exactly for the specified open finite-chain Hamiltonian using backward messages. Successive chains are independent; no equilibration or burn-in step is used.",
        "Each chain has n bonds, n+1 equal-mass sites and n−2 torsions. Torsions are trans=180°, gauche+=60°, gauche−=−60°. The internal bond angle is fixed.",
        "Coordinates share a fixed initial orientation. Distances are rotationally invariant; the plot is an arbitrary XY projection and does not establish the absence of chain crossings.",
        "No nonlocal excluded volume, solvent, electrostatics, intermolecular packing, or kinetics is present. This is an equilibrium local-conformation model.",
        "RIS does not use the MWD fields and does not assign a chemical repeat unit or molar mass to a skeletal site. It exports only its explicitly scoped parameters.",
        f"Sampling record: {values['n_chains']} independent chains; {values['n_bonds']} bonds; seed={values['seed']}; NumPy {np.__version__}. Standard errors describe sampling of this model, not experimental accuracy.",
    ]
    if values["n_bonds"] == 3:
        notes.append("Only one torsion is present, so the adjacent-pair fraction is undefined and reported as null.")
    rows = [[name, float(result["exact_state_fraction"][i]), float(result["sampled_state_fraction"][i])] for i, name in enumerate(["trans", "gauche+", "gauche−"])]
    rows.append(["opposite-gauche adjacent pair", result["exact_opposite_gauche_fraction"], result["sampled_opposite_gauche_fraction"]])
    parameters = [_parameter(quantity, values[key], unit, "synthetic_or_user_calibrated_uniform_three_state_finite_chain_RIS") for key, quantity, unit in [
        ("bond_length_nm", "bond_length", "nm"), ("bond_angle_deg", "bond_angle", "degree"),
        ("gauche_energy_kj_mol", "gauche_state_energy_relative_to_trans", "kJ/mol"),
        ("opposite_gauche_penalty_kj_mol", "adjacent_opposite_gauche_energy", "kJ/mol")]]
    return {"metrics": metrics, "plots": plots, "tables": [{"title": "Exact versus sampled local conformations", "columns": ["State or pair", "Exact fraction", "Sampled fraction"], "rows": rows}], "notes": notes, "warnings": [], "material_updates": {"parameters": parameters}, "handoffs": []}


def calculate(inputs: dict, material: dict | None = None) -> dict:
    """Calculate a bounded, finite, JSON-serializable research result."""
    values = _validated_inputs(inputs)
    if material is not None and not isinstance(material, dict):
        raise ValueError("Material must be an object when supplied.")
    material = material or {}
    chemistry = material.get("chemistry", {})
    if isinstance(chemistry, dict):
        architecture = str(chemistry.get("architecture", "")).strip().lower()
        if any(term in architecture for term in ("cyclic", "ring", "branch", "star", "brush", "network", "crosslink", "dendri")):
            raise ValueError("The selected chain models require open linear chains; the material declares an unsupported architecture.")
    if values["model"] == "petox_md":
        return md_output(values)
    if values["model"] == "intrinsic_viscosity_fit":
        return _viscosity_benchmark(values)
    if values["model"] == "intrinsic_viscosity_mwd":
        return _viscosity_mixture(values, material)
    if values["model"] == "ris":
        return _ris(values, material)
    return _population(values, material)


def _viscosity_benchmark(values):
    custom = values["calibration_csv"].strip()
    observations = read_observations(custom) if custom else literature_series()
    rows, predictions, k, a = benchmark(observations, values["viscosity_temp_C"])
    mass = np.asarray([r.mw for r in rows])
    observed = [r.eta for r in rows]
    held_out = [predictions[r.sample][0] for r in rows]
    errors = [predictions[r.sample][1] for r in rows]
    query = values["target_mw_g_mol"]
    in_domain = min(mass) <= query <= max(mass)
    target = k * query**a if in_domain else None
    metrics = [
        _metric("mh_a", "Sample-Mw scaling exponent", a),
        _metric("mh_k", "Sample-Mw prefactor K", k, "(mL/g)/(g/mol)^a"),
        _metric("mass_class_count", "Withheld mass classes", len({r.group for r in rows})),
        _metric("heldout_median_ape", "Median held-out sample error", np.median(errors), "%"),
        _metric("heldout_p90_ape", "90th percentile held-out sample error", np.percentile(errors, 90), "%"),
        _metric("heldout_mae", "Mean absolute held-out sample error", np.mean(np.abs(np.array(held_out) - observed)), "mL/g"),
        _metric("target_intrinsic_viscosity", "Interpolated intrinsic viscosity", target, "mL/g",
                f"Target Mw = {query:g} g/mol, same medium, temperature and comparable dispersity."),
    ]
    r0 = rows[0]
    notes = [
        f"Viscometry: {r0.medium} at {r0.temperature:g} °C. SEC–MALS: {r0.sec_eluent} at {r0.sec_temperature:g} °C.",
        "Each validation fold trains on all other mass classes at this temperature. Rows outside the remaining training mass range are identified in the last column.",
        "K and a come from unweighted log-space regression on sample Mw. Component-level constants for a broad distribution require their own calibration.",
    ]
    if not custom:
        notes += ["Gubarev et al., SI S1/S2: dilute-solution [η] from averaged Huggins/Kraemer extrapolations. Published dispersities of these samples span 1.04–1.18.",
                  "This cross-validation uses one experimental series. The authors' 37 °C scaling law used sedimentation–diffusion masses; this implementation uses SEC–MALS Mw."]
    warnings = [] if in_domain else ["The requested Mw falls outside the calibrated mass interval; the interpolation result is empty."]
    table_rows = [[r.sample, r.group, r.mw, r.eta, *predictions[r.sample], r.source] for r in rows]
    mass_grid = np.geomspace(min(mass), max(mass), 50)
    return {"metrics": metrics, "plots": [
        _plot("PEtOx intrinsic viscosity at fixed conditions", "Sample Mw (g/mol)", "[η] (mL/g)", [
            _series("Measured", mass, observed), _series("Held-out prediction", mass, held_out),
            _series("Fit using all classes", mass_grid, k * mass_grid**a)], "log", "log")],
        "tables": [{"title": "Leave-one-mass-class-out predictions",
                    "columns": ["Sample", "Mass class", "Mw (g/mol)", "Measured (mL/g)", "Held out (mL/g)", "Absolute error (%)", "Fold training interval", "Source"],
                    "rows": table_rows}],
        "notes": notes, "warnings": warnings, "material_updates": {}, "handoffs": []}


def _viscosity_mixture(values, material):
    condition = values["component_medium"].strip()
    source = values["component_source"].strip()
    if not condition or not source or max(len(condition), len(source)) > 500:
        raise ValueError("Supply a polymer/solvent/temperature description and a calibration source (1–500 characters each).")
    distribution = _distribution(values, material)
    k, a = values["component_k"], values["component_a"]
    eta, mv, components = mixture(distribution, k, a, values["component_low_g_mol"], values["component_high_g_mol"])
    contributions = components * distribution.mass_fraction
    metrics = [_metric("mixture_intrinsic_viscosity", "Intrinsic viscosity at infinite dilution", eta, "mL/g"),
               _metric("viscosity_average_mass", "Viscosity-average molar mass Mv", mv, "g/mol"),
               _metric("Mw", "Mass-average molar mass Mw", distribution.moments()["Mw"], "g/mol")]
    rows = [[float(m), float(w), float(e), float(c)] for m, w, e, c in zip(
        distribution.molar_mass_g_mol, distribution.mass_fraction, components, contributions)]
    return {"metrics": metrics, "plots": [
        _plot("Intrinsic viscosity contributions by mass", "Component M (g/mol)", "Contribution (mL/g)",
              [_series("wᵢ [η]ᵢ", distribution.molar_mass_g_mol, contributions)], "log")],
        "tables": [{"title": "Mass-weighted dilute-solution contributions",
                    "columns": ["M (g/mol)", "Mass fraction", "Component [η] (mL/g)", "Contribution (mL/g)"],
                    "rows": rows}],
        "notes": [f"Component calibration: {condition}. Source: {source}.",
                  "Weights are mass fractions after conversion of the supplied MWD. All components must share the calibration chemistry, linear architecture, solvent and temperature."],
        "warnings": [], "material_updates": {}, "handoffs": []}
