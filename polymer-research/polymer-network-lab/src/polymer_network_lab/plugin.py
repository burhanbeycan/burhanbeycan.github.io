"""JSON-friendly workbench adapter; the scientific model has no GUI dependency."""

from __future__ import annotations

import numpy as np

from .wbpu_ui import calculate_wbpu, extend_app
from .model import (
    MAX_SCAN_POINTS,
    PARAMETER_LIMITS,
    NetworkParameters,
    equilibrium_scan,
    finite_number,
    free_energy_components,
    integer_in_range,
    slab_release_fraction,
    solve_equilibrium,
    time_to_release_fraction,
)


REFERENCES = [
    {
        "id": "hong2008",
        "title": "Hong, Zhao, Zhou & Suo — A theory of coupled diffusion and large deformation in polymeric gels",
        "url": "https://doi.org/10.1016/j.jmps.2007.11.010",
        "year": 2008,
    },
    {
        "id": "hong2009",
        "title": "Hong, Liu & Suo — Inhomogeneous swelling of a gel in equilibrium with a solvent and mechanical load",
        "url": "https://doi.org/10.1016/j.ijsolstr.2009.04.022",
        "year": 2009,
    },
    {
        "id": "flory1943",
        "title": "Flory & Rehner — Statistical Mechanics of Cross-Linked Polymer Networks II. Swelling",
        "url": "https://doi.org/10.1063/1.1723792",
        "year": 1943,
    },
    {
        "id": "fick1855",
        "title": "Fick — On liquid diffusion",
        "url": "https://doi.org/10.1080/14786445508641925",
        "year": 1855,
    },
]

APP = {
    "id": "network",
    "name": "Polymer Network Lab",
    "version": "0.3.0",
    "tagline": "WBPU stoichiometry, film characterization and swelling",
    "description": (
        "Audit polyurethane feed, examine extraction and DMA, and compare withheld formulations under a common swelling protocol. "
        "Explore neutral-network equilibrium and preview solute release "
        "from the fully equilibrated slab."
    ),
    "accent": "#0f9d8a",
    "scope": (
        "Homogeneous, neutral, dry-reference Gaussian network in a pure solvent; constant χ; free isotropic swelling. "
        "The WBPU modes use their own mass balances and experimental protocols. Release uses fixed equilibrium geometry and a user-supplied D."
    ),
    "fields": [
        {"key": "chi", "label": "Flory–Huggins interaction parameter χ", "type": "number", "default": 0.35,
         "unit": "dimensionless", "min": -5.0, "max": 5.0, "step": 0.01, "group": "Network thermodynamics",
         "help": "Constant χ for this state. A demonstration value is not a calibrated material property; no χ(T,φ) law is inferred."},
        {"key": "nu_e_mol_m3", "label": "Effective elastic density νe", "type": "number", "default": 100.0,
         "unit": "mol/m³ dry polymer", "min": 1.0e-6, "max": 1.0e6, "step": "any", "group": "Network thermodynamics",
         "help": "Effective parameter multiplying the stated elastic free energy. Supply it independently; nominal crosslinker feed, loops and dangling chains do not determine it here."},
        {"key": "solvent_molar_volume_m3_mol", "label": "Solvent molar volume Vs", "type": "number", "default": 18.0e-6,
         "unit": "m³/mol", "min": 1.0e-8, "max": 1.0e-3, "step": "any", "group": "Network thermodynamics",
         "help": "Incompressible-solvent molar volume. Enter 18e-6 for 18 cm³/mol; the default is an illustrative input, not a solvent database."},
        {"key": "temperature_K", "label": "Temperature", "type": "number", "default": 298.15,
         "unit": "K", "min": 1.0, "max": 2000.0, "step": 0.1, "group": "Network thermodynamics",
         "help": "At fixed χ, νe and Vs, temperature rescales energy and modulus but does not change this ideal equilibrium Q."},
        {"key": "scan_variable", "label": "Equilibrium scan", "type": "select", "default": "chi", "unit": "", "group": "Parameter scan",
         "options": [{"value": "chi", "label": "Solvent interaction χ"}, {"value": "nu_e_mol_m3", "label": "Effective elastic density νe"},
                     {"value": "none", "label": "No scan"}],
         "help": "Each point is independently equilibrated; this is a parameter scan, not a swelling-time trajectory."},
        {"key": "chi_scan_min", "label": "Minimum χ", "type": "number", "default": 0.0,
         "unit": "dimensionless", "min": -5.0, "max": 5.0, "step": 0.01, "group": "Parameter scan", "show_if": {"key": "scan_variable", "value": "chi"},
         "help": "Lower endpoint of the linear χ scan."},
        {"key": "chi_scan_max", "label": "Maximum χ", "type": "number", "default": 0.8,
         "unit": "dimensionless", "min": -5.0, "max": 5.0, "step": 0.01, "group": "Parameter scan", "show_if": {"key": "scan_variable", "value": "chi"},
         "help": "Must exceed the minimum χ."},
        {"key": "nu_e_scan_min", "label": "Minimum νe", "type": "number", "default": 10.0,
         "unit": "mol/m³ dry polymer", "min": 1.0e-6, "max": 1.0e6, "step": "any", "group": "Parameter scan", "show_if": {"key": "scan_variable", "value": "nu_e_mol_m3"},
         "help": "Lower endpoint of the logarithmic elastic-density scan."},
        {"key": "nu_e_scan_max", "label": "Maximum νe", "type": "number", "default": 1000.0,
         "unit": "mol/m³ dry polymer", "min": 1.0e-6, "max": 1.0e6, "step": "any", "group": "Parameter scan", "show_if": {"key": "scan_variable", "value": "nu_e_mol_m3"},
         "help": "Must exceed the minimum νe."},
        {"key": "scan_points", "label": "Scan points", "type": "number", "default": 61,
         "unit": "", "min": 2, "max": MAX_SCAN_POINTS, "step": 1, "group": "Parameter scan",
         "help": "Bounded number of independent equilibrium calculations."},
        {"key": "dry_thickness_um", "label": "Full dry slab thickness", "type": "number", "default": 100.0,
         "unit": "µm", "min": 0.001, "max": 1.0e6, "step": "any", "group": "Equilibrated slab",
         "help": "For a freely swelling slab, full swollen thickness = full dry thickness × Q^(1/3). A substrate-constrained film is outside this model."},
        {"key": "release_preview", "label": "Fixed-geometry release preview", "type": "select", "default": "enabled", "unit": "", "group": "Equilibrated slab",
         "options": [{"value": "enabled", "label": "Show release after equilibration"}, {"value": "disabled", "label": "Equilibrium only"}],
         "help": "Initially uniform solute in an already swollen slab; perfect sinks on both faces. Dynamic swelling, relaxation and degradation are excluded."},
        {"key": "diffusivity_m2_s", "label": "Supplied effective solute diffusivity D", "type": "number", "default": 1.0e-12,
         "unit": "m²/s", "min": 1.0e-24, "max": 1.0e-4, "step": "any", "group": "Equilibrated slab",
         "show_if": {"key": "release_preview", "value": "enabled"},
         "help": "Independent effective solute diffusivity in the fully swollen material. This application does not infer D from νe, χ or a mesh size."},
        {"key": "duration_s", "label": "Release observation duration", "type": "number", "default": 3600.0,
         "unit": "s", "min": 1.0e-9, "max": 1.0e12, "step": "any", "group": "Equilibrated slab",
         "show_if": {"key": "release_preview", "value": "enabled"}, "help": "Time starts after swelling has fully equilibrated."},
        {"key": "n_points", "label": "Release curve points", "type": "number", "default": 240,
         "unit": "", "min": 20, "max": 600, "step": 1, "group": "Equilibrated slab",
         "show_if": {"key": "release_preview", "value": "enabled"}, "help": "Number of uniformly spaced times; t=0 is included exactly."},
        {"key": "input_origin", "label": "Parameter origin", "type": "select", "default": "demonstration", "unit": "", "group": "Provenance",
         "options": [{"value": "demonstration", "label": "Synthetic demonstration"}, {"value": "user_supplied", "label": "User supplied; not independently verified"}],
         "help": "Recording provenance does not experimentally validate the model or the parameters."},
        {"key": "parameter_provenance", "label": "Parameter sources / fitting context", "type": "textarea",
         "default": "Illustrative χ, νe, Vs and D. Replace with traceable values for a specific sample, temperature and solvent.",
         "unit": "", "group": "Provenance", "help": "Record measurements, fitting method, source identifiers and the elastic convention used."},
        {"key": "topology_notes", "label": "Loops, dangling chains and network preparation", "type": "textarea",
         "default": "Loops and dangling chains have not been quantified. No correction from nominal crosslinker feed has been applied.",
         "unit": "", "group": "Provenance", "help": "Preparation, loop/dangling-chain evidence and reference-state details are retained as metadata; they do not silently change νe."},
    ],
    "presets": [
        {"name": "Neutral-network demonstration", "description": "Synthetic neutral network with constant χ; no material-specific validation.",
         "inputs": {"chi": 0.35, "nu_e_mol_m3": 100.0, "scan_variable": "chi", "input_origin": "demonstration"}},
        {"name": "Elastic-density demonstration", "description": "Synthetic scan of effective elastic density at fixed solvent interaction.",
         "inputs": {"chi": 0.2, "nu_e_mol_m3": 100.0, "scan_variable": "nu_e_mol_m3", "nu_e_scan_min": 10.0,
                    "nu_e_scan_max": 1000.0, "input_origin": "demonstration"}},
        {"name": "Published-model numerical check", "description": "Hong et al. (2008) model parameters a=νeVs=0.001, χ=0.2 give λ≈3.215; a model benchmark, not experimental validation.",
         "inputs": {"chi": 0.2, "nu_e_mol_m3": 55.55555555555556, "solvent_molar_volume_m3_mol": 18.0e-6,
                    "temperature_K": 298.15, "input_origin": "demonstration",
                    "parameter_provenance": "Numerical model benchmark: Hong et al. (2008), a=νeVs=10^-3 and χ=0.2. The selected Vs only converts a to molar νe; it is not a measured material dataset."}},
    ],
    "models": [
        {"name": "Neutral free swelling: Hong/Suo elastic convention",
         "equation": "W/(RT/Vs) = (Q−1) ln(1−1/Q) + χ(1−1/Q) + (νeVs/2)[3Q^(2/3)−3−2ln Q]; g(φ)=ln(1−φ)+φ+χφ²+νeVs(φ^(1/3)−φ)=0",
         "assumptions": ["Dry polymer volume is the reference; Q=1/φ and λ=Q^(1/3).",
                         "Neutral, homogeneous, Gaussian network; constant χ; incompressible constituents; pure solvent reservoir; zero external load.",
                         "The same free energy supplies equilibrium and stability; the lowest-energy stable homogeneous root is selected.",
                         "The −2ln Q convention yields an elastic −φ term. Parameters from a historical −φ/2 equation require reconsideration."],
         "references": ["hong2008", "hong2009", "flory1943"]},
        {"name": "Small isochoric shear about the swollen state",
         "equation": "G_dry = νe RT; G_swollen = G_dry / λ",
         "assumptions": ["Ideal Gaussian elasticity with the same effective νe.", "Incremental shear at fixed solvent content; no additional phantom-network prefactor or topology correction."],
         "references": ["hong2008", "hong2009"]},
        {"name": "Solute release after full swelling equilibration",
         "equation": "L=L_dry Q^(1/3); Mt/M∞=1−(8/π²)Σ[n≥0] exp[−(2n+1)²π²Dt/L²]/(2n+1)²",
         "assumptions": ["L is FULL thickness; perfect sinks on both faces; initially uniform solute; fixed geometry and supplied constant D.",
                         "This analytical diffusion preview ignores dynamic swelling, binding, erosion, degradation and external mass-transfer resistance.",
                         "The slab series is derived from the stated diffusion boundary-value problem; the Fick reference establishes the constitutive law."],
         "references": ["fick1855"]},
    ],
    "references": REFERENCES,
}


extend_app(APP)

DEFAULTS = {field["key"]: field["default"] for field in APP["fields"]}
MODEL_SCOPE = "neutral_homogeneous_free_swelling_dry_reference_hong_neo_hookean_constant_chi"


def _validated_inputs(inputs: dict) -> dict:
    if not isinstance(inputs, dict):
        raise ValueError("Calculation inputs must be an object/dictionary.")
    unknown = sorted(set(inputs) - set(DEFAULTS), key=str)
    if unknown:
        raise ValueError(f"Unknown network input(s): {', '.join(str(key) for key in unknown)}.")
    values = {**DEFAULTS, **inputs}
    for field in APP["fields"]:
        key = field["key"]
        value = values[key]
        if field["type"] == "number":
            if key in ("scan_points", "n_points"):
                values[key] = integer_in_range(value, key, int(field["min"]), int(field["max"]))
            else:
                values[key] = finite_number(value, key, field["min"], field["max"])
        elif field["type"] == "select":
            if not isinstance(value, str) or value not in [option["value"] for option in field["options"]]:
                raise ValueError(f"{key} must be one of the available choices.")
        else:
            limit = 50000 if key.endswith('_csv') else 8000
            if not isinstance(value, str) or len(value) > limit:
                raise ValueError(f"{key} must be text of at most {limit} characters.")
    if values["scan_variable"] == "chi" and values["chi_scan_max"] <= values["chi_scan_min"]:
        raise ValueError("The maximum χ must exceed the minimum χ.")
    if values["scan_variable"] == "nu_e_mol_m3" and values["nu_e_scan_max"] <= values["nu_e_scan_min"]:
        raise ValueError("The maximum νe must exceed the minimum νe.")
    return values


def _metric(key: str, label: str, value: float, unit: str, note: str) -> dict:
    return {"key": key, "label": label, "value": float(value), "unit": unit, "note": note}


def _parameter(quantity: str, value: float, unit: str, origin: str, note: str, scope: str = MODEL_SCOPE) -> dict:
    return {"quantity": quantity, "value": float(value), "unit": unit, "origin": origin,
            "model_scope": scope, "notes": note}


def calculate(inputs: dict, material: dict | None = None) -> dict:
    values = _validated_inputs(inputs)
    if material is not None and not isinstance(material, dict):
        raise ValueError("The material record must be an object/dictionary.")
    material = material or {}
    if values["network_mode"] != "equilibrium":
        return calculate_wbpu(values)
    parameters = NetworkParameters(**{key: values[key] for key in PARAMETER_LIMITS})
    result = solve_equilibrium(parameters)
    root = result.selected
    q = root.volume_swelling_ratio
    stretch = root.isotropic_stretch
    full_thickness_um = values["dry_thickness_um"] * stretch
    full_thickness_m = full_thickness_um * 1.0e-6
    origin = values["input_origin"]
    notes = [
        "All equilibrium quantities use dry polymer volume as the reference. Q is a volume ratio, not a mass-uptake ratio.",
        "This is the Hong/Suo Gaussian elastic convention: −2 ln Q in the elastic energy and −φ in the equilibrium elastic term.",
        "The selected state minimizes the same free energy among the locally stable homogeneous roots found within 10⁻¹⁴ ≤ φ ≤ 1−10⁻¹⁴; this does not solve inhomogeneous phase coexistence.",
        "χ, νe and Vs are held constant. Temperature scales energy and modulus; a measured or fitted χ(T,φ) law is not supplied.",
        "Effective νe is an explicit input. Loops, dangling chains, nominal crosslinker percentage and preparation history are provenance, not inferred corrections.",
        "A network prepared in solvent can require a different reference-state model; the stated dry-reference constitutive law must be justified for the sample.",
        "The modulus is an ideal incremental shear modulus at fixed solvent content. It is not a drained Young's modulus or a topology measurement.",
        "These are model calculations and numerical verification targets; no experimental validation or material-specific fit has been performed by this application.",
    ]
    if material.get("chemistry"):
        notes.append("Imported chemistry, sequence, tacticity and end-group metadata are retained by the workbench but do not determine χ, νe or D in this model.")
    if material.get("molar_mass_distribution"):
        notes.append("Imported molar-mass distribution is not consumed by the network solver; νe and diffusivity are not inferred from it.")
    if material.get("parameters") or material.get("conditions"):
        notes.append("Imported material parameters and conditions are not automatically substituted. Explicit network-form inputs define the present calculation.")
    warnings = []
    if origin == "demonstration":
        warnings.append("The parameters are labelled as a synthetic demonstration; they are not measured properties of a named polymer.")
    if stretch > 5.0:
        warnings.append("The computed stretch exceeds 5. Gaussian chains have no finite-extension limit; assess strand extensibility before applying this result to a material.")
    if root.polymer_volume_fraction > 0.8:
        warnings.append("The selected state is polymer-rich. Rubbery-network, constant-χ and incompressibility assumptions require particular scrutiny near dry or glassy states.")
    metrics = [
        _metric("polymer_volume_fraction", "Equilibrium polymer volume fraction φ", root.polymer_volume_fraction, "", "Selected stable homogeneous equilibrium."),
        _metric("volume_swelling_ratio", "Equilibrium volume ratio Q", q, "V/Vdry", "Q=1/φ; distinct from solvent mass uptake."),
        _metric("isotropic_stretch", "Isotropic stretch λ", stretch, "", "λ=Q^(1/3); freely swelling geometry."),
        _metric("swollen_shear_modulus_pa", "Ideal swollen shear modulus", root.swollen_shear_modulus_pa, "Pa", "G=νeRT/λ; small isochoric shear at fixed solvent content."),
        _metric("dry_shear_modulus_pa", "Ideal dry shear parameter", parameters.dry_shear_modulus_pa, "Pa", "Gdry=νeRT in the selected elastic convention."),
        _metric("full_swollen_thickness_um", "Full equilibrium slab thickness", full_thickness_um, "µm", "Full dry thickness × λ; requires free isotropic swelling."),
        _metric("equilibrium_relative_residual", "Relative equilibrium residual", root.relative_residual, "", "|g| divided by the sum of absolute terms in the equilibrium equation; a numerical diagnostic."),
    ]
    plots = []
    tables = [
        {"title": "Stationary states and thermodynamic selection", "columns": ["State", "φ", "Q", "λ", "Local stability", "W (J/m³ dry)", "g=Δμ/RT", "Selected"],
         "rows": [[index + 1, state.polymer_volume_fraction, state.volume_swelling_ratio, state.isotropic_stretch,
                   state.stability, state.free_energy_j_m3_dry, state.chemical_potential_residual_over_rt,
                   "yes" if index == result.selected_index else "no"] for index, state in enumerate(result.roots)]},
        {"title": "Parameter and topology provenance", "columns": ["Item", "Value / explanation"],
         "rows": [["Input origin", origin], ["νe basis", "mol of effective elastic contribution per m³ of dry polymer"],
                  ["νeVs", parameters.elastic_to_mixing_ratio], ["Parameter sources", values["parameter_provenance"]],
                  ["Loops, dangling chains, preparation", values["topology_notes"]],
                  ["Topology correction", "None inferred from crosslinker feed, MWD or metadata; no additional phantom prefactor applied."],
                  ["Solvent reservoir", "Pure solvent, reference chemical potential zero; no applied mechanical load."]]},
    ]

    if values["scan_variable"] != "none":
        variable = values["scan_variable"]
        if variable == "chi":
            scan_values = np.linspace(values["chi_scan_min"], values["chi_scan_max"], values["scan_points"])
            x_label, x_scale = "Interaction parameter χ", "linear"
        else:
            scan_values = np.geomspace(values["nu_e_scan_min"], values["nu_e_scan_max"], values["scan_points"])
            x_label, x_scale = "Effective νe (mol/m³ dry polymer)", "log"
        scan = equilibrium_scan(parameters, variable, scan_values)
        plots.append({"title": "Equilibrium swelling across independently equilibrated states", "x_label": x_label,
                      "y_label": "Volume swelling ratio Q", "x_scale": x_scale, "y_scale": "log",
                      "series": [{"name": "Selected stable equilibrium", "x": [p.value for p in scan],
                                  "y": [p.volume_swelling_ratio for p in scan]}]})
        tables.append({"title": "Equilibrium scan", "columns": [x_label, "φ", "Q", "λ", "G swollen (Pa)", "Stable roots"],
                       "rows": [[p.value, p.polymer_volume_fraction, p.volume_swelling_ratio, p.isotropic_stretch,
                                 p.swollen_shear_modulus_pa, p.stable_root_count] for p in scan]})

    q_grid = np.unique(np.concatenate((np.geomspace(max(1.000001, q / 5.0), max(2.0, q * 5.0), 300), np.array([q]))))
    energy = free_energy_components(q_grid, parameters)
    plots.append({"title": "Consistent mixing and elastic free energy", "x_label": "Volume swelling ratio Q", "y_label": "Free energy (J/m³ dry polymer)",
                  "x_scale": "log", "y_scale": "linear", "series": [
                      {"name": label, "x": q_grid.tolist(), "y": np.asarray(energy[key]).tolist()}
                      for key, label in (("mixing", "Mixing"), ("elastic", "Elastic"), ("total", "Total"))
                  ]})

    output_parameters = [
        _parameter("chi", parameters.chi, "1", origin, values["parameter_provenance"]),
        _parameter("nu_e_mol_m3", parameters.nu_e_mol_m3, "mol/m3_dry", origin, values["parameter_provenance"] + " " + values["topology_notes"]),
        _parameter("solvent_molar_volume_m3_mol", parameters.solvent_molar_volume_m3_mol, "m3/mol", origin, values["parameter_provenance"]),
        _parameter("temperature_K", parameters.temperature_K, "K", origin, "Explicit calculation temperature; χ and Vs remain fixed."),
        _parameter("polymer_volume_fraction", root.polymer_volume_fraction, "1", "model_prediction", "Selected stable homogeneous equilibrium; no experimental validation."),
        _parameter("volume_swelling_ratio", q, "1", "model_prediction", "Vswollen/Vdry; free isotropic equilibrium."),
        _parameter("isotropic_stretch", stretch, "1", "model_prediction", "Q^(1/3); no substrate constraint."),
        _parameter("swollen_shear_modulus_pa", root.swollen_shear_modulus_pa, "Pa", "model_prediction", "Incremental isochoric shear with fixed solvent content."),
        _parameter("full_swollen_thickness_m", full_thickness_m, "m", "model_prediction", "Full dry thickness multiplied by equilibrium isotropic stretch."),
    ]
    handoffs = []
    if values["release_preview"] == "enabled":
        times = np.linspace(0.0, values["duration_s"], values["n_points"])
        released = np.asarray(slab_release_fraction(times, full_thickness_m, values["diffusivity_m2_s"]))
        t50 = time_to_release_fraction(0.5, full_thickness_m, values["diffusivity_m2_s"])
        t90 = time_to_release_fraction(0.9, full_thickness_m, values["diffusivity_m2_s"])
        metrics.extend([
            _metric("released_fraction_at_end", "Release at observation end", float(released[-1]), "fraction", "Fully swollen, fixed full thickness; constant supplied solute D; two perfect sinks."),
            _metric("release_t50_s", "Fixed-slab time to 50% release", t50, "s", "Calculated independently of the selected observation duration; ignores dynamic swelling."),
            _metric("release_t90_s", "Fixed-slab time to 90% release", t90, "s", "Calculated independently of the selected observation duration; ignores dynamic swelling."),
        ])
        plots.append({"title": "Release after equilibration — fixed geometry; dynamic swelling excluded", "x_label": "Time after equilibration (s)",
                      "y_label": "Released fraction Mt/M∞", "x_scale": "linear", "y_scale": "linear",
                      "series": [{"name": "Two-face perfect-sink slab", "x": times.tolist(), "y": released.tolist()}]})
        notes.append("Release starts with a uniformly loaded, fully equilibrated slab. D is an independent user input. Network geometry is handed to diffusion without solving coupled solvent uptake, relaxation or degradation.")
        if t50 > values["duration_s"]:
            notes.append("The predicted 50% release time lies beyond the displayed observation window.")
        output_parameters.append(_parameter("diffusivity_m2_s", values["diffusivity_m2_s"], "m2/s", origin,
                                            "Independently supplied effective solute diffusivity, never inferred from νe or mesh size.",
                                            "uniform_solute_fixed_fully_swollen_slab_constant_D_two_face_perfect_sink"))
        handoff_notes = [
            "The thickness_um value is the FULL slab thickness after free isotropic swelling, not a half-thickness.",
            "Fully equilibrated, fixed-geometry approximation; does not simulate dynamic swelling or solvent-uptake kinetics.",
            "diffusivity_m2_s is the independently supplied effective solute D; no mesh-size or crosslink-density conversion.",
            "Initially uniform solute; perfect sinks on BOTH faces; the material record is not experimentally verified.",
        ]
        if 0.001 <= full_thickness_um <= 1.0e6:
            handoffs.append({"target": "transport", "label": "Open fully swollen slab in Transport Lab",
                             "inputs": {"mode": "release", "thickness_um": full_thickness_um,
                                        "diffusivity_m2_s": values["diffusivity_m2_s"], "temperature_K": parameters.temperature_K,
                                        "duration_s": values["duration_s"], "n_points": values["n_points"]},
                             "notes": handoff_notes})
        else:
            warnings.append("The computed full swollen thickness exceeds Transport Lab's supported thickness range; no transport handoff was generated.")
    for parameter in output_parameters:
        parameter["source_input_origin"] = origin
    return {"metrics": metrics, "plots": plots, "tables": tables, "notes": notes,
            "warnings": warnings, "material_updates": {"parameters": output_parameters}, "handoffs": handoffs}
