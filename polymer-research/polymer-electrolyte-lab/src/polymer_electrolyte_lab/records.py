"""Editable GUI records, with source and specimen scope preserved for each measurement."""
import json

from .peo import analyze_peo


def literature_record(state_id, temperature, state, recipe, tracer, eis):
    common = {"state_id": state_id, "temperature_K": temperature}
    recipe = {"state_id": state_id, **recipe}
    if recipe.get("volume_at_state") is not None:
        recipe["volume_at_state"] = {"state_id": state_id, **recipe["volume_at_state"]}
    return {"state": {**common, **state}, "formulation": recipe,
            "tracer": {**common, **tracer}, "eis": {**common, **eis}}


PESKO = literature_record(
    "Pesko2017-5k-r008-90C", 363.15,
    {"peo_mass_kg_mol": 5, "mass_average": "Nominal supplier molar mass as described in Methods",
     "phase": "Homogeneous bulk PEO/LiTFSI; quantitative phase fractions unavailable",
     "thermal_history": "Dry handling under argon. EIS samples annealed 4 h at 90 °C.",
     "crystallinity": None, "source": "10.1149/2.0581711jes, Experimental, pp. E3570–E3571."},
    {"peo_g": 1, "eo_per_li": 12.5, "source": "1 g calculation basis; Table I reports r=0.08.",
     "volume_at_state": {"temperature_K": 363.15, "density_g_mL": 1.330,
                         "source": "Pesko et al. (2017), Table I. One density determination for this composition."}},
    {"D_Li_m2_s": 7.5e-12, "D_TFSI_m2_s": 3.0e-11, "relative_plot_reading_bound": 0.09,
     "source": "10.1149/2.0581711jes, Fig. 4, r=0.08, 90 °C. Approximate visual readings on a log axis."},
    {"sigma_mS_cm": 2.0, "source": "10.1149/2.0581711jes, Results for Fig. 2a, p. E3571. Reported rounded peak value."},
)

ZHENG = literature_record(
    "Zheng2018-100k-r008-90C", 363.15,
    {"peo_mass_kg_mol": 100, "mass_average": "Nominal supplier grade",
     "phase": "Bulk PEO/LiTFSI; crystalline fraction unavailable for this measurement",
     "thermal_history": "Vacuum-dried PEO and dry salt. Conductivity measured in sealed cells over 25–110 °C.",
     "crystallinity": None, "source": "10.1021/acs.macromol.7b02706, Electrolyte Preparation and Electrochemical Measurements."},
    {"peo_g": 1, "eo_per_li": 12.5, "source": "1 g calculation basis. Zheng et al., Table 3 uses r=0.08.",
     "volume_at_state": None},
    {"D_Li_m2_s": 6.4e-12, "D_TFSI_m2_s": 2.7e-11,
     "source": "10.1021/acs.macromol.7b02706, Table 3, PEO row, experiment (cm²/s converted to m²/s)."},
    {"sigma_mS_cm": 1.5, "source": "10.1021/acs.macromol.7b02706, Table 3, PEO row, experiment."},
)

SYNTHETIC = literature_record(
    "Example-film-1", 333.15,
    {"peo_mass_kg_mol": 100, "mass_average": "Mw (illustrative)", "phase": "Semicrystalline film (illustrative)",
     "thermal_history": "Illustrative heating branch", "crystallinity": {
         "fraction": 0.18, "basis": "PEO mass", "method": "Illustrative DSC-derived fraction"},
     "source": "Constructed example for checking units and specimen matching."},
    {"peo_g": 1, "salt_g": 0.3, "source": "Constructed recipe",
     "volume_at_state": {"temperature_K": 333.15, "final_volume_mL": 1.1, "source": "Constructed final mixture volume"}},
    {"D_Li_m2_s": 1e-12, "D_TFSI_m2_s": 3e-12, "source": "Synthetic PFG-NMR example"},
    {"thickness_um": 100, "area_cm2": 1, "bulk_resistance_ohm": 100,
     "source": "Synthetic blocking-electrode cell geometry and bulk resistance"},
)

RECORDS = [
    ("PEO/LiTFSI · Pesko 2017", "5 kg/mol PEO; EIS, mixture density and figure-read tracer diffusivities at 90 °C.", PESKO),
    ("PEO/LiTFSI · Zheng 2018", "100 kg/mol PEO; tabulated EIS and tracers at 90 °C; no borrowed mixture density.", ZHENG),
    ("PEO/LiTFSI · weighed film", "Constructed sample with mass-based composition, final volume, DSC metadata and cell geometry.", SYNTHETIC),
]
FORM_GROUPS = [
    ("state", "Material and thermal state", "A state ID identifies specimens matched in composition, thermal history and temperature. Supply crystallinity as null if unavailable."),
    ("formulation", "Dry formulation", "Use salt_g or eo_per_li, plus a dry PEO mass. Volume evidence can be mixture density, final mixture volume or independent molarity at the same temperature."),
    ("tracer", "Tracer self-diffusion", "Enter ⁷Li and ¹⁹F self-diffusivities in m²/s for this PEO/LiTFSI state and the measurement source."),
    ("eis", "Bulk impedance", "Use measured sigma_mS_cm or bulk_resistance_ohm, area_cm2 and thickness_um. Use the electrode overlap area and final film thickness at measurement."),
]


def gui_inputs(record):
    return {"model": "peo_matched", **{f"peo_{key}": json.dumps(value, ensure_ascii=False, indent=2)
                                       for key, value in record.items()}}


def peo_results(inputs):
    result = analyze_peo(inputs["peo_state"], inputs["peo_formulation"], inputs["peo_tracer"], inputs["peo_eis"])
    supplied = {key: json.loads(inputs[f"peo_{key}"]) if isinstance(inputs[f"peo_{key}"], str)
                else inputs[f"peo_{key}"] for key, *_ in FORM_GROUPS}
    if supplied == PESKO:
        result["notes"].append("Source check: the reported salt concentration is 1.59 mol/L after rounding. The σEIS/σNE comparison combines a rounded printed conductivity with approximate Fig. 4 marker readings. The source estimated a 2% density uncertainty.")
    elif supplied == ZHENG:
        result["notes"].append("Table 3 independently tabulates a tracer fraction of 0.19 and a steady-state polarization result of 0.08. The latter is a different measurement. The 5 kg/mol specimen's density has not been imported into this 100 kg/mol record.")
    elif supplied == SYNTHETIC:
        result["warnings"].append("This worked film example uses constructed measurements.")
    else:
        result["notes"].append("User-edited record. Verify each cited observation and specimen match before interpreting the comparison as literature evidence.")
    return result
