"""Dry binary PEO/LiTFSI composition and matched-state electrochemical observations."""
from __future__ import annotations

import json
import math

from .model import finite, nernst_einstein

# Precision used by Pesko et al., JES 164, E3569 (2017), equations 5 and 6.
EO_G_MOL = 44.05
LITFSI_G_MOL = 287.09


def object_record(raw, label, required, optional=()):
    if isinstance(raw, str):
        if len(raw) > 12000:
            raise ValueError(f"{label} is too large.")
        try:
            raw = json.loads(raw)
        except (ValueError, RecursionError):
            raise ValueError(f"{label} must contain valid JSON.") from None
    if not isinstance(raw, dict) or not set(required) <= raw.keys():
        raise ValueError(f"{label} requires these entries: {', '.join(required)}.")
    unknown = set(raw) - set(required) - set(optional)
    if unknown:
        raise ValueError(f"{label} contains unknown entries: {', '.join(map(str, unknown))}.")
    return raw


def label(value, field):
    if not isinstance(value, str) or not value.strip() or len(value) > 1000:
        raise ValueError(f"{field} requires non-empty text up to 1000 characters.")
    return value.strip()


def number(value, field):
    return finite(value, field, lower=1e-18, upper=1e14)


def match_temperature(value, state_T, method):
    temperature = finite(value, f"{method} temperature (K)", lower=1, upper=5000)
    if not math.isclose(temperature, state_T, rel_tol=0, abs_tol=1e-6):
        raise ValueError(f"{method} temperature differs from the material state. Use matched measurements.")


def formulation(recipe, temperature):
    recipe = object_record(
        recipe, "Dry binary formulation", ["peo_g", "source", "state_id"],
        ["salt_g", "eo_per_li", "eo_mass_fraction_in_peo", "volume_at_state"],
    )
    peo = number(recipe["peo_g"], "Dry PEO mass (g)")
    label(recipe["source"], "Formulation source")
    label(recipe["state_id"], "Formulation state_id")
    eo_fraction = finite(recipe.get("eo_mass_fraction_in_peo", 1.0),
                         "Repeat mass fraction of dry PEO", lower=0, strict=True, upper=1)
    eo_moles = peo * eo_fraction / EO_G_MOL
    if ("salt_g" in recipe) == ("eo_per_li" in recipe):
        raise ValueError("Supply exactly one of salt_g and eo_per_li.")
    if "salt_g" in recipe:
        salt = number(recipe["salt_g"], "Anhydrous LiTFSI mass (g)")
        moles = salt / LITFSI_G_MOL
        eo_per_li = eo_moles / moles
        basis = "Masses of dry components"
    else:
        eo_per_li = number(recipe["eo_per_li"], "EO units per lithium")
        moles = eo_moles / eo_per_li
        salt = moles * LITFSI_G_MOL
        basis = "Target EO:Li; LiTFSI mass calculated on the selected PEO basis"

    concentration = None
    volume_basis = "Unknown"
    volume_record = recipe.get("volume_at_state")
    if volume_record is not None:
        volume_record = object_record(
            volume_record, "Volume evidence", ["temperature_K", "source", "state_id"],
            ["density_g_mL", "final_volume_mL", "salt_mol_L"],
        )
        if volume_record["state_id"] != recipe["state_id"]:
            raise ValueError("Volume or concentration state_id differs from the dry formulation.")
        match_temperature(volume_record["temperature_K"], temperature, "Volume or concentration")
        label(volume_record["source"], "Volume source")
        quantities = {"density_g_mL", "final_volume_mL", "salt_mol_L"} & volume_record.keys()
        if len(quantities) != 1:
            raise ValueError("Supply one mixture density, final mixture volume or independent molarity.")
        key = quantities.pop()
        value = number(volume_record[key], key)
        if key == "density_g_mL":
            concentration = 1000 * moles * value / (peo + salt)
        elif key == "final_volume_mL":
            concentration = 1000 * moles / value
        else:
            concentration = value
        concentration = finite(concentration, "Salt molarity (mol/L)", lower=1e-9, upper=100)
        volume_basis = f"{key} = {value:g}; {volume_record['source']}"
    return {"state_id": recipe["state_id"], "peo_g": peo, "salt_g": salt, "salt_wt_pct": 100 * salt / (salt + peo),
            "molality": 1000 * moles / peo, "eo_per_li": eo_per_li, "concentration": concentration,
            "basis": basis, "volume_basis": volume_basis, "source": recipe["source"],
            "eo_fraction": eo_fraction}


def analyze_peo(state, recipe, tracer, eis):
    state = object_record(state, "Sample state", [
        "state_id", "temperature_K", "peo_mass_kg_mol", "mass_average", "phase",
        "thermal_history", "crystallinity", "source",
    ])
    for key in ("state_id", "mass_average", "phase", "thermal_history", "source"):
        label(state[key], key)
    mass = number(state["peo_mass_kg_mol"], "PEO molar mass (kg/mol)")
    temperature = finite(state["temperature_K"], "Sample temperature (K)", lower=1, upper=5000)
    xc = state["crystallinity"]
    xc_description = "Unreported"
    if xc is not None:
        xc = object_record(xc, "Crystallinity", ["fraction", "basis", "method"])
        xc_value = finite(xc["fraction"], "Measured crystalline mass fraction", lower=0, upper=1)
        if xc["basis"] not in ("PEO mass", "electrolyte mass"):
            raise ValueError("Set crystallinity basis to PEO mass or electrolyte mass.")
        xc_description = f"{100*xc_value:g}% of {xc['basis']}; {label(xc['method'], 'Crystallinity method')}"
    composition = formulation(recipe, temperature)
    if composition["state_id"] != state["state_id"]:
        raise ValueError("Formulation state_id differs from the material-state identifier.")

    observation_fields = ["state_id", "temperature_K", "source"]
    tracer = object_record(tracer, "Tracer record", observation_fields + ["D_Li_m2_s", "D_TFSI_m2_s"],
                           ["relative_plot_reading_bound"])
    eis = object_record(eis, "Impedance record", observation_fields,
                        ["sigma_mS_cm", "thickness_um", "area_cm2", "bulk_resistance_ohm"])
    for title, entry in (("PFG-NMR", tracer), ("EIS", eis)):
        if entry["state_id"] != state["state_id"]:
            raise ValueError(f"{title} state_id differs from the material-state identifier.")
        match_temperature(entry["temperature_K"], temperature, title)
        label(entry["source"], f"{title} source")
    dp = finite(tracer["D_Li_m2_s"], "Lithium tracer D", lower=1e-30, upper=1e-4)
    dm = finite(tracer["D_TFSI_m2_s"], "TFSI tracer D", lower=1e-30, upper=1e-4)
    geometry_keys = {"thickness_um", "area_cm2", "bulk_resistance_ohm"}
    if "sigma_mS_cm" in eis and not (geometry_keys & eis.keys()):
        sigma_eis = number(eis["sigma_mS_cm"], "EIS conductivity (mS/cm)")
        eis_basis = "Supplied bulk conductivity (mS/cm)"
    elif "sigma_mS_cm" not in eis and geometry_keys <= eis.keys():
        length_cm = number(eis["thickness_um"], "Thickness (µm)") * 1e-4
        area = number(eis["area_cm2"], "Active electrode area (cm²)")
        rb = number(eis["bulk_resistance_ohm"], "Bulk resistance (Ω)")
        sigma_eis = 1000 * length_cm / (area * rb)
        eis_basis = f"l / (A Rb); {eis['thickness_um']} µm, {area:g} cm², {rb:g} Ω"
    else:
        raise ValueError("Supply either sigma_mS_cm or all of thickness_um, area_cm2 and bulk_resistance_ohm.")
    c = composition["concentration"]
    sigma_ne = None if c is None else 10 * nernst_einstein(c * 1000, dp, dm, temperature)["conductivity_S_m"]
    ratio = None if sigma_ne is None else sigma_eis / sigma_ne

    def metric(key, title, value, unit, note):
        return dict(key=key, label=title, value=value, unit=unit, note=note)

    metrics = [
        metric("eo_per_li", "Ether oxygens per lithium", composition["eo_per_li"], "EO/Li", "Dry PEO repeat-unit basis"),
        metric("sigma_eis", "Bulk EIS conductivity", sigma_eis, "mS/cm", eis_basis),
        metric("t_tracer", "Tracer mobility fraction", dp / (dp + dm), "1", "⁷Li and ¹⁹F self-diffusion"),
    ]
    warnings = []
    if c is None:
        warnings.append("Supply the final mixture volume, mixture density or molarity at this temperature to compare conductivities. EO:Li alone does not determine molarity.")
    else:
        metrics += [
            metric("salt_molarity", "LiTFSI concentration", c, "mol/L", composition["volume_basis"]),
            metric("sigma_ne", "Tracer Nernst–Einstein baseline", sigma_ne, "mS/cm", "All LiTFSI formula units counted"),
            metric("sigma_ratio", "σEIS / σNE", ratio, "1", "Collective charge transport relative to a tracer baseline"),
        ]
    tables = [
        {"title":"Dry formulation ledger", "columns":["Quantity", "Value", "Unit or basis"], "rows":[
            ["PEO", composition["peo_g"], "g dry polymer"], ["LiTFSI", composition["salt_g"], "g dry salt"],
            ["Salt mass fraction", composition["salt_wt_pct"], "wt% of dry electrolyte"],
            ["Molality", composition["molality"], "mol/kg dry PEO"],
            ["EO repeat mass fraction", composition["eo_fraction"], "of dry PEO, excluding ends if specified"],
            ["Basis", composition["basis"], composition["source"]], ["Volume evidence", composition["volume_basis"], "at sample T"]]},
        {"title":"Matched material state", "columns":["Quantity","Value"], "rows":[
            ["State ID",state["state_id"]], ["PEO molar mass (kg/mol)",mass], ["Average as reported",state["mass_average"]],
            ["T (K)",temperature], ["Physical state",state["phase"]], ["Thermal history",state["thermal_history"]],
            ["Crystallinity",xc_description], ["State source",state["source"]]]},
        {"title":"Transport observations", "columns":["Quantity","Value","Unit","Source"], "rows":[
            ["⁷Li tracer D",dp,"m²/s",tracer["source"]], ["¹⁹F tracer D",dm,"m²/s",tracer["source"]],
            ["Bulk σ",sigma_eis,"mS/cm",eis["source"]]]},
    ]
    if "relative_plot_reading_bound" in tracer:
        delta = finite(tracer["relative_plot_reading_bound"], "Plot reading relative bound",
                       lower=0, upper=0.99)
        warnings.append(f"Tracer diffusivities were estimated from plotted markers (±{delta:.0%} reading allowance). Report the comparison at the precision of the source.")
        if ratio is not None:
            tables.append({"title": "Sensitivity to the stated tracer reading allowance",
                           "columns": ["Quantity", "Lower", "Upper", "Unit"],
                           "rows": [["σNE", sigma_ne * (1 - delta), sigma_ne * (1 + delta), "mS/cm"],
                                    ["σEIS / σNE", ratio / (1 + delta), ratio / (1 - delta), "1"]]})
            warnings.append("The sensitivity interval varies only the two tracer readings. It is not an experimental confidence interval.")
    return {"metrics": metrics, "plots": [], "tables": tables, "warnings": warnings,
            "notes": ["PFG-NMR measures self motion in the stationary bulk sample. DLi/(DLi+DTFSI) is the independent-ion mobility fraction. A polymer-fixed collective transference number requires independent transport information.",
                      "The conductivity ratio includes like- and unlike-ion correlations. It can exceed unity and cannot identify a free-ion fraction, an ion-pair population or a particular correlation term.",
                      "EO counts use dry PEO repeat mass. The default fraction of 1 uses the long-chain approximation. For oligomers with known ends, supply the independently established EO repeat mass fraction.",
                      "Compare the same composition, temperature, phase and thermal history (replicate specimens may share one state ID). PEO molar mass and crystalline fraction are recorded as evidence; the lab does not predict either D from these descriptors."],
            "material_updates": {}, "handoffs": []}
