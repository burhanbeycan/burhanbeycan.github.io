"""Mass, functional-group and assay accounting for waterborne polyurethane films."""

from __future__ import annotations

import csv
from io import StringIO
from math import isfinite

import numpy as np


RECIPE_COLUMNS = "component,stage,group,mass_g,basis,value,functionality"
DEFAULT_RECIPE = RECIPE_COLUMNS + """
PEG diol,prepolymer,OH,1,molar_mass,20000,2
PCL diol,prepolymer,OH,1,molar_mass,1000,2
IPDI,prepolymer,NCO,0.35,molar_mass,222.29,2
"1,4-butanediol",extension,OH,0.04506,molar_mass,90.12,2"""
EXTRACTION_COLUMNS = "sample,solvent,temperature_C,time_h,initial_dry_g,blotted_swollen_g,redried_g"
EXTRACTION_DEMO = EXTRACTION_COLUMNS + "\nexample film,toluene,23,168,0.100,0.220,0.080"
DMA_COLUMNS = "sample,frequency_Hz,temperature_C,E_storage_MPa"
DMA_DEMO = DMA_COLUMNS + "\nexample film,1,30,5.0\nexample film,1,45,3.0\nexample film,1,60,2.0"
SERIES_COLUMNS = ("formulation,specimen_family,solvent,temperature_C,time_h,protocol,hard_segment_wt_pct,"
                  "swell_initial_dry_pct,independent_E_MPa")
PAPER_SERIES = SERIES_COLUMNS + """
JPU 161,jatropha WBPU cast film,toluene,23,168,Saalah Sec 2.3.2,38,231.9,92.8
JPU 188,jatropha WBPU cast film,toluene,23,168,Saalah Sec 2.3.2,41,167.8,247.7
JPU 217,jatropha WBPU cast film,toluene,23,168,Saalah Sec 2.3.2,45,128.2,373.8"""
PAPER_DMA = "Saalah et al. §2.3.6 and Table 6: free film, tensile E′ at 25 °C, 1 Hz during 5 °C/min heating; strain amplitude unreported."
PAPER_DOI = "10.3390/polym13050795"
PAPER_SOURCE = f"Saalah et al., DOI {PAPER_DOI}, Tables 3, 4 and 6; methods §2.3.2 and §2.3.6."


def read_rows(data, columns, limit=200):
    if not isinstance(data, str) or len(data) > 50000:
        raise ValueError("CSV data must be text of at most 50000 characters.")
    try:
        reader = csv.DictReader(StringIO(data.strip()), skipinitialspace=True, strict=True)
        if reader.fieldnames != columns.split(","):
            raise ValueError(f"Expected CSV header: {columns}")
        rows = []
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"CSV row {number} has the wrong number of columns.")
            rows.append({key: value.strip() for key, value in row.items()})
            if len(rows) > limit:
                raise ValueError(f"Supply at most {limit} rows.")
    except csv.Error as error:
        raise ValueError(f"Invalid CSV: {error}") from error
    if not rows:
        raise ValueError("Supply at least one data row.")
    return rows


def numeric(value, name, low, high, *, include_low=False):
    try:
        value = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{name} must be numeric.") from error
    if not isfinite(value) or value > high or value < low or (value == low and not include_low):
        raise ValueError(f"{name} is outside its supported physical range.")
    return value


def label(row, key):
    value = row[key]
    if not value or len(value) > 160:
        raise ValueError(f"Provide {key} (1–160 characters).")
    return value


def recipe_balance(data):
    """Accumulate mmol of reactive equivalents. No conversion or topology model."""
    entries = []
    stages = {}
    groups = ("NCO", "OH", "amine", "reacted_water")
    totals = dict.fromkeys(groups, 0.0)
    for row in read_rows(data, RECIPE_COLUMNS):
        name = label(row, "component")
        stage = label(row, "stage")
        group = row["group"]
        if group not in groups:
            raise ValueError("group must be OH, NCO, amine or reacted_water.")
        mass = numeric(row["mass_g"], "mass_g", 0, 1e8)
        basis = row["basis"]
        functionality = row["functionality"]
        if basis == "molar_mass" and group != "reacted_water":
            molar_mass = numeric(row["value"], "molar mass (g/mol)", 0.1, 1e9)
            count = numeric(functionality, "reactive functionality", 0, 1e5)
            mmol = 1000 * mass * count / molar_mass
        elif functionality:
            raise ValueError("Only molar_mass basis uses the functionality column; assay and water bases include it.")
        elif basis == "OH_number" and group == "OH":
            oh_number = numeric(row["value"], "OH number (mg KOH/g)", 0, 4000)
            mmol = mass * oh_number / 56.1
        elif basis == "NCO_wt_pct" and group == "NCO":
            percentage = numeric(row["value"], "NCO wt%", 0, 100)
            mmol = 1000 * mass * percentage / (100 * 42.02)
        elif basis == "equivalent_weight" and group != "reacted_water":
            ew = numeric(row["value"], "equivalent weight (g/eq)", 1, 1e10)
            mmol = 1000 * mass / ew
        elif basis == "net_urea" and group == "reacted_water" and not row["value"]:
            mmol = 2000 * mass / 18.01528
        else:
            raise ValueError(f"Unsupported group/basis or unexpected value for {name}.")
        totals[group] += mmol
        stages.setdefault(stage, dict.fromkeys(groups, 0.0))[group] += mmol
        entries.append([name, stage, group, mass, basis, mmol])
    demand = totals["OH"] + totals["amine"] + totals["reacted_water"]
    if not totals["NCO"] or not demand:
        raise ValueError("Supply both an NCO source and at least one NCO-consuming group.")
    return {"entries": entries, "stages": stages, "totals": totals, "demand": demand,
            "balance": totals["NCO"] - demand, "index": 100 * totals["NCO"] / demand}


def extraction(data):
    rows = []
    for row in read_rows(data, EXTRACTION_COLUMNS):
        sample, solvent = label(row, "sample"), label(row, "solvent")
        temperature = numeric(row["temperature_C"], "extraction °C", -273.15, 600)
        time = numeric(row["time_h"], "extraction hours", 0, 1e6)
        initial = numeric(row["initial_dry_g"], "initial dry polymer mass (g)", 0, 1e8)
        wet = numeric(row["blotted_swollen_g"], "blotted swollen mass (g)", 0, 1e8, include_low=True)
        dry = numeric(row["redried_g"], "redried retained mass (g)", 0, 1e8, include_low=True)
        if dry > initial or dry > wet:
            raise ValueError("Redried residue cannot exceed either starting dry mass or blotted swollen mass in this balance.")
        gel = 100 * dry / initial
        retained_uptake = 100 * (wet - dry) / dry if dry else None
        rows.append([sample, solvent, temperature, time, gel, 100 - gel, 100 * (wet - initial) / initial, retained_uptake])
    return rows


def dma_sweeps(data):
    if not data.strip():
        return []
    sweeps = {}
    for row in read_rows(data, DMA_COLUMNS):
        sample = label(row, "sample")
        freq = numeric(row["frequency_Hz"], "DMA frequency (Hz)", 0, 1e8)
        temp = numeric(row["temperature_C"], "DMA °C", -273.15, 600)
        modulus = numeric(row["E_storage_MPa"], "tensile E′ (MPa)", 0, 1e6)
        points = sweeps.setdefault((sample, freq), {})
        if temp in points:
            raise ValueError("A DMA series has duplicate temperatures at the same frequency; separate the replicate specimen IDs.")
        points[temp] = modulus
    result = []
    for (sample, freq), points in sweeps.items():
        values = sorted(points.items())
        result.append({"name": f"{sample}, {freq:g} Hz", "x": [t for t, _ in values], "y": [e for _, e in values]})
    return result


def swelling_folds(data, independent_protocol=""):
    """Equal-weight holdout by formulation at a fixed protocol, using formulation means."""
    records = {}
    conditions = set()
    modulus_flags = set()
    for row in read_rows(data, SERIES_COLUMNS):
        sample = label(row, "formulation")
        conditions.add((label(row, "specimen_family"), label(row, "solvent"),
                        numeric(row["temperature_C"], "solvent °C", -273.15, 600),
                        numeric(row["time_h"], "exposure (h)", 0, 1e6), label(row, "protocol")))
        hard = numeric(row["hard_segment_wt_pct"], "hard segment wt%", 0, 100, include_low=True)
        swell = numeric(row["swell_initial_dry_pct"], "initial-dry-mass swell %", -100, 1e6, include_low=True)
        modulus_flags.add(bool(row["independent_E_MPa"]))
        modulus = numeric(row["independent_E_MPa"], "independent E′ (MPa)", 0, 1e6) if row["independent_E_MPa"] else None
        if sample in records and hard != records[sample]["hard"]:
            raise ValueError(f"Recipe hard-segment percentage changes within formulation {sample}.")
        point = records.setdefault(sample, {"hard": hard, "swell": [], "modulus": []})
        point["swell"].append(swell)
        if modulus is not None:
            point["modulus"].append(modulus)
    if len(conditions) != 1:
        raise ValueError("All formulations must share specimen family, solvent, temperature, exposure and protocol.")
    if len(modulus_flags) > 1:
        raise ValueError("For a fair paired comparison supply independent E′ for every observation or leave all blank.")
    augmented = True in modulus_flags
    if augmented and not independent_protocol.strip():
        raise ValueError("Specify the common pre-immersion DMA protocol for the independent E′ values.")
    if len(records) < 3:
        raise ValueError("At least three independent formulations are needed for recipe folds.")
    names = list(records)
    hard = np.array([v["hard"] for v in records.values()])
    target = np.array([np.mean(v["swell"]) for v in records.values()])
    design = np.column_stack([np.ones(len(records)), hard])
    if augmented:
        modulus = np.array([np.mean(v["modulus"]) for v in records.values()])
        design = np.column_stack([design, np.log(modulus)])
    models = [("Recipe (hard segment wt%)", [0, 1])]
    if augmented:
        models.append(("Independent dry-film E′", [0, 2]))
    if augmented and len(names) >= 4:
        models.append(("Recipe + independent dry-film E′", [0, 1, 2]))
    folds = []
    for index, name in enumerate(names):
        train = np.arange(len(names)) != index
        estimates = [float(np.mean(target[train]))]
        for _, columns in models:
            coefficients, _, rank, _ = np.linalg.lstsq(design[train][:, columns], target[train], rcond=None)
            if rank != len(columns):
                raise ValueError(f"Predictors lack independent variation in training fold excluding {name}.")
            estimates.append(float(design[index, columns] @ coefficients))
        in_range = np.all((design[index, 1:] >= design[train, 1:].min(axis=0)) &
                          (design[index, 1:] <= design[train, 1:].max(axis=0)))
        folds.append({"sample": name, "replicates_or_rows": len(records[name]["swell"]), "hard": hard[index],
                      "observed": target[index], "predictions": estimates, "outside_training_range": not in_range})
    predictions = np.array([v["predictions"] for v in folds])
    errors = abs(predictions - target[:, None])
    return {"folds": folds, "mae": np.mean(errors, axis=0).tolist(), "conditions": list(conditions.pop()),
            "augmented": augmented, "models": ["Training mean"] + [name for name, _ in models]}
