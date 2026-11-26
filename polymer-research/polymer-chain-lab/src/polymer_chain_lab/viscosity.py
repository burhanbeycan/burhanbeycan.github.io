"""Condition-specific dilute-solution viscosity for linear polymer chains."""

import csv
from dataclasses import dataclass
from importlib.resources import files
import io
import math

import numpy as np

from .model import finite_number

DOI = "10.1039/C8PY00255J"
CSV_FIELDS = ("sample_id", "mass_class", "Mw_g_mol", "eta_mL_g", "medium",
              "viscosity_temp_C", "sec_eluent", "sec_temp_C", "source")
PBS = "PBS (Sigma-Aldrich tablet per 200 mL)"
SEC = "DMAc + 50 mM LiCl"


@dataclass(frozen=True)
class Observation:
    sample: str
    group: str
    mw: float
    eta: float
    medium: str
    temperature: float
    sec_eluent: str
    sec_temperature: float
    source: str


def literature_series():
    """Experimental SI Tables S1/S2, joined on sample ID; absent cells stay absent."""
    path = files("polymer_chain_lab").joinpath("data").joinpath("petox_pbs.csv")
    observations = []
    with path.open(encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            for temperature in (15, 25, 37, 45, 55):
                eta = row[f"eta_{temperature}C"]
                if eta:
                    observations.append(Observation(
                        row["sample_id"], row["nominal_mass_kDa"] + " kDa nominal",
                        1000 * float(row["Mw_kDa"]), float(eta), PBS, temperature,
                        SEC, 50, DOI + ", SI Tables S1/S2"))
    return observations


def read_observations(raw):
    if not isinstance(raw, str) or len(raw) > 100_000:
        raise ValueError("The measurement CSV must be text under 100 kB.")
    reader = csv.DictReader(io.StringIO(raw.strip()), skipinitialspace=True, strict=True)
    observations = []
    try:
        if tuple(reader.fieldnames or ()) != CSV_FIELDS:
            raise ValueError("Measurement CSV header must be: " + ",".join(CSV_FIELDS))
        for row in reader:
            if len(observations) >= 256 or None in row or any(
                    v is None or not v.strip() or len(v) > 500 for v in row.values()):
                raise ValueError("Supply 1–256 complete measurement rows with bounded text.")
            row = {k: v.strip() for k, v in row.items()}
            observations.append(Observation(
                row["sample_id"], row["mass_class"],
                finite_number(row["Mw_g_mol"], "Mw (g/mol)", 1, 1e10),
                finite_number(row["eta_mL_g"], "Intrinsic viscosity (mL/g)", 1e-5, 1e6),
                row["medium"], finite_number(row["viscosity_temp_C"], "Viscosity temperature (°C)", -273.14, 600),
                row["sec_eluent"], finite_number(row["sec_temp_C"], "SEC temperature (°C)", -273.14, 600),
                row["source"]))
    except csv.Error as exc:
        raise ValueError("Invalid measurement CSV quoting.") from exc
    return observations


def power_fit(rows):
    """Unweighted least squares on ln(eta) versus ln(Mw) with centered abscissae."""
    x = np.log([r.mw for r in rows])
    y = np.log([r.eta for r in rows])
    centered = x - x.mean()
    denominator = np.dot(centered, centered)
    if denominator < 1e-12:
        raise ValueError("Each training fold needs at least two distinct, resolvable molar masses.")
    a = float(np.dot(centered, y - y.mean()) / denominator)
    if not 0 < a <= 2:
        raise ValueError("The measured series or a training fold has an unsupported exponent (expected 0 < a ≤ 2 for this model).")
    return float(y.mean() - a * x.mean()), a


def benchmark(observations, temperature):
    selected = sorted((row for row in observations if math.isclose(
        row.temperature, temperature, rel_tol=0, abs_tol=1e-6)), key=lambda row: row.mw)
    groups = {row.group for row in selected}
    if len(groups) < 3:
        raise ValueError("At least three complete mass classes are required at the chosen viscometry temperature.")
    if len({(r.medium, r.sec_eluent, r.sec_temperature) for r in selected}) != 1:
        raise ValueError("Use one viscometry medium, SEC eluent and SEC temperature per fit.")
    if len({r.sample for r in selected}) != len(selected):
        raise ValueError("Supply one measured value per unique sample ID at each temperature; aggregate technical replicates first.")
    mass_groups = {}
    for r in selected:
        if mass_groups.setdefault(r.mw, r.group) != r.group:
            raise ValueError("Identical measured masses cannot be split across mass classes.")
    index_predictions = {}
    for group in sorted(groups):
        training = [r for r in selected if r.group != group]
        log_k, a = power_fit(training)
        low, high = min(r.mw for r in training), max(r.mw for r in training)
        for row in selected:
            if row.group == group:
                pred = math.exp(log_k + a * math.log(row.mw))
                index_predictions[row.sample] = (pred, 100 * abs(pred / row.eta - 1),
                                                 "inside" if low <= row.mw <= high else "outside")
    log_k, a = power_fit(selected)
    return selected, index_predictions, math.exp(log_k), a


def mixture(distribution, k, a, lower, upper):
    """Intrinsic viscosity of a dilute noninteracting same-chemistry mixture."""
    k = finite_number(k, "Component-law K", 1e-12, 1e5)
    a = finite_number(a, "Component-law a", 0.001, 2)
    lower = finite_number(lower, "Lower calibration mass", 1, 1e12)
    upper = finite_number(upper, "Upper calibration mass", 1, 1e12)
    if upper <= lower:
        raise ValueError("Upper calibration mass must exceed lower mass.")
    mass, weight = distribution.molar_mass_g_mol, distribution.mass_fraction
    outside = ((mass < lower) | (mass > upper)) & (weight > 0)
    if any(outside):
        raise ValueError("A populated MWD bin lies outside the component-law calibration interval.")
    component_eta = k * mass**a
    intrinsic = float(np.dot(component_eta, weight))
    viscosity_mass = math.exp(math.log(intrinsic / k) / a)
    return intrinsic, viscosity_mass, component_eta
