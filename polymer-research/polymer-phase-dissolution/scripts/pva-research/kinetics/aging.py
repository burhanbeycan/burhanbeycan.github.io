"""Benchmarks of a published early-aging law for nearly fully hydrolyzed PVA/water.

Standard library only. Run as a script to regenerate the numerical results.
Time-zero falling-ball viscosity is a known conditioning measurement in the time-series test.
"""
import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_rows(filename, numeric):
    with (ROOT / filename).open() as file:
        rows = list(csv.DictReader(file))
    for row in rows:
        for k in numeric:
            row[k] = float(row[k])
    return rows


def group_log_coefficient(rows, exponent=2):
    """Estimate log(K) with recipe/batch groups weighted equally."""
    vals = {}
    for r in rows:
        value = math.log(r["alpha_per_hour"]) - exponent * math.log(r["wt_percent"])
        vals.setdefault(r["sample_group"], []).append(value)
    return math.exp(statistics.mean(statistics.mean(a) for a in vals.values()))


def apparent_viscosity(initial_eta_cp, age_h, wt_percent, K, *, storage_T_K=303.15,
                       residual_acetate_max_mol_percent=0.1, polymerization_T_K=333.15):
    """Conditional 0–16 h research prediction for water, same falling-ball protocol."""
    if not all(math.isfinite(v) for v in (initial_eta_cp, age_h, wt_percent, K,
                                            storage_T_K, residual_acetate_max_mol_percent, polymerization_T_K)):
        raise ValueError("A finite measured baseline and known protocol are required")
    valid_conditions = (abs(storage_T_K - 303.15) < 1e-7
                        and abs(polymerization_T_K - 333.15) < 1e-7
                        and 0 <= residual_acetate_max_mol_percent <= 0.1
                        and 14.63 <= wt_percent <= 16.04 and 0 <= age_h <= 16)
    if initial_eta_cp <= 0 or K <= 0 or not valid_conditions:
        raise ValueError("Outside early PVA-water aging benchmark protocol")
    return initial_eta_cp * (1 + K * wt_percent**2 * age_h)


def score(predictions, reference):
    errors = [100 * abs(p / y - 1) for p, y in zip(predictions, reference)]
    return dict(n=len(errors), mean_APE_pct=statistics.mean(errors),
                median_APE_pct=statistics.median(errors),
                within_10_pct_count=sum(e <= 10 for e in errors))


def main():
    coefficients = read_rows("early_rate_table.csv", ["wt_percent", "alpha_per_hour"])
    time_rows = read_rows("aging_fig1.csv", ["wt_percent", "age_h", "apparent_eta_cP"])
    # The blend has incompletely linked precursor fractions. Retain in the data but
    # exclude from both the fitting and the score to avoid shared-fraction leakage.
    records = [r for r in coefficients if r["sample_group"] != "low_blend"]
    groups = sorted({r["sample_group"] for r in records})
    rate_predictions = []
    for group in groups:
        train = [r for r in records if r["sample_group"] != group]
        K = group_log_coefficient(train)
        constant_rate = group_log_coefficient(train, 0)
        for r in records:
            if r["sample_group"] == group:
                rate_predictions.append(dict(group=group, material=r["material_id"],
                                            wt_percent=r["wt_percent"], observed=r["alpha_per_hour"],
                                            c_squared=K * r["wt_percent"]**2, constant=constant_rate))
    # Both illustrated samples withheld when estimating K; initial viscosity from each
    # illustrated curve remains necessary, since this law predicts subsequent change.
    training = [r for r in records if r["material_id"] not in {"L1490", "L1705"}]
    K = group_log_coefficient(training)
    initial = {p["material_id"]: p["apparent_eta_cP"] for p in time_rows if p["age_h"] == 0}
    pairs = []
    for p in time_rows:
        if p["age_h"] == 0:
            continue
        eta0 = initial[p["material_id"]]
        pairs.append(dict(material=p["material_id"], age_h=p["age_h"], wt_percent=p["wt_percent"],
                          observed_cP=p["apparent_eta_cP"], fresh_cP=eta0,
                          regime="early" if p["age_h"] <= 16 else "late_diagnostic",
                          early_law_cP=eta0 * (1 + K * p["wt_percent"]**2 * p["age_h"])))
    summary = {}
    for band in ("early", "late_diagnostic"):
        part = [p for p in pairs if p["regime"] == band]
        summary[band] = {column: score([p[column] for p in part], [p["observed_cP"] for p in part])
                         for column in ("fresh_cP", "early_law_cP")}
    group_scores = {}
    for column in ("constant","c_squared"):
        group_scores[column] = statistics.mean(score([p[column] for p in rate_predictions if p["group"] == g],
                                                      [p["observed"] for p in rate_predictions if p["group"] == g])
                                              ["mean_APE_pct"] for g in groups)
    results = dict(source_doi="10.1295/koron1944.15.348",
                   title="Early aging of aqueous PVA under the original falling-ball protocol",
                   n_digitized_fig1=len(time_rows), n_table1_derived_rates=len(coefficients),
                   n_parent_groups=len(groups), K_held_two_materials=K,
                   K_unit="h^-1 (numerical wt-percent)^-2", transfer_uses_time_zero_reference=True,
                   valid_early_range_h=[0,16], time_benchmark=summary,
                   rate_group_mean_MAPE_pct=group_scores, rate_predictions=rate_predictions,
                   time_predictions=pairs)
    (ROOT / "aging_results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({k:results[k] for k in ["n_digitized_fig1","n_table1_derived_rates","n_parent_groups",
                                             "K_held_two_materials","time_benchmark","rate_group_mean_MAPE_pct"]},indent=2))
    for r in pairs:
        if r["regime"] == "early":
            assert abs(apparent_viscosity(r["fresh_cP"],r["age_h"],r["wt_percent"],K)-r["early_law_cP"]) < 1e-7
    for t in [0, 10, 16]:
        assert apparent_viscosity(1250,t,15,K) >= 1250
    try:
        apparent_viscosity(1250,4,15,K,storage_T_K=298.)
    except ValueError:
        pass
    else:
        raise AssertionError("298 K lacks a calibrated PVA aging rate")

if __name__=="__main__":
    main()
