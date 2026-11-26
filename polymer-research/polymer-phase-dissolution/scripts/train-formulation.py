"""Reproduce the HPAM formulation experiment from the source CSV.

Run python train-formulation.py [input.csv] [formulation-model.js].
Requires numpy and scipy. Values for each distinct condition are averaged
before fitting and scoring; formulation group = (polymer ppm, NaCl g/L).
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

HERE = Path(__file__).resolve().parent
REFERENCE = np.array([2000., 12.5, 30.])  # ppm, g/L, 1/s
DESIGN_STEP = np.array([1000., 7.5, 20.])
RIDGE = 0.25  # fixed before evaluation; penalty on nonlinear terms, natural log scale
FEATURES = ("c_ppm", "salt_g_L", "shear_s")


def read_runs(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return [
        dict(id=str(r["id"]), x=[float(r[k]) for k in FEATURES], y=float(r["eta_mPa_s"]),
             partition=r["partition"]) for r in rows
    ]


def cluster(runs):
    groups = defaultdict(list)
    for run in runs:
        groups[tuple(run["x"])].append(run)
    return [dict(x=list(key), y=float(np.mean([r["y"] for r in rows])),
                 ids=[r["id"] for r in rows], n=len(rows),
                 sd=float(np.std([r["y"] for r in rows], ddof=1)) if len(rows)>1 else None)
            for key, rows in sorted(groups.items())]


def features(x, nonlinear):
    z = np.log(np.asarray(x) / REFERENCE)
    u, v, w = z.T
    cols = [np.ones(len(z)), u, v, w]
    if nonlinear:
        cols.extend([u*u, v*v, w*w, u*v, u*w, v*w])
    return np.column_stack(cols)


def fit(rows, key):
    nonlinear = key != "baseline"
    x = [r["x"] for r in rows]
    design = features(x, nonlinear)
    target = np.log([r["y"] for r in rows])
    # Log-linear coefficients have no penalty. Six nonlinear features receive L2 shrinkage.
    if key == "ridge":
        penalty = np.diag([0.]*4 + [RIDGE**0.5]*6)
        design = np.vstack([design, penalty])
        target = np.r_[target, np.zeros(len(penalty))]
    return np.linalg.lstsq(design, target, rcond=None)[0]


def predict(rows, coef):
    return np.exp(features([r["x"] for r in rows], len(coef) > 4) @ coef)


def scores(points, key):
    err = np.asarray([p[key] - p["y"] for p in points])
    return {"mae": float(np.mean(abs(err))), "rmse": float(np.sqrt(np.mean(err**2)))}


MODEL_KEYS = ("baseline", "quadratic", "ridge")


def validate(conditions):
    held_points = []
    folds = []
    for group in sorted({tuple(r["x"][:2]) for r in conditions}):
        train = [r for r in conditions if tuple(r["x"][:2]) != group]
        test = [r for r in conditions if tuple(r["x"][:2]) == group]
        if not set(sum([r["ids"] for r in train], [])).isdisjoint(sum([r["ids"] for r in test], [])):
            raise ValueError("A run appears in both partitions.")
        by_model = {key: predict(test, fit(train, key)) for key in MODEL_KEYS}
        records = [dict(r, **{key: float(pred[i]) for key,pred in by_model.items()})
                   for i,r in enumerate(test)]
        held_points.extend(records)
        quadratic_rank = int(np.linalg.matrix_rank(features([r["x"] for r in train], True)))
        folds.append(dict(composition=list(group), trainConditions=len(train),
                          testConditions=len(test), quadraticRank=quadratic_rank,
                          scores={key:scores(records,key) for key in MODEL_KEYS}))
    return held_points, folds


def main():
    runs = read_runs(Path(sys.argv[1]) if len(sys.argv)>1 else HERE/"hpam-akbari.csv")
    if any(min(r["x"]) <= 0 or r["y"]<=0 for r in runs):
        raise ValueError("Log features require positive quantities.")
    conditions = cluster([r for r in runs if r["partition"]=="calibration"])
    if not conditions:
        raise ValueError("No calibration rows.")
    x = np.array([r["x"] for r in conditions])
    planes = ConvexHull((x-REFERENCE)/DESIGN_STEP).equations
    coefficients = {k: fit(conditions, k) for k in MODEL_KEYS}
    held, folds = validate(conditions)
    challenge = cluster([r for r in runs if r["partition"]=="challenge"])
    for key, coef in coefficients.items():
        if challenge:
            for r,p in zip(challenge, predict(challenge, coef)):
                r[key]=float(p)
    data = dict(reference=REFERENCE.tolist(), step=DESIGN_STEP.tolist(), lambdaRidge=RIDGE,
                axes=list(FEATURES), lower=x.min(axis=0).tolist(), upper=x.max(axis=0).tolist(),
                planes=np.unique(np.round(planes,12), axis=0).tolist(),
                coef={k:v.tolist() for k,v in coefficients.items()}, test=held, folds=folds,
                rankFull=int(np.linalg.matrix_rank(features(x,True))),
                scores={k:scores(held,k) for k in coefficients}, challenge=challenge,
                challengeScores={k:scores(challenge,k) for k in coefficients} if challenge else {},
                nGroups=len({tuple(r["x"][:2]) for r in conditions}), nRuns=sum(r["n"] for r in conditions))
    output = Path(sys.argv[2]) if len(sys.argv)>2 else HERE/"formulation-model.js"
    output.write_text("/* Generated by train-formulation.py */\n"
        "(function(r){\n  const data = " + json.dumps(data, ensure_ascii=False,indent=2,allow_nan=False) +
        ";\n  r.HPAM_MODEL=data;\n  if(typeof module!=='undefined' && module.exports) module.exports=data;\n"
        "})(typeof window!=='undefined' ? window : globalThis);\n", encoding="utf-8")
    print(len(held), "conditions,", data["nGroups"], "formulations", data["scores"], data["challengeScores"])


if __name__ == "__main__":
    main()
