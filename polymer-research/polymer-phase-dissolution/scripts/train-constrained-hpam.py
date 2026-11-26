"""Fit a shape-constrained HPAM response surface and recompute each composition holdout."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, minimize
from scipy.spatial import ConvexHull

REFERENCE = np.array([2000.0, 12.5, 30.0])
NAMES = ('c_ppm', 'salt_g_L', 'shear_s')
DOI = 'https://doi.org/10.1007/s13202-017-0375-6'


def design(z):
    u, v, w = z.T
    return np.column_stack([np.ones(len(z)), u, v, w, u*u, v*v, w*w, u*v, u*w, v*w])


def derivative_constraints(points):
    rows = []
    for u, v, w in points:
        rows.extend([[0, 1, 0, 0, 2*u, 0, 0, v, w, 0],
                     [0, 0, -1, 0, 0, -2*v, 0, -u, 0, -w],
                     [0, 0, 0, -1, 0, 0, -2*w, 0, -u, -v]])
    return np.asarray(rows)


def fit(x, y):
    z = np.log(x / REFERENCE)
    matrix = design(z)
    target = np.log(y)
    if np.linalg.matrix_rank(matrix) != matrix.shape[1]:
        raise ValueError('The training design is rank deficient.')
    coef = np.linalg.lstsq(matrix, target, rcond=None)[0]
    hull = ConvexHull(z)
    limits = derivative_constraints(z[hull.vertices])
    def loss(beta):
        return 0.5 * np.sum((matrix @ beta - target)**2)
    def grad(beta):
        return matrix.T @ (matrix @ beta - target)
    solution = minimize(loss, coef, jac=grad, method='SLSQP',
                        constraints=LinearConstraint(limits, 0, np.inf),
                        options={'ftol': 1e-13, 'maxiter': 600})
    constrained = solution.x
    if not solution.success or np.min(limits @ constrained) < -1e-7:
        raise RuntimeError(solution.message)
    return coef, constrained, hull


def load(path):
    buckets = defaultdict(list)
    with path.open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            if row['partition'] == 'calibration':
                buckets[tuple(float(row[k]) for k in NAMES)].append(float(row['eta_mPa_s']))
    return np.array(sorted(buckets)), np.array([np.mean(buckets[key]) for key in sorted(buckets)])


def benchmark(x, y):
    rows = []
    for composition in sorted(set(map(tuple, x[:, :2]))):
        hold = np.all(x[:, :2] == composition, axis=1)
        plain, constrained, hull = fit(x[~hold], y[~hold])
        for target, measured in zip(x[hold], y[hold]):
            z = np.log(target / REFERENCE)
            a = design(z.reshape(1,3))[0]
            in_hull = bool(np.all(hull.equations[:,:-1] @ z + hull.equations[:,-1] <= 1e-8))
            c, s, g = target
            rows.append(dict(id='/'.join(f'{v:g}' for v in target), group=f'{c:g} ppm · {s:g} g/L NaCl',
                             label=f'{g:g} s⁻¹', observed=float(measured), baseline=float(np.exp(a @ plain)),
                             candidate=float(np.exp(a @ constrained)), sourceUrl=DOI, sourceLabel='Akbari ve ark.',
                             condition=f'24 °C · {g:g} s⁻¹ · '+('eğitim koşullarının içinde' if in_hull else 'eğitim koşullarının dışında')))
    _, beta, hull = fit(x,y)
    return dict(rows=rows,coef=beta.tolist(),reference=REFERENCE.tolist(),logPlanes=hull.equations.tolist())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    x,y=load(args.csv)
    result=benchmark(x,y)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    for key in ('baseline','candidate'):
        error=np.array([r[key]-r['observed'] for r in result['rows']])
        print(key,len(error),np.sqrt(np.mean(error**2)),np.mean(abs(error)))
if __name__=='__main__':
    main()
