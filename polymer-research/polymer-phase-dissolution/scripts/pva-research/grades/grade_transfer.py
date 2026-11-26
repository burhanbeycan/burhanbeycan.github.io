"""Grade-held-out benchmark for the digitized PVA/water flow curves of Parisi et al.

Fits a regularized empirical response in logarithmic c, shear rate, and (if known)
Mn and hydrolysis. These measurements are at 298.15 K. No virtual grade is a label.
"""
from pathlib import Path
import csv
import json
from collections import defaultdict
import argparse
import numpy as np

ROOT = Path(__file__).resolve().parent
KNOT = 60.0
FORMS = ('conditions_only', 'mass', 'mass_hydrolysis')

def features(row, form):
    c = np.log(float(row['concentration_wt_pct']) / 10)
    g = np.log(float(row['shear_rate_s_inv']) / KNOT)
    values = [1., c, g, c*c, c*g, g*g]
    if form != 'conditions_only':
        m = np.log(float(row['mn_g_mol']) / 50000)
        values.extend([m, c*m, g*m])
    if form == 'mass_hydrolysis':
        h = (float(row['dh_mol_pct']) - 80) / 10
        values.extend([h, c*h])
    return values

def fit(rows, form):
    key = lambda r: (r['grade_id'],r['concentration_wt_pct'])
    counts = defaultdict(int)
    for row in rows: counts[key(row)] += 1
    # Equal total loss per flow curve. Dimensionless features have fixed references.
    A = np.array([features(r,form) for r in rows])
    y = np.log([float(r['viscosity_cP']) for r in rows])
    w = np.array([1/counts[key(r)] for r in rows])
    # Fixed, never selected on held-out grades.
    return np.linalg.solve(A.T@(w[:,None]*A)+np.diag([0]+[0.15]*(A.shape[1]-1)), A.T@(w*y))

def stats(actual, predicted):
    errors = np.abs(np.array(predicted)/np.array(actual)-1)
    return dict(n=len(errors),median_APE_pct=100*float(np.median(errors)),
                MAPE_pct=100*float(errors.mean()),within10_pct=100*float(np.mean(errors<=.10)),
                MAE_log10=float(np.abs(np.log10(predicted)-np.log10(actual)).mean()))

def run(src, dest):
    records=list(csv.DictReader(src.open()))
    keys=('grade_id','mn_g_mol','dh_mol_pct','concentration_wt_pct','shear_rate_s_inv','viscosity_cP')
    data=[r for r in records if all(r.get(k) for k in keys)]
    assert all(float(r['temperature_K'])==298.15 and float(r['viscosity_cP'])>0 for r in data)
    grades=sorted({r['grade_id'] for r in data})
    if len(grades)<3:raise ValueError('At least three identified grades required.')
    held=[]; metrics={}
    for grade in grades:
        train=[r for r in data if r['grade_id']!=grade]
        test=[r for r in data if r['grade_id']==grade]
        predictions={m:np.exp(np.array([features(r,m) for r in test])@fit(train,m)) for m in FORMS}
        for i,r in enumerate(test):
            held.append({k:r[k] for k in keys}|{m:float(p[i]) for m,p in predictions.items()})
    for m in FORMS:
        groups={g:stats([float(r['viscosity_cP']) for r in held if r['grade_id']==g],
                         [r[m] for r in held if r['grade_id']==g]) for g in grades}
        metrics[m]={'groups':groups,'macro_by_grade':{k:float(np.mean([v[k] for v in groups.values()])) for k in ('MAPE_pct','within10_pct','MAE_log10')}}
    by_grade={}
    for g in grades:
        gr=[r for r in data if r['grade_id']==g]
        curves={}
        for wt in sorted({float(r['concentration_wt_pct']) for r in gr}):
            xy=sorted([(float(r['shear_rate_s_inv']),float(r['viscosity_cP'])) for r in gr if float(r['concentration_wt_pct'])==wt])
            curves[str(float(wt))]=xy
        by_grade[g]={k:float(gr[0][k]) for k in ('mn_g_mol','mw_g_mol','dh_mol_pct','dispersity')}|{'curves':curves}
    dest.mkdir(parents=True,exist_ok=True)
    with (dest/'heldout_grades.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=held[0]);w.writeheader();w.writerows(held)
    out=dict(source='10.1122/8.0000435',T_K=298.15,points=len(data),grades=by_grade,
             metrics=metrics,fitted={m:fit(data,m).tolist() for m in FORMS},regularization=0.15,
             validation='Each grade, with all of its concentrations and rates, is excluded in turn. Equal grade weights in macro metrics. Single publication.')
    (dest/'grade_state.json').write_text(json.dumps(out,indent=2))
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--csv',type=Path,default=ROOT.parent/'parisi'/'parisi_water.csv')
    parser.add_argument('--out',type=Path,default=ROOT)
    a=parser.parse_args()
    result=run(a.csv,a.out)
    print(result['points'],{k:v['macro_by_grade'] for k,v in result['metrics'].items()})
