"""Invert PVA steady shear curves and predict shear viscosity for the documented grade.

Only columns concentration (wt%), T (K) and shear rate (1/s) enter prediction.
Requirements: numpy, scipy. CLI examples are in README.md.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
REFERENCE_RATE = 30.0  # s^-1; fixed basis knot inside the observed shear windows
REFERENCE_T = 298.15  # K; centering convention, not the user-facing target T
REFERENCE_C = 7.0  # reference numeric wt% (not mass fraction)
C_SCALE = 3.0  # wt%-points; fixed before evaluation
THETA_SCALE = 0.08  # inverse-temperature scale for local correction
LAMBDA_GRID = 0.15  # diagonal regularization of curve-parameter kernel correction
SHAPE_PENALTY = 0.025  # regularization between left/right exponents per curve
PROTOCOL = "pva-water-merck-141351-table2-v1"
METRIC_NAMES = ("mae_cp","mae_log10","mape_pct","median_ape_pct","within10_pct")
MODELS = ("direct_ridge","shared_power_law","inverse_curve_ml","inverse_parametric")


def load_data(path):
    with open(path, newline="", encoding="utf-8") as stream:
        raw = list(csv.DictReader(stream))
    xy=[]
    for r in raw:
        c, t, g, tau = (float(r[k]) for k in (
            "wt_percent", "temperature_C", "shear_rate_s_inv", "shear_stress_Pa"))
        if not all(math.isfinite(v) for v in (c,t,g,tau)) or min(c,g,tau)<=0 or t<=-273.15:
            raise ValueError("Encountered an invalid concentration, rate, temperature or stress.")
        xy.append([c,t+273.15,g,1000*tau/g])
    if len(xy)==0: raise ValueError("The measurement file has no rows.")
    return np.array(xy,dtype=float)


def theta(t):
    return 1000*(1/np.asarray(t)-1/REFERENCE_T)


def features(X):
    c,t,g=X.T
    return np.log(c/REFERENCE_C),theta(t),np.log(g/REFERENCE_RATE)


def design(X):
    # Direct empirical polynomial in c, T and gamma, with shear curvature and interactions.
    # This comparator does not constrain the sign or magnitude of local rheological slopes.
    c,t,g=features(X)
    t=t/THETA_SCALE
    return np.column_stack([np.ones(len(c)),c,t,g,c*c,g*g,c*t,c*g,t*g])


def parameter_design(C):
    c=np.log(C[:,0]/REFERENCE_C)
    t=theta(C[:,1])/THETA_SCALE
    return np.column_stack([np.ones_like(c), c, t, c*c, c*t])


def kernel(X,Y):
    # Two independent length scales; no response information in the kernel features.
    distance_c = (X[:,None,0]-Y[None,:,0])/C_SCALE
    distance_t = (theta(X[:,None,1])-theta(Y[None,:,1]))/THETA_SCALE
    return np.exp(-0.5*(distance_c**2+distance_t**2))


def bracket(values,value):
    a=sorted(set(float(i) for i in values))
    below=[x for x in a if x<=value+1e-8]
    above=[x for x in a if x>=value-1e-8]
    if not below or not above:return []
    return sorted(set((below[-1],above[0])))


def coverage(windows,c,t,g):
    cs=bracket((w[0] for w in windows),c); ts=bracket((w[1] for w in windows),t)
    by_ct={(w[0],w[1]):(w[2],w[3]) for w in windows}
    if not cs or not ts:return 'outside_concentration_or_temperature_range'
    bounds=[]
    for ci in cs:
        for ti in ts:
            if (ci,ti) not in by_ct: return 'incomplete_supporting_grid'
            bounds.append(by_ct[(ci,ti)])
    if not all(lo<=g<=hi for lo,hi in bounds):return 'outside_shared_measured_shear_window'
    return 'supported_by_neighboring_curve_windows'


def fit(X,Y):
    ln=np.log(Y)
    groups,inv,counts=np.unique(X[:,:2],axis=0,return_inverse=True,return_counts=True)
    # All measured flow curves receive equal weight in shared fits.
    w=1/np.sqrt(counts[inv]);A=design(X)
    direct=np.linalg.solve((w[:,None]*A).T@(w[:,None]*A)+np.diag([0]+[0.015]*(A.shape[1]-1)),A.T@(w*w*ln))
    cf,tt,l=features(X)
    base_design=np.column_stack([np.ones_like(l),cf,tt,-l])
    sol=lsq_linear(base_design*w[:,None],ln*w,bounds=([-np.inf,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,1]),tol=1e-11)
    base=sol.x
    params=[];rows=[];windows=[]
    for i,(ci,ti) in enumerate(groups):
        use=inv==i;li=l[use]
        # ln(eta) = a - b_low*min(ln(g/g_ref),0) - b_high*max(ln(g/g_ref),0).
        # eta>0; -1 <= d ln(eta)/d ln(g) <= 0; hence stress is nondecreasing.
        curve=np.column_stack([np.ones_like(li),-np.minimum(li,0),-np.maximum(li,0)])
        if np.linalg.matrix_rank(curve)<3:
            raise ValueError(f"Rate range does not identify both shear exponents at {ci},{ti}.")
        # Equal point weighting within a curve after division by sqrt(sample size).
        fit_design=np.vstack([curve/np.sqrt(len(li)),[0,np.sqrt(SHAPE_PENALTY),-np.sqrt(SHAPE_PENALTY)]])
        obs=np.r_[ln[use]/np.sqrt(len(li)),0]
        local=lsq_linear(fit_design,obs,bounds=([-np.inf,0,0],[np.inf,1,1]),tol=1e-11)
        a,blo,bhi=local.x
        pred_log=curve@local.x
        rows.append(dict(concentration_wt_percent=float(ci),temperature_K=float(ti),n_points=int(sum(use)),
                         eta30_cp=float(np.exp(a)),low_shear_n=float(1-blo),high_shear_n=float(1-bhi),
                         design_condition=float(np.linalg.cond(curve)),mae_log10_within_curve=float(np.mean(abs(ln[use]-pred_log))/np.log(10))))
        params.append(local.x)
        windows.append([float(ci),float(ti),float(np.min(X[use,2])),float(np.max(X[use,2]))])
    prior=np.column_stack([base[0]+base[1]*np.log(groups[:,0]/REFERENCE_C)+base[2]*theta(groups[:,1]),
                           np.full(len(groups),base[3]),np.full(len(groups),base[3])])
    # Smooth conditional inverse model: map three fitted constitutive parameters to c/T.
    # Each of the 16 independently inverted curves contributes one vector, with equal weight.
    B=parameter_design(groups)
    beta=np.linalg.solve(B.T@B+np.diag([0,0.015,0.015,0.015,0.015]),B.T@np.array(params))
    delta=np.array(params)-prior
    coeff=np.linalg.solve(kernel(groups,groups)+LAMBDA_GRID*np.eye(len(groups)),delta)
    return dict(version=1,protocol=PROTOCOL,unit='cP',columns=['wt_percent','temperature_K','rate_s_inv'],
                reference_rate_s=REFERENCE_RATE,reference_temperature_K=REFERENCE_T,
                reference_concentration_wt_percent=REFERENCE_C,concentration_kernel_scale=C_SCALE,
                inverse_temperature_kernel_scale=THETA_SCALE,kernel_regularization=LAMBDA_GRID,
                within_curve_slope_difference_regularization=SHAPE_PENALTY,
                direct_ridge=direct.tolist(),shared_power_law=base.tolist(),
                curve_conditions=groups.tolist(),correction_kernel_weights=coeff.tolist(),conditional_parameter_coefficients=beta.tolist(),
                shear_windows=windows,curve_fits=rows)


def predictions(state,X):
    X=np.atleast_2d(np.asarray(X,dtype=float))
    if X.shape[1]!=3 or not np.isfinite(X).all() or not np.all(X>0):
        raise ValueError("Prediction requires positive finite c (wt%), T (K) and shear rate (s^-1).")
    cf,tt,l=features(X);base=np.array(state['shared_power_law'])
    a=base[0]+base[1]*cf+base[2]*tt
    prior=a-base[3]*l
    corr=kernel(X[:,:2],np.array(state['curve_conditions']))@np.array(state['correction_kernel_weights'])
    low=np.clip(base[3]+corr[:,1],0,1);high=np.clip(base[3]+corr[:,2],0,1)
    final=a+corr[:,0]-low*np.minimum(l,0)-high*np.maximum(l,0)
    pars=parameter_design(X[:,:2])@np.array(state['conditional_parameter_coefficients'])
    blo=np.clip(pars[:,1],0,1);bhi=np.clip(pars[:,2],0,1)
    structured=pars[:,0]-blo*np.minimum(l,0)-bhi*np.maximum(l,0)
    return dict(direct_ridge=np.exp(design(X)@np.array(state['direct_ridge'])),
                shared_power_law=np.exp(prior),inverse_curve_ml=np.exp(final),inverse_parametric=np.exp(structured))


def predict_one(state,c,t,g,anchor_rate=None,anchor_cp=None):
    if min(c,t,g)<=0 or not all(math.isfinite(v) for v in (c,t,g)):raise ValueError('All inputs must be positive and finite.')
    domain=coverage(state['shear_windows'],c,t,g)
    if (anchor_rate is None)!=(anchor_cp is None):
        raise ValueError('Supply both anchor shear rate and anchor viscosity, or neither.')
    if anchor_rate is not None and (not all(math.isfinite(x) and x>0 for x in [anchor_rate,anchor_cp])):
        raise ValueError('Anchor shear rate and viscosity must be positive and finite.')
    eta=None
    if domain.startswith('supported'):
        eta=float(predictions(state,[[c,t,g]])['inverse_parametric'][0])
        if anchor_rate is not None:
            if not coverage(state['shear_windows'],c,t,anchor_rate).startswith('supported'):
                raise ValueError('The anchor lies outside the local measured rate window.')
            baseline_at_anchor=float(predictions(state,[[c,t,anchor_rate]])['inverse_parametric'][0])
            eta *= anchor_cp/baseline_at_anchor
    return dict(protocol=state['protocol'],model='inverse_parametric',temperature_K=t,concentration_wt_percent=c,
                rate_s_inv=g,domain=domain,eta_app_cp=eta,one_reference_used=(anchor_rate is not None))


def metrics(y,p):
    y,p=np.asarray(y),np.asarray(p)
    pct=100*np.abs(p/y-1)
    return dict(mae_cp=float(np.mean(abs(p-y))),mae_log10=float(np.mean(abs(np.log10(p/y)))),
                mape_pct=float(np.mean(pct)),median_ape_pct=float(np.median(pct)),within10_pct=float(np.mean(pct<=10)*100))


def summary(rows):
    y=[r['observed_cp'] for r in rows];p=[r['predicted_cp'] for r in rows]
    key=lambda r:(r['concentration_wt_percent'],r['temperature_K'])
    curves=sorted(set(map(key,rows)))
    means=[]
    for curve in curves:
        cur=[r for r in rows if key(r)==curve]
        means.append(metrics([r['observed_cp'] for r in cur],[r['predicted_cp'] for r in cur]))
    return dict(n_points=len(rows),n_curves=len(curves),point_weighted=metrics(y,p),
                equal_curve_mean={k:float(np.mean([m[k] for m in means])) for k in METRIC_NAMES})


def splits(X,mode):
    if mode=='whole_curve':keys=[tuple(x) for x in X[:,:2]]
    elif mode=='whole_concentration':keys=list(X[:,0])
    elif mode=='whole_temperature':keys=list(X[:,1])
    else:raise ValueError(mode)
    for label in sorted(set(keys)):
        test=np.array([v==label for v in keys])
        yield str(label).replace('np.float64(','').replace('),',',').replace('))',')'),~test,test


def evaluate(data):
    X,Y=data[:,:3],data[:,3];rows=[];folds=[];overall={}
    for mode in ('whole_concentration','whole_temperature','whole_curve'):
        overall[mode]={}
        for fold,train,test in splits(X,mode):
            state=fit(X[train],Y[train]);p=predictions(state,X[test]);xt=X[test];yt=Y[test]
            foldrows=[]
            for i,x in enumerate(xt):
                for model in MODELS:
                    foldrows.append(dict(split=mode,fold=fold,model=model,concentration_wt_percent=x[0],temperature_K=x[1],
                                     rate_s_inv=x[2],observed_cp=yt[i],predicted_cp=float(p[model][i]),
                                     domain=coverage(state['shear_windows'],*x)))
            for model in MODELS:
                tmp=[r for r in foldrows if r['model']==model]
                folds.append(dict(split=mode,fold=fold,model=model,n_points=len(tmp),
                                  **metrics([r['observed_cp'] for r in tmp],[r['predicted_cp'] for r in tmp])))
            rows.extend(foldrows)
        for model in MODELS:
            subset=[r for r in rows if r['split']==mode and r['model']==model]
            overall[mode][model]=summary(subset)
    return rows,folds,overall


def anchoring(rows):
    # The withheld observation nearest 30 s^-1 is used only for calibration.
    # Evaluate both zero-shot and one-reference models on exactly the remaining rows.
    groups={}
    for row in rows:
        key=(row['split'],row['fold'],row['model'],row['concentration_wt_percent'],row['temperature_K'])
        groups.setdefault(key,[]).append(row)
    results=[];summary_results={}
    for items in groups.values():
        anchor=min(items,key=lambda row:(abs(math.log(row['rate_s_inv']/REFERENCE_RATE)),row['rate_s_inv']))
        correction=anchor['observed_cp']/anchor['predicted_cp']
        for row in items:
            if row is anchor:continue
            results.append(dict(**row,anchor_rate_s_inv=anchor['rate_s_inv'],anchor_cp=anchor['observed_cp'],
                                predicted_cp_with_one_reference=row['predicted_cp']*correction))
    for split in sorted({r['split'] for r in results}):
        summary_results[split]={}
        for model in MODELS:
            matched=[r for r in results if r['split']==split and r['model']==model]
            updated=[dict(r,predicted_cp=r['predicted_cp_with_one_reference']) for r in matched]
            summary_results[split][model]=dict(without_reference=summary(matched),with_one_reference=summary(updated))
    return results,summary_results


def writecsv(path,rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='task',required=True)
    b=sub.add_parser('build');b.add_argument('--data',type=Path,default=HERE/'pva_table2.csv');b.add_argument('--out',type=Path,default=HERE)
    p=sub.add_parser('predict');p.add_argument('--state',type=Path,default=HERE/'model_state.json')
    p.add_argument('--wt',type=float,required=True);p.add_argument('--kelvin',type=float,default=298.0);p.add_argument('--rate',type=float,required=True)
    p.add_argument('--anchor-rate',type=float);p.add_argument('--anchor-cp',type=float)
    args=parser.parse_args()
    if args.task=='predict':
        state=json.loads(args.state.read_text())
        print(json.dumps(predict_one(state,args.wt,args.kelvin,args.rate,args.anchor_rate,args.anchor_cp),indent=2));return
    data=load_data(args.data);state=fit(data[:,:3],data[:,3]);rows,folds,overall=evaluate(data)
    anchored,reference_metrics=anchoring(rows)
    args.out.mkdir(parents=True,exist_ok=True)
    for name,obj in [('model_state.json',state),('benchmark.json',overall),('one_reference_benchmark.json',reference_metrics)]:
        (args.out/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    for name,vals in [('heldout_predictions.csv',rows),('fold_metrics.csv',folds),('one_reference_predictions.csv',anchored),('inverse_parameters.csv',state['curve_fits'])]:
        writecsv(args.out/name,vals)

if __name__=='__main__':main()
