"""Refit one binary Flory–Huggins potential to cloud and dilute osmometry records.

Usage: python scripts/fit_transfer.py --cloud dist/ps-cloud-points.csv --osmotic
dist/transfer-inputs/osmotic.json --coexist dist/transfer-inputs/coexist.json
--out dist/transfer-results.js. Requires SciPy and NumPy. Output is deterministic.
The existing fit_ps.binodal_chi0 solves chi_0 for a measured cloud composition.
"""
import argparse
import csv
import json
from pathlib import Path
import math

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar

from fit_ps import binodal_chi0, VBAR_PS, VREF

G = VBAR_PS**2 / VREF
SCALE_T = 1.0
SCALE_A2 = 1e-5
# Each observable contributes one group mean. Objective scales are algorithm settings.
WEIGHTS = [0, 0.25, 1, 4]
TREF = 300.0
CONSTANTS = dict(vbar_ml_g=VBAR_PS, V_s_ml_mol=VREF,
                 scale_T_K=SCALE_T, scale_A2_mol_ml_g2=SCALE_A2)

def predicted_a2(T, a,b,c):
    return G*(0.5 - a - b/np.asarray(T) + c)

def read_cloud(path):
    with path.open() as fh:
        rows=list(csv.DictReader(fh))
    return [dict(T_K=float(r['T_cloud_observed_K']),phi_p=float(r['phi_polymer']),
                 N=float(r['N_lattice_sites']),Mw_g_mol=float(r['Mw_g_mol']),
                 source=r['source_doi'],series=r['series'],split=r['split']) for r in rows]

def fit(cloud, osm, weight):
    """Cloud temperatures with only selected osmometry as fitting observations."""
    pts=[p for p in cloud if p['split']=='calibration']
    if not pts or (weight and not osm):
        raise ValueError('Missing calibration data')
    To=np.array([p['T_K'] for p in pts])
    Ta=np.array([p['T_K'] for p in osm])
    Ya=np.array([p['A2_mol_ml_g2'] for p in osm])
    def profile(c):
        roots=np.array([binodal_chi0(p['phi_p'],p['N'],c) for p in pts])
        b0,a0=np.polyfit(1/To,roots,1)
        def residual(v):
            alpha,beta=v
            tcp=1/((roots-alpha)/beta+1/TREF)
            eT=(tcp-To)/(math.sqrt(len(To))*SCALE_T)
            if not weight: return eT
            est=G*(.5-alpha-beta*(1/Ta-1/TREF)+c)
            eA2=math.sqrt(weight/len(Ya))*(est-Ya)/SCALE_A2
            return np.r_[eT,eA2]
        rr=least_squares(residual,[a0+b0/TREF,b0],ftol=1e-13,gtol=1e-13,xtol=1e-13)
        return (float(rr.x[0]-rr.x[1]/TREF),float(rr.x[1]),float(np.sum(rr.fun**2)))
    optimum=minimize_scalar(lambda c:profile(c)[2],bounds=(-.10,.24),method='bounded',options={'xatol':1e-10})
    a,b,_=profile(optimum.x)
    return (a,b,float(optimum.x))

def merge_tangent(n,t,a,b,c):
    """Coexistence from equal slope and equal solvent chemical potential.

    Solve on each stable branch, using log(phi) on the dilute side. A missed
    measured pair is returned as one phase and retained in coverage/width loss.
    """
    c0=a+b/t
    def curv(p):
        return 1/(n*p)+1/(1-p)-2*c0+2*c*(1-3*p)
    def third(p):
        return -1/(n*p*p)+1/(1-p)**2-6*c
    pc=brentq(third,1e-12,1-1e-12)
    if curv(pc)>=-1e-13:
        return dict(phases=1,phi_lean=None,phi_rich=None)
    left=brentq(curv,1e-14,pc,xtol=1e-15)
    right=brentq(curv,pc,1-1e-13,xtol=1e-15)
    def slope(y):
        p=math.exp(y)
        return (y+1)/n-math.log1p(-p)-1+c0*(1-2*p)+c*p*(2-3*p)
    def solvent(p):
        dlog = -p*p*(.5+p*(1/3+p*(.25+p/5))) if p<1e-3 else math.log1p(-p)+p
        return dlog-p/n+p*p*(c0-c+2*c*p)
    yr, yl=math.log(right),math.log(left)
    m0,m1=slope(yr),slope(yl)
    def sides(m):
        if m>=m1: ly=yl
        else: ly=brentq(lambda y:slope(y)-m,min(-42.,n*(m+1-c0)-7.,yl-5),yl,xtol=1e-14)
        if m<=m0: ry=yr
        else: ry=brentq(lambda y:slope(y)-m,yr,-1e-14,xtol=1e-14)
        return math.exp(ly),math.exp(ry)
    def residual(m):
        p,q=sides(m)
        return solvent(p)-solvent(q)
    p,q=sides(brentq(residual,m0,m1,xtol=5e-15))
    return dict(phases=2,phi_lean=p,phi_rich=q)

def error_summary(residuals, unit):
    d=np.asarray(residuals)
    if not len(d): return dict(n=0)
    return dict(n=len(d),**{'rmse_'+unit:float(np.mean(d**2)**.5),
                            'mae_'+unit:float(np.mean(abs(d))),
                            'bias_'+unit:float(np.mean(d))})
def endpoints_summary(records, predictions):
    if not records: return dict(n_total=0,n_two_phase=0,n_missed=0,
                                mae_lean=None,mae_rich=None,rmse_pooled=None,mae_gap=None)
    lean=[]; rich=[]; widths=[]
    for r,p in zip(records,predictions,strict=True):
        gap=0.
        if p['phases']==2:
            lean.append(p['phi_lean']-r['phi_lean']); rich.append(p['phi_rich']-r['phi_rich'])
            gap=p['phi_rich']-p['phi_lean']
        widths.append(abs(gap-(r['phi_rich']-r['phi_lean'])))
    return dict(n_total=len(records),n_two_phase=len(lean),n_missed=len(records)-len(lean),
                mae_lean=float(np.mean(np.abs(lean))) if lean else None,
                mae_rich=float(np.mean(np.abs(rich))) if rich else None,
                rmse_pooled=float(np.mean(np.square(lean+rich))**.5) if lean else None,
                mae_gap=float(np.mean(widths)))

def evaluate(cloud, osm, coex, weight):
    a,b,c=fit(cloud,osm,weight)
    theta=b/(.5-a+c) if .5-a+c>0 else None
    cp_pred=np.array([b/(binodal_chi0(r['phi_p'],r['N'],c)-a) for r in cloud])
    errors=cp_pred-np.array([r['T_K'] for r in cloud])
    train=np.array([r['split']=='calibration' for r in cloud])
    a2_pred=predicted_a2([r['T_K'] for r in osm],a,b,c)
    a2err=a2_pred-np.array([r['A2_mol_ml_g2'] for r in osm])
    pairs=[]
    for r in coex:
        n=r['Mn_used_g_mol']*VBAR_PS/VREF
        pred=merge_tangent(n,r['T_K'],a,b,c)
        phi0=r['phi_total']
        if pred['phases']==2 and not (pred['phi_lean']<phi0<pred['phi_rich']):
            pred=dict(phases=1,phi_lean=None,phi_rich=None,
                      global_binodal=[pred['phi_lean'],pred['phi_rich']],
                      reason='outside binodal at sample composition')
        pairs.append(pred)
    metrics=endpoints_summary(coex,pairs)
    model=dict(A=a,B=b,C=c,weight=weight,theta_K=theta,train=error_summary(errors[train],'K'),
               cloudHoldout=error_summary(errors[~train],'K'),osmoticFit=error_summary(a2err,'mol_ml_g2'),
               endpointMetrics=metrics)
    model['criticalBySeries']={}
    for r in coex:
        if r['series'] in model['criticalBySeries']: continue
        n=r['Mn_used_g_mol']*VBAR_PS/VREF
        x=brentq(lambda p:-1/(n*p*p)+1/(1-p)**2-6*c,1e-12,1-1e-12)
        crit=.5*(1/(n*x)+1/(1-x))+c*(1-3*x)
        model['criticalBySeries'][r['series']]=dict(phi_c=x,Tc_K=b/(crit-a),observed_Tc_K=r.get('Tc_measured_K'))
    return model,(cp_pred,a2_pred,pairs)

def run(cloud, osmotic, coex):
    op= osmotic['points']
    points=coex['points']
    models=[]; sensitivity=[]
    for w in WEIGHTS:
        m,preds=evaluate(cloud,op,points,w)
        held = osmotic.get('holdoutPoints',[])
        hp = predicted_a2([r['T_K'] for r in held],m['A'],m['B'],m['C'])
        m['osmoticHoldout'] = error_summary(hp - np.array([r['A2_mol_ml_g2'] for r in held]),'mol_ml_g2')
        m['thetaDelta_K'] = m['theta_K']-osmotic['thetaHoldout']['T_K']
        m['id']={0:'cloud',1:'combined'}.get(w,'lambda_'+str(w))
        m['label']={0:'Yalnız bulutlanma',1:'Bulutlanma + ozmometri'}.get(w,'λ = '+str(w))
        sensitivity.append(m.copy())
        if w not in [0,1]: continue
        models.append(m)
        for row,h in zip(held,hp,strict=True):
            row.setdefault('predictions',{})[m['id']]=float(h)
        for rows, y in zip([cloud,op,points],preds,strict=True):
            for row,p in zip(rows,y,strict=True):
                row.setdefault('predictions',{})[m['id']]=float(p) if not isinstance(p,dict) else p
    shared=[r for r in points if all(r['predictions'][m['id']]['phases']==2 for m in models)]
    for m in models:
        m['sharedEndpointMetrics']=endpoints_summary(shared,[r['predictions'][m['id']] for r in shared])
    return dict(version=1,title='Bulutlanma ve ozmometriden faz bileşimlerine',
                constants=CONSTANTS, models=models, osmotic=osmotic,coexistence=coex,
                cloudPoints=cloud,sensitivity=sensitivity,
                method=dict(groupWeight=1,sensitivityWeights=WEIGHTS,cloudUnit='K',massConcentrationUnit='g/mL',
                    osmoticUnit='mL mol g⁻²',temperatureCloudScale=SCALE_T,osmoticScale=SCALE_A2,
                    sampleWeights='Each cloud record equally weighted within its block; each osmometry record equally weighted within its block.',
                    conditionalEndpointErrors='Composition errors conditional on two predicted phases. Gap MAE includes missed pairs with predicted gap 0.'))
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('cloud','osmotic','coexist','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    result=run(read_cloud(args.cloud),json.loads(args.osmotic.read_text()),json.loads(args.coexist.read_text()))
    args.out.write_text('window.TRANSFER_RESULTS='+json.dumps(result,ensure_ascii=False,allow_nan=False,separators=(',',':'))+';\n')
    print(json.dumps(result['models'],indent=2,ensure_ascii=False))
if __name__=='__main__':
    main()
