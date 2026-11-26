"""Post-hoc single-parameter alternative: chi_E = A + B/T + C*phi + E*phi^3.
Usage: python scripts/probe-e.py [dist-or-flat-web-folder] [output-json]
Reports a separate exploratory analysis; leaves the Studio's active model unchanged.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import xlogy
from scipy.spatial import ConvexHull

here = Path(__file__).resolve().parent.parent
web = Path(sys.argv[1]) if len(sys.argv)>1 else (here/'dist' if (here/'dist').is_dir() else here)
raw = (web/'dense-results.js').read_text().split('=',1)[1].rstrip().rstrip(';')
report = json.loads(raw)
fit = report['base']
Mn,v,Vs=(report['sample'][k] for k in ('Mn_g_mol','vbar_ml_g','Vs_ml_mol'))
N=Mn*v/Vs
R=8.31446261815324
def state_grid(T,n):
    p=np.unique(np.r_[0,np.geomspace(1e-30,1e-3,160),np.linspace(0,1,n)])
    g=xlogy(p,p)/N+xlogy(1-p,1-p)+(fit['A']+fit['B']/T+fit['C']*p)*p*(1-p)
    return p,g,p**4*(1-p)
def grid_prediction(T,c,E,base):
    p,g,h=base
    f=g+E*h
    # The convex hull in increasing phi from leftmost to rightmost is the lower boundary.
    # Parabolic composition gaps span non-consecutive vertices; exact homogeneous analytic
    # mu is used on each resolved part of a convex branch.
    verts=ConvexHull(np.c_[p,f]).vertices
    i0=np.flatnonzero(verts==0)[0]
    cyclic=np.r_[verts[i0:],verts[:i0]]
    lower=cyclic[:np.flatnonzero(cyclic==len(p)-1)[0]+1]
    x=c*v/1000
    c0=fit['A']+fit['B']/T
    mu=np.log1p(-x)+x-x/N+x*x*(c0-fit['C']+2*fit['C']*x)+E*x**4*(4*x-3)
    j=np.searchsorted(p[lower],x)
    a,b=lower[j-1],lower[j]
    gap=(b-a)>1
    intercept=(p[b]*f[a]-p[a]*f[b])/(p[b]-p[a])
    return -1000*R*T/Vs*np.where(gap,intercept,mu)
def evaluate(rows,E,grids):
    values=np.array([x['c_g_L'] for x in rows])
    temps=np.array([x['T_K'] for x in rows])
    calc=np.empty(len(rows))
    for T in np.unique(temps):
        mask=temps==T
        calc[mask]=grid_prediction(T,values[mask],E,grids[T])
    return calc
def test_fold(f,n):
    rows=f['rows']
    train=[r for r in rows if r['id'] in f['train_ids']]
    test=[r for r in rows if r['id'] in f['score_ids']]
    grids={T:state_grid(T,n) for T in sorted({r['T_K'] for r in rows})}
    def obj(x):
        residual=evaluate(train,x,grids)-[r['observed_kPa'] for r in train]
        return np.mean(residual**2)
    xs=np.linspace(-1,1,65)
    ys=np.array([obj(t) for t in xs])
    candidates=[(-1.,ys[0]),(0.,ys[len(xs)//2]),(1.,ys[-1])]
    for i in range(1,len(xs)-1):
        if ys[i]<=min(ys[i-1],ys[i+1]):
            opt=minimize_scalar(obj,bounds=(xs[i-1],xs[i+1]),method='bounded',
                                options={'xatol':1e-9})
            candidates.append((opt.x,opt.fun))
    E=min(candidates,key=lambda item:item[1])[0]
    held=evaluate(test,E,grids)
    diff=held-np.array([r['observed_kPa'] for r in test])
    return dict(E=E,mae_kPa=np.mean(np.abs(diff)),rmse_kPa=np.sqrt(np.mean(diff**2)),
                held_kPa=held.tolist())
out={'status':'post-hoc structural sensitivity; pressure and third virial coefficient examined before this model family was chosen',
     'chi_term':'E*phi^3', 'A3':fit['A3'], 'units_A3':'mL² mol g⁻³',
     'algorithm':'Global convex envelope on 4001 and 8001 uniform nodes with an additional log grid at low phi; 65 E grid samples on [-1,1] and local loss refinement.',
     'folds':[]}
for f in report['folds']:
    coarse,fine=[test_fold(f,n) for n in [4001,8001]]
    out['folds'].append(dict(id=f['id'],**fine,
      max_pressure_difference_from_coarse_kPa=max(abs(x-y) for x,y in zip(fine['held_kPa'],coarse['held_kPa'])),
      E_difference_from_coarse=abs(fine['E']-coarse['E'])))
doc=json.dumps(out,indent=2,ensure_ascii=False)+'\n'
if len(sys.argv)>2: Path(sys.argv[2]).write_text(doc)
else: print(doc)
