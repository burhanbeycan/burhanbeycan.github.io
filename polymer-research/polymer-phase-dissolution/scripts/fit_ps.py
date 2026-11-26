"""Rebuild the PS evidence and compare two Flory–Huggins models.

Run python fit_ps.py [local 3PDb HTML] [output JS]. Requires NumPy + SciPy.
Only narrow-distribution, hydrogenated PS and cyclohexane, upper branch, at
0.10 MPa. Two publications are assigned to distinct subsets before fitting.
"""
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import urlopen
import json
import csv
import math
import sys

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar

DB = 'https://pppdb.uchicago.edu/cloud_points'
PUBS = {
    'sip': ('10.1039/B301190A', 'Siporska, Szydłowski & Rebelo (2003)', 'calibration'),
    'sae': ('10.1021/ma60032a020', 'Saeki, Kuwahara, Konno & Kaneko (1973)', 'holdout')
}
M_SOLVENT = 84.16   # g/mol
RHO_SOLVENT = 0.779 # g/mL, 20 °C reference
VBAR_PS = 0.934     # mL/g; Table 4.1
VREF = M_SOLVENT / RHO_SOLVENT
VOLUME_REF = 'https://mro.massey.ac.nz/server/api/core/bitstreams/2daff9ba-d17d-4471-afca-f57b06741a7f/content'

class Rows(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell = [], [], None
    def handle_starttag(self, name, attrs):
        if name == 'tr': self.row = []
        if name == 'td': self.cell = []
    def handle_endtag(self, name):
        if name == 'td' and self.cell is not None:
            self.row.append(''.join(self.cell).strip())
            self.cell = None
        if name == 'tr' and self.row: self.rows.append(self.row)
    def handle_data(self, data):
        if self.cell is not None: self.cell.append(data)

def select_rows(data):
    parser = Rows()
    parser.feed(data)
    samples, points = {}, []
    for row in parser.rows:
        if len(row) != 14: continue
        polymer, pcas, _, solvent, scas, _, mw, disp, phi, w, pressure, cp, direction, ref = row
        source = next((key for key, record in PUBS.items() if ref.lower() == 'doi.org/' + record[0].lower()), None)
        if source is None or polymer != 'polystyrene' or solvent != 'cyclohexane': continue
        if pcas != '9003-53-6' or scas != '110-82-7' or float(pressure) != 0.1 or direction != 'positive': continue
        id = 'ps-%s-%d' % (source, float(mw))
        mw, disp = float(mw), float(disp)
        if not (1.0 <= disp < 1.07): continue
        mn = mw / disp
        n = mn * VBAR_PS / VREF
        samples[id] = dict(id=id,label=f'{mw/1000:g} kg/mol · {PUBS[source][1]}',Mw_g_mol=mw,
                           Mn_g_mol=mn,dispersity=disp,N=n,source=source,split=PUBS[source][2])
        w = float(w)
        points.append(dict(series=id,T_K=round(float(cp)+273.15, 5),w_p=w if math.isfinite(w) else None,
                           phi_p=float(phi),Mw_g_mol=mw,dispersity=disp,kind='cloud_point',source=source,
                           pressure_MPa=0.1,split=PUBS[source][2]))
    for id, sample in samples.items():
        included = [p for p in points if p['series']==id]
        sample['T_range_K'] = [min(p['T_K'] for p in included),max(p['T_K'] for p in included)]
        sample['phi_range'] = [min(p['phi_p'] for p in included),max(p['phi_p'] for p in included)]
    return samples, sorted(points,key=lambda x:(x['series'],x['phi_p'],x['T_K']))

def entropy(x,n):
    return x*math.log(x)/n+(1-x)*math.log1p(-x)
def d_entropy(x,n):
    return (math.log(x)+1)/n - math.log1p(-x)-1

def tie_error(x,y,n,c):
    """Remove the quadratic A+B/T term from equality of secant and tangent."""
    d = y-x
    mid = (x+y)/2
    if abs(d) < 0.02*min(mid,1-mid):
        e3 = -1/(n*mid**2) + 1/(1-mid)**2
        e5 = -6/(n*mid**4) + 6/(1-mid)**4
        return c - e3/6 - d*d*e5/240
    return c+2*((entropy(y,n)-entropy(x,n))/d-(d_entropy(x,n)+d_entropy(y,n))/2)/(d*d)

def binodal_chi0(x,n,c):
    """χ0 at which x is a binodal endpoint, not the value at the spinodal."""
    at_x = tie_error(x,x,n,c)
    if abs(at_x) < 1e-9:
        return (1/(n*x)+1/(1-x))/2+c*(1-3*x)
    if at_x < 0:  # x on rich side, use log y to capture tiny lean fraction
        y=math.exp(brentq(lambda z:tie_error(x,math.exp(z),n,c),-720,math.log(x),xtol=3e-13))
    else:
        y=brentq(lambda u:tie_error(x,u,n,c),x,1-1e-12,xtol=5e-14)
    d = y-x
    slope_diff = (math.log(y/x)/n + math.log((1-x)/(1-y)))/d
    return slope_diff/2+c*(1-1.5*(x+y))

def metrics(err):
    return {'n':len(err),'rmse_K':float(np.mean(err**2)**0.5),
            'mae_K':float(np.mean(abs(err))), 'bias_K':float(np.mean(err))}

def fit(training, points, samples):
    nlist = [samples[p['series']]['N'] for p in training]
    train_T = np.array([p['T_K'] for p in training])
    def profile(c):
        roots = np.array([binodal_chi0(p['phi_p'],n,c) for p,n in zip(training,nlist)])
        b, a = np.polyfit(1/train_T, roots, 1)
        offset = 293.15
        # Centered parameterization reduces A/B numerical correlation.
        init = [a+b/offset,b]
        def err(v):
            return 1/((roots-v[0])/v[1]+1/offset) - train_T
        fit = least_squares(err,init,ftol=2e-12,xtol=2e-12,gtol=2e-12)
        return (float(fit.x[0]-fit.x[1]/offset),float(fit.x[1]),float(np.mean(fit.fun**2)))
    bounds=(-0.10, 0.24)
    # Use the calibration source only for this 1D profile over concentration dependence.
    optimum=minimize_scalar(lambda c:profile(c)[2],bounds=bounds,method='bounded',options={'xatol':1e-10})
    models = []
    for id, label, c in [('fh0','χ(T)',0.),('fh1','χ(φ,T)',optimum.x)]:
        a,b,_=profile(c)
        pred = np.array([b/(binodal_chi0(p['phi_p'],samples[p['series']]['N'],c)-a) for p in points])
        obs = np.array([p['T_K'] for p in points])
        err = pred-obs
        for point,value in zip(points,pred):
            point.setdefault('predicted_T_K',{})[id] = round(float(value),6)
        grouped = {'trainMetric': metrics(err[[p['split']=='calibration' for p in points]]),
                   'evalMetric': metrics(err[[p['split']=='holdout' for p in points]])}
        train_masses = {p['Mw_g_mol'] for p in training}
        novel = [p['split']=='holdout' and p['Mw_g_mol'] not in train_masses for p in points]
        grouped['newMassEvalMetric'] = metrics(err[novel])
        by = {k:metrics(err[[p['series']==k for p in points]]) for k in samples}
        models.append(dict(id=id,label=label,A=a,B=b,C=float(c),**grouped,bySeries=by,
            fitSeries=sorted({p['series'] for p in training}),
            holdoutSeries=sorted({p['series'] for p in points if p['split']=='holdout'}),
            N=samples['ps-sip-123000']['N']))
    return models

def main():
    data = Path(sys.argv[1]).read_text() if len(sys.argv)>1 else urlopen(DB,timeout=30).read().decode()
    samples, points = select_rows(data)
    train = [p for p in points if p['split']=='calibration']
    if len(train)!=73 or len(points)-len(train)!=24: raise ValueError('Snapshot or source selection changed.')
    models = fit(train,points,samples)
    sources = [dict(id='3pdb',name='University of Chicago · Binary Polymer Solution Cloud Point Database',
                    url=DB),dict(id='volumes',name='Polymer dynamics studied by NMR and light scattering methods, Table 4.1',
                    url=VOLUME_REF)]+[dict(id=k,doi=v[0],name=v[1],url='https://doi.org/'+v[0]) for k,v in PUBS.items()]
    ps=dict(label='Polistiren / siklohekzan',sources=sources,samples=sorted(samples.values(),key=lambda x:(x['split'],x['Mw_g_mol'])),
            measurements=points,fit=models[1],fits=models,defaultSeries='ps-sip-123000',pressure_MPa=0.1,
            phi_basis='3PDb volume fraction φ (as reported in curated table); no w→φ reconstruction.',
            N_basis='N=(Mw/Đ) v̄_PS / V_s; v̄_PS=0.934 mL/g, V_s=84.16/0.779 mL/mol, fixed near room temperature. This is a narrow-distribution pseudocomponent.',
            model_note='Effective integral interaction parameter χ=A+B/T+Cφ fitted by unweighted least squares in cloud temperature. Experimental cloud points approximate coexistence for narrow distributions. Only low-temperature UCST branch at 0.10 MPa; no independent experimental spinodal calibration.',
            volume_model=dict(polymer_specific_volume_mL_g=VBAR_PS,solvent_molar_volume_mL_mol=VREF,solvent_reference_C=20),
            calibration_range=dict(Mw_g_mol=[25000,400000],phi=[min(p['phi_p'] for p in train),max(p['phi_p'] for p in train)],T_K=[min(p['T_K'] for p in train),max(p['T_K'] for p in train)]))
    root = Path(__file__).resolve().parents[1]
    web_root = root / 'dist' if (root / 'dist' / 'index.html').is_file() else root
    out = Path(sys.argv[2]) if len(sys.argv)>2 else web_root / 'measurements.js'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text('window.EVIDENCE=window.EVIDENCE||{version:1,systems:{}};\nwindow.EVIDENCE.systems.ps='+json.dumps(ps,ensure_ascii=False,separators=(',',':'),allow_nan=False)+';\n')
    csv_path = out.parent / 'ps-cloud-points.csv'
    fields = ['source_doi','series','split','Mw_g_mol','Mn_g_mol','dispersity','N_lattice_sites','phi_polymer','w_polymer','T_cloud_observed_K','pressure_MPa','T_cloud_FH0_K','T_cloud_FH1_K']
    with csv_path.open('w', encoding='utf8', newline='') as fp:
        writer=csv.writer(fp)
        writer.writerow(fields)
        for row in points:
            sample=samples[row['series']]
            writer.writerow([PUBS[row['source']][0], row['series'],row['split'],row['Mw_g_mol'],round(sample['Mn_g_mol'],6),row['dispersity'],round(sample['N'],6),row['phi_p'],row['w_p'],row['T_K'],row['pressure_MPa'],row['predicted_T_K']['fh0'],row['predicted_T_K']['fh1']])
    for m in models: print(m['id'],m['A'],m['B'],m['C'],m['trainMetric'],m['evalMetric'])
if __name__ == '__main__': main()
