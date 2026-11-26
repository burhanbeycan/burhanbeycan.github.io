"""Build the browser data and a reproducible, downloadable package. Standard library."""
from pathlib import Path
from collections import defaultdict
import csv,json,zipfile,shutil,sys
BASE=Path(__file__).resolve().parent
def load(path):return json.loads((BASE/path).read_text())
def rows(path):
    with (BASE/path).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def number(x):return f"{x:.2f}".replace('.',',')
def research_data():
    flow=load('model/model_state.json'); grades=load('grades/grade_state.json'); age=load('kinetics/aging_results.json')
    stats=load('model/benchmark.json');external=load('model/external_product_benchmark.json')['source_scores']['PerezPuyana2018']['point_weighted']
    observations={r['id']:r for r in rows('transfer/transfer_observations.csv')}
    grouped={}
    for r in rows('transfer/heldout_predictions.csv'):
        if (r['split'],r['mode'],r['model'])!=('group','zero_reference','explicit_plus_polybert'):continue
        ref=observations[r['id']]; units={'wt_percent':'wt%','ppm_reported':'ppm','g_per_100mL':'g/100 mL'}
        label=f"{ref['grade']} · {ref['conc_value']} {units[ref['conc_basis']]} · {r['T_K']} K"
        if float(ref['NaCl_g_L']):label+=f" · NaCl {ref['NaCl_g_L']} g/L"
        d=grouped.setdefault(r['curve'],dict(material=r['material'],label=label,t=float(r['T_K']),doi=ref['source_doi'],rows=[]))
        d['rows'].append([float(r['shear_s_inv']),float(r['eta_cP']),float(r['estimate_cP'])])
    for group in grouped.values():group['rows'].sort()
    transfer=list(grouped.values())
    macro=next(r for r in rows('transfer/macro_scores.csv') if (r['model'],r['mode'],r['split'])==('explicit_plus_polybert','zero_reference','group'))
    evidence=[]
    for choice,label in [('whole_curve','Ayrılmış akış eğrisi'),('whole_temperature','Ayrılmış sıcaklık'),('whole_concentration','Ayrılmış derişim')]:
        v=stats[choice];evidence.append([label+' · aynı PVA ürünü','Ortak kuvvet yasası / ters fiziksel model','MAPE %'+number(v['shared_power_law']['equal_curve_mean']['mape_pct'])+' / %'+number(v['inverse_parametric']['equal_curve_mean']['mape_pct'])])
    evidence.append(['Başka PVA ürünü ve laboratuvarı','Pérez-Puyana · desteklenen üç koşul','MAPE %'+number(external['mape_pct'])+' · ±%10 içinde 1/3'])
    evidence.append(['Yeni Mₙ · dört ticari numune','Numune dışı regresyon · Mₙ terimi var','Numune ağırlıklı MAPE %'+number(grades['metrics']['mass']['macro_by_grade']['MAPE_pct'])])
    evidence.append(['Yaşlanma · 0–16 saat, 303,15 K','Başlangıç sabit / kinetik düzeltme','MAPE %'+number(age['time_benchmark']['early']['fresh_cP']['mean_APE_pct'])+' / %'+number(age['time_benchmark']['early']['early_law_cP']['mean_APE_pct'])])
    evidence.append(['Üç kimyasal grup sırayla ayrıldı','Yapı tanımlayıcıları + polyBERT · ölçümsüz','286 noktada aile ağırlıklı ±%10 oranı %'+number(float(macro['macro_within10_pct']))])
    return dict(flow=flow,mass=grades,aging=age,evidence=evidence,transfer=transfer)
def main():
    app=BASE/'app'; state=research_data()
    (app/'pva-research-data.js').write_text('globalThis.PvaResearch = '+json.dumps(state,separators=(',',':'),ensure_ascii=False,allow_nan=False)+';\n')
    shutil.copyfile(BASE/'model/inverse-engine.js',app/'inverse-engine.js')
    if len(sys.argv)<2:return
    target=Path(sys.argv[1]);target.mkdir(exist_ok=True,parents=True)
    for name in ['pva-inverse.html','pva-inverse.css','pva-inverse.js','pva-research-data.js','inverse-engine.js']:
        shutil.copyfile(app/name,target/name)
    shutil.copyfile(BASE/'data/pva-inverse-data.csv',target/'pva-inverse-data.csv')
    with zipfile.ZipFile(target/'pva-research.zip','w',zipfile.ZIP_DEFLATED,compresslevel=8) as z:
        for file in sorted(BASE.rglob('*')):
            rel=file.relative_to(BASE)
            if any(t.startswith('.') or t in {'__pycache__','qa','publish'} for t in rel.parts) or not file.is_file() or file.suffix not in {'.md','.csv','.json','.py','.js','.cjs','.zip','.npz','.html','.css','.txt'} or file.suffix=='.zip': continue
            if rel.as_posix()=='app/pva-inverse.html':
                s=file.read_text().replace('href="index.html"','href="../README.md"').replace('href="viscosity-improvement.html"','href="../model/README.md"').replace('href="pva-inverse-data.csv"','href="../data/pva-inverse-data.csv"').replace('href="pva-research.zip" download>Analiz kodu, tanımlayıcılar ve testler (ZIP)','href="../README.md">Kod ve testlerin kullanım kılavuzu')
                z.writestr(rel.as_posix(),s)
            else:z.write(file,rel.as_posix())
if __name__=='__main__':main()
