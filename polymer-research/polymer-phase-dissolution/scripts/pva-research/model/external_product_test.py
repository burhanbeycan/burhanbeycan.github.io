"""Stress-test the frozen Carrasco product model on independent PVA/water products.

Output values for independent grades are retrospective challenge predictions.
The production model is scoped to Merck-Peru 141351, studied by Carrasco et al.
Run from any directory after regenerating the literatures' primary observations.
"""
import argparse
import csv
import json
from pathlib import Path
from pva_model import metrics, predict_one, writecsv

HERE = Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--measurements',type=Path,default=HERE.parent/'literature'/'pva_measurements.csv')
    parser.add_argument('--parisi',type=Path,default=HERE.parent/'parisi'/'parisi_water.csv')
    parser.add_argument('--out',type=Path,default=HERE)
    args=parser.parse_args()
    state=json.loads((HERE/'model_state.json').read_text())
    with args.measurements.open(encoding='utf-8',newline='') as f:
        records=list(csv.DictReader(f))
    selected=[r for r in records if r['measurement_type']=='Apparent steady shear viscosity'
              and r['solvent']=='water' and r['concentration_unit']=='wt%' and r['temperature_K'] and r['shear_rate_s_inv']]
    rows=[]
    for r in selected:
        c=float(r['concentration']);t=float(r['temperature_K']);g=float(r['shear_rate_s_inv'])
        observed=float(r['eta_cP'])
        estimate=predict_one(state,c,t,g)
        pred=estimate['eta_app_cp']
        err=100*abs(pred/observed-1) if pred is not None else None
        rows.append(dict(study_id=r['study_id'],source_url=r['source_url'],source_location=r['source_location'],
                         external_grade=r['grade_id'],grade_seen_in_training=False,concentration_wt_percent=c,
                         temperature_K=t,rate_s_inv=g,observed_cp=observed,reported_plusminus_cp=r['plusminus_cP'],
                         numeric_support=estimate['domain'],challenge_prediction_cp=pred,abs_percentage_error=err))
    matched=[r for r in rows if r['challenge_prediction_cp'] is not None]
    report=dict(training_grade='Merck Peru 141351 (Carrasco et al., 2023)',training_source='https://doi.org/10.3390/polym15071743',
                aim='Stress-test for independent product and laboratory given matched nominal concentration, temperature and shear rate; this is outside the chemical grade scope of the deployed Carrasco model.',
                specified_points=len(rows),numerically_supported_points=len(matched),source_scores={})
    for study in sorted({r['study_id'] for r in rows}):
        group=[r for r in rows if r['study_id']==study]
        subset=[r for r in matched if r['study_id']==study]
        score=metrics([r['observed_cp'] for r in subset],[r['challenge_prediction_cp'] for r in subset]) if subset else None
        report['source_scores'][study]=dict(attempted=len(group),numeric_matches=len(subset),point_weighted=score)
    args.out.mkdir(exist_ok=True,parents=True)
    writecsv(args.out/'external_product_predictions.csv',rows)
    (args.out/'external_product_benchmark.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    # Separate challenge register: measured symbols digitized by the companion workflow.
    if args.parisi.exists():
        with args.parisi.open(encoding='utf-8',newline='') as f:
            digitized=list(csv.DictReader(f))
        comparisons=[]
        for r in digitized:
            c=float(r['concentration_wt_pct']); t=float(r['temperature_K']);g=float(r['shear_rate_s_inv'])
            y=float(r['viscosity_cP']); result=predict_one(state,c,t,g); value=result['eta_app_cp']
            comparisons.append(dict(study_id=r['study_id'],doi=r['doi'],figure=r['source_figure'],panel=r['panel'],
                                    grade_id=r['grade_id'],temperature_K=t,concentration_wt_percent=c,rate_s_inv=g,
                                    observed_cp=y,challenge_prediction_cp=value,availability=result['domain'],
                                    absolute_pct_error=100*abs(value/y-1) if value is not None else None))
        scored=[r for r in comparisons if r['challenge_prediction_cp'] is not None]
        bygrade={}
        for grade in sorted({r['grade_id'] for r in scored}):
            subset=[r for r in scored if r['grade_id']==grade]
            bygrade[grade]=dict(n=len(subset),**metrics([r['observed_cp'] for r in subset],
                                                       [r['challenge_prediction_cp'] for r in subset]))
        writecsv(args.out/'external_parisi_product_predictions.csv',comparisons)
        (args.out/'external_parisi_product_benchmark.json').write_text(json.dumps(
            dict(training_source='10.3390/polym15071743',test_source='10.1122/8.0000435',
                 source_figure='4 (water)',attempted=len(comparisons),numeric_matches=len(scored),
                 distinct_concentrations_wt_percent=sorted({r['concentration_wt_percent'] for r in scored}),
                 by_grade=bygrade),indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
