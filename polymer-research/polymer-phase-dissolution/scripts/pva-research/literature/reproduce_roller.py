"""Transcribe Table 1 of 10.1002/pen.21599; no assumed rheometry temperature."""
import csv
from pathlib import Path
from math import log, exp
from statistics import mean

root=Path(__file__).resolve().parent
table=[('PVA11',67000,12.5,.2),('PVA12',67000,14.5,.48),('PVA13',67000,17,.95),
       ('PVA21',80000,10,.21),('PVA22',80000,12,.48),('PVA23',80000,14,.98),
       ('PVA31',150000,7.8,.25),('PVA32',150000,9,.48),('PVA33',150000,10.5,.98)]
with (root/'pva_measurements.csv').open() as baseline:
    columns=next(csv.reader(baseline))
with (root/'pva_roller.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=columns)
    w.writeheader()
    for sample,mw,c,eta in table:
        w.writerow(dict(study_id='Roller2010',source_url='https://doi.org/10.1002/pen.21599',
                        source_location='Table 1, p.937, '+sample,grade_id='roller-atactic-'+str(mw),Mw_g_mol=mw,
                        DH_mol_pct=88,solvent='water',concentration=c,concentration_unit='wt%',
                        shear_rate_s_inv=100,eta_cP=round(eta*1000,8),measurement_type='Apparent steady shear viscosity',
                        instrument='Haake RotoVisco1',eta_source_value=eta,eta_source_unit='Pa s',source_to_cP=1000,
                        notes='Atactic PVA; same 88% hydrolysis for all three grades. Table 4 temperature is for surface tension, not rheometry. Formulations chosen for approximately matched viscosities across molar masses.'))

# Exploratory leave-one-weight-class-out, not a new-polymer or independent-laboratory test.
# Estimate the common within-grade concentration slope on the two training grades.
# Then relate their intercepts to log Mw (two unique values, no tunable hyperparameters).
def fit(train, with_mass):
    if not with_mass:
        x=[log(c) for _,m,c,y in train]
        y=[log(y) for _,m,c,y in train]
        slope=sum((u-mean(x))*(v-mean(y)) for u,v in zip(x,y))/sum((u-mean(x))**2 for u in x)
        return lambda mass,c: exp(mean(y)+slope*(log(c)-mean(x)))
    groups={m:[(log(c),log(y)) for _,mw,c,y in train if mw==m] for _,m,_,_ in train}
    group_means={m:(mean(x for x,y in xy),mean(y for x,y in xy)) for m,xy in groups.items()}
    numer=denom=0
    for m,xy in groups.items():
        mx,my=group_means[m]
        numer+=sum((u-mx)*(v-my) for u,v in xy)
        denom+=sum((u-mx)**2 for u,v in xy)
    slope=numer/denom
    (m1,m2)=sorted(groups)
    b1=group_means[m1][1]-slope*group_means[m1][0]
    b2=group_means[m2][1]-slope*group_means[m2][0]
    mass_slope=(b2-b1)/log(m2/m1)
    return lambda mass,c: exp(b1+mass_slope*log(mass/m1)+slope*log(c))

checks=[]
for held in sorted(set(r[1] for r in table)):
    train=[r for r in table if r[1]!=held]
    for model,flag in [('logC',False),('logC+logMw',True)]:
        predict=fit(train,flag)
        for sample,mw,c,observed in table:
            if mw==held:
                pred=predict(mw,c)
                checks.append(dict(model=model,sample=sample,heldout_Mw_g_mol=held,concentration_wt_pct=c,
                                   observed_cP=round(1000*observed,8),predicted_cP=1000*pred,APE_pct=100*abs(pred/observed-1)))
with (root/'roller_holdout.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(checks[0]));w.writeheader();w.writerows(checks)
for model in ['logC','logC+logMw']:
    subset=[r for r in checks if r['model']==model]
    print(model,'MAPE =',round(mean(r['APE_pct'] for r in subset),2),'±10:',sum(r['APE_pct']<=10 for r in subset),'/9; per grade',[(m,round(mean(r['APE_pct'] for r in subset if r['heldout_Mw_g_mol']==m),1)) for m in [67000,80000,150000]])
