"""Compile independent PVA/water observations from explicit primary-source numbers.

Each source tuple is a direct observation, not a fit, simulated target, or reading from
a smoothed figure. Field definitions and protocol restrictions are in SOURCES.md.
Run with python build_literature.py, using Python 3.9 or newer.
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = (
    'study_id source_url source_location grade_id Mn_g_mol Mw_g_mol '
    'untyped_molar_mass_min_g_mol untyped_molar_mass_max_g_mol '
    'DH_mol_pct DH_lower_pct DH_upper_pct solvent concentration concentration_unit '
    'temperature_K temperature_C_reported shear_rate_s_inv eta_cP plusminus_cP '
    'measurement_type instrument eta_source_value eta_source_unit '
    'source_to_cP geometry_or_speed notes'
).split()
rows = []


def observation(study, *, concentration, concentration_unit, value, eta_unit,
                t_C=None, rate=None, spread=None, **fields):
    factors = {'Pa s': 1000, 'cP': 1, 'mPa s': 1, '10^-1 Pa s': 100,
               '10^-3 Pa s': 1}
    factor = factors[eta_unit]
    item = {**study, **fields, 'solvent': 'water', 'concentration': concentration,
            'concentration_unit': concentration_unit, 'eta_source_value': value,
            'eta_source_unit': eta_unit, 'source_to_cP': factor,
            'eta_cP': round(value*factor, 8), 'shear_rate_s_inv': rate,
            'plusminus_cP': round(spread * factor, 8) if spread is not None else None,
            'temperature_C_reported': t_C,
            'temperature_K': round(t_C+273.15, 2) if t_C is not None else None}
    rows.append(item)


modarress = dict(study_id='Modarress2005',
    source_url='https://citeseerx.ist.psu.edu/document?doi=9f1c27f6955ee64d785ae0bd225e605396d158b3&repid=rep1&type=pdf',
    source_location='Experimental, p.182; Table 1, p.183', grade_id='Merck-Mn49000-DH98',
    Mn_g_mol=49000, DH_mol_pct=98, measurement_type='Newtonian dynamic viscosity',
    instrument='Haake Rotovisco RV100/CV100 concentric cylinders',
    notes='Published concentration denominator retained; authors identify the linear stress versus rate regime. Individual shear rates are unreported.')
for c, values in (('0.6107',(3.91,3.00,2.31)),('0.5381',(3.20,2.40,1.80)),
                   ('0.2467',(1.41,1.11,0.90)),('0.1180',(1.11,0.90,0.74))):
    for t,v in zip((15,25,35),values):
        observation(modarress, concentration=c, concentration_unit='g/20 mL (source)', t_C=t, value=v, eta_unit='10^-3 Pa s')

perez = dict(study_id='PerezPuyana2018',source_url='https://doi.org/10.1088/2053-1591/aab164',
    source_location='Sections 2.1–2.3, 2.6; Table 1, four pure PVA rows',grade_id='Sigma-Mw130000-DH86.7-88.7',
    Mw_g_mol=130000,DH_lower_pct=86.7,DH_upper_pct=88.7,
    measurement_type='Apparent steady shear viscosity',instrument='TA Instruments AR2000',
    geometry_or_speed='40 mm parallel plates; controlled stress 1–200 Pa',
    notes='Aqueous PVA, stirred 4 h at 40 °C. The authors report at least three replicates. Pure 10 wt% appears again as 10/0 mixture in Table 1; counted once.')
for c,v,s in ((2.5,.0103,.00255),(5,.126,.0372),(7.5,.405,.150),(10,1.49,.18)):
    observation(perez,concentration=c,concentration_unit='wt%',t_C=25,rate=10,value=v,spread=s,eta_unit='Pa s')

turan = dict(study_id='Turan2021',source_url='https://doi.org/10.35378/gujs.794638',
    source_location='Section 2.1; Table 3, five additive-free PVA rows',grade_id='Acros-average88000-DH88',
    untyped_molar_mass_min_g_mol=88000,untyped_molar_mass_max_g_mol=88000,DH_mol_pct=88,
    measurement_type='Apparent steady shear viscosity',instrument='Lamy Rheology First Plus',
    notes='Paper identifies average M.W.; moment not defined. Temperature taken from each viscosity row, not electrospinning. N=5. The ± statistic is preserved as printed.')
for c,t,g,v,s in ((5,20.8,100,60.313,.373),(8,22,100,674.91,9.48),
                    (10,22.4,100,1247.,10.13),(12,21.8,100,3522.67,41.26),(15,22.1,25,10642.11,265.57)):
    observation(turan,concentration=c,concentration_unit='wt%',t_C=t,rate=g,value=v,spread=s,eta_unit='mPa s')

saxena = dict(study_id='Saxena2018',source_url='https://doi.org/10.1155/2018/1738612',
    source_location='Section 2; Table 2',grade_id='PVA-approx140000-untyped',
    untyped_molar_mass_min_g_mol=140000,untyped_molar_mass_max_g_mol=140000,
    measurement_type='Capillary viscometry, unspecified effective shear rate',instrument='Ostwald capillary',
    notes='Approximate molar mass with unidentified moment. Source table caption gives 10^-1 Pa s. 1 MHz is the separate acoustic measurement frequency. Ostwald result is not tagged as zero-shear plateau.')
for c,v in ((.3,.072),(.4,.073),(.5,.074),(.6,.077),(.8,.086),(1.,.089)):
    observation(saxena,concentration=c,concentration_unit='% m/v',t_C=35,value=v,eta_unit='10^-1 Pa s')

lip = dict(study_id='Lipophilic2021',source_url='https://doi.org/10.3390/nano11102514',
    source_location='Sections 2.1–2.3; Table 1',grade_id='Hefei-Mw80000',Mw_g_mol=80000,
    measurement_type='Spindle viscosity, rate unresolved',instrument='Fungilab Viscolead One',
    geometry_or_speed='L3 rotor; 10–100 rpm method range',
    notes='Temperature of viscosity measurements 20 °C. No per-row equivalent shear rate is supplied for Table 1. ± value preserved in reported units.')
for c,v,s in ((8,241.34,3.33),(9,457.24,8.17),(10,655.4,6.50),(11,1101.84,16.82),
              (12,1602.94,31.19),(13,2056.52,33.73),(14,3615.84,27.77),(15,5845.16,69.42)):
    observation(lip,concentration=c,concentration_unit='wt%',t_C=20,value=v,spread=s,eta_unit='cP')

electrode = dict(study_id='Electrode2022',source_url='https://doi.org/10.3390/membranes12060609',
    source_location='Sections 2.1–2.2, 3.1; exact numerical values accompanying Fig. 1',
    measurement_type='Spindle viscosity, fixed instrument setting',instrument='Fungilab Viscolead adv',
    geometry_or_speed='L6 tip, 100 rpm',
    notes='Molecular-weight range as reported without explicit moment. Section 2.2 sample temperature 20 ±0.15 °C. Only values printed in text were included, not the graph line.')
for c,v,grade,lo,hi,dl,dh in ((8,15,'Sigma-30-70k',30000,70000,87,90),
                              (6,30,'Sigma-89-98k',89000,98000,99,99),
                              (8,90,'Sigma-89-98k',89000,98000,99,99)):
    observation(electrode,grade_id=grade,concentration=c,concentration_unit='wt%',t_C=20,value=v,eta_unit='mPa s',
                untyped_molar_mass_min_g_mol=lo,untyped_molar_mass_max_g_mol=hi,
                DH_mol_pct=dl if dl==dh else None,DH_lower_pct=dl if dl!=dh else None,DH_upper_pct=dh if dl!=dh else None)

observation(dict(study_id='Kan2022',source_url='https://doi.org/10.3390/nano12060998',
    source_location='Sections 2.1, 2.3, 2.5, 3.1; Table 2 row PVA 10%',grade_id='AppliChem-72000',
    untyped_molar_mass_min_g_mol=72000,untyped_molar_mass_max_g_mol=72000,DH_lower_pct=85,DH_upper_pct=89,
    measurement_type='Apparent steady shear viscosity',instrument='EMS-1000 rotating sphere',
    geometry_or_speed='2 mm aluminum sphere; 500 rpm for 5 s',
    notes='Room temperature as reported. Solution was stirred at 250 rpm for 2 h at 80 °C; these are preparation conditions. Per-row shear rate comes from Table 2.'),
    concentration=10,concentration_unit='wt%',value=1058.9,rate=3.1224,eta_unit='mPa s')

observation(dict(study_id='FexOy2024',source_url='https://doi.org/10.3390/membranes14090189',
    source_location='Sections 2.1–2.2; Table 1 pure PVA row',grade_id='Sigma-Mowiol10-98',
    untyped_molar_mass_min_g_mol=61000,untyped_molar_mass_max_g_mol=61000,DH_mol_pct=98,
    measurement_type='Spindle viscosity, fixed instrument setting',instrument='Brookfield LV DV-E',
    geometry_or_speed='S31 spindle, 50 rpm',
    notes='25 °C is specified for electrospinning, not for the viscosity determination. Pure PVA row only. No conversion from rpm to s^-1.'),
    concentration=12,concentration_unit='wt%',value=309.33,spread=.65,eta_unit='cP')

with (HERE/'pva_measurements.csv').open('w',newline='',encoding='utf-8') as handle:
    writer=csv.DictWriter(handle,fieldnames=COLUMNS)
    writer.writeheader()
    writer.writerows(rows)

assert len(rows)==40
assert len({(r['study_id'],r['grade_id'],r['concentration'],r['temperature_K'],r['shear_rate_s_inv']) for r in rows})==40
assert len({r['study_id'] for r in rows})==8
assert sum(r['shear_rate_s_inv'] is not None and r['temperature_K'] is not None for r in rows)==9
assert all(r['eta_cP']>0 for r in rows)
print('40 verified observations, 8 sources, 9 with both numerical shear rate and temperature.')
