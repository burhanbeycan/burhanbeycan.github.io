"""Unify observations without conflating methods, molar-mass averages or concentration bases."""
from pathlib import Path
from collections import Counter
import csv
import json
ROOT=Path(__file__).resolve().parent.parent
def load(name):
    with (ROOT/name).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def build():
    rows=[]
    def add(study,source,location,sample,series,**values):
        rows.append(dict(observation_id=f'PVA-{len(rows)+1:04}',study_id=study,source=source,source_location=location,grade_id=sample,protocol=series,**values))
    for r in load('model/pva_table2.csv'):
        t,g,tau=map(float,[r['temperature_C'],r['shear_rate_s_inv'],r['shear_stress_Pa']])
        add('Carrasco2023','https://doi.org/10.3390/polym15071743','Table 2','Merck-141351','rotational_steady',solvent='water',
            concentration=r['wt_percent'],concentration_unit='wt%',T_K=t+273.15,rate_s_inv=g,shear_stress_Pa=tau,eta_cP=round(1000*tau/g,8),origin='tabulated stress / rate',scope='grade Mn and DH unreported')
    for filename,solvent,resolved in [('parisi_water.csv','water',True),('parisi_unresolved.csv','water',False),('parisi_dmso.csv','DMSO',True)]:
        for r in load('parisi/'+filename):
            add('Parisi2022','https://doi.org/'+r['doi'],'Figure '+r['source_figure']+r['panel'],r['grade_id'],'rotational_steady',
                solvent=solvent,concentration=r['concentration_wt_pct'],concentration_unit='wt%',T_K=r['temperature_K'],
                rate_s_inv=r['shear_rate_s_inv'],eta_cP=r['viscosity_cP'],Mn_g_mol=r['mn_g_mol'],Mw_g_mol=r['mw_g_mol'],
                dispersity=r['dispersity'],DH_mol_pct=r['dh_mol_pct'],pixel_sensitivity_px=r['pixel_bound_px'],origin='digitized experimental marker',
                scope='grade matched' if resolved else 'unresolved legend/table grade; excluded from all grade models')
    for r in load('literature/pva_measurements.csv') + load('literature/pva_roller.csv'):
        add(r['study_id'],r['source_url'],r['source_location'],r['grade_id'],r['measurement_type'],
            solvent=r['solvent'],concentration=r['concentration'],concentration_unit=r['concentration_unit'],T_K=r['temperature_K'],
            rate_s_inv=r['shear_rate_s_inv'],eta_cP=r['eta_cP'],plusminus_cP=r['plusminus_cP'],Mn_g_mol=r['Mn_g_mol'],Mw_g_mol=r['Mw_g_mol'],
            DH_mol_pct=r['DH_mol_pct'],DH_low=r['DH_lower_pct'],DH_high=r['DH_upper_pct'],
            nominal_mass_lower=r['untyped_molar_mass_min_g_mol'],nominal_mass_upper=r['untyped_molar_mass_max_g_mol'],
            origin='published numerical value',scope=r['notes'])
    for r in load('kinetics/aging_fig1.csv'):
        add('Ohyanagi1958','https://doi.org/10.1295/koron1944.15.348','Figure 1',r['material_id'],'falling_ball_aging',
            solvent='water',concentration=r['wt_percent'],concentration_unit='wt%',T_K=r['storage_T_K'],time_h=r['age_h'],
            eta_cP=r['apparent_eta_cP'],DP_reported=r['reported_dp_index'],DH_low=99.9,origin='digitized experimental marker',scope='falling ball; time-dependent')
    for r in load('kinetics/terminal_relaxation.csv'):
        add('Wang2025','https://doi.org/'+r['study_doi'],r['source_location'],r['material_id'],'SAOS_terminal_extrapolation',
            solvent=r['solvent'],concentration=r['wt_percent'],concentration_unit='wt%',T_K=r['temperature_K'],terminal_eta0_cP=r['terminal_eta0_cP'],
            terminal_tau_ms=r['tau_d_ms'],nominal_mass_lower=r['nominal_mass_g_mol'],nominal_mass_upper=r['nominal_mass_g_mol'],DH_mol_pct=r['hydrolysis_percent'],
            origin='published terminal estimate from SAOS',scope='independent measurement target; no steady-shear rate')
    cols=['observation_id','study_id','source','source_location','grade_id','protocol','solvent','concentration','concentration_unit','T_K','rate_s_inv','time_h','eta_cP','plusminus_cP','terminal_eta0_cP','terminal_tau_ms','Mn_g_mol','Mw_g_mol','nominal_mass_lower','nominal_mass_upper','dispersity','DH_mol_pct','DH_low','DH_high','DP_reported','shear_stress_Pa','pixel_sensitivity_px','origin','scope']
    with (ROOT/'data'/'pva-inverse-data.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
    summary=dict(total=len(rows),studies=len(set(r['study_id'] for r in rows)),per_study=Counter(r['study_id'] for r in rows),nonempty_Mn=sum(bool(r.get('Mn_g_mol')) for r in rows),unresolved=sum('unresolved' in r.get('scope','') for r in rows),original_concentration_units=sorted({r['concentration_unit'] for r in rows}))
    (ROOT/'data'/'manifest.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    return summary
if __name__=='__main__':print(build())
