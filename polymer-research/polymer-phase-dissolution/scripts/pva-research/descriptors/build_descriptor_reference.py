"""Write structure-only chemical features and primary PVA grade reference values."""
import csv
from structure import HERE,CATALOG,pva_composition,repeat_descriptors
from rdkit import rdBase

# Parisi et al., Tables I and III. Nominal DH and supplier molar masses.
PVA_GRADES=[
    ('171k-80',112300,171000,1.52,80,78,170),
    ('144k-88',90900,144500,1.59,88,89,163),
    ('49k-75',28000,49500,1.77,75,24,58),
    ('48k-88',30000,48100,1.61,88,42,68),
    ('24k-80',15100,24500,1.63,80,25,40),
]


def main():
    with (HERE/'pva_grade_descriptors.csv').open('w',newline='') as out:
        w=csv.writer(out)
        w.writerow(['source_doi','grade','Mn_g_mol','Mw_g_mol','reported_dispersity','DH_nominal_mol_pct',
                    'weighted_repeat_g_mol','alcohol_OH_mmol_g','acetate_ester_mmol_g','CH2_mmol_g',
                    'intrinsic_water_298_15K_mL_g','intrinsic_DMSO_298_15K_mL_g'])
        for grade,mn,mw,dd,dh,ew,ed in PVA_GRADES:
            r=pva_composition(dh)
            b=r['sites_mmol_per_g_on_selected_basis']
            w.writerow(['10.1122/8.0000435',grade,mn,mw,dd,dh,round(r['repeat_mass_g_mol_on_selected_basis'],6),
                        round(b['alcoholic_or_phenolic_OH'],6),round(b['acetate_ester'],6),round(b['CH2'],6),ew,ed])
    with (HERE/'motif_descriptors.csv').open('w',newline='') as out:
        w=csv.writer(out)
        w.writerow(['motif','psmiles','covalent_repeat_g_mol','CH2_per_repeat','alcoholic_or_phenolic_OH_per_repeat',
                    'acetate_ester_per_repeat','HBD_per_repeat','HBA_per_repeat','formal_charge','scope','rdkit'])
        for key,s in CATALOG.items():
            r=repeat_descriptors(s['psmiles'])
            w.writerow([key,s['psmiles'],round(r['covalent_repeat_mass_g_mol'],6),*r['sites_per_repeat'].values(),s['semantic'],rdBase.rdkitVersion])


if __name__=='__main__': main()
