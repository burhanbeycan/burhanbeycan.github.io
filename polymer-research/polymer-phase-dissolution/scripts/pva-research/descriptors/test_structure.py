import json
import unittest

from rdkit import Chem
from structure import CATALOG, HERE, describe, mass_concentration, molar_masses, pva_composition, repeat_descriptors, solvent_info, statistical_repeat
from encode_polybert import normalized_input


class TestDescriptors(unittest.TestCase):
    def test_pva_composition(self):
        self.assertAlmostEqual(pva_composition(100)['covalent_repeat_mass_g_mol'],44.053)
        self.assertAlmostEqual(pva_composition(0)['covalent_repeat_mass_g_mol'],86.090)
        specimen=pva_composition(88)
        self.assertAlmostEqual(specimen['repeat_mass_g_mol_on_selected_basis'],49.09744)
        self.assertAlmostEqual(specimen['sites_mmol_per_g_on_selected_basis']['alcoholic_or_phenolic_OH'],17.9235,places=4)
        self.assertAlmostEqual(specimen['sites_per_repeat']['acetate_ester'],.12)

    def test_invalid_graphs(self):
        for text in ('CCO','[*]CCO','*C(*)C*','*CO.*CO','[*]=CC[*]','*[*]'):
            with self.subTest(text=text),self.assertRaises(ValueError):
                repeat_descriptors(text)
        for row in CATALOG.values():
            self.assertGreater(repeat_descriptors(row['psmiles'])['covalent_repeat_mass_g_mol'],40)

    def test_registered_isomer_and_epimers(self):
        def key(name):
            graph=CATALOG[name]['psmiles'].replace('[*]O[','O[').replace('([*])','(O)').replace('[O-]','O')
            return Chem.MolToInchiKey(Chem.MolFromSmiles(graph))
        self.assertEqual(key('CMC'),'WQZGKKKJIJFFOK-VFUOTHLCSA-N')
        self.assertEqual(key('ALG_M'),'AEMOLEFTQBMNLQ-SYJWYVCOSA-N')
        self.assertEqual(key('ALG_G'),'AEMOLEFTQBMNLQ-AZLKCVHYSA-N')
        self.assertEqual(normalized_input('[*:5]CC([*:2])O', CATALOG),CATALOG['PVA']['psmiles'])

    def test_units(self):
        sample=json.loads((HERE/'sample.json').read_text())
        self.assertAlmostEqual(describe(sample)['c_times_intrinsic'],.105)
        sample['concentration']={'unit':'wt%','value':4}
        self.assertIsNone(describe(sample)['c_times_intrinsic'])
        sample['concentration']['solution_density_g_mL']=1.02
        self.assertAlmostEqual(describe(sample)['c_times_intrinsic'],1.7136)
        with self.assertRaises(ValueError):
            mass_concentration(0.4,'g per solvent volume')
        sample['temperature_K']=298
        with self.assertRaises(ValueError):
            describe(sample)

    def test_unknown_moments_and_dh(self):
        moments=molar_masses(mn=10000)
        self.assertIsNone(moments['Mw_g_mol'])
        self.assertIsNone(moments['Mv_g_mol'])
        self.assertEqual(molar_masses(mn=10000,dispersity=1.5)['Mw_g_mol'],15000)
        s={'family':'PVA','solvent':'water','temperature_K':298,'Mn_g_mol':10000,'Mn_origin':'hypothetical'}
        self.assertIsNone(describe(s)['structure'])
        with self.assertRaises(ValueError):molar_masses(mn=10000,mw=8000)
        self.assertEqual(solvent_info('water')['donor_sites_per_molecule'],1)
        self.assertEqual(solvent_info('dmso')['donor_sites_per_molecule'],0)

    def test_end_groups(self):
        s={'family':'PEO','solvent':'water','temperature_K':298,'Mn_g_mol':458.545,
           'Mn_origin':'theoretical neutral PEG','combined_end_group_mass_g_mol':18.015,
           'end_group_basis':'H and OH relative to C2H4O repeats'}
        self.assertAlmostEqual(describe(s)['number_average_repeat_count'],10)
        del s['combined_end_group_mass_g_mol']
        self.assertIsNone(describe(s)['number_average_repeat_count'])

    def test_counterion_mass(self):
        comp=[{'mole_fraction':.95,'psmiles':CATALOG['HPAM']['psmiles']},
              {'mole_fraction':.05,'psmiles':CATALOG['AA_anion']['psmiles']}]
        r=statistical_repeat(comp,sodium_per_repeat=.05)
        self.assertAlmostEqual(r['repeat_mass_g_mol_on_selected_basis']-r['covalent_repeat_mass_g_mol'],1.1495)
        with self.assertRaises(ValueError):statistical_repeat(comp,sodium_per_repeat=.5)


if __name__=='__main__': unittest.main()
