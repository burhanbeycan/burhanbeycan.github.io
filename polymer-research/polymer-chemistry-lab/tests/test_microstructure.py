import unittest

from rdkit import Chem
from polymer_chemistry_lab.plugin import APP, calculate
from polymer_chemistry_lab.microstructure import held_out_predictions,SAMPLES
from polymer_chemistry_lab.end_groups import petox_candidates

class NMRHoldoutTests(unittest.TestCase):
    def test_published_rows_and_no_test_label_leakage(self):
        folds=held_out_predictions(SAMPLES)
        self.assertEqual(len(folds),6)
        self.assertTrue(all(row[2] == 2 for row in folds))
        self.assertAlmostEqual(folds[1][5],(.72+.95)/2)
        changed=list(SAMPLES)
        name,route,da,_,oh,ac=changed[1]
        changed[1]=(name,route,da,1.2,oh*2,ac*3)
        self.assertEqual(held_out_predictions(changed)[1][5],folds[1][5])

    def test_score_units_and_controlled_composition(self):
        r=calculate({"chemistry_mode":"pva_benchmark"})
        self.assertAlmostEqual(r["metrics"][0]["value"],.11)
        self.assertAlmostEqual(r["metrics"][1]["value"],.4466666666666667)
        first,second=r["pva_benchmark"]["models"]
        self.assertEqual((first[1],second[1]),(.15,.15))
        self.assertAlmostEqual(first[3],.24225)
        self.assertAlmostEqual(second[3],.0765)
        self.assertEqual(len(r["molecular_views"]),2)
        self.assertEqual(r["handoffs"],[])
    
    def test_boundary_and_finite_sampling(self):
        for x in (0,1,-1,True,float("nan")):
            with self.subTest(x=x), self.assertRaises(ValueError):
                calculate({"chemistry_mode":"pva_benchmark","matched_acetate":x})
        for r in calculate({"chemistry_mode":"pva_benchmark","ensemble_count":2000})["pva_benchmark"]["models"]:
            self.assertAlmostEqual(r[3],r[4],delta=.02)

class HydrolysisTests(unittest.TestCase):
    def simulate(self,pattern,anion="chloride"):
        return calculate({"chemistry_mode":"petox_hydrolysis","ei_pattern":pattern,"ei_anion":anion})

    def test_intact_agrees_with_existing_end_graph(self):
        for dp in (1,6):
            o=self.simulate("E"*dp)
            self.assertEqual(len(o["molecular_views"]),1)
            mol=Chem.MolFromSmiles(o["molecular_views"][0]["smiles"])
            self.assertEqual(Chem.MolToSmiles(mol),petox_candidates(dp)[1]["smiles"])
            self.assertIsNone(o["metrics"][1]["value"])
    
    def test_all_amine_and_protonated_one_repeat(self):
        i=self.simulate("I"); p=self.simulate("P")
        for r,expected in [(i,"CNCCO"),(p,"C[NH2+]CCO")]:
            last=r["tables"][0]["rows"][-1]
            self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(last[5])),Chem.MolToSmiles(Chem.MolFromSmiles(expected)))
        self.assertEqual(p["hydrolysis"]["anion_count"],1)
    
    def test_sequence_not_only_composition(self):
        a=self.simulate("EIIP"); b=self.simulate("PEII")
        self.assertNotEqual(a["tables"][0]["rows"][-1][5],b["tables"][0]["rows"][-1][5])
        self.assertEqual(a["tables"][0]["rows"][-1][2],b["tables"][0]["rows"][-1][2])
    
    def test_salt_and_site_ledger(self):
        for anion in ("chloride","acetate","bromide"):
            r=self.simulate("EEIPEEIP",anion)
            covalent=r["tables"][0]["rows"][-1]
            self.assertEqual(covalent[2],"C29H62N8O5+2")
            self.assertEqual(r["hydrolysis"]["neutral_ei"],2)
            self.assertEqual(r["hydrolysis"]["protonated_ei"],2)
            salt=Chem.MolFromSmiles(r["tables"][1]["rows"][-1][2])
            self.assertEqual(Chem.GetFormalCharge(salt),0)
            self.assertEqual(len(Chem.GetMolFrags(salt)),3)
            self.assertIn("Counterion",r["molecular_views"][-1]["title"])
            self.assertEqual(r["material_updates"],{})

    def test_invalid_pattern_and_remainder(self):
        for x in ("","EX","P"*65,None,3):
            with self.assertRaises(ValueError):self.simulate(x)
        with self.assertRaises(ValueError):self.simulate("PI","divalent")
        self.assertEqual(self.simulate("e i p")["hydrolysis"]["dp"],3)
