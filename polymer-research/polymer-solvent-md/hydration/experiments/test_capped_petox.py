import unittest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from capped_petox import EXPLICIT, identify_sites


class TestCappedTrimer(unittest.TestCase):
    def test_exact_formula_endpoints_and_aligned_sites(self):
        mol=Chem.AddHs(Chem.MolFromSmiles(EXPLICIT))
        sites=identify_sites(mol)
        self.assertEqual((sites["formula"],sites["explicit_atoms"],sites["charge_e"]),("C16H31N3O4",54,0))
        refs=sites["ordered_sites_alpha_to_omega"]
        self.assertEqual([r["site"] for r in refs],["alpha","middle","omega"])
        for i,r in enumerate(refs):
            n=r["amide_n"]
            self.assertEqual(len(Chem.GetShortestPath(mol,sites["alpha_methyl_c"],n))-1,1+3*i)
            self.assertTrue(all(a.GetSymbol()!="H" for a in mol.GetAtomWithIdx(n).GetNeighbors()))
        donor=sites["omega_hydroxyl_o"]
        self.assertEqual(len(sites["omega_hydroxyl_h"]),1)
        self.assertEqual(mol.GetBondBetweenAtoms(donor,sites["omega_hydroxyl_h"][0]).GetBondType(),Chem.BondType.SINGLE)

    def test_equal_formula_but_misplaced_alcohol_fails(self):
        positional_isomer=Chem.AddHs(Chem.MolFromSmiles(EXPLICIT[:-3]+"C(O)C"))
        self.assertEqual(rdMolDescriptors.CalcMolFormula(positional_isomer),"C16H31N3O4")
        with self.assertRaisesRegex(ValueError,"end-group"):
            identify_sites(positional_isomer)


if __name__=="__main__":
    unittest.main()
