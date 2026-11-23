import math
import unittest

from rdkit import Chem
from rdkit.Chem import Descriptors

from polymer_chemistry_lab.model import (
    H_EXACT_MASS, H_MASS, MAX_ATOMS, TEMPLATES, build_oligomer,
    molecular_view, parse_template, repeat_information, sequence_statistics,
)
from polymer_chemistry_lab.plugin import APP, calculate


class ChemistryTests(unittest.TestCase):
    def test_repeat_formulas_are_polymerized_contributions(self):
        for key, item in TEMPLATES.items():
            with self.subTest(key=key):
                self.assertEqual(repeat_information(parse_template(item["smiles"]))["formula"], item["formula"])

    def test_pe_h_caps_and_mass_balance(self):
        for n in (1, 2, 11):
            with self.subTest(n=n):
                oligomer = build_oligomer(TEMPLATES["pe"]["smiles"], "A" * n)
                self.assertEqual(oligomer.summary["formula"], f"C{2*n}H{4*n+2}")
                self.assertEqual(oligomer.molecule.GetNumAtoms(), 2*n)
                self.assertEqual(oligomer.molecule.GetNumBonds(), 2*n-1)
                repeat = oligomer.repeats["A"]
                self.assertAlmostEqual(oligomer.summary["average_molar_mass_g_mol"], n*repeat["average_mass_g_mol"]+2*H_MASS)
                self.assertAlmostEqual(oligomer.summary["exact_mass_Da"], n*repeat["exact_mass_Da"]+2*H_EXACT_MASS)

    def test_peo_and_pva_have_same_composition_different_connectivity(self):
        peo = build_oligomer(TEMPLATES["peo"]["smiles"], "AAAA")
        pva = build_oligomer(TEMPLATES["pva"]["smiles"], "AAAA")
        self.assertEqual(peo.summary["formula"], "C8H18O4")
        self.assertEqual(peo.summary["formula"], pva.summary["formula"])
        self.assertNotEqual(peo.summary["canonical_isomeric_smiles"], pva.summary["canonical_isomeric_smiles"])
        # Only one terminal OH in the H-capped PEO model; four OH in PVA.
        hydroxyl = Chem.MolFromSmarts("[O;H1]")
        self.assertEqual(len(peo.molecule.GetSubstructMatches(hydroxyl)), 1)
        self.assertEqual(len(pva.molecule.GetSubstructMatches(hydroxyl)), 4)

    def test_petox_is_n_propionyl_not_n_ethyl(self):
        oligomer = build_oligomer(TEMPLATES["petox"]["smiles"], "AAAA")
        self.assertEqual(oligomer.summary["formula"], "C20H38N4O4")
        motif = Chem.MolFromSmarts("N-C(=O)-C-C")
        self.assertEqual(len(oligomer.molecule.GetSubstructMatches(motif)), 4)
        secondary_amide = Chem.MolFromSmarts("[N;H1]-C(=O)")
        self.assertEqual(len(oligomer.molecule.GetSubstructMatches(secondary_amide)), 1)

    def test_block_and_alternating_same_mass_different_sequence(self):
        block = build_oligomer(TEMPLATES["petox"]["smiles"], "AAABBB", TEMPLATES["pmetox"]["smiles"])
        alternate = build_oligomer(TEMPLATES["petox"]["smiles"], "ABABAB", TEMPLATES["pmetox"]["smiles"])
        self.assertAlmostEqual(block.summary["average_molar_mass_g_mol"], alternate.summary["average_molar_mass_g_mol"])
        self.assertNotEqual(block.summary["canonical_isomeric_smiles"], alternate.summary["canonical_isomeric_smiles"])
        self.assertEqual(sequence_statistics("AAABBB")["dyad_counts"], {"AA": 2, "AB": 1, "BA": 0, "BB": 2})
        self.assertEqual(sequence_statistics("ABABAB")["dyad_counts"], {"AA": 0, "AB": 3, "BA": 2, "BB": 0})

    def test_internal_explicit_stereo_retained_and_caps_may_remove_stereocentre(self):
        oligomer = build_oligomer("[*:1]C[C@H](F)[*:2]", "AAAA")
        reference = Chem.MolFromSmiles("C[C@H](F)C[C@H](F)C[C@H](F)CCF")
        self.assertEqual(oligomer.summary["canonical_isomeric_smiles"], Chem.MolToSmiles(reference))
        self.assertEqual(oligomer.summary["assigned_tetrahedral_centres"], 3)
        # Replacing terminal port by H makes its CHF site CH2F, correctly achiral.
        single = build_oligomer("[*:1]C[C@H](F)[*:2]", "A")
        self.assertEqual(single.summary["assigned_tetrahedral_centres"], 0)

    def test_stereo_is_not_invented(self):
        oligomer = build_oligomer(TEMPLATES["pva"]["smiles"], "AAAA")
        self.assertEqual(oligomer.summary["assigned_tetrahedral_centres"], 0)
        self.assertGreater(oligomer.summary["unassigned_tetrahedral_centres"], 0)
        self.assertTrue(calculate({"template_a": TEMPLATES["pva"]["smiles"], "sequence": "AAAA"})["warnings"])

    def test_isotopic_mass_balance(self):
        oligomer = build_oligomer("[*:1][13CH2]CO[*:2]", "AAA")
        self.assertIn("[13C]3", oligomer.summary["formula"])
        self.assertAlmostEqual(oligomer.summary["average_molar_mass_g_mol"], 3*oligomer.repeats["A"]["average_mass_g_mol"]+2*H_MASS)

    def test_reject_invalid_ports_or_templates(self):
        invalid = ["CC", "[*:1]CC*", "[*:1]CC[*:1]", "[*:1]CC([*:2])[*:3]", "[*:1]=CC[*:2]", "[*:1]CC[*:2].O", "[*:1][CH]C[*:2]", "[*:1][CH2:7]C[*:2]", "[*:1]CC[*:2] name", "[*:1][*:2]", "CC([*:1])C([*:2]) invalid", "[*:1]CO(C)C[*:2]"]
        for smiles in invalid:
            with self.subTest(smiles=smiles), self.assertRaises(ValueError):
                parse_template(smiles)

    def test_limits_and_sequence_validation(self):
        for sequence in ("", "AC", "A"*65, None):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                build_oligomer(TEMPLATES["pe"]["smiles"], sequence)
        with self.assertRaises(ValueError):
            build_oligomer("[*:1]" + "C"*61 + "[*:2]", "A"*10)
        self.assertEqual(build_oligomer(TEMPLATES["pe"]["smiles"], "aa aa").sequence, "AAAA")
        self.assertIsNone(sequence_statistics("A")["dyad_fractions"]["AA"])

    def test_two_dimensional_graph_has_valid_lineage_and_finite_coordinates(self):
        oligomer = build_oligomer(TEMPLATES["petox"]["smiles"], "AAB", TEMPLATES["pmetox"]["smiles"])
        view = molecular_view(oligomer)
        self.assertEqual(len(view["atoms"]), oligomer.molecule.GetNumAtoms())
        for atom in view["atoms"]:
            self.assertTrue(math.isfinite(atom["x"]) and math.isfinite(atom["y"]))
            self.assertIn(atom["unit_index"], range(3))
        joins = []
        for bond in view["bonds"]:
            a, b = view["atoms"][bond["begin"]], view["atoms"][bond["end"]]
            if a["unit_index"] != b["unit_index"]:
                joins.append(bond)
        self.assertEqual(len(joins), 2)

    def test_handoff_transfers_repeat_contribution_and_cap_correction_only(self):
        result = calculate({})
        handoff = result["handoffs"][0]
        self.assertEqual(handoff["target"], "chain")
        self.assertEqual(set(handoff["inputs"]), {"repeat_mass_g_mol", "end_group_mass_g_mol"})
        self.assertAlmostEqual(handoff["inputs"]["repeat_mass_g_mol"], 99.133)
        self.assertAlmostEqual(handoff["inputs"]["end_group_mass_g_mol"], 2.016)
        self.assertEqual(calculate({"sequence": "AAABBB"})["handoffs"], [])

    def test_all_presets_calculate(self):
        for preset in APP["presets"]:
            with self.subTest(name=preset["name"]):
                output = calculate(preset["inputs"])
                self.assertEqual(output["material_updates"]["chemistry"]["architecture"], "linear")
                self.assertTrue(output["molecular_views"][0]["atoms"])

    def test_runtime_output_contract_when_available(self):
        try:
            from polymer_workbench.exchange import default_material, make_result
        except ImportError:
            self.skipTest("Optional workbench runtime is not installed.")
        result = make_result(APP, {}, calculate({}), default_material())
        self.assertEqual(result["application"]["id"], "chemistry")


if __name__ == "__main__":
    unittest.main()
