"""Chemical regression checks at junctions, stereocenters and H-capped ends."""
import unittest

from rdkit import Chem
from polymer_identity import build_structure


PE = "[*:1]CC[*:2]"
PEO = "[*:1]CCO[*:2]"
PP = "[*:1]CC(C)[*:2]"


def build(repeat, dp=6, tacticity="as_drawn", **kw):
    return build_structure(dict(mode="repeat", smiles=repeat, dp=dp, tacticity=tacticity, **kw))


class IdentityTests(unittest.TestCase):
    def test_petox_mass_and_caps(self):
        a = build("[*:1]CCN(C(=O)CC)[*:2]", 3)
        self.assertEqual(a["formula"], "C15H29N3O3")
        self.assertAlmostEqual(a["molar_mass"], 299.415, places=3)
        self.assertEqual(a["h_caps"], 2)
        self.assertTrue(a["stereo_resolved"])

    def test_finite_graph_and_repeat_can_share_identity(self):
        a = build(PE, 3)
        b = build_structure(dict(mode="finite",smiles="CCCCCC"))
        self.assertEqual(a["structure_id"], b["structure_id"])
        self.assertNotEqual(a["recipe_id"], b["recipe_id"])
        self.assertEqual(b["composition_class"], "unspecified")

    def test_stars_count_arms_atoms_hydrogens(self):
        a = build_structure(dict(mode="star",smiles=PEO,core_smiles="C([*:1])([*:2])([*:3])[*:4]",dp=2))
        self.assertEqual(a["formula"], "C17H36O8")
        self.assertEqual((a["arm_count"],a["repeat_counts"],a["h_caps"]), (4, {"A":8}, 4))
        self.assertEqual(a["architecture"], "star")
        self.assertTrue(a["stereo_resolved"])

    def test_finite_custom_branch(self):
        self.assertEqual(build_structure(dict(mode="finite",smiles="CC(C)(CC)CC"))["formula"], "C7H16")

    def test_sequence_graph_matter_and_types(self):
        spec = dict(mode="sequence",smiles=PE,repeat_b=PEO)
        a = build_structure(dict(**spec,sequence="AABB"))
        b = build_structure(dict(**spec,sequence="ABAB"))
        self.assertEqual(a["formula"], b["formula"])
        self.assertNotEqual(a["structure_id"], b["structure_id"])
        self.assertEqual(a["composition_class"], "copolymer")
        r = build_structure(dict(mode="sequence",smiles="[*:1]CCO[*:2]",repeat_b="[*:2]CCO[*:1]",sequence="AB"))
        self.assertEqual(r["composition_class"], "homopolymer_repeats")

    def test_relative_parity_against_independent_explicit_finite_structures(self):
        for mode,dyads,expected in [
            ("syndiotactic","rrr","CCC[C@@H](C)C[C@H](C)C[C@@H](C)C[C@H](C)CC(C)C"),
            ("isotactic","mmm","CCC[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)CC(C)C"),
        ]:
            x = build(PP, tacticity=mode)
            self.assertEqual(x["stereo_pattern"]["dyads"], dyads)
            self.assertEqual(x["stereo_pattern"]["retained_units"], [2, 3, 4, 5])
            self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(expected)), x["canonical_smiles"])
            # Same oriented repeat expressed with different atom traversal.
            self.assertEqual(build("CC([*:2])C[*:1]",tacticity=mode)["structure_id"], x["structure_id"])

    def test_atactic_recipe_seed_and_identity(self):
        def get(seed):
            return build(PP, dp=12, tacticity="atactic", stereo_seed=seed)
        self.assertEqual(get(51)["recipe_id"],get(51)["recipe_id"])
        self.assertNotEqual(get(51)["structure_id"],get(52)["structure_id"])
        self.assertEqual(len(get(51)["stereo_pattern"]["dyads"]), 9)

    def test_tacticity_requires_eligible_vinyl_and_retained_dyads(self):
        for s, n in [(PEO,6),(PE,6),(PP,2),("[*:1]CCN(C(=O)CC)[*:2]",8)]:
            with self.subTest(s=s, n=n),self.assertRaises(ValueError):
                build(s, dp=n, tacticity="syndiotactic")

    def test_h_cap_loss_depends_on_side_group(self):
        pattern = build("[*:1]CC(O)[*:2]",6,"syndiotactic")["stereo_pattern"]
        self.assertEqual(pattern["dyads"],"rrrr")
        self.assertEqual(pattern["retained_units"],[1,2,3,4,5])

    def test_unknown_stereo_gate(self):
        self.assertEqual(build(PP)["unspecified_stereo_count"],4)
        for s in ("FC(Cl)Br","CC=CC"):
            x=build_structure(dict(mode="finite",smiles=s))
            self.assertEqual(x["unspecified_stereo_count"],1)
            self.assertFalse(x["stereo_resolved"])

    def test_directional_stereo_at_junctions(self):
        for s, tag in [(r"[*:1]/C=C/[*:2]",Chem.StereoDescriptor.Bond_Trans),
                       (r"[*:1]/C=C\[*:2]",Chem.StereoDescriptor.Bond_Cis)]:
            x=build(s,dp=5)
            actual=[d.descriptor for d in Chem.FindPotentialStereo(Chem.MolFromSmiles(x["canonical_smiles"]))]
            self.assertEqual(actual, [tag]*3)
            self.assertTrue(x["stereo_resolved"])

    def test_as_drawn_tetrahedral_and_isotopes(self):
        l=build("[*:1]O[C@H](C)C(=O)[*:2]",dp=4)
        r=build("[*:1]O[C@@H](C)C(=O)[*:2]",dp=4)
        self.assertTrue(l["stereo_resolved"])
        self.assertNotEqual(l["structure_id"],r["structure_id"])
        self.assertIn("13C",build("[*:1][13CH2]C[*:2]",dp=2)["formula"])

    def test_bad_input_and_domain_limits(self):
        for spec in [
                dict(mode="repeat",smiles=PE,dp=True),dict(mode="repeat",smiles=PE,dp=33),
                dict(mode="repeat",smiles="C=C",dp=4),dict(mode="finite",smiles=PE),
                dict(mode="finite",smiles="CC.O"), dict(mode="finite",smiles="C"*151),
                dict(mode="finite",smiles="[CH3]"),dict(mode="sequence",smiles=PE,repeat_b=PEO,sequence="ABX"),
                dict(mode="star",smiles=PE,core_smiles="C([*:1])([*:2])[*:2]",dp=3),
                dict(mode="repeat",smiles="[*:1]CCCCCC[*:2]",dp=32)]:
            with self.subTest(spec=spec),self.assertRaises(ValueError):
                build_structure(spec)

if __name__ == "__main__":
    unittest.main()
