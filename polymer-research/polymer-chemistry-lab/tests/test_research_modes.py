import unittest

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from polymer_chemistry_lab.end_groups import petox_candidates
from polymer_chemistry_lab.plugin import APP, calculate
from polymer_chemistry_lab.sequences import pva_ensemble


class EndGroupTests(unittest.TestCase):
    def test_dp_one_matches_independently_written_covalent_structures(self):
        reference = ["CCNC(=O)CC", "CN(CCO)C(=O)CC", "CN(CCN=[N+]=[N-])C(=O)CC"]
        for expected, candidate in zip(reference, petox_candidates(1)):
            mol = Chem.MolFromSmiles(expected)
            self.assertEqual(candidate["smiles"], Chem.MolToSmiles(mol))
            self.assertAlmostEqual(candidate["neutral_Da"], rdMolDescriptors.CalcExactMolWt(mol))

    def test_same_dp_sites_formula_and_charge(self):
        hh, oh, n3 = petox_candidates(6)
        self.assertEqual([c["formula"] for c in (hh, oh, n3)], ["C30H56N6O6","C31H58N6O7","C31H57N9O6"])
        for c in (hh, oh, n3):
            self.assertEqual(Chem.GetFormalCharge(Chem.MolFromSmiles(c["smiles"])),0)
            self.assertAlmostEqual(c["mz_H"] - c["neutral_Da"], 1.007276, places=7)
            self.assertAlmostEqual(c["mz_Na"] - c["neutral_Da"], 22.989218, places=7)
            self.assertAlmostEqual(c["mz_2H"]*2-c["neutral_Da"],2.014552, places=7)
        self.assertAlmostEqual(oh["end_mass_g_mol"]-hh["end_mass_g_mol"],30.026,places=3)
        smarts = Chem.MolFromSmarts("N-C-C-N=[N+]=[N-]")
        self.assertTrue(Chem.MolFromSmiles(n3["smiles"]).HasSubstructMatch(smarts))

    def test_alignment_and_colour_roles(self):
        baseline, *others = petox_candidates(4)
        atoms = baseline["view"]["atoms"]
        for candidate in others:
            view = candidate["view"]
            for ref, actual in zip(atoms, view["atoms"]):
                self.assertAlmostEqual(ref["x"], actual["x"])
                self.assertAlmostEqual(ref["y"], actual["y"])
            for role in ("alpha","omega"):
                self.assertTrue(any(a.get("end_role") == role for a in view["atoms"]))
                self.assertTrue(any(b.get("end_role") == role for b in view["bonds"]))

    def test_explicit_selection_for_handoff_and_evidence(self):
        self.assertFalse(calculate({"chemistry_mode":"petox"})["material_updates"])
        result = calculate({"chemistry_mode":"petox","selected_ends":"me_oh","end_evidence":"Illustrative candidate; no measurement."})
        self.assertAlmostEqual(result["handoffs"][0]["inputs"]["end_group_mass_g_mol"],32.042,places=3)
        self.assertIn("hydroxyl", result["material_updates"]["chemistry"]["end_groups"])
        self.assertTrue(any(t["title"] == "User-supplied evidence" for t in result["tables"]))

    def test_input_bounds(self):
        for bad in (-1, 0, 1.5, 65, True, None, float("nan")):
            with self.subTest(dp=bad), self.assertRaises(ValueError):
                petox_candidates(bad)
        self.assertEqual(len(petox_candidates(64)),3)
        with self.assertRaises(ValueError):
            calculate({"chemistry_mode":"petox","selected_ends":"Et/H"})


class SequenceEnsembleTests(unittest.TestCase):
    def test_primary_nmr_benchmark(self):
        r = calculate({"chemistry_mode":"pva_nmr"})["sequence_ensemble"]
        # Budhlall et al., original Poval 217 C-NMR example, equations 1–8.
        self.assertAlmostEqual(r["oh_fraction"], .881)
        self.assertAlmostEqual(r["eta"], .42, places=2)
        self.assertAlmostEqual(r["mean_run_oh"], 20.02, places=2)
        self.assertAlmostEqual(r["mean_run_oac"], 2.71, delta=.01)

    def test_chain_composition_histogram_matches_binomial_bernoulli(self):
        h = .7
        r = pva_ensemble("pva_random",h,[],16,4096,1)
        self.assertAlmostEqual(r["eta"],1)
        var = sum(n*(k-16*h)**2 for k,n in r["histogram"])/4096
        self.assertAlmostEqual(var,16*h*(1-h),delta=.21)
        self.assertAlmostEqual(r["sample_oh_fraction"],h,delta=.012)

    def test_markov_samples_reproduce_target_and_differ_from_random(self):
        r = pva_ensemble("pva_nmr",None,[.837,.088,.075],40,2000,4)
        for target, observed in zip(r["target"],r["observed_dyads"]):
            self.assertAlmostEqual(target, observed,delta=.012)
        random = pva_ensemble("pva_random",r["oh_fraction"],[],40,2000,4)
        self.assertGreater(random["target"][1],2*r["target"][1])

    def test_reproducibility_and_chemical_repeat_labels(self):
        a = calculate({"chemistry_mode":"pva_nmr"})
        b = calculate({"chemistry_mode":"pva_nmr"})
        self.assertEqual(a["sequence_ensemble"], b["sequence_ensemble"])
        seq = a["sequence_ensemble"]["examples"][0]
        drawing = Chem.MolFromSmiles(a["molecular_views"][0]["smiles"])
        self.assertEqual(len(drawing.GetSubstructMatches(Chem.MolFromSmarts("OC(=O)C"))), seq.count("B"))
        self.assertFalse(a["material_updates"])
        self.assertFalse(a["handoffs"])

    def test_boundary_models(self):
        for h, state in ((0,"B"),(1,"A")):
            r = pva_ensemble("pva_random",h,[],20,8,98)
            self.assertIsNone(r["eta"])
            self.assertEqual(r["examples"][0],state*20)
        r = pva_ensemble("pva_nmr",None,[0,1,0],20,8,98)
        self.assertEqual(r["eta"],2)
        self.assertEqual(r["observed_dyads"],[0,1,0])
        r = pva_ensemble("pva_nmr",None,[.8,0,.2],20,8,98)
        self.assertEqual(r["observed_dyads"][1],0)

    def test_invalid_or_unattributed_dyads(self):
        for dyads in ([.2,.3,.4],[-.2,.8,.4],[float("nan"),0,1]):
            with self.subTest(dyads=dyads), self.assertRaises(ValueError):
                pva_ensemble("pva_nmr",0,dyads,4,4,0)
        with self.assertRaisesRegex(ValueError, "source"):
            calculate({"chemistry_mode":"pva_nmr","nmr_oo":.5,"nmr_mixed":.4,"nmr_aa":.1})
        self.assertTrue(calculate({"chemistry_mode":"pva_nmr","nmr_oo":.5,"nmr_mixed":.4,"nmr_aa":.1,"dyad_source":"Hypothetical example."})["plots"])

    def test_optional_runtime_new_presets(self):
        try:
            from polymer_workbench.exchange import make_result
        except ImportError:
            self.skipTest("Workbench is optional for scientific tests.")
        for preset in APP["presets"]:
            with self.subTest(preset=preset["name"]):
                doc = make_result(APP,preset["inputs"],calculate(preset["inputs"]),None)
                self.assertTrue(doc["checksum_sha256"])
