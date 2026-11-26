"""Regression checks for real computed values, cache isolation and unavailable results."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import math
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from rdkit import Chem
from rdkit.Chem import AllChem

from energy_store import analyze, expected_model_id


class EnergyScreenTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.db = Path(temp.name, "calculations.db")

    def test_coordinates_reproduce_energy_and_radius(self):
        result = analyze("CCO", self.db, count=4)
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["converged"], 4)
        self.assertEqual(result["conformers"][0]["delta_kcal_mol"], 0)
        self.assertEqual(result["model_id"], expected_model_id("OCC", count=4))
        rows = result["conformers"]
        self.assertEqual(sorted(r["energy_kcal_mol"] for r in rows), [r["energy_kcal_mol"] for r in rows])
        geom = result["lowest_geometry"]
        mol = Chem.AddHs(Chem.MolFromSmiles(result["canonical_smiles"]))
        conf = Chem.Conformer(mol.GetNumAtoms())
        self.assertEqual(len(geom["atoms"]), 9)
        self.assertEqual(len(geom["bonds"]), 8)
        for i, atom in enumerate(geom["atoms"]):
            self.assertEqual(mol.GetAtomWithIdx(i).GetSymbol(), atom["element"])
            conf.SetAtomPosition(i, (atom["x"], atom["y"], atom["z"]))
        mol.AddConformer(conf, assignId=True)
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        energy = AllChem.MMFFGetMoleculeForceField(mol, props).CalcEnergy()
        self.assertAlmostEqual(rows[0]["energy_kcal_mol"], energy, places=7)
        masses = [atom.GetMass() for atom in mol.GetAtoms()]
        weighted = lambda values: sum(m*v for m, v in zip(masses, values)) / sum(masses)
        com = [weighted(a[axis] for a in geom["atoms"]) for axis in "xyz"]
        rg = math.sqrt(weighted(sum((a[k]-c)**2 for k, c in zip("xyz", com)) for a in geom["atoms"]))
        self.assertAlmostEqual(rows[0]["rg_angstrom"], rg, places=7)

    def test_canonical_graph_cache_and_settings(self):
        r = analyze("COC(C)C", self.db, count=3)
        with patch("energy_store._screen", side_effect=AssertionError("cache hit must skip optimizer")):
            cached = analyze("CC(OC)C", self.db, count=3)
        self.assertTrue(cached.pop("cache_hit"))
        self.assertFalse(r.pop("cache_hit"))
        self.assertEqual(cached, r)
        for n, seed in [(2, 2403), (3, 1)]:
            response = analyze("COC(C)C", self.db, count=n, seed=seed)
            self.assertFalse(response["cache_hit"])
            self.assertNotEqual(response["model_id"], r["model_id"])

    def test_sample_is_reproducible_in_separate_databases(self):
        first = analyze("C1CCCCC1", self.db, count=4, seed=42)
        second = analyze("C1CCCCC1", ":memory:", count=4, seed=42)
        self.assertEqual(first, second)

    def test_specified_and_unspecified_stereochemistry(self):
        ccw = analyze("O[C@H](F)C", self.db, count=2)
        cw = analyze("O[C@@H](F)C", self.db, count=2)
        for result in (ccw, cw, analyze("F/C=C/F", self.db, count=2)):
            self.assertEqual(result["status"], "ok")
        self.assertNotEqual(cw["model_id"], ccw["model_id"])
        for smiles in ("OC(F)C", "FC=CF", "[Na+]", "[CH3]"):
            unavailable = analyze(smiles, self.db, count=2)
            self.assertEqual(unavailable["status"], "unavailable")
            self.assertFalse(unavailable["conformers"])
            self.assertIsNone(unavailable["lowest_geometry"])
        with self.assertRaises(ValueError):
            analyze("C.C", self.db)

    def test_large_valid_structure_returns_cached_unavailable(self):
        smiles = "C" * 101
        with patch.object(AllChem, "EmbedMultipleConfs", side_effect=AssertionError("must not embed")):
            result = analyze(smiles, self.db, count=2)
            self.assertEqual(result["status"], "unavailable")
            self.assertEqual(result["canonical_smiles"], smiles)
            self.assertEqual(result["model_id"], expected_model_id(smiles, count=2))
            self.assertTrue(analyze(smiles, self.db, count=2)["cache_hit"])
        with self.assertRaises(ValueError):
            expected_model_id("C" * 151)

    def test_nonconverged_candidates_are_excluded(self):
        original = AllChem.MMFFOptimizeMoleculeConfs

        def one_failure(*args, **kwargs):
            results = list(original(*args, **kwargs))
            results[0] = (1, -10000)
            return results

        with patch.object(AllChem, "MMFFOptimizeMoleculeConfs", side_effect=one_failure):
            result = analyze("CCCO", self.db, count=3)
        self.assertEqual(result["status"], "partial")
        self.assertEqual(result["converged"], 2)
        self.assertNotIn(0, [r["id"] for r in result["conformers"]])
        self.assertEqual(result["conformers"][0]["delta_kcal_mol"], 0.0)

    def test_concurrent_requests_leave_one_complete_row(self):
        with ThreadPoolExecutor(max_workers=3) as pool:
            results = list(pool.map(lambda _: analyze("CCO", self.db, count=3), range(3)))
        self.assertEqual(sum(not r["cache_hit"] for r in results), 1)
        self.assertEqual(len(set(r["model_id"] for r in results)), 1)
        with sqlite3.connect(self.db) as db:
            self.assertEqual(db.execute("SELECT count(*) FROM conformer_runs").fetchone()[0], 1)
        self.assertTrue(all(r["converged"] == 3 for r in results))


if __name__ == "__main__":
    unittest.main()
