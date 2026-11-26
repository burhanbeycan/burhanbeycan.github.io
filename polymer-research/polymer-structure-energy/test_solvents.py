from pathlib import Path
import tempfile
import unittest

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

from energy_store import analyze, expected_model_id
from solvent_registry import (
    COUNT, SEED, analyze_with_registry, build_registry, get_component,
    load_registry, lookup_energy,
)


class SolventRegistryTests(unittest.TestCase):
    def test_all_exported_geometries_have_traceable_model_and_stoichiometry(self):
        rows = load_registry()["components"]
        formulas = {
            "water": "H2O", "ethanol": "C2H6O", "acetonitrile": "C2H3N",
            "acetone": "C3H6O", "thf": "C4H8O", "dmf": "C3H7NO", "dmso": "C2H6OS",
        }
        self.assertEqual({s["id"]: s["formula"] for s in rows}, formulas)
        self.assertEqual(len(rows), len({s["model_id"] for s in rows}))
        for record in rows:
            with self.subTest(solvent=record["id"]):
                mol = Chem.MolFromSmiles(record["smiles"])
                raw = record["energy_record"]
                self.assertEqual(rdMolDescriptors.CalcMolFormula(mol), record["formula"])
                self.assertAlmostEqual(Descriptors.MolWt(mol), record["molar_mass_g_mol"], places=3)
                self.assertEqual(len(Chem.AddHs(mol).GetAtoms()), len(record["geometry"]["atoms"]))
                self.assertEqual(raw["model_id"], record["model_id"])
                self.assertEqual(raw["lowest_geometry"], record["geometry"])
                self.assertEqual(raw["protocol"], record["protocol"])
                self.assertGreater(record["converged"], 0)
                self.assertEqual(record["converged"], len(record["conformers"]))
                self.assertEqual(min(c["delta_kcal_mol"] for c in record["conformers"]), 0.0)
                self.assertEqual(record["protocol"]["force_field"], "MMFF94s")
                self.assertEqual(record["protocol"]["energy_unit"], "kcal/mol")
                self.assertNotIn("cache_hit", raw)
        with self.assertRaises(KeyError):
            get_component("absent")

    def test_portable_hit_requires_matching_graph_and_method(self):
        with tempfile.TemporaryDirectory() as directory:
            registry = build_registry(Path(directory) / "generated.sqlite3")
            trial = lookup_energy("OCC", registry=registry)
            self.assertTrue(trial["cache_hit"])
            self.assertEqual(trial["cache_source"], "bundled_solvent_registry")
            self.assertEqual(trial["model_id"], expected_model_id("CCO", count=COUNT, seed=SEED))
            self.assertIsNone(lookup_energy("CCO", seed=SEED+1, registry=registry))
            self.assertIsNone(lookup_energy("CCO", count=COUNT-1, registry=registry))
            self.assertIsNone(lookup_energy("CCCO", registry=registry))
            self.assertIsNone(lookup_energy("[13CH3]CO", registry=registry))
            # The archived geometry is not mutated when a caller modifies its copy.
            trial["conformers"][0]["delta_kcal_mol"] = 100
            self.assertEqual(lookup_energy("CCO", registry=registry)["conformers"][0]["delta_kcal_mol"], 0)

    def test_regular_entrypoint_uses_local_sqlite_on_protocol_miss(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fallback.sqlite3"
            baseline = get_component("ethanol")
            if baseline["protocol"]["rdkit_version"] == rdBase.rdkitVersion and baseline["model_id"] == expected_model_id("OCC"):
                match = analyze_with_registry("OCC", path)
                self.assertFalse(path.exists())
                self.assertEqual(match["cache_source"], "bundled_solvent_registry")
            first = analyze_with_registry("CCO", path, seed=117)
            second = analyze_with_registry("OCC", path, seed=117)
            self.assertTrue(path.exists())
            self.assertEqual(first["cache_source"], "computed")
            self.assertEqual(second["cache_source"], "sqlite")
            self.assertEqual(first["model_id"], second["model_id"])
            self.assertEqual(first["conformers"], second["conformers"])

    def test_fresh_result_matches_archived_ethanol_with_compatible_runtime(self):
        entry = get_component("ethanol")
        if entry["model_id"] != expected_model_id(entry["smiles"]):
            self.skipTest("Kurulu RDKit/protokol ile kayıtlı model sürümleri farklı.")
        with tempfile.TemporaryDirectory() as directory:
            result = analyze(entry["smiles"], Path(directory) / "fresh.sqlite3")
        self.assertEqual(result["converged"], entry["converged"])
        self.assertAlmostEqual(result["conformers"][0]["energy_kcal_mol"], entry["lowest_energy_kcal_mol"], places=5)


if __name__ == "__main__":
    unittest.main()
