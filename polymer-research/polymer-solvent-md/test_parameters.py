"""Regression tests for protocol identity, scientific gates and published parameters."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import parameters as p


class ParameterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = json.loads((p.BUNDLE.parent / "source_model.json").read_text())["energy"]
        cls.smiles = cls.source["canonical_smiles"]

    def test_exact_bundled_reuse_and_input(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(p, "_amber_home", return_value=None):
            record = p.parameterize(self.smiles, tmp, source_record=self.source)
            self.assertEqual(record["cache_source"], "bundled")
            self.assertEqual(record["native_validation"]["num_atoms"], 50)
            self.assertLess(abs(record["native_validation"]["native_charge_e"]), 1e-7)
            self.assertAlmostEqual(record["charge_balance"]["raw_net_e"], -0.001998)
            self.assertAlmostEqual(record["charge_balance"]["uniform_offset_per_atom_e"], 0.00003996)
            _, sdf = p._source_input(self.smiles, self.source)
            self.assertEqual((p.BUNDLE / "solute.sdf").read_text(), sdf)

    def test_other_coordinates_cannot_consume_pilot(self):
        changed = copy.deepcopy(self.source)
        changed["lowest_geometry"]["atoms"][0]["x"] += .001
        with tempfile.TemporaryDirectory() as tmp, patch.object(p, "_amber_home", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "AmberTools"):
                p.parameterize(self.smiles, tmp, source_record=changed)

    def test_unsupported_input_is_rejected_before_external_programs(self):
        with tempfile.TemporaryDirectory() as tmp:
            for smiles in ("[NH4+]", "CC.[Cl-]", "CS", "CC(O)CC", "C[CH2]", "[13CH3]C"):
                with self.subTest(smiles=smiles), self.assertRaises(ValueError):
                    p.parameterize(smiles, tmp)

    def test_known_analogies_and_missing_parameter(self):
        data = (p.BUNDLE / "solute.frcmod").read_text()
        audit = p._supplements(data)
        self.assertEqual(audit["max_reported_analogy_penalty"], 6)
        self.assertEqual(audit["supplementary_rows"], {"DIHE": 6, "IMPROPER": 4})
        with self.assertRaises(ValueError):
            p._supplements(data + "ATTN, need revision")
        with self.assertRaises(ValueError):
            p._supplements(data + "penalty score= 14.0")

    def test_modified_native_file_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            b = shutil.copytree(p.BUNDLE, Path(tmp) / "parameters")
            with (b / "solute.prmtop").open("a") as output:
                output.write("changed")
            d, _ = p._source_input(self.smiles, self.source)
            with self.assertRaisesRegex(ValueError, "bütünlük"):
                p._cached(b, d)


if __name__ == "__main__":
    unittest.main()
