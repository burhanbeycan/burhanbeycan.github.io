"""End-to-end identity, provenance, and artifact reuse across callers."""
import tempfile
import unittest
from pathlib import Path

from app import calculate, find_record


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.store = tempfile.TemporaryDirectory()
        self.addCleanup(self.store.cleanup)
        self.db = Path(self.store.name) / "cache" / "models.db"

    def test_same_calculation_reused_across_different_recipes(self):
        linear = calculate({"mode": "repeat", "smiles": "[*:1]CC[*:2]", "dp": 3}, self.db)
        explicit = calculate({"mode": "finite", "smiles": "CCCCCC"}, self.db)
        self.assertEqual(linear["structure"]["structure_id"], explicit["structure"]["structure_id"])
        self.assertNotEqual(linear["structure"]["recipe_id"], explicit["structure"]["recipe_id"])
        self.assertEqual(linear["energy"]["model_id"], explicit["energy"]["model_id"])
        self.assertTrue(explicit["energy"]["cache_hit"])
        self.assertNotEqual(linear["artifact_id"], explicit["artifact_id"])
        replay = calculate({"mode": "repeat", "smiles": "[*:1]CC[*:2]", "dp": 3}, self.db)
        self.assertEqual(replay["artifact_id"], linear["artifact_id"])
        portable = find_record(self.db, replay["artifact_id"])
        self.assertEqual(portable["structure"]["formula"], "C6H14")
        self.assertNotIn("cache_hit", portable["energy"])

    def test_solvent_uses_same_protocol_and_unknown_stereo_has_no_energy(self):
        cached = calculate({"mode": "finite", "smiles": "OCC"}, self.db)
        self.assertEqual(cached["energy"]["cache_source"], "bundled_solvent_registry")
        blocked = calculate({"mode": "finite", "smiles": "CC(O)Cl"}, self.db)
        self.assertGreater(blocked["structure"]["unspecified_stereo_count"], 0)
        self.assertEqual(blocked["energy"]["status"], "unavailable")
        self.assertEqual(blocked["energy"]["conformers"], [])
