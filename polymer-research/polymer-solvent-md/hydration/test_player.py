"""Check displayed donations against independent, published per-frame counts."""
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
if not (BASE / "pilot").exists():
    BASE = BASE / "mix08"  # separate development checkout
sys.path.insert(0,str(BASE))
from hydration.player import trajectory_frames


class ProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads((BASE/"pilot/study.json").read_bytes())

    def test_archived_frame_identities_and_metrics(self):
        seed = self.record["method"]["seeds"][0]
        movie = trajectory_frames(BASE/"pilot",self.record,seed)
        rep = self.record["analysis"]["replicas"][0]
        refs = self.record["analysis"]["meta"]["carbonyl_oxygens"]
        self.assertEqual(len(movie["frames"]),rep["frames"])
        for n, frame in enumerate(movie["frames"]):
            self.assertAlmostEqual(len(frame["hydrogen_links"])/refs, rep["hydrogen_bonds_per_carbonyl"][n])
            self.assertAlmostEqual(frame["production_time_ps"],rep["time_ps"][n])
            ids = {a["topology_atom_index"]:a for a in frame["atoms"]}
            for o,h,acc in frame["hydrogen_links"]:
                self.assertEqual((ids[o]["group"],ids[h]["element"],ids[acc]["group"]),("water","H","polymer"))

    def test_unknown_seed_rejected(self):
        with self.assertRaises(ValueError):
            trajectory_frames(BASE/"pilot",self.record,2)


if __name__=="__main__":
    unittest.main()
