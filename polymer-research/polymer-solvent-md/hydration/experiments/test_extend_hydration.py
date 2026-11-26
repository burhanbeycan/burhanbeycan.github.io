"""Regression checks for lineage integrity and exact graph reconstruction."""
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
from openmm import app
import extend_hydration as e


def base_path():
    paths=[Path(__file__).resolve().parent.parent.parent,Path("mix08")]
    if "POLYMER_MD_BASE" in os.environ:
        paths.insert(0,Path(os.environ["POLYMER_MD_BASE"]))
    return next(p for p in paths if (p/"amber_bridge.py").is_file())


class TestRestart(unittest.TestCase):
    def test_seed_mapping_and_integration_units(self):
        streams=[i for s in (13031,27011,49033) for i in e.seed_stream(s)]
        self.assertEqual(len(set(streams)),9)
        self.assertTrue(all(0 < i < 2**31 for i in streams))
        self.assertEqual(e.steps(1,.002),500)
        for t in (-1,.0027,0):
            with self.assertRaises(ValueError):
                e.steps(t,.002)

    def test_solute_and_water_graph_reconstructed_in_exact_order(self):
        base=base_path()
        params=base/"pilot/parameters"
        topo=json.loads((base/"pilot/run/topology.json").read_bytes())
        amber=app.AmberPrmtopFile(str(params/"solute.prmtop"))
        system=e.rebuild_topology(amber,topo,np.diag([2.5]*3))
        self.assertEqual(system.getNumAtoms(),1550)
        self.assertEqual(system.getNumBonds(),1049)
        mutated=dict(topo,atom_names=list(topo["atom_names"]))
        mutated["atom_names"][0]="WRONG"
        with self.assertRaisesRegex(ValueError,"atom kimliği"):
            e.rebuild_topology(app.AmberPrmtopFile(str(params/"solute.prmtop")),mutated,np.diag([2.5]*3))

    def test_incomplete_or_changed_checkpoint_chain_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            e.save_npz(d/"00000.npz",run_id="original",stop_frame=12)
            self.assertEqual(e.ready_chunks(d,"original"),[d/"00000.npz"])
            with self.assertRaisesRegex(ValueError,"kimliği"):
                e.ready_chunks(d,"changed")
            shutil.copyfile(d/"00000.npz",d/"00002.npz")
            with self.assertRaisesRegex(ValueError,"sırası"):
                e.ready_chunks(d,"original")


if __name__=="__main__":
    unittest.main()
