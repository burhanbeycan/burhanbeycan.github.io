"""The local process adapter can be exercised without installed MD engines."""
import json
import os
from pathlib import Path
import tempfile
import unittest

from polymer_network_lab.bulk import build_plan, execute

STUB="""#!/usr/bin/env python3
import sys
from pathlib import Path
a=sys.argv[1:]
if a[0]=='gen_coords': Path(a[a.index('-o')+1]).write_text('stub packing')
elif a[0]=='grompp':
 if '-t' in a: assert Path(a[a.index('-t')+1]).is_file()
 for flag in ('-o','-po','-pp'): Path(a[a.index(flag)+1]).write_text('stub preprocessing')
elif a[0]=='mdrun':
 for suffix in ('.gro','.cpt'): Path(a[a.index('-deffnm')+1]+suffix).write_text('stub coordinates/checkpoint')
else: sys.exit(9)
"""


class BulkAdapter(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        system=self.root/'bundle';system.mkdir()
        for f in ('topol.top','polymer.itp','nvt.mdp','prod.mdp','a.gro'):
            (system/f).write_text('test input')
        self.config=self.root/'config.json'
        self.spec={"evidence":"Synthetic adapter test only",
                   "stages":[{"id":"nvt","mdp":"nvt.mdp"},{"id":"prod","mdp":"prod.mdp","continue_from_previous_checkpoint":True}],
                   "replicas":[{"id":"prebuilt","system_dir":"bundle","topology":"topol.top","coordinates":"a.gro"},
                               {"id":"pack","system_dir":"bundle","topology":"topol.top","polyply_box_nm":[5,5,5]}]}
        self.out=self.root/'runs'

    def plan(self,tool="gmx"):
        self.config.write_text(json.dumps(self.spec))
        return build_plan(self.config,self.out,tool,tool)

    def test_dry_plan_is_side_effect_free(self):
        p=self.plan()
        self.assertFalse(self.out.exists())
        self.assertEqual(len(p["replicas"]),2)
        self.assertIn("polymer.itp",p["replicas"][0]["inputs_sha256"])
        self.assertIn("-t",p["replicas"][0]["steps"][2]["argv"])
        self.assertEqual(p["replicas"][1]["steps"][0]["argv"][1],"gen_coords")

    @unittest.skipIf(os.name=="nt","POSIX test executable")
    def test_stub_roundtrip_and_fresh_directory(self):
        tool=self.root/'stub';tool.write_text(STUB);tool.chmod(0o755)
        p=self.plan(str(tool));execute(p,self.out)
        m=json.loads((self.out/'manifest.json').read_text())
        self.assertTrue(all(x['status']=='completed' for x in m['replicas']))
        self.assertTrue((self.out/'pack'/'packed.gro').exists())
        self.assertTrue((self.out/'prebuilt'/'stage-prod.cpt').exists())
        with self.assertRaises(FileExistsError):execute(p,self.out)
        self.assertEqual(m['replicas'][1]['steps'][0]['status'],'completed')

    def test_missing_executable_records_failed_stage(self):
        with self.assertRaises(FileNotFoundError):execute(self.plan(str(self.root/'absent')),self.out)
        m=json.loads((self.out/'manifest.json').read_text())
        self.assertEqual(m['status'],'failed')
        self.assertEqual(m['replicas'][1]['status'],'planned')

    def test_conflicting_coordinate_sources_and_path_traversal(self):
        self.spec['replicas'][0]['polyply_box_nm']=[2,2,2]
        with self.assertRaises(ValueError):self.plan()
        del self.spec['replicas'][0]['polyply_box_nm']
        self.spec['replicas'][0]['topology']='../config.json'
        with self.assertRaises(ValueError):self.plan()
