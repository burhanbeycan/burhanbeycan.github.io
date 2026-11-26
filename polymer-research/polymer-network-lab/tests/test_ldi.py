"""Reference values and input boundaries for the LDI polyurethane workflows."""
import copy
import io
import json
import unittest

import numpy as np

from polymer_network_lab.ldi import (
    DEMO_SERIES, SERIES_NOTE, contact_analysis, contact_folds, feed_scenario, fixture,
)
from polymer_network_lab.plugin import APP, calculate
from polymer_network_lab.trajectory import gro_frames


class LdiChemistry(unittest.TestCase):
    def test_urethane_feed_and_bounded_lysine_reaction(self):
        b=feed_scenario(.05,1,1.65,.5,"BDO","",.2)
        g=feed_scenario(.05,1,1.65,.5,"glycerol",.3,.2)
        self.assertAlmostEqual(b["index"],106.45161290)
        self.assertAlmostEqual(g["index"],91.66666667)
        self.assertIsNone(b["urea_upper"])
        self.assertEqual(g["urea_upper"],.3)
        self.assertEqual(g["triol"],.5)
        self.assertAlmostEqual(feed_scenario(.05,1,1.65,.5,"PDO",.5,.2)["urea_upper"],.4)
        with self.assertRaises(ValueError):
            feed_scenario(.05,1,1.65,.5,"BDO",4,.2)

    def test_contact_counts_normalization_and_intermolecular_exclusion(self):
        m=fixture()
        r=contact_analysis(m,m["frames"])
        self.assertEqual(r["total"],5)
        self.assertAlmostEqual(r["lag_survival"],2/3)
        self.assertEqual(r["rows"][0],["PEG_ether",2,2,3,.75,.5,1.2])
        self.assertEqual(r["rows"][1],["PCL_ester",1,1,2,.5,.25,1.6])
        self.assertEqual(r["rows"][2][1:4],[2,1,0])
        for frame in m["frames"]:
            frame["xyz_nm"]=((np.asarray(frame["xyz_nm"])+.5)%2).tolist()
        self.assertEqual(contact_analysis(m,m["frames"])["rows"],r["rows"])

    def test_contacts_without_events_and_angular_filter(self):
        m=fixture()
        for frame in m["frames"]:
            frame["xyz_nm"][2]=[1.85,.57,.5]  # Close to N, but points behind H.
        r=contact_analysis(m,m["frames"])
        self.assertEqual(r["rows"][0][3],0)
        no=contact_analysis(m,m["frames"],.05)
        self.assertEqual(no["total"],0)
        self.assertIsNone(no["rows"][0][-1])
        self.assertIsNone(no["lag_survival"])

    def test_contact_input_contract(self):
        for change in ["same chain","bad frame","nonuniform","duplicate oxygen","bad bond","no chains","no frames"]:
            with self.subTest(change=change),self.assertRaises(ValueError):
                m=fixture()
                if change=="same chain": m["atom_chain"]=["A"]*7
                if change=="bad frame": m["frames"][0]["box_nm"]=[.5]*3
                if change=="nonuniform": m["frames"][1]["time_ps"]=.5
                if change=="duplicate oxygen":m["acceptors"].append(m["acceptors"][0])
                if change=="bad bond":m["frames"][0]["xyz_nm"][1][0]=1.0
                if change=="no chains":m["atom_chain"]=None
                if change=="no frames":m["frames"]=[]
                contact_analysis(m,m["frames"])

    def test_gro_nm_and_wrapped_indices(self):
        m=fixture()
        lines=[]
        for f in m["frames"]:
            lines.extend(["test",str(len(m["atom_chain"]))])
            for i,(x,y,z) in enumerate(f["xyz_nm"]):
                lines.append(f"{1:5d}{'MON':<5}{'A':>5}{i+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
            lines.append("2.00 2.00 2.00")
        g="\n".join(lines)
        self.assertEqual(contact_analysis(m,gro_frames(io.StringIO(g),1.0))["rows"],
                         contact_analysis(m,m["frames"])["rows"])
        with self.assertRaises(ValueError):list(gro_frames(io.StringIO(g.replace("MON", "NEW", 1)),1))
        with self.assertRaises(ValueError):list(gro_frames(io.StringIO(g)))
        with self.assertRaises(ValueError):list(gro_frames(io.StringIO(g.replace("2.00 2.00 2.00","2 2 2 0 0 .3 0 0 0")),1))

    def test_holdout_reference_and_blinding(self):
        r=contact_folds(DEMO_SERIES)
        np.testing.assert_allclose(r["mae_pp"],[5.5891679759,1.4534656635])
        modified=DEMO_SERIES.replace("0.65,0.35,sim1","0.65,0.99,sim1")
        np.testing.assert_allclose(contact_folds(modified)["rows"][0][3:],r["rows"][0][3:])
        for bad in [DEMO_SERIES.replace("example8","example7"),DEMO_SERIES.replace("sim8","sim7"),
                    DEMO_SERIES.replace("F1-a","F2-a"),
                    DEMO_SERIES.replace("illustrative dry films","different protocol",1)]:
            with self.assertRaises(ValueError):
                contact_folds(bad)
        repeated=DEMO_SERIES+"\n"+DEMO_SERIES.splitlines()[1].replace("F1-a","F1-b")
        self.assertEqual(contact_folds(repeated)["rows"][0][1],2)
        self.assertEqual(contact_folds(repeated)["mae_pp"],r["mae_pp"])

    def test_plugin_routes_and_sources(self):
        presets=[p for p in APP["presets"] if p["inputs"].get("network_mode","").startswith("ldi_")]
        self.assertEqual(len(presets),3)
        for p in presets:
            result=calculate(p["inputs"])
            json.dumps(result,allow_nan=False)
            self.assertTrue(result["metrics"] and result["tables"])
        with self.assertRaises(ValueError):
            calculate({"network_mode":"ldi_folds","ldi_ftirs":DEMO_SERIES+"\n"+DEMO_SERIES.splitlines()[1],
                       "ldi_ftir_source":SERIES_NOTE})
        with self.assertRaises(ValueError):
            calculate({"network_mode":"ldi_contacts","ldi_traj":"{\"frames\":[],\"atom_chain\":null}","ldi_traj_source":"manual"})
