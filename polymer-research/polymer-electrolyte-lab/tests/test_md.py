import unittest

import numpy as np

from polymer_electrolyte_lab.displacements import displacement_result, ion_moments
from polymer_electrolyte_lab.md_views import contact_example, moment_example, published_shell
from polymer_electrolyte_lab.shell import orthorhombic_contacts, shell_runs, summarize_shell, survival


class OxygenCoordination(unittest.TestCase):
    def test_two_primary_reference_compositions(self):
        for ratio, refs in [(10,(5.40,5.33,1.17,-1.24)),(18,(5.74,5.42,.79,-1.11))]:
            report=published_shell(ratio)
            rows=report["tables"][0]["rows"]
            np.testing.assert_allclose([rows[0][3],rows[1][3],*[v["value"] for v in report["metrics"][:2]]],refs)

    def test_reentry_breaks_continuous_identity(self):
        a,b=(0,1),(0,2)
        # Only first and second observations are consecutive; a returns at frame 3.
        self.assertEqual(survival([{a},{a},{b},{a}],2),[1,1/3,0])
        self.assertEqual(survival([set(),set(),set()],2),[None,None,None])

    def test_replicas_are_independent_units(self):
        record=contact_example()
        record["runs"]=[record["runs"][0]]
        self.assertIsNone(summarize_shell(record)["tables"][0]["rows"][0][4])
        record["runs"]*=2
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            shell_runs(record)

    def test_contacts_have_stable_atom_indices(self):
        rec=contact_example()
        rec["runs"][0]["frames"][0]["tfsi"]=[[0,0],[0,0]]
        with self.assertRaisesRegex(ValueError,"Duplicate pair"):
            shell_runs(rec)
        rec["runs"][0]["frames"][0]["tfsi"]=[[0,4]]
        with self.assertRaisesRegex(ValueError,"oxygen index"):
            shell_runs(rec)

    def test_missing_silent_frame_and_mixed_state_are_rejected(self):
        rec=contact_example()
        del rec["runs"][0]["frames"][0]["ether"]
        with self.assertRaisesRegex(ValueError,"each frame"):
            shell_runs(rec)
        rec=contact_example()
        rec["runs"][1]["state_id"]="other temperature"
        with self.assertRaisesRegex(ValueError,"same material"):
            shell_runs(rec)

    def test_genuine_minimum_image_with_distinct_oxygens(self):
        li=[[[0,0,0]]]; ether=[[[.99,0,0],[.4,0,0]]]; tfsi=[[[.1,0,0],[.99,0,0],[.49,0,0],[.47,0,0]]]
        self.assertEqual(orthorhombic_contacts(li,ether,tfsi,[1,1,1],.15),
                         [{"ether":[[0,0]],"tfsi":[[0,0],[0,1]]}])
        with self.assertRaises(ValueError):
            orthorhombic_contacts(li,ether,tfsi,[.1,1,1],.15)


class ChargeMotion(unittest.TestCase):
    def test_half_collective_for_analytic_fixture(self):
        rows=displacement_result(moment_example())["tables"][1]["rows"]
        np.testing.assert_allclose(rows[0][2:],(1e-11,2e-11,.18592487783177,.37184975566354,.5))

    def test_ballistic_curves_cannot_produce_conductivity(self):
        rec=moment_example()
        for row in rec["runs"][0]["moments"]:
            row[3]*=row[0]/80
        row=displacement_result(rec)["tables"][1]["rows"][0]
        self.assertIsNone(row[4])
        self.assertIsNotNone(row[5])

    def test_nonstationary_halves_are_withheld(self):
        rec=moment_example()
        r=rec["runs"][0]["moments"]
        for row in r:
            row[3]*=min(1, 40/row[0])
        self.assertIsNone(displacement_result(rec)["tables"][1]["rows"][0][4])

    def test_matched_transport_frame_required(self):
        rec=moment_example()
        rec["transport_charge"]="force_field_partial"
        with self.assertRaisesRegex(ValueError,"formal"):
            displacement_result(rec)

    def test_ion_positions_and_frame_invariance(self):
        t=np.arange(8)
        r=np.zeros((len(t),2,3)); r[:,0,0]=t; r[:,1,0]=-t/2
        c=np.zeros((len(t),3)); c[:,2]=10*t
        values=ion_moments(r+c[:,None,:],c,1,1,[1,3])
        np.testing.assert_allclose(values,[[1,1,.25,2.25],[3,9,2.25,20.25]])
        with self.assertRaises(ValueError):
            ion_moments(r,c,1,2,[1])

    def test_moment_bound_and_lag_order_checked(self):
        rec=moment_example();rec["runs"][0]["moments"][1][3]=10
        with self.assertRaisesRegex(ValueError,"bound"):
            displacement_result(rec)
        rec=moment_example();rec["runs"][0]["moments"][1][0]=10
        with self.assertRaisesRegex(ValueError,"increasing"):
            displacement_result(rec)
