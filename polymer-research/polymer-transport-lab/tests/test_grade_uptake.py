import csv
import io
import unittest
from polymer_transport_lab import grade_uptake as g
from polymer_transport_lab.plugin import calculate


def records(filename):
    return g.DATA.joinpath(filename).read_text()


def alter(body, func):
    reader=csv.DictReader(io.StringIO(body)); stream=io.StringIO()
    writer=csv.DictWriter(stream,reader.fieldnames,lineterminator="\n")
    writer.writeheader()
    for row in reader:
        func(row)
        writer.writerow(row)
    return stream.getvalue()


class GradeUptakeTests(unittest.TestCase):
    def setUp(self):
        self.spec = records("pva_grade_spec.csv")
        self.obs = records("pva_grade_obs.csv")
        self.xyz = records("pva_site_fixture.csv")

    def test_comparable_three_models_and_grade_holdout(self):
        cov, fold, pred, skip = g.compare(self.spec,self.obs,"grade","synthetic")
        self.assertEqual((len(cov),len(fold),len(pred),skip),(12,36,108,[]))
        for row in fold:
            self.assertEqual(row[2:4],[11,1])
        self.assertEqual({r[2] for r in pred if r[0]=="fixture-G01"},{"fixture-F01"})

    def test_holdout_predictions_do_not_see_targets(self):
        changed = alter(self.obs,lambda row:row.update(water_g_g=2) if row["study"]=="fixture-paper1" else None)
        def results(body):
            return [r[6] for r in g.compare(self.spec,body,"study","synthetic")[2] if r[0]=="fixture-paper1"]
        self.assertEqual(results(self.obs),results(changed))

    def test_unseen_study_reused_grade_purged_only_under_joint_holdout(self):
        changed=alter(self.spec,lambda row:row.update(grade="fixture-G01") if row["specimen"]=="fixture-F02" else None)
        def train(policy):
            return g.compare(changed,self.obs,policy,"synthetic")[1][0][2]
        self.assertEqual(train("study"),8)
        self.assertEqual(train("study_new_grade"),7)

    def test_missing_md_excluded_from_all_models(self):
        spec = alter(self.spec,lambda row:row.update(water_per_oh="",water_per_oac="") if row["specimen"]=="fixture-F01" else None)
        cov,fold,pred,skip = g.compare(spec,self.obs,"grade","synthetic")
        self.assertEqual(cov[0][-1],"excluded")
        self.assertFalse(any(r[2]=="fixture-F01" for r in pred))
        self.assertEqual(len(fold),33)

    def test_invalid_sequences_and_outcome_conditioned_md(self):
        for edit in [dict(mixed_dyad=.95),dict(md_origin="conditioned_by_measured_uptake"),
                     dict(md_wref_g_g=.15),dict(nmr_source="")]:
            spec=alter(self.spec,lambda row:row.update(edit) if row["specimen"]=="fixture-F01" else None)
            with self.assertRaises(ValueError):g.compare(spec,self.obs,"grade","synthetic")
        with self.assertRaises(ValueError):g.compare(self.spec,self.obs,"grade","user_supplied")

    def test_missing_groups_and_composition(self):
        spec=alter(self.spec,lambda r:r.update(grade="one-grade"))
        cov,fold,pred,skips=g.compare(spec,self.obs,"grade","synthetic")
        self.assertEqual((fold,pred,skips),([],[],["one-grade"]))
        with self.assertRaises(ValueError):g.read_observations(self.obs+"\n"+self.obs.splitlines()[1],g.read_specimens(self.spec)[0])

    def test_measured_evidence_has_no_artificial_scores(self):
        result=calculate({"mode":g.MODE})
        self.assertEqual(result["metrics"][1]["value"],0)
        self.assertEqual(result["plots"],[])
        self.assertTrue(any("Nakajima" in line for line in result["notes"]))

    def test_site_hydration_replica_weighting_and_periodic_wrap(self):
        per,mean = g.hydration_frames(self.xyz,.35)
        self.assertEqual(mean,[2,.5])
        self.assertEqual(per,[["packing-A",2,2,1],["packing-B",1,2,0]])
        shifted=alter(self.xyz,lambda r:r.update(x_nm=float(r["x_nm"])+4) if r["atom"]=="W1" else None)
        self.assertEqual(g.hydration_frames(shifted,.35)[1],mean)

    def test_invalid_trajectory_geometry_and_roles(self):
        for edit in [dict(Lx_nm=.6),dict(role="N"),dict(atom="W1")]:
            xyz=alter(self.xyz,lambda row:row.update(edit) if row["atom"]=="P1" and row["frame"]=="1" else None)
            with self.assertRaises(ValueError):g.hydration_frames(xyz,.35)

    def test_full_exchange_models(self):
        for dataset in ("csv","coordinates","papers"):
            out=calculate({"mode":g.MODE, "grade_data":dataset})
            self.assertIn("material_updates",out)
            self.assertTrue(out["metrics"])


if __name__=="__main__":unittest.main()
