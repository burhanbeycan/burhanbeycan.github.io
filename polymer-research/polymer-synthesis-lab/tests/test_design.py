"""Prospective acquisition, support gates and experimental provenance."""
import csv
import io
import math
import unittest

import numpy as np

from polymer_synthesis_lab.crop import HEADER, Recipe, load_cohort, parse_cohort
from polymer_synthesis_lab.design import POOL, candidates, design, features, fit, folds, series_id
from polymer_synthesis_lab.plugin import calculate


def synthetic(change_paper=None):
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(HEADER)
    for p in ("mock-A", "mock-B"):
        for m in (1, 2, 4):
            for feed in (20, 40, 80):
                r = Recipe(m, m / feed, 100, 10, "Me/OH")
                y = r.predict()[1] * math.exp(features(m, feed) @ [0.1, 0.2, -0.3])
                writer.writerow([p, f"M{m}-F{feed}", "1", "mock-batch", "mock-SEC", "Me/OH",
                                 m, m/feed, 100, 10, y*2 if p==change_paper else y, 0.57])
    return output.getvalue()


class DesignTests(unittest.TestCase):
    def test_real_cohort_score_and_support(self):
        rows, ranked, model, why = design(load_cohort("flow_feed"), POOL, 160, 5)
        self.assertEqual(len(rows),4)
        self.assertIsNone(model)
        self.assertIn("varied concentration",why)
        self.assertEqual(len(ranked),6)
        recipe, gain, empirical = ranked[0]
        self.assertEqual((recipe.M,recipe.feed),(2,120))
        self.assertAlmostEqual(recipe.I,1/60)
        self.assertAlmostEqual(gain,.3939193183944)
        self.assertIsNone(empirical)
        report=calculate({"model":"petox_design"})
        scores = next(t for t in report["tables"] if t["title"].startswith("Blind validation"))["rows"]
        self.assertAlmostEqual(scores[0][3],49.00744783, places=7)
        self.assertEqual(scores[2][2:], [0, None, None])

    def test_inputs_at_unsupported_condition_and_empty_pool(self):
        with self.assertRaisesRegex(ValueError,"No measured"):
            design(load_cohort("flow_feed"),POOL,140,5)
        rows,ranked,_,_=design(load_cohort("flow_feed"), "M0_M,feed\n4,30",160,5)
        self.assertEqual(ranked,[])
        for raw in ("M0_M,feed\n1,0","M0_M,feed\n2,nan", "M0_M,feed\n2,30\n2,30","different,header"):
            with self.assertRaises(ValueError):
                candidates(raw,100,10,"Me/OH")

    def test_completed_measurement_changes_next_choice(self):
        raw=load_cohort("flow_feed")
        sample = list(csv.DictReader(io.StringIO(raw)))[-1]
        sample.update(paper="mock-design-update", series="illustrative", M0_M="2", I0_M=str(2/120), Mn_g_mol="8100")
        output=io.StringIO()
        csv.DictWriter(output,HEADER,lineterminator="\n").writerow(sample)
        rows,ranked,_,_=design(raw.strip()+"\n"+output.getvalue(),POOL,160,5)
        self.assertEqual(len(rows),5)
        self.assertNotIn((2,120),[(r.M,r.feed) for r,*_ in ranked])

    def test_identifiable_synthetic_surface_and_inside_hull(self):
        rows=parse_cohort(synthetic())
        model,_=fit(rows)
        np.testing.assert_allclose(model.coefficients,[.1,.2,-.3])
        self.assertTrue(model.interpolates(Recipe(2,.04,100,10,"Me/OH")))
        self.assertFalse(model.interpolates(Recipe(3,.02,100,10,"Me/OH")))
        points=[r for r in rows if r.recipe.M!=4 or r.recipe.feed==20]
        triangle,_=fit(points)
        self.assertFalse(triangle.interpolates(Recipe(3,.04,100,10,"Me/OH")))
        self.assertTrue(triangle.interpolates(Recipe(3,.1,100,10,"Me/OH")))

    def test_paper_holdout_labels_never_leak(self):
        def estimates(raw):
            groups=list(folds(parse_cohort(raw),"paper"))
            return [v[2] for group,_,tested,_ in groups if group=="mock-a" for v in tested]
        old=estimates(synthetic())
        self.assertTrue(all(x is not None for x in old))
        self.assertEqual(old,estimates(synthetic(change_paper="mock-A")))

    def test_joint_holdout_protects_series_and_coordinate(self):
        rows=parse_cohort(synthetic())
        for group,train_n,tested,_ in folds(rows,"joint"):
            self.assertEqual(train_n,8)
            for r,_,predicted in tested:
                self.assertAlmostEqual(predicted/r.Mn,1.0)
        # Collapse series IDs across all recipes in one source to model repeated correlated runs.
        from dataclasses import replace
        related=[replace(r,series="common") for r in rows]
        for _,n,tested,_ in folds(related,"joint"):
            self.assertEqual(n,0)
            self.assertTrue(all(x[2] is None for x in tested))

    def test_negative_control_ignores_conversion_and_candidate_order(self):
        raw=synthetic()
        pool="M0_M,feed\n1.5,30\n5,90"
        def report(obs,requested):
            return [(r.M,r.feed,g,p) for r,g,p in design(obs,requested,100,10)[1]]
        self.assertEqual(report(raw,pool),report(raw.replace(',0.57\n',',0.01\n'),"M0_M,feed\n5,90\n1.5,30"))
        self.assertIsNone(report(raw,pool)[0][-1] if report(raw,pool)[0][0]==5 else report(raw,pool)[1][-1])
