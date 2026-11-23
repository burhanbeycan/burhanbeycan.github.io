from dataclasses import replace
import unittest

import numpy as np

from polymer_transport_lab import pva
from polymer_transport_lab.plugin import calculate
from polymer_transport_lab.pva_workflow import EQ_DEMO, KIN_DEMO


class FixedFilmWaterTests(unittest.TestCase):
    def setUp(self):
        self.iso = pva.read_isotherm(EQ_DEMO)
        self.curves = pva.read_curves(KIN_DEMO, self.iso)
        self.curve = self.curves[0]

    def test_one_face_fourfold_time_rule(self):
        for beta in (0, 2):
            one = replace(self.curve, faces=1)
            two = replace(self.curve, faces=2, times=self.curve.times / 4)
            np.testing.assert_allclose(pva.loading_curve(self.iso, one, 1e-12, beta),
                                       pva.loading_curve(self.iso, two, 1e-12, beta), atol=1e-12)

    def test_constant_limit_refines_toward_exact_solution(self):
        ref = pva.loading_curve(self.iso, self.curve, 1e-12, 0)
        errors = [np.max(np.abs(ref - pva.loading_curve(self.iso, self.curve, 1e-12, 1e-8, count))) for count in (24,48,96)]
        self.assertGreater(errors[0], 3 * errors[1])
        self.assertGreater(errors[1], 3 * errors[2])

    def test_desorption_branch(self):
        iso = pva.read_isotherm(EQ_DEMO.replace("sorption", "desorption"))
        c = replace(self.curve,start=.85,end=.4)
        y = pva.loading_curve(iso, c, 1e-12, 2)
        self.assertAlmostEqual(y[0], .16)
        self.assertTrue((np.diff(y) < 0).all())
        self.assertTrue((y > .04).all())
        with self.assertRaisesRegex(ValueError,"branch"):
            pva.loading_curve(iso,self.curve,1e-12,2)

    def test_holdout_all_curves_no_target_leakage(self):
        m, train, test = pva.compare_models(self.iso,self.curves,90)
        self.assertEqual(len(test),2)
        self.assertTrue(all(c.thickness == 90 for c in test))
        self.assertTrue(all(c.thickness == 50 for c in train))
        changed = [replace(c,observed=c.observed+.05) if c.thickness == 90 else c for c in self.curves]
        m2,_,_ = pva.compare_models(self.iso,changed,90)
        for a,b in zip(m,m2):
            self.assertEqual((a["dref"],a["contrast"]),(b["dref"],b["contrast"]))
            self.assertGreater(b["test_rmse"],5*a["test_rmse"])
        self.assertAlmostEqual(m[1]["dref"]/1.5e-12,1,places=2)
        self.assertAlmostEqual(m[1]["contrast"],2,places=2)

    def test_requires_condition_matched_independent_equilibrium(self):
        for body in [KIN_DEMO.replace("synthetic-a-0.4","example-equilibrium"),
                     KIN_DEMO.replace("298.15","303.15",1), KIN_DEMO.replace("synthetic-PVA-A","unrelated",1),
                     KIN_DEMO.replace(",1,0,0.4",",3,0,0.4",1)]:
            with self.assertRaises(ValueError):
                pva.read_curves(body,self.iso)
        with self.assertRaises(ValueError):
            self.iso.step(0,.9)
        with self.assertRaises(ValueError):
            pva.read_isotherm(EQ_DEMO.replace("0.7,0.10","0.7,0.01"))

    def test_holdout_coverage_and_unresolved_design(self):
        with self.assertRaisesRegex(ValueError,"holdout"):
            pva.compare_models(self.iso,self.curves,80)
        with self.assertRaisesRegex(ValueError,"two independent humidity steps"):
            pva.compare_models(self.iso,[c for c in self.curves if c.start==0],90)

    def test_publication_units_isotherm_and_distinct_observables(self):
        out = calculate({"mode":"pva_published"})
        table = out["tables"][0]["rows"]
        np.testing.assert_allclose([row[3] for row in table],[9.4e-12,1.4e-11,2.4e-11],rtol=1e-13,atol=0)
        self.assertAlmostEqual(table[2][1], .1173465659, places=9)
        self.assertAlmostEqual(out["metrics"][0]["value"],2.4/.94)

    def test_user_record_origin_and_no_extrapolation(self):
        with self.assertRaisesRegex(ValueError,"synthetic"):
            calculate({"mode":"pva_holdout","pva_origin":"user_supplied_unverified"})
        with self.assertRaises(ValueError):
            calculate({"mode":"pva_forward","pva_end":1})

    def test_early_time_similarity_and_experimental_humidity_interpolation(self):
        self.assertAlmostEqual(pva.apparent_d_ratio(1e-8), 1, places=6)
        for b in (0.3,1.5,3.5):
            self.assertAlmostEqual(pva.apparent_d_ratio(b)/pva.apparent_d_ratio(b,10,1e-6),1,places=5)
        result = pva.paper_humidity_check()
        low, mid, high = pva.literature_rows()
        for row in (low,high):
            self.assertAlmostEqual(result["dry_d"]*pva.apparent_d_ratio(result["contrast"]*row[1]/high[1])/row[3],1,places=6)
        self.assertAlmostEqual(result["predicted"] / 1e-11,1.509464,places=4)
        self.assertEqual(result["observed"],mid[3])
