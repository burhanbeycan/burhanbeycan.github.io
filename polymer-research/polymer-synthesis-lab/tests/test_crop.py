"""Numerical, provenance and leakage regressions for EtOx CROP."""
import csv
import io
import math
import unittest

from polymer_synthesis_lab.crop import Recipe, audit, load_cohort, parse_cohort
from polymer_synthesis_lab.plugin import calculate


def modified(raw, update):
    rows = list(csv.DictReader(io.StringIO(raw)))
    output = io.StringIO()
    writer = csv.DictWriter(output, rows[0].keys(), lineterminator="\n")
    writer.writeheader()
    for i, row in enumerate(rows):
        update(i, row)
        writer.writerow(row)
    return output.getvalue()


class CROPTests(unittest.TestCase):
    def test_reference_kinetics_and_end_mass(self):
        r = Recipe(3, .03, 100, 30, "Me/OH")
        x, m, full = r.predict()
        self.assertAlmostEqual(r.kp, .01006287124946)  # Table 2 Arrhenius regression at 373.15 K
        self.assertAlmostEqual(x, 1 - math.exp(-r.kp * .03 * 1800))
        self.assertAlmostEqual(m, x * 99.133 * 100 + 32.042)
        self.assertAlmostEqual(full, 9945.342)
        azide = Recipe(3, .03, 100, 30, "Me/N3").predict()
        self.assertAlmostEqual(azide[0], x)
        self.assertAlmostEqual(azide[1] - m, 57.056 - 32.042)

    def test_time_and_active_site_limits(self):
        r = Recipe(4, .04, 100, 0, "Me/OH")
        self.assertEqual(r.predict()[0], 0)
        a, b = (Recipe(4*d, .04*d, 100, 10, "Me/OH") for d in (.5, 1))
        self.assertEqual(a.feed, b.feed)
        self.assertGreater(b.predict()[0], a.predict()[0])
        self.assertAlmostEqual(Recipe(1, 1, 200, 10000, "Me/OH").predict()[0], 1)

    def test_guard_unknown_chemistry_and_invalid_inputs(self):
        for args in [(4, 0, 100, 5, "Me/OH"),(4, .04, float("nan"), 1, "Me/OH"),(4,.04,100,-1,"Me/OH"),(4,.04,100,1,"MeOTs")]:
            with self.assertRaises(ValueError):
                Recipe(*args)

    def test_bundled_figure_and_table(self):
        batch = audit(load_cohort("batch_azide"), "paper")
        self.assertEqual([row.Mn for row in batch[0]], [3600,5900,9600,11100])
        self.assertAlmostEqual(batch[-1]["kinetic"][1], 10.30259, places=4)
        self.assertEqual(batch[-1]["scaled"], (0, None, None))
        flow = audit(load_cohort("flow_feed"), "feed")
        self.assertEqual(len(flow[1]), 4)
        self.assertTrue(all(f[1:3] == [3, 1] for f in flow[1]))
        self.assertAlmostEqual(flow[-1]["scaled"][1], 28.54213, places=4)
        self.assertAlmostEqual(flow[-1]["scaled"][2], 46.65953, places=4)

    def test_heldout_Mn_and_conversion_never_enter_its_own_predictions(self):
        raw = load_cohort("flow_feed")
        changed = modified(raw, lambda i, r: r.update(Mn_g_mol="100000", X_observed="0.01") if i==0 else None)
        before, after = (audit(s, "feed")[2][0] for s in [raw, changed])
        for key in ("kinetic","scaled","X","full"):
            self.assertEqual(before[key], after[key])
        self.assertNotEqual(before["ape_scaled"], after["ape_scaled"])

    def test_new_feed_alone_is_not_a_new_concentration(self):
        raw = load_cohort("flow_feed")
        for method in ("concentration", "concentration_and_feed", "paper"):
            result = audit(raw, method)
            self.assertTrue(all(f[1] == 0 for f in result[1]))
            self.assertEqual(result[-1]["scaled"], (0,None,None))
        for method in ("series","feed"):
            self.assertTrue(all(f[1] == 3 for f in audit(raw, method)[1]))

    def test_joint_split_drops_either_matching_coordinate(self):
        raw = modified(load_cohort("flow_feed"), lambda i,r: r.update(
            M0_M=[2,2,4,4][i], I0_M=[.1,.2,.2,.4][i],paper=["lab1","lab1","lab2","lab2"][i]))
        folds = audit(raw,"concentration_and_feed")[1]
        self.assertTrue(all(f[1:3]==[1,1] for f in folds))
        self.assertTrue(all(f[1:3]==[2,2] for f in audit(raw,"paper")[1]))

    def test_reject_mixed_assay_and_duplicate_measurement_or_coefficient_source(self):
        raw = load_cohort("batch_azide")
        for change in [lambda i,r: r.update(assay="SEC-MALS") if i==1 else None,
                       lambda i,r: r.update(sample="same"),
                       lambda i,r: r.update(paper="https://doi.org/10.1021/acsmacrolett.5b00392"),
                       lambda i,r: r.update(M0_M="inf")]:
            with self.assertRaises(ValueError):
                parse_cohort(modified(raw, change))

    def test_missing_conversions_do_not_block_validation(self):
        raw = modified(load_cohort("flow_feed"), lambda i,r: r.update(X_observed=""))
        rows, _, results, _ = audit(raw, "feed")
        self.assertTrue(all(r.conversion is None for r in rows))
        self.assertTrue(all(r["scaled"] is not None for r in results))
        output = calculate({"model":"petox_benchmark","crop_csv":raw,"crop_split":"feed"})
        self.assertTrue(any("Cross-concentration" in msg for msg in output["warnings"]))

    def test_invalid_group_and_csv_report_useful_errors(self):
        with self.assertRaisesRegex(ValueError, "grouping"):
            audit(load_cohort("flow_feed"), "random")
        with self.assertRaisesRegex(ValueError, "columns"):
            parse_cohort("source,time,Mn\na,1,10")
