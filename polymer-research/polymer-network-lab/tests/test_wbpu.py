"""WBPU functional groups, extraction mass bases and group-preserving comparisons."""

import json
import unittest

from polymer_network_lab.plugin import calculate
from polymer_network_lab.wbpu import (
    DEFAULT_RECIPE, DMA_COLUMNS, EXTRACTION_COLUMNS, PAPER_DMA, PAPER_SERIES,
    RECIPE_COLUMNS, SERIES_COLUMNS, dma_sweeps, extraction, recipe_balance, swelling_folds,
)


class WbpuMeasurements(unittest.TestCase):
    def test_mass_and_functionality_for_peg_pcl_bdo(self):
        r = recipe_balance(DEFAULT_RECIPE)
        self.assertAlmostEqual(r["entries"][0][-1], 0.1)
        self.assertAlmostEqual(r["entries"][1][-1], 2)
        self.assertEqual(r["stages"]["extension"]["OH"], 1)
        self.assertAlmostEqual(r["index"], 101.58192074)

    def test_measured_assays_and_extender_amine_sites(self):
        recipe = RECIPE_COLUMNS + """
polyol,pre,OH,1,OH_number,112.2,
blend,pre,NCO,1,NCO_wt_pct,42.02,
diamine,extend,amine,0.03,molar_mass,60,2
water,extend,reacted_water,0.01801528,net_urea,,"""
        r = recipe_balance(recipe)
        self.assertAlmostEqual(r["totals"]["OH"], 2)
        self.assertAlmostEqual(r["totals"]["NCO"], 10)
        self.assertAlmostEqual(r["totals"]["amine"], 1)
        self.assertAlmostEqual(r["totals"]["reacted_water"], 2)
        self.assertAlmostEqual(r["index"], 200)
        self.assertAlmostEqual(r["balance"], 5)

    def test_equivalent_weight_takes_no_second_functionality_factor(self):
        row = RECIPE_COLUMNS + "\npolyol,one,OH,1,equivalent_weight,500,\niso,one,NCO,1,equivalent_weight,250,"
        self.assertAlmostEqual(recipe_balance(row)["index"], 200)
        for bad in [row.replace("500,", "500,2"), row.replace("polyol,one,OH", "polyol,one,water"),
                    row.replace(",1,equivalent_weight,", ",nan,equivalent_weight,"), row.replace(",one,OH", ",,OH")]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                recipe_balance(bad)

    def test_solvent_extraction_three_mass_bases_and_zero_residue(self):
        data = EXTRACTION_COLUMNS + "\na,THF,23,24,1,0.9,0.5\nb,THF,23,24,1,0,0"
        a, b = extraction(data)
        self.assertAlmostEqual(a[4], 50)
        self.assertAlmostEqual(a[5], 50)
        self.assertAlmostEqual(a[6], -10)
        self.assertAlmostEqual(a[7], 80)
        self.assertIsNone(b[-1])
        with self.assertRaises(ValueError):
            extraction(data.replace(",0.9,0.5", ",0.9,0.95"))

    def test_dma_keeps_frequencies_separate(self):
        series = dma_sweeps(DMA_COLUMNS + "\na,1,40,2\na,1,25,4\na,10,25,5")
        self.assertEqual(series[0]["x"], [25, 40])
        self.assertEqual(series[0]["y"], [4, 2])
        self.assertEqual(len(series), 2)
        self.assertEqual(dma_sweeps(""), [])
        with self.assertRaises(ValueError):
            dma_sweeps(DMA_COLUMNS + "\na,1,20,2\na,1,20,4")

    def test_published_whole_formulation_folds(self):
        data = swelling_folds(PAPER_SERIES, PAPER_DMA)
        self.assertEqual(len(data["folds"]), 3)
        self.assertEqual(len(data["models"]), 3)
        for result, expected in zip(data["mae"], [55.933333, 33.307937, 17.361253]):
            self.assertAlmostEqual(result, expected, places=5)
        self.assertEqual([f["outside_training_range"] for f in data["folds"]], [True, False, True])
        # The middle formulation's held-out recipe prediction lies on the straight
        # line through the two endpoint observations.
        expected = 231.9 + (128.2 - 231.9) * 3 / 7
        self.assertAlmostEqual(expected, data["folds"][1]["predictions"][1])

    def test_held_out_replicates_never_change_own_predictions(self):
        rows = PAPER_SERIES.splitlines()
        duplicate = rows[1].replace("231.9", "900")
        original = swelling_folds(PAPER_SERIES, PAPER_DMA)["folds"][0]
        new = swelling_folds(PAPER_SERIES + "\n" + duplicate, PAPER_DMA)["folds"][0]
        self.assertEqual(new["predictions"], original["predictions"])
        self.assertAlmostEqual(new["observed"], (231.9 + 900) / 2)
        self.assertEqual(new["replicates_or_rows"], 2)

    def test_fixed_conditions_and_measurement_independence(self):
        malformed = [
            PAPER_SERIES.replace(",23,168", ",22,168", 1),
            PAPER_SERIES.replace("92.8", ""),
            PAPER_SERIES.replace("247.7", "nan"),
            PAPER_SERIES.replace(",41,167", ",38,167"),
            PAPER_SERIES.replace("jatropha WBPU cast film", "mat", 1),
        ]
        for sample in malformed:
            with self.subTest(sample=sample), self.assertRaises(ValueError):
                swelling_folds(sample, PAPER_DMA)
        with self.assertRaises(ValueError):
            swelling_folds(PAPER_SERIES)

    def test_four_groups_enable_joint_model_and_recipe_only_needs_no_dma(self):
        rows = SERIES_COLUMNS + "\na,film,water,23,2,x,10,2,3\nb,film,water,23,2,x,20,4,8\nc,film,water,23,2,x,30,9,5\nd,film,water,23,2,x,50,1,2"
        full = swelling_folds(rows, "synthetic pre-soak 1 Hz 20 °C")
        self.assertEqual(len(full["models"]), 4)
        stripped = "\n".join([SERIES_COLUMNS] + [r.rsplit(",", 1)[0] + "," for r in rows.splitlines()[1:]])
        self.assertEqual(len(swelling_folds(stripped)["models"]), 2)

    def test_wbpu_output_does_not_transfer_unmeasured_network_density(self):
        for mode in ("recipe", "evidence", "holdout"):
            res = calculate({"network_mode": mode})
            self.assertEqual(res["material_updates"], {})
            self.assertEqual(res["handoffs"], [])
            json.dumps(res, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
