import csv
from dataclasses import replace
import io
import unittest

import numpy as np

from polymer_chain_lab.model import make_distribution
from polymer_chain_lab.plugin import calculate
from polymer_chain_lab.viscosity import (CSV_FIELDS, benchmark, literature_series, mixture,
                                        read_observations)


class ViscometryTests(unittest.TestCase):
    def setUp(self):
        self.data = literature_series()
        self.pbs37 = [r for r in self.data if r.temperature == 37]

    def test_reference_transcription_and_missing_cells(self):
        self.assertEqual(len(self.data), 26)
        self.assertEqual([r.mw for r in self.pbs37], [11000, 18800, 43400, 62500, 107600, 216500, 330500])
        self.assertEqual([r.eta for r in self.pbs37], [13.8, 19.0, 31.4, 39.6, 55.7, 84.4, 99])
        self.assertFalse(any(r.sample == "P2" and r.temperature == 55 for r in self.data))
        self.assertEqual(len({r.group for r in self.data}), 7)

    def test_independent_linear_algebra_reference_and_fold_boundaries(self):
        rows, predictions, k, a = benchmark(self.data, 37)
        x = np.log([r.mw for r in rows])
        y = np.log([r.eta for r in rows])
        expected_a, intercept = np.polyfit(x, y, 1)
        self.assertAlmostEqual(a, expected_a)
        self.assertAlmostEqual(k, np.exp(intercept))
        for i, row in enumerate(rows):
            mask = np.arange(len(x)) != i
            slope, intercept = np.polyfit(x[mask], y[mask], 1)
            self.assertAlmostEqual(predictions[row.sample][0], np.exp(intercept + slope * x[i]))
        self.assertEqual(sum(v[2] == "outside" for v in predictions.values()), 2)

    def test_group_holdout_removes_every_member(self):
        r = self.pbs37[2]
        extra = replace(r, sample="P3b", mw=44000, eta=30.5)
        _, p1, _, _ = benchmark([*self.pbs37, extra], 37)
        _, p2, _, _ = benchmark([*self.pbs37, replace(extra, eta=48)], 37)
        self.assertEqual(p1["P3"][0], p2["P3"][0])
        self.assertEqual(p1["P3b"][0], p2["P3b"][0])

    def test_protocol_mixing_and_duplicate_samples_rejected(self):
        r = self.pbs37[0]
        for other in (replace(r, sample="extra", medium="water"),
                      replace(r, sample="extra", sec_temperature=60),
                      replace(r, sample="extra", sec_eluent="THF"),
                      replace(r, sample="extra", group="different"),
                      r):
            with self.subTest(other=other), self.assertRaises(ValueError):
                benchmark([*self.pbs37, other], 37)
        with self.assertRaises(ValueError):
            benchmark(self.data, 28)

    def test_reference_metrics_and_extrapolation(self):
        out = calculate({"model": "intrinsic_viscosity_fit", "target_mw_g_mol": 1000})
        metrics = {m["key"]: m["value"] for m in out["metrics"]}
        self.assertAlmostEqual(metrics["heldout_median_ape"], 2.33253, places=4)
        self.assertAlmostEqual(metrics["heldout_p90_ape"], 7.11989, places=4)
        self.assertIsNone(metrics["target_intrinsic_viscosity"])
        self.assertTrue(out["warnings"])
        self.assertFalse(out["material_updates"])

    def test_custom_csv_matches_fixed_temperature_reference(self):
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(CSV_FIELDS)
        for r in self.pbs37:
            writer.writerow((r.sample, r.group, r.mw, r.eta, r.medium, r.temperature,
                             r.sec_eluent, r.sec_temperature, r.source))
        default = calculate({"model": "intrinsic_viscosity_fit"})
        custom = calculate({"model": "intrinsic_viscosity_fit", "calibration_csv": buf.getvalue()})
        self.assertEqual(default["metrics"], custom["metrics"])
        self.assertEqual(default["tables"], custom["tables"])

    def test_reject_csv_with_missing_source_or_wrong_units_header(self):
        for raw in ("Mw,eta\n18,1", ",".join(CSV_FIELDS) + "\np1,low,12000,0.8,PBS,37,DMAc,50,"):
            with self.assertRaises(ValueError):
                read_observations(raw)

    def test_mwd_uses_mass_fractions_and_component_law(self):
        d = make_distribution([1e4, 4e4], [1, 1], "number")
        eta, mv, components = mixture(d, 0.1, 0.5, 1e4, 4e4)
        self.assertAlmostEqual(eta, 18)
        self.assertAlmostEqual(mv, 180**2)
        np.testing.assert_allclose(components, [10, 20])

    def test_mixture_monodisperse_and_bounds(self):
        d = make_distribution([3e4], [1], "number")
        for a in (0.05, 0.5, 1, 1.5):
            self.assertAlmostEqual(mixture(d, .02, a, 1e4, 4e4)[1], 3e4)
        with self.assertRaises(ValueError):
            mixture(d, .02, .5, 3.1e4, 4e4)

    def test_imported_mwd_and_basis_consistency(self):
        params = {"model":"intrinsic_viscosity_mwd", "mwd_source":"material"}
        mat = {"molar_mass_distribution":{"basis":"number_fraction", "molar_mass_g_mol":[2e4,8e4],
                                          "number_fraction":[.5,.5]}}
        via_record = calculate(params, mat)
        via_mass = calculate({"model":"intrinsic_viscosity_mwd", "distribution_basis":"mass",
                              "distribution_csv":"20000,0.2\n80000,0.8"})
        self.assertEqual(via_mass["metrics"], via_record["metrics"])


if __name__ == "__main__":
    unittest.main()
