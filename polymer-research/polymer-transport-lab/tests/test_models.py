"""Numerical verification, not experimental or material validation."""

import json
import unittest

import numpy as np

from polymer_transport_lab.model import (
    characteristic_times, finite_volume_slab, fit_diffusivity,
    fraction_transferred, fraction_transferred_fourier, parse_fraction_data,
    steady_permeation,
)
from polymer_transport_lab.plugin import APP, calculate


class AnalyticSlabTests(unittest.TestCase):
    def test_initial_and_long_time_limits(self):
        fractions = fraction_transferred_fourier([0, 1e-20, 1e-12, 0.01, 0.1, 1, 10])
        self.assertEqual(fractions[0], 0)
        self.assertEqual(fractions[-1], 1)
        self.assertTrue(np.all(np.diff(fractions) >= 0))
        self.assertTrue(np.all((fractions >= 0) & (fractions <= 1)))

    def test_short_time_half_space_limit(self):
        fo = np.array([1e-20, 1e-12, 1e-8, 1e-4])
        np.testing.assert_allclose(fraction_transferred_fourier(fo), 4*np.sqrt(fo/np.pi), rtol=2e-14)

    def test_image_and_independent_fourier_sum_agree(self):
        fo = np.array([0.001, 0.01, 0.0499999, 0.05, 0.0500001, 0.2])
        # Independent direct sum with ample terms, including across the branch.
        odd = np.arange(1, 1000, 2)
        reference = 1 - 8/np.pi**2 * np.sum(np.exp(-fo[:, None]*np.pi**2*odd**2)/odd**2, axis=1)
        np.testing.assert_allclose(fraction_transferred_fourier(fo), reference, atol=4e-15, rtol=0)

    def test_full_thickness_scaling(self):
        times = np.linspace(0, 5000, 50)
        np.testing.assert_allclose(
            fraction_transferred(times, 1e-12, 1e-4),
            fraction_transferred(4*times, 1e-12, 2e-4), atol=1e-14,
        )
        first = characteristic_times(1e-12, 1e-4)
        second = characteristic_times(1e-12, 2e-4)
        self.assertAlmostEqual(second["t50_s"] / first["t50_s"], 4.0, places=12)

    def test_characteristic_fraction_times(self):
        scales = characteristic_times(1e-12, 1e-4)
        self.assertAlmostEqual(scales["first_mode_time_s"], 1e4/np.pi**2, places=10)
        for key, target in [("t50_s", 0.5), ("t90_s", 0.9), ("t99_s", 0.99)]:
            actual = fraction_transferred([scales[key]], 1e-12, 1e-4)[0]
            self.assertAlmostEqual(actual, target, places=11)

    def test_invalid_physical_inputs(self):
        for d, length in [(0, 1e-4), (-1e-12, 1e-4), (1e-12, 0), (np.inf, 1e-4), (1e-12, np.nan)]:
            with self.subTest(d=d, length=length), self.assertRaises(ValueError):
                fraction_transferred([0, 1], d, length)
        with self.assertRaises(ValueError):
            fraction_transferred([-1], 1e-12, 1e-4)


class FiniteVolumeTests(unittest.TestCase):
    def test_independent_solver_against_analytic(self):
        times = np.r_[0, np.geomspace(10, 40000, 90)]
        numerical = finite_volume_slab(times, 1e-12, 1e-4, 120)
        analytic = fraction_transferred(times, 1e-12, 1e-4)
        self.assertLess(np.max(np.abs(numerical.fraction_transferred - analytic)), 4e-4)
        self.assertTrue(np.all(np.diff(numerical.fraction_transferred) >= -2e-9))

    def test_independently_integrated_flux_conserves_mass(self):
        result = finite_volume_slab(np.r_[0, np.geomspace(0.01, 1e5, 100)], 1e-12, 1e-4, 80)
        balance = result.remaining_profiles.mean(axis=1) + result.integrated_boundary_transfer
        np.testing.assert_allclose(balance, np.ones_like(balance), atol=2e-9, rtol=0)
        self.assertLess(result.mass_balance_max_error, 1e-10)
        self.assertEqual(result.fraction_transferred[0], 0)
        self.assertGreater(result.fraction_transferred[-1], 0.99999)

    def test_grid_convergence(self):
        times = np.r_[0, np.geomspace(10, 10000, 70)]
        exact = fraction_transferred(times, 1e-12, 1e-4)
        errors = []
        for cells in [20, 40, 80]:
            result = finite_volume_slab(times, 1e-12, 1e-4, cells)
            errors.append(float(np.max(np.abs(result.fraction_transferred - exact))))
        self.assertLess(errors[1], errors[0] / 2)
        self.assertLess(errors[2], errors[1] / 2)

    def test_symmetric_bounded_profiles(self):
        result = finite_volume_slab([0, 20, 100, 1000], 1e-12, 1e-4, 48)
        np.testing.assert_allclose(result.remaining_profiles, result.remaining_profiles[:, ::-1], atol=1e-12)
        self.assertTrue(np.all(result.remaining_profiles >= 0))
        self.assertTrue(np.all(result.remaining_profiles <= 1))
        self.assertTrue(np.all(result.remaining_profiles[0] == 1))

    def test_zero_only_and_large_time_are_bounded(self):
        zero = finite_volume_slab([0, 0], 1e-12, 1e-4, 20)
        np.testing.assert_array_equal(zero.fraction_transferred, [0, 0])
        tail = finite_volume_slab([0, 1e12], 1e-4, 1e-9, 20)
        self.assertTrue(tail.time_tail_held_at_equilibrium)
        self.assertGreater(tail.fraction_transferred[-1], 1 - 1e-8)

    def test_grid_and_time_validation(self):
        for cells in [7, 161, 20.5, True]:
            with self.subTest(cells=cells), self.assertRaises(ValueError):
                finite_volume_slab([0, 1], 1e-12, 1e-4, cells)
        with self.assertRaises(ValueError):
            finite_volume_slab([0, 2, 1], 1e-12, 1e-4, 20)


class FitTests(unittest.TestCase):
    def test_synthetic_parameter_recovery(self):
        times = np.array([0, 30, 60, 120, 240, 480, 900, 1800, 3600])
        truth = 2e-12
        fractions = fraction_transferred(times, truth, 1e-4)
        result = fit_diffusivity(times, fractions, 1e-4, 1e-16, 1e-8)
        self.assertAlmostEqual(result.diffusivity_m2_s / truth, 1, places=7)
        self.assertLess(result.rmse, 1e-9)
        self.assertFalse(result.at_bound)
        self.assertIsNone(result.log10_d_standard_error)

    def test_parameter_recovery_from_independent_numeric_data(self):
        times = np.array([0, 60, 120, 240, 480, 900, 1800, 3600])
        numerical = finite_volume_slab(times, 2e-12, 1e-4, 160)
        result = fit_diffusivity(times, numerical.fraction_transferred, 1e-4, 1e-15, 1e-9)
        self.assertLess(abs(result.diffusivity_m2_s / 2e-12 - 1), 0.002)

    def test_diffusivity_thickness_nonidentifiability(self):
        times = np.array([0, 60, 120, 240, 480, 900, 1800])
        data = fraction_transferred(times, 2e-12, 1e-4)
        fit = fit_diffusivity(times, data, 2e-4, 1e-15, 1e-9)
        self.assertAlmostEqual(fit.diffusivity_m2_s / 2e-12, 4, places=6)

    def test_bound_and_plateau_warning(self):
        result = fit_diffusivity([0, 1, 2, 3, 4], [0, 1, 1, 1, 1], 1e-4, 1e-14, 1e-10)
        self.assertTrue(result.at_bound)
        self.assertTrue(result.warnings)
        self.assertIsNone(result.log10_d_standard_error)

    def test_wide_bounds_include_initialization_margin_in_boundary_diagnostic(self):
        # Zero transfer favors the lower bound. Across twenty decades the
        # feasible initialization margin exceeds the old fixed 1e-5 threshold,
        # while the tiny gradient lets least_squares stop at that initial point.
        lower = 1e-24
        result = fit_diffusivity([0, 1, 2, 3, 4], [0, 0, 0, 0, 0], 1e-4, lower, 1e-4)
        self.assertLess(result.diffusivity_m2_s / lower, 1.001)
        self.assertTrue(result.at_bound)
        self.assertTrue(any("bound" in warning for warning in result.warnings))
        self.assertIsNone(result.log10_d_standard_error)

    def test_parser_and_bad_data(self):
        times, fractions = parse_fraction_data("time_s,fraction\n0,0\n4,.4\n1,.1\n2,.2\n3,.3")
        np.testing.assert_array_equal(times, [0, 1, 2, 3, 4])
        np.testing.assert_allclose(fractions, [0, .1, .2, .3, .4])
        invalid = ["", "0,0\n1,.1", "0,0\n1,.1\n2,nan\n3,.3\n4,.4", "0,0\n1,.1\n2,1.2\n3,.3\n4,.4", "0,0\n0,.1\n0,.2\n1,.3\n1,.4"]
        for text in invalid:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_fraction_data(text)


class PluginAndPermeationTests(unittest.TestCase):
    def test_all_presets_strict_json(self):
        json.dumps(APP, allow_nan=False)
        for preset in APP["presets"]:
            with self.subTest(preset=preset["name"]):
                result = calculate(preset["inputs"])
                json.dumps(result, allow_nan=False)
                for plot in result["plots"]:
                    for series in plot["series"]:
                        self.assertEqual(len(series["x"]), len(series["y"]))
                        self.assertLessEqual(len(series["x"]), 1200)
                for parameter in result["material_updates"]["parameters"]:
                    self.assertTrue({"quantity", "unit", "origin", "model_scope"}.issubset(parameter))

    def test_release_uptake_profiles_are_complementary(self):
        release = calculate({"mode": "release", "n_points": 20, "n_cells": 20})
        uptake = calculate({"mode": "uptake", "n_points": 20, "n_cells": 20})
        np.testing.assert_array_equal(release["plots"][0]["series"][0]["y"], uptake["plots"][0]["series"][0]["y"])
        for left, right in zip(release["plots"][-1]["series"], uptake["plots"][-1]["series"]):
            np.testing.assert_allclose(np.array(left["y"])+right["y"], 1, atol=1e-14)

    def test_units_and_signed_permeation(self):
        result = steady_permeation(2e-12, 1e-4, 3e-5, 2e5)
        self.assertAlmostEqual(result["permeability_mol_m_s_Pa"] / 6e-17, 1)
        self.assertAlmostEqual(result["steady_flux_mol_m2_s"] / 1.2e-7, 1)
        negative = steady_permeation(2e-12, 1e-4, 3e-5, -2e5)
        self.assertEqual(negative["steady_flux_mol_m2_s"], -result["steady_flux_mol_m2_s"])

    def test_invalid_backend_fields(self):
        invalid = [{"mode": "sphere"}, {"thickness_um": 0}, {"diffusivity_m2_s": float("nan")}, {"duration_s": -1}, {"n_points": 1001}, {"n_points": 20.5}, {"n_cells": True}, {"temperature_K": 0}, {"fit_lower_d": 1e-8, "fit_upper_d": 1e-10}, {"permeation_mode": "unknown"}, {"fit_data": None}]
        for values in invalid:
            with self.subTest(values=values), self.assertRaises(ValueError):
                calculate(values)

    def test_unknown_input_key_is_rejected_by_direct_api(self):
        with self.assertRaisesRegex(ValueError, "diffusvity_m2_s"):
            calculate({"diffusvity_m2_s": 1e-9})

    def test_metadata_does_not_invent_diffusivity(self):
        material = {"name": "Unverified sample", "chemistry": {"repeat_unit": "C"}, "molar_mass_distribution": {"basis": "number_fraction", "molar_mass_g_mol": [1000, 2000], "number_fraction": [0.4, 0.6]}}
        snapshot = json.dumps(material, sort_keys=True)
        result = calculate({"n_points": 20, "n_cells": 20}, material)
        d = next(row["value"] for row in result["metrics"] if row["key"] == "diffusivity_m2_s")
        self.assertEqual(d, 1e-12)
        self.assertEqual(snapshot, json.dumps(material, sort_keys=True))
        self.assertTrue(any("not used to infer" in note for note in result["notes"]))


if __name__ == "__main__":
    unittest.main()
