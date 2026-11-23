"""Analytical, enumeration, invariant, interface, and invalid-input checks.

These verify the implementation; they do not validate any polymer chemistry.
Run: PYTHONPATH=src python -m unittest discover -s tests -v
"""

import itertools
import json
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad

from polymer_chain_lab.model import (
    R_KJ_MOL_K, debye_form_factor, degrees_of_polymerization, gaussian_dimensions, make_distribution,
    parse_distribution, population_form_factor, ris_statistics, sample_ris,
    wormlike_dimensions,
)
from polymer_chain_lab.plugin import APP, DEFAULT_CSV, TAIL_CSV, calculate


class DistributionTests(unittest.TestCase):
    def test_monodisperse_moments_and_weights(self):
        distribution = make_distribution([50000], [7])
        moments = distribution.moments()
        self.assertEqual(moments, {"Mn": 50000., "Mw": 50000., "Mz": 50000., "dispersity": 1.})
        assert_allclose(distribution.mass_fraction, [1])
        assert_allclose(distribution.scattering_fraction, [1])

    def test_discrete_reference_population(self):
        distribution = make_distribution([1, 3], [1, 1])
        assert_allclose(list(distribution.moments().values()), [2, 2.5, 2.8, 1.25])
        assert_allclose(distribution.mass_fraction, [.25, .75])
        assert_allclose(distribution.scattering_fraction, [.1, .9])

    def test_mass_input_recovers_number_distribution(self):
        distribution = make_distribution([20, 80], [.2, .8], "mass")
        assert_allclose(distribution.number_fraction, [.5, .5])
        self.assertAlmostEqual(distribution.moments()["Mn"], 50)

    def test_duplicate_discrete_masses_combined_and_zero_rows_removed(self):
        distribution = make_distribution([20, 10, 10, 30], [2, 1, 1, 0])
        assert_allclose(distribution.molar_mass_g_mol, [10, 20])
        assert_allclose(distribution.number_fraction, [.5, .5])

    def test_number_log_density_uses_unequal_grid_widths(self):
        distribution = make_distribution([10, 100, 10000], [1, 1, 1], "number_log_density")
        assert_allclose(distribution.number_fraction, [1 / 6, 1 / 2, 1 / 3])

    def test_mass_log_density_converts_after_quadrature(self):
        distribution = make_distribution([10, 100, 10000], [1, 1, 1], "mass_log_density")
        assert_allclose(distribution.mass_fraction, [1 / 6, 1 / 2, 1 / 3])
        weights = np.array([.5 / 10, 1.5 / 100, 1 / 10000])
        assert_allclose(distribution.number_fraction, weights / weights.sum())

    def test_same_mn_and_dispersity_do_not_fix_mz(self):
        first, second = parse_distribution(DEFAULT_CSV), parse_distribution(TAIL_CSV)
        for key in ("Mn", "Mw", "dispersity"):
            self.assertAlmostEqual(first.moments()[key], second.moments()[key])
        self.assertAlmostEqual(first.moments()["Mn"], 50000)
        self.assertAlmostEqual(first.moments()["Mw"], 68000)
        self.assertAlmostEqual(first.moments()["Mz"], 76470.58823529412)
        self.assertAlmostEqual(second.moments()["Mz"], 118823.5294117647)

    def test_csv_header_comments_and_tabs(self):
        distribution = parse_distribution("# example\nmolar_mass_g_mol\tamount\n20\t1\n80\t1\n")
        self.assertAlmostEqual(distribution.moments()["Mn"], 50)

    def test_weight_scale_does_not_change_moments(self):
        first = make_distribution([20, 80], [1e-60, 2e-60])
        second = make_distribution([20, 80], [1e90, 2e90])
        assert_allclose(list(first.moments().values()), list(second.moments().values()))

    def test_moment_ordering(self):
        rng = np.random.default_rng(2026)
        distribution = make_distribution(np.exp(rng.uniform(2, 13, 100)), rng.random(100))
        moments = distribution.moments()
        self.assertLessEqual(moments["Mn"], moments["Mw"])
        self.assertLessEqual(moments["Mw"], moments["Mz"])
        self.assertGreaterEqual(moments["dispersity"], 1)

    def test_invalid_distribution_arrays(self):
        for mass, amount, basis in [([], [], "number"), ([0], [1], "number"), ([1], [-1], "number"),
                                    ([1], [0], "number"), ([1], [np.inf], "number"),
                                    ([np.nan], [1], "number"), ([1, 2], [1], "number"),
                                    ([1], [1], "unknown"), ([1], [1], "number_log_density"),
                                    ([1, 1], [1, 1], "mass_log_density")]:
            with self.subTest(mass=mass, amount=amount, basis=basis), self.assertRaises(ValueError):
                make_distribution(mass, amount, basis)

    def test_invalid_csv(self):
        for value in ("", "# no data", "bad,header", "1,2,3", "1,nan", "nan,1", "molar_mass_g_mol,amount", "1,1\n" * 513):
            with self.subTest(value=value[:30]), self.assertRaises(ValueError):
                parse_distribution(value)


class ChainDimensionTests(unittest.TestCase):
    def test_gaussian_dimensions_and_six_ratio(self):
        re2, rg2 = gaussian_dimensions([6, 60, 600], 2)
        assert_allclose(re2, [12, 120, 1200])
        assert_allclose(rg2, [2, 20, 200])
        assert_allclose(re2 / rg2, [6, 6, 6])

    def test_wlc_rod_limit_without_cancellation(self):
        length = np.array([1e-7, 1e-5, 1e-3])
        re2, rg2 = wormlike_dimensions(length, 1e5)
        assert_allclose(re2 / length**2, 1, rtol=4e-9)
        assert_allclose(rg2 / length**2, 1 / 12, rtol=3e-9)

    def test_wlc_long_chain_gaussian_limit(self):
        length = np.array([1e7, 1e8])
        re2, rg2 = wormlike_dimensions(length, 2)
        assert_allclose(re2, 4 * length, rtol=3e-7)
        assert_allclose(rg2, 2 * length / 3, rtol=7e-7)

    def test_wlc_gyration_agrees_with_pair_distance_integral(self):
        persistence = 2.7
        for ratio in (.02, .1, 1, 7, 100):
            length = ratio * persistence
            _, rg2 = wormlike_dimensions([length], persistence)
            # Rg² = L^-2 ∫_0^L (L-s)<|r(s)-r(0)|²> ds.
            numerical = quad(lambda t: (1 - t) * 2 * persistence**2 *
                             (ratio * t + np.expm1(-ratio * t)), 0, 1,
                             epsabs=1e-13, epsrel=1e-11)[0]
            self.assertAlmostEqual(float(rg2[0]) / numerical, 1, places=10)

    def test_wlc_dimensions_do_not_exceed_contour_limits(self):
        length = np.geomspace(1e-7, 1e7, 200)
        re2, rg2 = wormlike_dimensions(length, 1)
        self.assertTrue(np.all(re2 > 0) and np.all(rg2 > 0))
        self.assertTrue(np.all(re2 <= length**2 * (1 + 1e-13)))
        self.assertTrue(np.all(rg2 <= length**2 / 12 * (1 + 1e-13)))

    def test_wlc_series_boundary_continuity(self):
        re2, rg2 = wormlike_dimensions([.1 * (1 - 1e-8), .1 * (1 + 1e-8)], 1)
        self.assertLess(abs(re2[1] / re2[0] - 1), 5e-8)
        self.assertLess(abs(rg2[1] / rg2[0] - 1), 5e-8)

    def test_invalid_dimension_parameters(self):
        for function in (gaussian_dimensions, wormlike_dimensions):
            for lengths, stiffness in [([0], 1), ([1], 0), ([-1], 1), ([np.inf], 1), ([1], np.nan)]:
                with self.subTest(function=function.__name__, lengths=lengths, stiffness=stiffness), self.assertRaises(ValueError):
                    function(lengths, stiffness)


class ScatteringTests(unittest.TestCase):
    def test_zero_q_normalization(self):
        assert_allclose(debye_form_factor([0], [1, 10, 100]), [[1, 1, 1]], atol=0)

    def test_debye_reference_at_unit_argument(self):
        self.assertAlmostEqual(debye_form_factor([1], [1])[0, 0], 2 / math.e)

    def test_debye_guinier_limit(self):
        q, rg2 = 1e-5, 15
        slope = (1 - debye_form_factor([q], [rg2])[0, 0]) / q**2
        self.assertAlmostEqual(slope / (rg2 / 3), 1, places=6)

    def test_intensity_mixture_uses_squared_mass_weights(self):
        distribution = make_distribution([1, 3], [1, 1])
        components = debye_form_factor([1], [1, 9])[0]
        expected = .1 * components[0] + .9 * components[1]
        actual = population_form_factor([1], [1, 9], distribution)[0]
        self.assertAlmostEqual(actual, expected)
        self.assertNotAlmostEqual(actual, .5 * components.sum())

    def test_population_guinier_radius_is_intensity_weighted(self):
        distribution = make_distribution([1, 3], [1, 1])
        q = 1e-5
        estimate = 3 * (1 - population_form_factor([q], [1, 9], distribution)[0]) / q**2
        self.assertAlmostEqual(estimate, 8.2, places=5)

    def test_form_factor_positive_bounded_decreasing(self):
        values = debye_form_factor(np.geomspace(1e-8, 1e4, 500), [1, 100])
        self.assertTrue(np.all(values > 0))
        self.assertTrue(np.all(values <= 1))
        self.assertTrue(np.all(np.diff(values, axis=0) <= 1e-15))

    def test_invalid_scattering_inputs(self):
        for q, rg2 in [([-1], [1]), ([1], [0]), ([np.nan], [1])]:
            with self.assertRaises(ValueError):
                debye_form_factor(q, rg2)


class RISTests(unittest.TestCase):
    def test_zero_energy_states_have_exact_uniform_probabilities(self):
        result = ris_statistics(12, 300, 0, 0)
        assert_allclose(result["marginals"], 1 / 3)
        self.assertAlmostEqual(result["log_partition"], 12 * math.log(3))
        self.assertAlmostEqual(result["opposite_gauche_fraction"], 2 / 9)

    def test_no_pair_energy_recovers_independent_boltzmann_states(self):
        result = ris_statistics(7, 330, 2, 0)
        expected = np.exp(-np.array([0, 2, 2]) / (R_KJ_MOL_K * 330))
        expected /= expected.sum()
        assert_allclose(result["marginals"], np.tile(expected, (7, 1)))

    def test_transfer_matrix_matches_full_finite_chain_enumeration(self):
        n, temperature, gauche, penalty = 4, 280, 1.3, 3.8
        combinations = np.array(list(itertools.product(range(3), repeat=n)))
        weights = []
        for states in combinations:
            energy = gauche * np.count_nonzero(states)
            energy += penalty * sum({int(a), int(b)} == {1, 2} for a, b in zip(states[:-1], states[1:]))
            weights.append(math.exp(-energy / (R_KJ_MOL_K * temperature)))
        weights = np.asarray(weights)
        partition = weights.sum()
        weights /= partition
        exact_marginals = np.array([[weights[combinations[:, i] == state].sum() for state in range(3)] for i in range(n)])
        result = ris_statistics(n, temperature, gauche, penalty)
        self.assertAlmostEqual(result["log_partition"], math.log(partition))
        assert_allclose(result["marginals"], exact_marginals, atol=1e-14)
        expected_pairs = np.zeros((n - 1, 3, 3))
        for i in range(n - 1):
            for a, b in itertools.product(range(3), repeat=2):
                expected_pairs[i, a, b] = weights[(combinations[:, i] == a) & (combinations[:, i + 1] == b)].sum()
        assert_allclose(result["pair_marginals"], expected_pairs, atol=1e-14)

    def test_reproducible_independent_sampling(self):
        first = sample_ris(n_bonds=20, n_chains=20, seed=201)
        second = sample_ris(n_bonds=20, n_chains=20, seed=201)
        self.assertTrue(np.array_equal(first["states"], second["states"]))
        self.assertTrue(np.array_equal(first["coordinates_nm"], second["coordinates_nm"]))
        third = sample_ris(n_bonds=20, n_chains=20, seed=202)
        self.assertFalse(np.array_equal(first["states"], third["states"]))

    def test_reconstructed_bond_lengths_and_internal_angles(self):
        result = sample_ris(n_bonds=30, n_chains=10)
        vectors = np.diff(result["coordinates_nm"], axis=1)
        assert_allclose(np.linalg.norm(vectors, axis=2), .154, atol=1e-14)
        cosine = -np.sum(vectors[:, :-1] * vectors[:, 1:], axis=2) / .154**2
        assert_allclose(cosine, math.cos(math.radians(112)), atol=1e-13)

    def test_all_trans_chain_has_planar_zigzag_end_to_end_limit(self):
        for bonds in (10, 11):
            result = sample_ris(n_bonds=bonds, n_chains=2, temperature_K=1, gauche_energy_kj_mol=50)
            theta = math.radians(112)
            expected = .154**2 * (bonds**2 * math.sin(theta / 2)**2 + bonds % 2 * math.cos(theta / 2)**2)
            assert_allclose(result["re2_nm2"], expected, rtol=1e-13)
            assert_allclose(result["states"], 0)

    def test_gyration_radius_matches_all_pair_distance_identity(self):
        result = sample_ris(n_bonds=10, n_chains=3)
        for coordinates, actual in zip(result["coordinates_nm"], result["rg2_nm2"]):
            pairs = coordinates[:, None, :] - coordinates[None, :, :]
            expected = np.sum(pairs**2) / (2 * len(coordinates)**2)
            self.assertAlmostEqual(actual, expected)

    def test_sampled_fractions_agree_with_exact_marginals(self):
        result = sample_ris(n_bonds=25, n_chains=1000, seed=819)
        assert_allclose(result["sampled_state_fraction"], result["exact_state_fraction"], atol=.015)
        self.assertLess(abs(result["sampled_opposite_gauche_fraction"] - result["exact_opposite_gauche_fraction"]), .01)

    def test_single_torsion_has_no_pair_fraction(self):
        result = sample_ris(n_bonds=3, n_chains=4)
        self.assertIsNone(result["exact_opposite_gauche_fraction"])
        self.assertIsNone(result["sampled_opposite_gauche_fraction"])

    def test_extreme_allowed_energies_remain_finite(self):
        result = ris_statistics(498, 1, -50, -50)
        self.assertTrue(np.isfinite(result["marginals"]).all())
        self.assertTrue(math.isfinite(result["log_partition"]))
        assert_allclose(result["marginals"].sum(axis=1), 1)

    def test_ris_invalid_and_bounded_inputs(self):
        for inputs in ({"n_bonds": 2}, {"n_bonds": 3.5}, {"n_chains": 1}, {"n_bonds": 500, "n_chains": 2000},
                       {"temperature_K": 0}, {"bond_angle_deg": 180}, {"seed": -1}, {"gauche_energy_kj_mol": np.inf}):
            with self.subTest(inputs=inputs), self.assertRaises(ValueError):
                sample_ris(**inputs)


class ChemicalMassMappingTests(unittest.TestCase):
    def test_terminal_atoms_do_not_count_as_repeats(self):
        # Independent exact mass balance, not the mass of capped monomer × DP.
        assert_allclose(degrees_of_polymerization([282.556, 563.096], 28.054, 2.016), [10, 20])

    def test_zero_terminal_correction_preserves_original_mapping(self):
        assert_allclose(degrees_of_polymerization([1000, 10000], 100), [10, 100])

    def test_fractional_effective_dp_is_retained(self):
        assert_allclose(degrees_of_polymerization([112], 100, 2), [1.1])

    def test_impossible_terminal_corrected_masses_are_rejected(self):
        for masses, repeat, terminal in [([2, 10], 28, 2), ([1], 28, 2), ([10], 0, 2), ([10], 28, -2), ([10], 28, float("nan"))]:
            with self.subTest(masses=masses, repeat=repeat, terminal=terminal), self.assertRaises(ValueError):
                degrees_of_polymerization(masses, repeat, terminal)

    def test_terminal_correction_changes_geometry_but_not_total_mass_moments(self):
        plain = calculate({"distribution_csv": "1000,1\n2000,3", "repeat_mass_g_mol": 100, "contour_increment_nm": 1, "kuhn_length_nm": 1})
        capped = calculate({"distribution_csv": "1000,1\n2000,3", "repeat_mass_g_mol": 100, "end_group_mass_g_mol": 100, "contour_increment_nm": 1, "kuhn_length_nm": 1})
        a = {x["key"]: x["value"] for x in plain["metrics"]}
        b = {x["key"]: x["value"] for x in capped["metrics"]}
        for key in ("Mn", "Mw", "Mz", "dispersity", "tail_number", "tail_mass", "tail_scattering"):
            self.assertEqual(a[key], b[key])
        self.assertAlmostEqual(b["dp_number"], 16.5)
        self.assertAlmostEqual(a["rg_number_nm"]**2 - b["rg_number_nm"]**2, 1 / 6)
        self.assertAlmostEqual(a["rg_mass_nm"]**2 - b["rg_mass_nm"]**2, 1 / 6)
        self.assertAlmostEqual(a["rg_scattering_nm"]**2 - b["rg_scattering_nm"]**2, 1 / 6)

    def test_end_correction_applies_to_wlc_before_finite_chain_geometry(self):
        result = calculate({"model": "wlc", "distribution_csv": "1100,1", "repeat_mass_g_mol": 100, "end_group_mass_g_mol": 100, "contour_increment_nm": .25, "persistence_length_nm": 1})
        metric = {x["key"]: x["value"] for x in result["metrics"]}
        expected_re2, expected_rg2 = wormlike_dimensions([2.5], 1)
        self.assertAlmostEqual(metric["dp_number"], 10)
        self.assertAlmostEqual(metric["rg_number_nm"]**2, expected_rg2[0])
        self.assertAlmostEqual(metric["re_number_nm"]**2, expected_re2[0])

    def test_declared_non_linear_architecture_requires_another_model(self):
        for architecture in ("ring", "cyclic", "branched", "star", "brush", "crosslinked network", "dendrimer"):
            for mode in ("gaussian", "wlc", "ris"):
                with self.subTest(architecture=architecture, mode=mode), self.assertRaisesRegex(ValueError, "open linear"):
                    calculate({"model": mode}, {"chemistry": {"architecture": architecture}})
        calculate({}, {"chemistry": {"architecture": "linear"}})


class PluginTests(unittest.TestCase):
    def test_default_and_presets_are_finite_json_with_bounded_plots(self):
        json.dumps(APP, allow_nan=False)
        for inputs in [{}] + [preset["inputs"] for preset in APP["presets"]]:
            result = calculate(inputs)
            json.dumps(result, allow_nan=False)
            for plot in result["plots"]:
                for series in plot["series"]:
                    self.assertEqual(len(series["x"]), len(series["y"]))
                    self.assertLessEqual(len(series["x"]), 1200)

    def test_wlc_mode_does_not_claim_gaussian_scattering(self):
        result = calculate({"model": "wlc"})
        self.assertFalse(any("form factor" in plot["title"].lower() for plot in result["plots"]))
        self.assertTrue(any("not implemented" in note for note in result["notes"]))

    def test_material_number_distribution_used_explicitly(self):
        material = {"molar_mass_distribution": {"molar_mass_g_mol": [10000, 50000], "number_fraction": [.25, .75], "basis": "number_fraction"}}
        result = calculate({"mwd_source": "material"}, material)
        metrics = {metric["key"]: metric["value"] for metric in result["metrics"]}
        self.assertAlmostEqual(metrics["Mn"], 40000)
        assert_allclose(result["material_updates"]["molar_mass_distribution"]["number_fraction"], [.25, .75])

    def test_missing_or_wrong_basis_material_does_not_fall_back(self):
        for material in ({}, {"molar_mass_distribution": {"basis": "mass_fraction", "molar_mass_g_mol": [100], "number_fraction": [1]}}):
            with self.assertRaises(ValueError):
                calculate({"mwd_source": "material"}, material)

    def test_invalid_plugin_enums_numbers_and_grid(self):
        for inputs in ({"model": "bad"}, {"distribution_basis": "bad"}, {"temperature_K": True},
                       {"n_points": 18.2}, {"q_max_nm_inverse": .0001}, {"kuhn_length_nm": "nan"}):
            with self.subTest(inputs=inputs), self.assertRaises(ValueError):
                calculate(inputs)

    def test_unknown_input_key_is_rejected_by_direct_api(self):
        with self.assertRaisesRegex(ValueError, "kuhn_lenght_nm"):
            calculate({"kuhn_lenght_nm": 100})

    def test_ris_does_not_export_invented_molar_mass(self):
        result = calculate({"model": "ris", "n_bonds": 3, "n_chains": 2,
                            "temperature_K": 1, "gauche_energy_kj_mol": 50})
        self.assertNotIn("molar_mass_distribution", result["material_updates"])
        json.dumps(result, allow_nan=False)
        self.assertEqual(result["plots"][0]["series"][0]["y"], [1.0])
        self.assertGreater(result["plots"][0]["series"][0]["x"][0], 0)
        for parameter in result["material_updates"]["parameters"]:
            self.assertTrue({"quantity", "unit", "origin", "model_scope"} <= parameter.keys())


if __name__ == "__main__":
    unittest.main()
