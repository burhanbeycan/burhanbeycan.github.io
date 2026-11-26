"""Numerical, analytical and contract checks; not experimental validation."""

from copy import deepcopy
from dataclasses import replace
import json
import math
import unittest

import numpy as np
from scipy.special import erfc

from polymer_network_lab.model import (
    GAS_CONSTANT_J_MOL_K,
    NetworkParameters,
    chemical_potential_residual,
    equilibrium_scan,
    free_energy_components,
    free_energy_j_m3_dry,
    residual_derivative_phi,
    slab_release_fraction,
    solve_equilibrium,
    time_to_release_fraction,
)
from polymer_network_lab.plugin import APP, calculate


class ThermodynamicsTests(unittest.TestCase):
    def test_published_model_equilibrium_not_experimental_data(self):
        p = NetworkParameters(chi=0.2, nu_e_mol_m3=0.001 / 18.0e-6)
        root = solve_equilibrium(p).selected
        # Hong et al. (2008), free swelling before the layer is constrained:
        # their Nv=0.001, chi=0.2 give the reported stretch 3.215.
        self.assertAlmostEqual(root.isotropic_stretch, 3.215, delta=0.0005)
        self.assertEqual(root.stability, "stable")
        self.assertAlmostEqual(root.volume_swelling_ratio, root.isotropic_stretch**3, places=11)
        self.assertAlmostEqual(root.polymer_volume_fraction * root.volume_swelling_ratio, 1.0, places=14)

    def test_residual_and_energy_derivative_are_consistent(self):
        p = NetworkParameters(chi=0.27, nu_e_mol_m3=234.0)
        scale = GAS_CONSTANT_J_MOL_K * p.temperature_K / p.solvent_molar_volume_m3_mol
        for q in [1.2, 2.0, 7.5, 40.0]:
            h = q * 1.0e-5
            derivative = (free_energy_j_m3_dry(q + h, p) - free_energy_j_m3_dry(q - h, p)) / (2.0 * h * scale)
            phi = 1.0 / q
            independent_g = math.log1p(-phi) + phi + p.chi * phi**2 + p.nu_e_mol_m3 * p.solvent_molar_volume_m3_mol * (phi ** (1.0 / 3.0) - phi)
            self.assertAlmostEqual(derivative, independent_g, delta=max(1.0e-10, abs(independent_g) * 2.0e-8))

    def test_curvature_matches_second_derivative_and_minimum(self):
        p = NetworkParameters()
        root = solve_equilibrium(p).selected
        q = root.volume_swelling_ratio
        h = q * 2.0e-4
        numerical_curvature = (free_energy_j_m3_dry(q + h, p) - 2.0 * free_energy_j_m3_dry(q, p)
                               + free_energy_j_m3_dry(q - h, p)) / (h**2 * p.energy_scale_j_m3_dry)
        self.assertAlmostEqual(numerical_curvature / root.curvature_dimensionless, 1.0, delta=2.0e-5)
        self.assertGreater(root.curvature_dimensionless, 0.0)
        grid = np.geomspace(1.000001, q * 100.0, 1000)
        self.assertLessEqual(root.free_energy_j_m3_dry, float(np.min(free_energy_j_m3_dry(grid, p))) + 1.0e-7)
        self.assertLess(root.free_energy_j_m3_dry, 0.0)

    def test_analytic_residual_derivative(self):
        p = NetworkParameters(chi=0.4)
        for phi in [0.001, 0.1, 0.7]:
            h = phi * 1.0e-5
            numeric = (chemical_potential_residual(phi + h, p) - chemical_potential_residual(phi - h, p)) / (2.0 * h)
            self.assertAlmostEqual(numeric, residual_derivative_phi(phi, p), delta=max(1.0e-9, abs(numeric) * 1.0e-7))

    def test_dry_energy_reference_and_nonnegative_elastic_energy(self):
        p = NetworkParameters()
        self.assertEqual(free_energy_components(1.0, p), {"mixing": 0.0, "elastic": 0.0, "total": 0.0})
        q = np.array([1.0, 1.0 + 1.0e-12, 1.00001, 1.1, 3.0, 100.0])
        self.assertTrue(np.all(free_energy_components(q, p)["elastic"] >= 0.0))
        tiny_strain = 1.0e-6
        energy = free_energy_components(math.exp(tiny_strain), p)["elastic"]
        expected_leading = p.dry_shear_modulus_pa * tiny_strain**2 / 3.0
        self.assertAlmostEqual(energy / expected_leading, 1.0, delta=1.0e-6)

    def test_more_crosslinking_or_worse_solvent_reduces_equilibrium_volume(self):
        p = NetworkParameters(chi=0.2)
        crosslink_scan = equilibrium_scan(p, "nu_e_mol_m3", [1.0, 10.0, 100.0, 1000.0])
        chi_scan = equilibrium_scan(p, "chi", [-0.2, 0.2, 0.5, 0.8, 1.2])
        self.assertTrue(np.all(np.diff([point.volume_swelling_ratio for point in crosslink_scan]) < 0.0))
        self.assertTrue(np.all(np.diff([point.volume_swelling_ratio for point in chi_scan]) < 0.0))

    def test_temperature_rescales_energy_and_modulus_at_fixed_chi(self):
        p = NetworkParameters(temperature_K=250.0)
        low = solve_equilibrium(p).selected
        high = solve_equilibrium(replace(p, temperature_K=500.0)).selected
        self.assertEqual(low.polymer_volume_fraction, high.polymer_volume_fraction)
        self.assertAlmostEqual(high.free_energy_j_m3_dry / low.free_energy_j_m3_dry, 2.0, places=13)
        self.assertAlmostEqual(high.swollen_shear_modulus_pa / low.swollen_shear_modulus_pa, 2.0, places=13)

    def test_dimensionless_group_collapse(self):
        p1 = NetworkParameters(nu_e_mol_m3=100.0, solvent_molar_volume_m3_mol=20.0e-6)
        p2 = NetworkParameters(nu_e_mol_m3=50.0, solvent_molar_volume_m3_mol=40.0e-6)
        self.assertAlmostEqual(solve_equilibrium(p1).selected.volume_swelling_ratio,
                               solve_equilibrium(p2).selected.volume_swelling_ratio, places=13)

    def test_dilute_good_solvent_asymptote(self):
        p = NetworkParameters(chi=0.2, nu_e_mol_m3=1.0e-6, solvent_molar_volume_m3_mol=1.0e-4)
        phi = solve_equilibrium(p).selected.polymer_volume_fraction
        expected = (p.elastic_to_mixing_ratio / (0.5 - p.chi)) ** (3.0 / 5.0)
        self.assertAlmostEqual(phi / expected, 1.0, delta=0.0005)

    def test_parameter_boundaries_are_finite_and_root_selection_is_stable(self):
        for chi in [-5.0, 0.5, 5.0]:
            for nu_e in [1.0e-6, 1.0e6]:
                for volume in [1.0e-8, 1.0e-3]:
                    with self.subTest(chi=chi, nu_e=nu_e, volume=volume):
                        p = NetworkParameters(chi=chi, nu_e_mol_m3=nu_e, solvent_molar_volume_m3_mol=volume)
                        equilibrium = solve_equilibrium(p)
                        root = equilibrium.selected
                        self.assertTrue(0.0 < root.polymer_volume_fraction < 1.0)
                        self.assertTrue(math.isfinite(root.free_energy_j_m3_dry))
                        self.assertGreater(root.volume_swelling_ratio, 1.0)
                        self.assertLess(root.relative_residual, 2.0e-8)
                        self.assertEqual(root.stability, "stable")
                        stable_energies = [r.free_energy_j_m3_dry for r in equilibrium.roots if r.stability == "stable"]
                        self.assertEqual(root.free_energy_j_m3_dry, min(stable_energies))

    def test_invalid_model_inputs(self):
        for field, value in [("chi", float("nan")), ("temperature_K", 0), ("temperature_K", float("inf")),
                             ("nu_e_mol_m3", 0), ("nu_e_mol_m3", -1), ("chi", "0.2"), ("chi", True),
                             ("solvent_molar_volume_m3_mol", 18.0)]:
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                NetworkParameters(**{field: value})
        for phi in [0.0, 1.0, -0.2, float("nan"), [True], [], [0.1] * 1201]:
            with self.subTest(phi=str(phi)[:30]), self.assertRaises(ValueError):
                chemical_potential_residual(phi, NetworkParameters())
        with self.assertRaises(ValueError):
            free_energy_j_m3_dry(0.99, NetworkParameters())
        for variable, values in [("temperature", [1, 2]), ("chi", [0.2, 0.1]), ("chi", [0.2]),
                                 ("chi", np.linspace(0.0, 0.8, 161))]:
            with self.assertRaises(ValueError):
                equilibrium_scan(NetworkParameters(), variable, values)


class ReleaseTests(unittest.TestCase):
    def test_zero_monotonic_bounded_and_long_time_limits(self):
        times = np.concatenate(([0.0], np.geomspace(1.0e-12, 1.0e7, 500)))
        released = slab_release_fraction(times, 1.0e-4, 1.0e-12)
        self.assertEqual(released[0], 0.0)
        self.assertEqual(released[-1], 1.0)
        self.assertTrue(np.all(np.diff(released) >= -1.0e-14))
        self.assertTrue(np.all((released >= 0.0) & (released <= 1.0)))

    def test_short_time_two_face_solution_and_full_thickness(self):
        thickness = 1.0e-4
        diffusion = 1.0e-12
        time = 1.0e-6
        expected = 4.0 * math.sqrt(diffusion * time / math.pi) / thickness
        self.assertAlmostEqual(slab_release_fraction(time, thickness, diffusion), expected, places=15)

    def test_fourier_solution_against_independent_image_representation(self):
        # Integral of the Poisson-transformed surface flux; this is a separate
        # representation from the odd-mode series used at these times.
        for fo in [0.00500001, 0.01, 0.03, 0.08]:
            image_sum = 0.0
            erfc_sum = 0.0
            for k in range(1, 30):
                image_sum += (-1.0)**k * math.exp(-k**2 / (4.0 * fo))
                erfc_sum += (-1.0)**k * k * erfc(k / (2.0 * math.sqrt(fo)))
            expected = 4.0 * math.sqrt(fo / math.pi) * (1.0 + 2.0 * image_sum) - 4.0 * erfc_sum
            observed = slab_release_fraction(fo, 1.0e-2, 1.0e-4)
            self.assertAlmostEqual(observed, expected, places=13)

    def test_diffusion_and_geometry_similarity_and_release_times(self):
        l = 2.0e-4
        d = 1.0e-12
        t = np.linspace(0.0, 3600.0, 100)
        np.testing.assert_allclose(slab_release_fraction(t, l, d), slab_release_fraction(4.0 * t, 2.0 * l, d), atol=1.0e-14)
        np.testing.assert_allclose(slab_release_fraction(t, l, d), slab_release_fraction(t / 3.0, l, 3.0 * d), atol=1.0e-14)
        for fraction in [0.1, 0.5, 0.9, 0.99]:
            time = time_to_release_fraction(fraction, l, d)
            self.assertAlmostEqual(slab_release_fraction(time, l, d), fraction, places=12)
            self.assertAlmostEqual(time_to_release_fraction(fraction, 2.0 * l, d) / time, 4.0, places=12)

    def test_invalid_release_inputs(self):
        for time, thickness, diffusion in [(-1.0, 1.0e-4, 1.0e-12), (1.0, 0.0, 1.0e-12),
                                           (1.0, 1.0e-4, 0.0), (float("nan"), 1.0e-4, 1.0e-12),
                                           (1.0, 1.0e-4, float("inf"))]:
            with self.assertRaises(ValueError):
                slab_release_fraction(time, thickness, diffusion)
        with self.assertRaises(ValueError):
            time_to_release_fraction(1.0, 1.0e-4, 1.0e-12)


class PluginContractTests(unittest.TestCase):
    def test_default_and_all_presets_serialize_without_nonfinite_numbers(self):
        json.dumps(APP, allow_nan=False)
        for inputs in [{}] + [preset["inputs"] for preset in APP["presets"]]:
            result = calculate(inputs)
            json.dumps(result, allow_nan=False)
            for key in ["metrics", "plots", "tables", "notes", "warnings", "material_updates", "handoffs"]:
                self.assertIn(key, result)
            for plot in result["plots"]:
                for series in plot["series"]:
                    self.assertEqual(len(series["x"]), len(series["y"]))
                    self.assertLessEqual(len(series["x"]), 1200)

    def test_full_thickness_handoff_and_supplied_diffusivity(self):
        result = calculate({"scan_variable": "none", "dry_thickness_um": 80.0, "diffusivity_m2_s": 4.0e-13})
        metrics = {metric["key"]: metric["value"] for metric in result["metrics"]}
        handoff = result["handoffs"][0]
        self.assertEqual(handoff["target"], "transport")
        self.assertEqual(handoff["inputs"]["mode"], "release")
        self.assertEqual(handoff["inputs"]["diffusivity_m2_s"], 4.0e-13)
        self.assertAlmostEqual(handoff["inputs"]["thickness_um"], 80.0 * metrics["isotropic_stretch"], places=12)
        self.assertIn("FULL", " ".join(handoff["notes"]))
        self.assertIn("dynamic swelling", " ".join(handoff["notes"]))

    def test_material_metadata_is_not_silently_used_or_mutated(self):
        material = {"chemistry": {"repeat_unit": "unspecified", "tacticity": "unspecified"},
                    "molar_mass_distribution": {"molar_mass_g_mol": [10000], "number_fraction": [1.0]},
                    "conditions": {"temperature_K": 500}, "parameters": [{"quantity": "nu_e_mol_m3", "value": 999}]}
        before = deepcopy(material)
        result = calculate({"scan_variable": "none", "release_preview": "disabled"}, material)
        self.assertEqual(before, material)
        self.assertTrue(any("molar-mass" in note for note in result["notes"]))
        self.assertTrue(any("not automatically" in note for note in result["notes"]))
        for parameter in result["material_updates"]["parameters"]:
            for key in ["quantity", "unit", "origin", "model_scope", "source_input_origin"]:
                self.assertIn(key, parameter)
        self.assertEqual(result["handoffs"], [])

    def test_scan_off_release_off_and_out_of_range_handoff(self):
        result = calculate({"scan_variable": "none", "release_preview": "disabled"})
        self.assertEqual(len(result["plots"]), 1)
        self.assertEqual(result["handoffs"], [])
        large = calculate({"scan_variable": "none", "dry_thickness_um": 1.0e6})
        self.assertEqual(large["handoffs"], [])
        self.assertTrue(any("thickness range" in warning for warning in large["warnings"]))

    def test_backend_rejects_invalid_inputs_without_gui(self):
        invalid = [{"n_points": 20.5}, {"n_points": 1000000}, {"scan_points": 161}, {"chi": float("nan")},
                   {"scan_variable": "ionic"}, {"release_preview": True}, {"nu_e_mol_m3": -1},
                   {"diffusivity_m2_s": "1e-12"}, {"temperature_K": 0}, {"unknown_key": 1},
                   {"chi_scan_min": 1, "chi_scan_max": 0}, {"parameter_provenance": "x" * 8001}]
        for inputs in invalid:
            with self.subTest(inputs=str(inputs)[:100]), self.assertRaises(ValueError):
                calculate(inputs)
        with self.assertRaises(ValueError):
            calculate([])
        with self.assertRaises(ValueError):
            calculate({}, material=[])


if __name__ == "__main__":
    unittest.main()
