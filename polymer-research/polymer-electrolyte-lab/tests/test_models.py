import json
import math
import unittest
import numpy as np
from polymer_electrolyte_lab.model import FARADAY, GAS_CONSTANT, nernst_einstein, onsager
from polymer_electrolyte_lab.plugin import APP, calculate


class TransportTheoryTests(unittest.TestCase):
    def test_ne_independent_constant_reference(self):
        # Direct elementary-charge number-density expression, independent units.
        c, dp, dm, T = 500, 1e-12, 2e-12, 298.15
        expected = (1.602176634e-19**2)*(c*6.02214076e23)*(dp+dm)/(1.380649e-23*T)
        actual = nernst_einstein(c,dp,dm,T)
        self.assertAlmostEqual(actual["conductivity_S_m"], expected, places=14)
        self.assertAlmostEqual(actual["cation_current_fraction_ideal"], 1/3)

    def test_collective_recovers_ideal(self):
        ne=nernst_einstein(1000,2e-11,3e-11,320)
        actual=onsager(ne["Lpp_ideal"],ne["Lmm_ideal"],0)
        self.assertAlmostEqual(actual["conductivity_S_m"],ne["conductivity_S_m"])
        self.assertAlmostEqual(actual["t_plus_polymer_fixed"],0.4)

    def test_negative_current_fraction_is_admissible(self):
        result=onsager(2e-13,8e-13,3e-13)
        self.assertGreater(result["minimum_eigenvalue"],0)
        self.assertAlmostEqual(result["t_plus_polymer_fixed"],-0.25)
        self.assertGreater(result["conductivity_S_m"],0)

    def test_positive_dissipation(self):
        L=np.array([[2e-13,3e-13],[3e-13,8e-13]])
        for x in np.random.default_rng(10).normal(size=(100,2)):
            self.assertGreaterEqual(float(x@L@x),0)

    def test_singular_charge_mode(self):
        result=onsager(1e-13,1e-13,1e-13)
        self.assertEqual(result["conductivity_S_m"],0)
        self.assertIsNone(result["t_plus_polymer_fixed"])

    def test_symmetric_fraction(self):
        for pm in (-1e-13,0,2e-13):
            self.assertAlmostEqual(onsager(4e-13,4e-13,pm)["t_plus_polymer_fixed"],0.5)

    def test_invalid_inputs(self):
        for args in ((-1,1e-12,1e-12,300),(1,0,1e-12,300),(1,1e-12,1e-12,0),(1,float("nan"),1e-12,300)):
            with self.assertRaises(ValueError): nernst_einstein(*args)
        for args in ((1e-13,1e-13,2e-13),(-1,1,0),(1,1,float("inf"))):
            with self.assertRaises(ValueError): onsager(*args)

    def test_presets_json_and_handoff_safety(self):
        for preset in APP["presets"]:
            output=calculate(preset["inputs"])
            json.dumps(output,allow_nan=False)
            self.assertEqual(output["handoffs"],[])

    def test_molar_concentration_conversion(self):
        output=calculate({"concentration_mol_L":1})
        direct=nernst_einstein(1000,1e-12,2e-12,298.15)
        self.assertAlmostEqual(output["metrics"][0]["value"],direct["conductivity_S_m"])


class BoundaryRegressionTests(unittest.TestCase):
    def test_rank_one_psd_at_multiple_scales(self):
        # L = s [3,2]^T [3,2] is PSD with charge response F²*s.
        for scale in (1e-300, 1e-13, 1e-3):
            with self.subTest(scale=scale):
                response = onsager(9*scale, 4*scale, 6*scale)
                self.assertAlmostEqual(response["conductivity_S_m"] / (FARADAY**2*scale), 1.0, places=12)
                self.assertAlmostEqual(response["t_plus_polymer_fixed"], 3.0, places=12)
                self.assertAlmostEqual(response["cross_correlation_ratio"], 1.0, places=14)
        for scale in (1e-300, 1e-14, 1e-13, 0.01):
            response = onsager(scale, scale, scale)
            self.assertEqual(response["conductivity_S_m"], 0)
            self.assertIsNone(response["t_plus_polymer_fixed"])

    def test_real_psd_violation_is_rejected_at_small_scale(self):
        # A fixed dimensional tolerance would incorrectly admit these matrices.
        for scale in (1e-300, 1e-30, 1e-13, 1e-3):
            for sign in (-1, 1):
                with self.subTest(scale=scale, sign=sign), self.assertRaises(ValueError):
                    onsager(scale, scale, sign*scale*(1+1e-10))

    def test_zero_diagonal_requires_exactly_zero_cross_coefficient(self):
        for args in ((0, 1e-12, 1e-300), (1e-12, 0, -1e-300),
                     (0, 0, np.nextafter(0.0, 1.0))):
            with self.subTest(args=args), self.assertRaises(ValueError):
                onsager(*args)

    def test_zero_and_single_species_collective_responses(self):
        zero = onsager(0, 0, 0)
        self.assertEqual(zero["conductivity_S_m"], 0)
        self.assertIsNone(zero["t_plus_polymer_fixed"])
        self.assertIsNone(zero["cross_correlation_ratio"])
        anion_only = onsager(0, 2e-13, 0)
        cation_only = onsager(2e-13, 0, 0)
        self.assertEqual(anion_only["t_plus_polymer_fixed"], 0)
        self.assertEqual(cation_only["t_plus_polymer_fixed"], 1)

    def test_current_fraction_is_invariant_under_matrix_rescaling(self):
        # The dimensionless current fraction must not depend on the units-scale
        # of a well-conditioned matrix, including representable tiny entries.
        for scale in (np.nextafter(0.0, 1.0), 1e-300, 1e-13, 1e-3):
            with self.subTest(scale=scale):
                response = onsager(2*scale, 8*scale, 3*scale)
                self.assertAlmostEqual(response["t_plus_polymer_fixed"], -0.25, places=12)
                self.assertGreater(response["conductivity_S_m"], 0)

    def test_current_fraction_matches_independent_flux_accounting(self):
        matrix = np.array([[2e-13, 3e-13], [3e-13, 8e-13]])
        electric_field = 7.0
        charge_numbers = np.array([1.0, -1.0])
        electrochemical_gradient = -FARADAY * charge_numbers * electric_field
        molar_flux = -matrix @ electrochemical_gradient
        ionic_currents = FARADAY * charge_numbers * molar_flux
        response = onsager(matrix[0, 0], matrix[1, 1], matrix[0, 1])
        self.assertAlmostEqual(response["conductivity_S_m"], sum(ionic_currents)/electric_field, places=14)
        self.assertAlmostEqual(response["t_plus_polymer_fixed"], ionic_currents[0]/sum(ionic_currents), places=14)

    def test_backend_rejects_values_below_declared_ranges(self):
        invalid = [{"temperature_K": 1e-320}, {"temperature_K": 0.5},
                   {"concentration_mol_L": 1e-10}, {"d_plus_m2_s": 1e-31},
                   {"d_minus_m2_s": 1e-31}, {"Lpp": -1e-30},
                   {"Lmm": 1.01}, {"Lpm": 1.01}]
        for values in invalid:
            with self.subTest(values=values), self.assertRaises(ValueError):
                calculate(values)
        for args in ((500, 1e-12, 2e-12, 1e-320), (1e-12, 1e-12, 2e-12, 300),
                     (500, 1e-40, 2e-12, 300)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                nernst_einstein(*args)

    def test_bad_payload_and_numeric_overflow_raise_value_error(self):
        for payload in (None, [], "not a dictionary", {"temperature_K": 10**400},
                        {"temperature_K": True}, {"Lpp": np.bool_(True)}, {"model": None}):
            with self.subTest(payload_type=type(payload).__name__), self.assertRaises(ValueError):
                calculate(payload)
        with self.assertRaises(ValueError):
            calculate({}, material=[])

    def test_hidden_numeric_fields_must_also_be_finite(self):
        for key, value in (("Lpp", float("nan")), ("Lmm", float("inf")), ("Lpm", -float("inf"))):
            with self.subTest(key=key), self.assertRaises(ValueError):
                calculate({"model": "nernst_einstein", key: value})

    def test_unknown_keys_do_not_silently_use_defaults(self):
        for values in ({"temperature_k": 320}, {"d_pluss_m2_s": 1e-10}, {7: 1}):
            with self.subTest(values=values), self.assertRaisesRegex(ValueError, "Unknown input fields"):
                calculate(values)

    def test_supported_extreme_inputs_produce_finite_json(self):
        for values in (
            {"temperature_K": 1, "concentration_mol_L": 100, "d_plus_m2_s": 1e-4, "d_minus_m2_s": 1e-4},
            {"temperature_K": 5000, "concentration_mol_L": 1e-9, "d_plus_m2_s": 1e-30, "d_minus_m2_s": 1e-30},
            {"model": "onsager", "Lpp": 1, "Lmm": 1, "Lpm": -1},
        ):
            with self.subTest(values=values):
                output = calculate(values)
                json.dumps(output, allow_nan=False)
                if values.get("model", "nernst_einstein") == "nernst_einstein":
                    concentrations = output["plots"][0]["series"][0]["x"]
                    self.assertGreaterEqual(min(concentrations), 1e-9)
                    self.assertLessEqual(max(concentrations), 100)

    def test_undefined_cross_ratio_is_not_plotted(self):
        for pp, mm in ((0, 0), (0, 1e-12), (1e-12, 0)):
            with self.subTest(pp=pp, mm=mm):
                output = calculate({"model": "onsager", "Lpp": pp, "Lmm": mm, "Lpm": 0})
                self.assertEqual(output["plots"], [])
                self.assertTrue(any("cross-correlation ratio is undefined" in note for note in output["notes"]))
                json.dumps(output, allow_nan=False)

    def test_plot_and_table_payload_contract(self):
        for preset in APP["presets"]:
            output = calculate(preset["inputs"])
            for plot in output["plots"]:
                self.assertIn(plot["x_scale"], ("linear", "log"))
                self.assertIn(plot["y_scale"], ("linear", "log"))
                for series in plot["series"]:
                    self.assertEqual(len(series["x"]), len(series["y"]))
                    self.assertLessEqual(len(series["x"]), 1200)
                    self.assertTrue(all(math.isfinite(value) for value in series["x"] + series["y"]))
            for table in output["tables"]:
                self.assertTrue(all(len(row) == len(table["columns"]) for row in table["rows"]))

if __name__ == "__main__": unittest.main()
