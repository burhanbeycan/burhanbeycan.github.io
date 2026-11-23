import json
import math
import unittest

from polymer_synthesis_lab.model import step_growth, living_chain_growth, terminal_copolymer
from polymer_synthesis_lab.plugin import APP, calculate


class StepGrowthTests(unittest.TestCase):
    def test_balanced_carothers(self):
        v = step_growth(1, 1, .99, 146.14, 116.21, 18.015)
        self.assertAlmostEqual(v["DPn_original_monomer_units"], 100)
        self.assertAlmostEqual(v["remaining_A_mol"]+v["remaining_B_mol"], 2*v["remaining_molecules_mol"])
        self.assertAlmostEqual(v["initial_mass_g"], v["byproduct_mass_g"]+v["retained_population_mass_g"])
        self.assertAlmostEqual(v["Mn_population_g_mol"]*v["remaining_molecules_mol"], v["retained_population_mass_g"])

    def test_published_99_100_full_conversion_case(self):
        # Meaurio 2024 Table 1, system 2: 99 mol AA,100 mol BB -> Xn199.
        v = step_growth(99, 100, 1, 100, 100)
        self.assertEqual(v["DPn_original_monomer_units"], 199)
        self.assertEqual(v["remaining_A_mol"], 0)
        self.assertEqual(v["remaining_B_mol"], 2)
        self.assertEqual(v["balanced_Flory_length_dispersity"], None)

    def test_labels_can_be_interchanged(self):
        v = step_growth(2, 1, .5, 110, 90)
        w = step_growth(1, 2, .5, 90, 110)
        self.assertEqual(v["DPn_original_monomer_units"], w["DPn_original_monomer_units"])
        self.assertEqual(v["Mn_population_g_mol"], w["Mn_population_g_mol"])
        self.assertEqual(v["conversion_A"], .25)
        self.assertEqual(v["conversion_B"], .5)

    def test_zero_conversion_and_singular_limit(self):
        v = step_growth(1, 1, 0, 100, 120)
        self.assertEqual(v["DPn_original_monomer_units"], 1)
        self.assertEqual(v["Mn_population_g_mol"], 110)
        v = step_growth(1, 1, 1, 100, 120)
        self.assertIsNone(v["Mn_population_g_mol"])
        self.assertIsNone(v["DPn_original_monomer_units"])
        json.dumps(v, allow_nan=False)

    def test_impossible_loss_and_invalid_inputs(self):
        for p in [-.1, 1.1, float("nan"), True]:
            with self.assertRaises(ValueError):
                step_growth(1, 1, p, 100, 120)
        with self.assertRaises(ValueError):
            step_growth(1, 1, .9, 100, 120, 110)


class ChainGrowthTests(unittest.TestCase):
    def test_chain_count_and_mass_conservation(self):
        v = living_chain_growth(.01, .0001, .9, .5, 1, 100, 30)
        self.assertAlmostEqual(v["DPn"], 180)
        self.assertAlmostEqual(v["Mn_g_mol"], 18030)
        self.assertAlmostEqual(v["chain_population_mass_g"], v["Mn_g_mol"]*v["initiated_chains_mol"])
        self.assertAlmostEqual(v["reacted_monomer_mol"]+v["unreacted_monomer_mol"], .01)

    def test_zero_addition_and_integer_chain_count(self):
        v = living_chain_growth(1, .1, 0, 1, 1, 100, 30)
        self.assertEqual(v["Mn_g_mol"], 30)
        self.assertEqual(v["end_mass_fraction"], 1)
        with self.assertRaises(ValueError):
            living_chain_growth(1, .1, .5, 1, 1.5, 100, 30)
        with self.assertRaises(ValueError):
            living_chain_growth(1, .1, .5, 0, 1, 100, 30)


class CopolymerTests(unittest.TestCase):
    def test_bernoulli_limit(self):
        v = terminal_copolymer(.3, 1, 1)
        self.assertAlmostEqual(v["F_A"], .3)
        self.assertAlmostEqual(v["dyad_AA"], .09)
        self.assertAlmostEqual(v["dyad_BB"], .49)

    def test_alternating_limit(self):
        v = terminal_copolymer(.25, 0, 0)
        self.assertAlmostEqual(v["F_A"], .5)
        self.assertEqual(v["dyad_AA"], 0)
        self.assertEqual(v["dyad_BB"], 0)
        self.assertEqual(v["mean_A_run"], 1)

    def test_stationarity_normalization_and_label_swap(self):
        v = terminal_copolymer(.21, 4, .1)
        w = terminal_copolymer(.79, .1, 4)
        self.assertAlmostEqual(v["F_A"], 1-w["F_A"])
        self.assertAlmostEqual(v["dyad_AB"], v["dyad_BA"])
        self.assertAlmostEqual(sum(v["dyad_"+s] for s in ["AA", "AB", "BA", "BB"]), 1)
        self.assertAlmostEqual(v["P_AA"]+v["P_AB"], 1)
        self.assertAlmostEqual(v["P_BA"]+v["P_BB"], 1)
        self.assertAlmostEqual(v["F_A"]*v["P_AA"]+v["F_B"]*v["P_BA"], v["F_A"])

    def test_absent_monomers_and_negative_ratios_rejected(self):
        for args in [(0, 1, 1), (1, 1, 1), (.3, -1, 1)]:
            with self.assertRaises(ValueError):
                terminal_copolymer(*args)


class AdapterTests(unittest.TestCase):
    def test_presets_and_default_export_finite_json(self):
        for inputs in [{}]+[x["inputs"] for x in APP["presets"]]:
            v = calculate(inputs)
            json.dumps(v, allow_nan=False)
            self.assertTrue(v["metrics"])
            self.assertTrue(v["tables"])

    def test_unknown_inputs_and_invalid_hidden_fields(self):
        with self.assertRaises(ValueError):
            calculate({"invented_rate": 3})
        with self.assertRaises(ValueError):
            calculate({"model": "step_growth", "initiation_efficiency": 0})

    def test_singular_export_explanation(self):
        v = calculate({"p_groups": 1})
        self.assertTrue(any("infinite-chain" in s for s in v["warnings"]))


if __name__ == "__main__":
    unittest.main()
