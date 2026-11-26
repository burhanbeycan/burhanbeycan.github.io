import copy
import json
import unittest

from polymer_electrolyte_lab.peo import formulation
from polymer_electrolyte_lab.plugin import calculate
from polymer_electrolyte_lab.records import PESKO, ZHENG, SYNTHETIC, gui_inputs


def run(record):
    return calculate(gui_inputs(record))


def metrics(result):
    return {metric["key"]: metric["value"] for metric in result["metrics"]}


class MatchedPEOTests(unittest.TestCase):
    def test_primary_data_against_printed_composition_and_conductivity(self):
        p = metrics(run(PESKO))
        self.assertEqual(round(p["salt_molarity"], 2), 1.59)
        self.assertEqual(p["sigma_eis"], 2.0)
        self.assertAlmostEqual(p["sigma_ne"], 1.84, delta=0.01)
        self.assertAlmostEqual(p["sigma_ratio"], 1.09, delta=0.01)
        self.assertGreater(p["sigma_ratio"], 1)  # the program must not clip the ratio
        allowances = run(PESKO)["tables"][-1]["rows"][1]
        self.assertLess(allowances[1], 1.0)
        self.assertGreater(allowances[2], 1.1)

    def test_second_specimen_matches_reported_fraction_without_density(self):
        p = metrics(run(ZHENG))
        self.assertEqual(round(p["t_tracer"], 2), 0.19)
        self.assertEqual(p["sigma_eis"], 1.5)
        self.assertNotIn("sigma_ne", p)
        self.assertNotIn("sigma_ratio", p)
        self.assertNotIn("salt_molarity", p)
        self.assertIn("1.59", run(PESKO)["notes"][-1])

    def test_weighed_film_geometry_unit_conversion(self):
        p = metrics(run(SYNTHETIC))
        self.assertAlmostEqual(p["salt_molarity"], (0.3 / 287.09) / 0.0011)
        self.assertAlmostEqual(p["sigma_eis"], 0.1)
        self.assertAlmostEqual(p["eo_per_li"], (1/44.05)/(0.3/287.09))

    def test_salt_mass_and_target_ratio_modes_agree(self):
        original = formulation(PESKO["formulation"], 363.15)
        weighed = dict(PESKO["formulation"])
        weighed.pop("eo_per_li")
        weighed["salt_g"] = original["salt_g"]
        measured = formulation(weighed, 363.15)
        for key in ("eo_per_li", "molality", "salt_g", "salt_wt_pct", "concentration"):
            self.assertAlmostEqual(original[key], measured[key])

    def test_unknown_or_mismatched_volume_is_never_silently_used(self):
        for evidence in [
            {"density_g_mL": 1, "temperature_K": 298.15, "source": "Lab"},
            {"final_volume_mL": 2, "density_g_mL": 1, "temperature_K": 363.15, "source": "Lab"},
            {"pure_PEO_density": 1, "temperature_K": 363.15, "source": "Lab"},
        ]:
            r = copy.deepcopy(PESKO)
            r["formulation"]["volume_at_state"] = {"state_id": r["state"]["state_id"], **evidence}
            with self.subTest(evidence=evidence), self.assertRaises(ValueError):
                run(r)
        r = copy.deepcopy(ZHENG)
        r["formulation"]["volume_at_state"] = {"state_id": r["state"]["state_id"], "salt_mol_L": 0.8, "temperature_K": 363.15, "source": "Test only"}
        self.assertAlmostEqual(metrics(run(r))["salt_molarity"], 0.8)

    def test_mismatched_tracer_eis_specimens_or_temperatures_are_rejected(self):
        for part in ("tracer", "eis"):
            for key, value in (("temperature_K", 343.15), ("state_id", "other")):
                r = copy.deepcopy(PESKO)
                r[part][key] = value
                with self.subTest(part=part, key=key), self.assertRaises(ValueError):
                    run(r)

    def test_cross_paper_formulation_or_density_at_same_temperature_is_rejected(self):
        r = copy.deepcopy(ZHENG)
        r["formulation"] = copy.deepcopy(PESKO["formulation"])
        with self.assertRaisesRegex(ValueError, "Formulation state_id"):
            run(r)
        r = copy.deepcopy(ZHENG)
        r["formulation"]["volume_at_state"] = copy.deepcopy(PESKO["formulation"]["volume_at_state"])
        with self.assertRaisesRegex(ValueError, "Volume or concentration state_id"):
            run(r)

    def test_ambiguous_eis_and_ambiguous_composition_are_rejected(self):
        for entry, extras in (("formulation", {"salt_g": 0.4}), ("eis", {"area_cm2": 0.1})):
            r = copy.deepcopy(PESKO)
            r[entry].update(extras)
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                run(r)

    def test_oligomer_end_group_mass_fraction(self):
        r = dict(PESKO["formulation"], eo_mass_fraction_in_peo=0.9)
        self.assertAlmostEqual(formulation(r, 363.15)["salt_g"], formulation(PESKO["formulation"], 363.15)["salt_g"] * 0.9)
        for fraction in (0, 1.1, True):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                formulation(dict(r, eo_mass_fraction_in_peo=fraction), 363.15)

    def test_missing_crystallinity_is_distinct_from_zero(self):
        r = copy.deepcopy(SYNTHETIC)
        r["state"]["crystallinity"]["fraction"] = 0
        known = str(run(r)["tables"][1])
        self.assertIn("0% of PEO mass", known)
        r["state"]["crystallinity"] = None
        self.assertIn("Unreported", str(run(r)["tables"][1]))
        r["state"]["crystallinity"] = {"fraction": 0.2, "basis": "unknown", "method": "DSC"}
        with self.assertRaises(ValueError):
            run(r)

    def test_edited_citation_is_not_represented_as_verified_preset(self):
        r = copy.deepcopy(PESKO)
        r["eis"]["sigma_mS_cm"] = 2.2
        self.assertIn("User-edited", run(r)["notes"][-1])

    def test_reading_bound_and_json_payload_validation(self):
        r = copy.deepcopy(PESKO)
        r["tracer"]["relative_plot_reading_bound"] = 1
        with self.assertRaises(ValueError):
            run(r)
        args = gui_inputs(PESKO)
        args["peo_state"] = '{"state":'
        with self.assertRaises(ValueError):
            calculate(args)
        json.dumps(run(PESKO), allow_nan=False)
