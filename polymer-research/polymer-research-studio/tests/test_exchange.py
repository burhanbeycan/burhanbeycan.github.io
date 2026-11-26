from copy import deepcopy
import unittest
from polymer_workbench.exchange import default_material, validate_material, make_result, validate_exchange, handoff, digest

APP={"id":"test","name":"Test Model","version":"0.1.0","references":[]}

class ExchangeTests(unittest.TestCase):
    def result(self):
        return make_result(APP,{"D":1e-12},{"metrics":[],"plots":[],"tables":[],"notes":[],"warnings":[],"material_updates":{},"handoffs":[{"target":"transport","inputs":{"thickness_um":120},"notes":["Equilibrium geometry only"]}]},default_material())

    def test_record_round_trip(self):
        doc=self.result()
        self.assertEqual(validate_exchange(doc),doc)
        self.assertEqual(doc["material"]["history"][-1]["run_id"],doc["run_id"])
        self.assertEqual(doc["checksum_sha256"],digest(doc))

    def test_tampered_result_rejected(self):
        doc=self.result();doc["inputs"]["D"]=1e-3
        with self.assertRaisesRegex(ValueError,"checksum"):validate_exchange(doc)

    def test_schema_rejected(self):
        with self.assertRaisesRegex(ValueError,"schema"):validate_exchange({"schema":"polymer-material/2.0"})

    def test_nonfinite_rejected(self):
        m=default_material();m["conditions"]["temperature_K"]=float('nan')
        with self.assertRaises(ValueError):validate_material(m)

    def test_mwd_semantics(self):
        m=default_material();m["molar_mass_distribution"]={"molar_mass_g_mol":[100,300],"number_fraction":[.5,.5],"basis":"number_fraction"}
        self.assertEqual(validate_material(m),m)
        m["molar_mass_distribution"]["basis"]="mass_fraction"
        with self.assertRaises(ValueError):validate_material(m)

    def test_no_implicit_coefficient_mapping(self):
        link=handoff(self.result(),"electrolyte",[])
        self.assertEqual(link["inputs"],{})
        self.assertEqual(link["material"]["sample_id"],"demonstration-01")

    def test_explicit_compatible_mapping(self):
        link=handoff(self.result(),"transport",[{"key":"thickness_um"}])
        self.assertEqual(link["inputs"],{"thickness_um":120})
        with self.assertRaises(ValueError):handoff(self.result(),"transport",[{"key":"unrelated"}])

    def test_identity_cannot_be_overwritten(self):
        with self.assertRaises(ValueError):make_result(APP,{}, {"material_updates":{"sample_id":"changed"}},default_material())

    def test_session_round_trip_and_missing_material(self):
        doc={"schema":"polymer-session/1.0","material":default_material(),"runs":[self.result()]}
        self.assertEqual(validate_exchange(doc),doc)
        del doc["material"]
        with self.assertRaises(ValueError):validate_exchange(doc)

    def test_bad_distribution_weights(self):
        m=default_material();m["molar_mass_distribution"]={"molar_mass_g_mol":[100,300],"number_fraction":[.6,.5]}
        with self.assertRaises(ValueError):validate_material(m)

    def test_bad_history_cannot_break_ui(self):
        m=default_material();m["history"]=["not an object"]
        with self.assertRaises(ValueError):validate_material(m)

    def test_structurally_bad_signed_result_rejected(self):
        # A checksum is not authentication: a self-consistent malformed file
        # must still be rejected by the schema validator.
        for change in ("application", "plot", "table", "handoff"):
            with self.subTest(change=change):
                doc=self.result()
                if change=="application":doc["application"]["name"]={"bad":"type"}
                if change=="plot":doc["results"]["plots"]=[{"title":"curve","x_label":"t","y_label":"F","series":[{"name":"f","x":[1,2],"y":[1]}]}]
                if change=="table":doc["results"]["tables"]=[{"title":"data","columns":["A"],"rows":[[1,2]]}]
                if change=="handoff":doc["results"]["handoffs"][0]["notes"]="not a list"
                doc["checksum_sha256"]=digest(doc)
                with self.assertRaises(ValueError):validate_exchange(doc)

if __name__=='__main__':unittest.main()
