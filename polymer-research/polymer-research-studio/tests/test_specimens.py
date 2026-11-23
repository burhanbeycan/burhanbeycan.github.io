import unittest
from copy import deepcopy
from fastapi.testclient import TestClient
from polymer_workbench.exchange import default_material, handoff, make_result, validate_exchange
from polymer_workbench.lineage import derive_material
from polymer_workbench.server import create_app


class PhysicalSpecimens(unittest.TestCase):
    def setUp(self):
        self.batch = default_material()
        self.batch.update(sample_id="B1", lineage={"form": "batch", "process": "synthesis", "ancestors": []})
        self.batch["chemistry"]["repeat_unit"] = "PVA"
        self.batch["parameters"] = [{"quantity": "test observation", "unit": "1"}]
        self.film = derive_material(self.batch, "F1", "Film", "cast_film", "solution casting")
        self.mat = derive_material(self.batch, "M1", "Mat", "electrospun_mat", "electrospinning")

    def test_siblings_share_only_chemical_identity(self):
        self.assertEqual(self.film["lineage"]["ancestors"], ["B1"])
        self.assertEqual(self.mat["chemistry"]["repeat_unit"], "PVA")
        for item in (self.film, self.mat):
            self.assertEqual(item["parameters"], [])
            self.assertEqual(item["evidence"], [])
            self.assertIsNone(item["conditions"]["temperature_K"])
            self.assertEqual(item["history"], [])
        self.film["chemistry"]["repeat_unit"] = "changed"
        self.assertEqual(self.batch["chemistry"]["repeat_unit"], "PVA")

    def test_transfer_between_two_specimens(self):
        app = {"id": "network", "name": "Network", "version": "0.3.0"}
        doc = make_result(app, {}, {"handoffs": [{"target": "transport", "inputs": {"thickness_um": 200}}]}, self.film)
        fields = [{"key": "thickness_um"}]
        self.assertEqual(handoff(doc, "transport", fields, self.mat)["inputs"], {})
        self.assertEqual(handoff(doc, "transport", fields, self.mat)["material"], self.mat)
        self.assertEqual(handoff(doc, "transport", fields, self.film)["inputs"], {"thickness_um": 200})

    def test_source_article_attribution_roundtrip(self):
        obs = dict(origin="secondary_citation", quantity="diffusivity", value="illustrative", unit="m²/s",
                   conditions="unit test", reported_in="article A", attributed_to="article B", locator="table 3")
        self.mat["evidence"] = [obs]
        self.assertEqual(validate_exchange(self.mat)["evidence"], [obs])
        del obs["attributed_to"]
        with self.assertRaises(ValueError):
            validate_exchange(self.mat)

    def test_cycle_and_duplicate_registry_fail(self):
        session = dict(schema="polymer-session/1.0", material=self.mat, runs=[], samples=[self.batch, self.mat, self.film])
        validate_exchange(session)
        broken=deepcopy(session); broken["samples"].append(self.film)
        with self.assertRaises(ValueError):
            validate_exchange(broken)
        self.film["lineage"]["ancestors"]=["M1"]
        with self.assertRaises(ValueError):
            validate_exchange(session)
        with self.assertRaises(ValueError):
            derive_material(self.mat, "B1", "impossible", "cast_film", "casting")

    def test_older_exchanges_and_api(self):
        client=TestClient(create_app())
        legacy=default_material()
        self.assertEqual(validate_exchange(legacy), legacy)
        result=client.post("/api/material/derive", json=dict(parent=self.batch, sample_id="F2", name="Film",
                            form="cast_film", process="solution casting"))
        self.assertEqual(result.status_code,200,result.text)
        self.assertEqual(result.json()["lineage"]["ancestors"],["B1"])
