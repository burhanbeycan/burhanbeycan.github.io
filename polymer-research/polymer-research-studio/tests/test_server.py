"""HTTP contract checks that also run without any laboratory installed."""
from types import SimpleNamespace
from unittest import TestCase, main
from unittest.mock import patch

from fastapi.testclient import TestClient
from polymer_workbench.exchange import validate_exchange
from polymer_workbench.server import create_app


class ServerTests(TestCase):
    def setUp(self):
        plugin=SimpleNamespace(
            APP={"id":"demo","name":"Test equation","version":"0.1.0",
                 "fields":[{"key":"x","default":2}],"references":[]},
            calculate=lambda inputs,material:{"metrics":[{"label":"x squared","value":inputs["x"]**2}],"material_updates":{}})
        with patch("polymer_workbench.server.discover", return_value={"demo":plugin}):
            self.client=TestClient(create_app())

    def test_calculation_is_an_exchange_record(self):
        response=self.client.post("/api/calculate/demo",json={"inputs":{"x":3}})
        self.assertEqual(response.status_code,200)
        doc=validate_exchange(response.json())
        self.assertEqual(doc["results"]["metrics"][0]["value"],9)
        self.assertEqual(doc["inputs"],{"x":3})

    def test_unknown_inputs_and_unknown_application(self):
        self.assertEqual(self.client.post("/api/calculate/demo",json={"inputs":{"typo":3}}).status_code,422)
        self.assertEqual(self.client.post("/api/calculate/missing",json={}).status_code,404)

    def test_nonfinite_request_is_readable_error(self):
        response=self.client.post("/api/calculate/demo",content='{"inputs":{"x":NaN}}',headers={"Content-Type":"application/json"})
        self.assertEqual(response.status_code,422)
        self.assertIn("NaN",response.json()["detail"])

    def test_local_origin_and_host_boundary(self):
        self.assertEqual(self.client.get("/api/catalog",headers={"Origin":"https://unrelated.example"}).status_code,403)
        self.assertEqual(self.client.get("/api/catalog",headers={"Origin":"http://testserver"}).status_code,200)
        self.assertEqual(self.client.get("/api/catalog",headers={"Host":"unrelated.example"}).status_code,400)

    def test_oversized_request_and_static_assets(self):
        self.assertEqual(self.client.post("/api/exchange/validate",content=' ' * 4_000_001).status_code,413)
        for path in ("/","/static/app.js","/static/styles.css"):
            response=self.client.get(path)
            self.assertEqual(response.status_code,200)
            self.assertIn("default-src 'self'",response.headers["Content-Security-Policy"])


if __name__=="__main__":main()
