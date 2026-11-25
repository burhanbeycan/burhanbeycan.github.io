"""Regression checks for reusable chain records and the public API."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from threading import Thread
import unittest
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from http.server import ThreadingHTTPServer
from app import Handler, Studies, archived
from study import PILOT, check_export, inputs, read, request_key

REPEAT = "[*:1]CCN(C(=O)CC)[*:2]"


class TestIdentity(unittest.TestCase):
    def setUp(self):
        self.source = read(PILOT/"source_model.json")
        self.payload = {"repeat_smiles":REPEAT,"dp":3,"temperature_K":298.15}

    def test_transfer_ignores_acquisition_fields_but_checks_geometry(self):
        record = deepcopy(self.source)
        record["energy"].update(cache_source="local",cache_hit=True)
        self.assertEqual(check_export(record),self.source)
        record["energy"]["lowest_geometry"]["atoms"][0]["x"] += 0.001
        with self.assertRaisesRegex(ValueError,"özeti"):
            check_export(record)

    def test_structure_and_conditions_are_part_of_key(self):
        s,_,t = inputs(self.payload)
        self.assertEqual(s["canonical_smiles"],self.source["structure"]["canonical_smiles"])
        key = request_key(s,self.source,t)
        self.assertNotEqual(key,request_key(s,self.source,310.0))
        different,_,_ = inputs({**self.payload,"dp":2})
        self.assertNotEqual(key,request_key(different,None,t))

    def test_unsupported_inputs(self):
        for patch in ({"solvent":"ethanol"},{"temperature_K":float("nan")},{"repeat_smiles":"[*:1]CC[*:2]"},
                      {"dp":0},{"temperature_K":True}):
            with self.subTest(patch=patch), self.assertRaises(ValueError):
                inputs({**self.payload,**patch})


class TestDashboard(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.service = Studies(cls.tmp.name)
        cls.srv = ThreadingHTTPServer(("127.0.0.1",0),Handler)
        cls.srv.studies = cls.service
        cls.thread = Thread(target=cls.srv.serve_forever,daemon=True)
        cls.thread.start()
        cls.base = "http://127.0.0.1:"+str(cls.srv.server_port)

    @classmethod
    def tearDownClass(cls):
        cls.srv.shutdown()
        cls.srv.server_close()
        cls.thread.join()
        cls.service.executor.shutdown(wait=True)
        cls.tmp.cleanup()

    def post(self,obj):
        r=Request(self.base+"/api/studies",data=json.dumps(obj).encode(),headers={"Content-Type":"application/json"})
        with urlopen(r,timeout=10) as resp:
            return json.load(resp)

    def test_manual_import_and_pilot_resolve_same_archived_model(self):
        with urlopen(self.base+"/api/pilot") as resp:
            study=json.load(resp)
        src=read(PILOT/"source_model.json")
        for form in ({"repeat_smiles":REPEAT,"dp":3},{"material":src}):
            self.assertEqual(self.post(form)["study_id"],study["study_id"])
        self.assertEqual(study["analysis"]["meta"]["n_frames"],360)
        self.assertEqual(study["system"]["water_count"],500)
        self.assertEqual(study["system"]["source_artifact_id"],src["artifact_id"])
        self.assertEqual(archived(PILOT)["study_id"],study["study_id"])
        url="/api/data/"+study["study_id"]+"/analysis/rg_trace.csv"
        with urlopen(self.base+url) as resp:
            self.assertEqual(resp.read(),(PILOT/"analysis/rg_trace.csv").read_bytes())

    def test_rejects_wrong_inputs_and_unlisted_download(self):
        with self.assertRaises(HTTPError) as bad:
            self.post({"dp":2,"repeat_smiles":"not-a-structure"})
        self.assertEqual(bad.exception.code,400)
        with self.assertRaises(HTTPError) as disallowed:
            urlopen(self.base+"/api/data/"+read(PILOT/"study.json")["study_id"]+"/../study.py")
        self.assertEqual(disallowed.exception.code,404)


if __name__ == "__main__":
    unittest.main()
