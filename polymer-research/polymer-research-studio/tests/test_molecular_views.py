"""Validate graph payloads before an imported result reaches the molecular renderer."""
from copy import deepcopy
import unittest
from polymer_workbench.exchange import validate_output, make_result, validate_exchange, default_material, digest

GRAPH = {"title":"Ethane", "caption":"Two-dimensional molecular connectivity", "smiles":"CC",
         "atoms":[{"index":0,"element":"C","x":0,"y":0,"label":"C","unit_index":0},
                  {"index":1,"element":"C","x":1.5,"y":0,"label":"C","unit_index":0}],
         "bonds":[{"begin":0,"end":1,"order":1,"direction":"NONE","stereo":"STEREONONE"}]}

class MolecularViews(unittest.TestCase):
    def test_valid_graph_round_trip(self):
        doc=make_result({"id":"chemistry","name":"Chemistry","version":"0.2.0"}, {},
                        {"molecular_views":[deepcopy(GRAPH)]},default_material())
        self.assertEqual(validate_exchange(doc), doc)

    def test_invalid_graphs_rejected(self):
        for mutation in ("duplicate_atom", "missing_endpoint", "nonfinite", "duplicate_bond", "bad_order", "bad_membership", "boolean_endpoint", "text_label"):
            with self.subTest(mutation=mutation):
                view=deepcopy(GRAPH)
                if mutation=="duplicate_atom":view["atoms"][1]["index"]=0
                if mutation=="missing_endpoint":view["bonds"][0]["end"]=2
                if mutation=="nonfinite":view["atoms"][0]["x"]=float("inf")
                if mutation=="duplicate_bond":view["bonds"].append(deepcopy(view["bonds"][0]))
                if mutation=="bad_order":view["bonds"][0]["order"]=4
                if mutation=="bad_membership":view["atoms"][0]["unit_index"]="red"
                if mutation=="boolean_endpoint":view["bonds"][0]["begin"]=False
                if mutation=="text_label":view["atoms"][0]["label"]={"bad":"label"}
                with self.assertRaises(ValueError):validate_output({"molecular_views":[view]})

    def test_recomputed_checksum_does_not_bypass_graph_validation(self):
        doc=make_result({"id":"chemistry","name":"Chemistry","version":"0.2.0"}, {},
                        {"molecular_views":[deepcopy(GRAPH)]},default_material())
        doc["results"]["molecular_views"][0]["bonds"][0]["end"]=999
        doc["checksum_sha256"]=digest(doc)
        with self.assertRaises(ValueError):validate_exchange(doc)

    def test_microstructure_record_fields_preserved(self):
        material=default_material();material["chemistry"].update({"repeat_smiles":"[*:1]CC[*:2]", "stereochemistry_source":"Unspecified"})
        self.assertEqual(validate_exchange(material),material)

if __name__=="__main__":unittest.main()
