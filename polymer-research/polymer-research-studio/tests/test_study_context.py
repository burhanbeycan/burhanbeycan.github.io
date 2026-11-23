import unittest
from copy import deepcopy
from polymer_workbench.exchange import default_material, validate_material, make_result, validate_exchange, digest
from polymer_workbench.lineage import derive_material
from polymer_workbench.study import study_audit,validate_study


class StudyContext(unittest.TestCase):
    def test_partial_declaration_records_missing_evidence(self):
        val=validate_study(dict(method='planned_md',population='EtOx hexamers',replicas=1))
        notes=study_audit(val)
        self.assertTrue(any('counterions' in n for n in notes))
        self.assertTrue(any('trajectory' in n for n in notes))
        self.assertTrue(any('variability' in n for n in notes))

    def test_bad_replica_values(self):
        for n in [True,0,-1,1.5,'2']:
            with self.assertRaises(ValueError):
                validate_study(dict(method='executed_md',replicas=n))

    def test_run_snapshot_and_derived_specimen(self):
        m=default_material(); m['study_context']=dict(method='executed_md',replicas=3,force_field='test force field')
        run=make_result(dict(id='test',name='test',version='0.4.0'),{},{},m)
        self.assertEqual(run['study']['declaration'],m['study_context'])
        validate_exchange(run)
        copy=deepcopy(run);copy['study']['specimen']='other';copy['checksum_sha256']=digest(copy)
        with self.assertRaises(ValueError):validate_exchange(copy)
        self.assertNotIn('study_context', derive_material(m,'mat','mat','electrospun_mat','electrospinning'))
        self.assertNotIn('study_context',validate_material(default_material()))
