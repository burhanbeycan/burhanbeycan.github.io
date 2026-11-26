"""Actual application interoperability; skipped if the full suite is absent."""
import importlib.util
import json
import unittest

from fastapi.testclient import TestClient
from polymer_workbench.exchange import validate_exchange
from polymer_workbench.server import APPS, create_app

INSTALLED=all(importlib.util.find_spec(module.split('.')[0]) for module in APPS.values())


@unittest.skipUnless(INSTALLED,"Install all six laboratories to run suite integration checks.")
class SuiteIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.client=TestClient(create_app())
        cls.catalog=cls.client.get("/api/catalog").json()

    def calculate(self,app_id,inputs=None,material=None):
        response=self.client.post(f"/api/calculate/{app_id}",json={"inputs":inputs or {},"material":material})
        self.assertEqual(response.status_code,200,response.text[:1200])
        result=response.json()
        validate_exchange(result)
        json.dumps(result,allow_nan=False)
        return result

    def test_all_defaults_and_reference_scenarios(self):
        self.assertEqual({a['id'] for a in self.catalog['apps']},set(APPS))
        for app in self.catalog['apps']:
            for name,inputs in [('default',{})]+[(p['name'],p['inputs']) for p in app['presets']]:
                with self.subTest(application=app['id'],scenario=name):
                    self.calculate(app['id'],inputs)

    def test_standalone_launchers_expose_only_their_application(self):
        for app_id in APPS:
            with self.subTest(app=app_id):
                client=TestClient(create_app(app_id))
                catalog=client.get('/api/catalog').json()
                self.assertEqual(catalog['mode'],'standalone')
                self.assertEqual([a['id'] for a in catalog['apps']],[app_id])
                self.assertEqual(client.post('/api/calculate/missing',json={}).status_code,404)

    def test_network_geometry_and_transport_agree(self):
        import numpy as np
        from polymer_network_lab.model import slab_release_fraction
        from polymer_transport_lab.model import fraction_transferred
        network=self.calculate('network')
        response=self.client.post('/api/exchange/transfer',json={'document':network,'target':'transport'})
        self.assertEqual(response.status_code,200,response.text)
        transfer=response.json()
        mapped=transfer['inputs']
        self.assertEqual(mapped['mode'],'release')
        self.assertEqual(mapped['diffusivity_m2_s'],network['inputs']['diffusivity_m2_s'])
        thickness=mapped['thickness_um']*1e-6
        diffusion=mapped['diffusivity_m2_s']
        times=np.r_[0.0,np.geomspace(1e-10,3,120)]*thickness**2/diffusion
        np.testing.assert_allclose(slab_release_fraction(times,thickness,diffusion),fraction_transferred(times,diffusion,thickness),rtol=1e-10,atol=2e-14)
        transported=self.calculate('transport',mapped,transfer['material'])
        self.assertEqual(transported['material']['sample_id'],network['material']['sample_id'])
        self.assertEqual(transported['material']['history'][-2]['run_id'],network['run_id'])

    def test_chemistry_repeat_mass_handoff_is_explicit(self):
        chemistry=self.calculate('chemistry')
        views=chemistry['results']['molecular_views']
        self.assertTrue(views and views[0]['atoms'] and views[0]['bonds'])
        response=self.client.post('/api/exchange/transfer',json={'document':chemistry,'target':'chain'})
        self.assertEqual(response.status_code,200,response.text)
        mapped=response.json()['inputs']
        self.assertEqual(set(mapped),{'repeat_mass_g_mol','end_group_mass_g_mol'})
        self.assertGreater(mapped['repeat_mass_g_mol'],0)
        self.assertAlmostEqual(mapped['end_group_mass_g_mol'],2.016,places=3)
        self.assertNotIn('contour_increment_nm',mapped)
        transferred=self.calculate('chain',mapped,response.json()['material'])
        self.assertEqual(transferred['material']['sample_id'],chemistry['material']['sample_id'])

    def test_population_round_trip_and_no_implicit_diffusivity(self):
        first=self.calculate('chain')
        second=self.calculate('chain',{'mwd_source':'material'},first['material'])
        self.assertEqual(first['material']['molar_mass_distribution'],second['material']['molar_mass_distribution'])
        electrolyte=self.calculate('electrolyte',material=first['material'])
        response=self.client.post('/api/exchange/transfer',json={'document':electrolyte,'target':'transport'})
        self.assertEqual(response.status_code,200)
        self.assertEqual(response.json()['inputs'],{})
        session={'schema':'polymer-session/1.0','material':electrolyte['material'],'runs':[first,second,electrolyte]}
        self.assertEqual(self.client.post('/api/exchange/validate',json=session).status_code,200)


if __name__=='__main__':unittest.main()
