import csv
import io
import math
from pathlib import Path
import tempfile
import unittest
from importlib.util import find_spec

import numpy as np
from polymer_chain_lab.classical_md import (ACCELERATION, baoab_step, frame_observables,
                                            petox_cluster, run_mmff)
from polymer_chain_lab.md_workflow import analyze_gro, md_summary
from polymer_chain_lab.plugin import calculate


class Integrator(unittest.TestCase):
    def test_acceleration_and_verlet_energy(self):
        self.assertAlmostEqual(ACCELERATION, 4184e20 / 1e30 * 1000, places=12)
        mass=np.array([1.]); x=np.array([[0.05,0,0]]); v=np.zeros_like(x); spring=400.
        force=lambda p:-spring*p
        f=force(x); e0=spring*0.05**2/2
        for _ in range(1000):
            x,v,f=baoab_step(x,v,f,mass,0.25,1.,0.,None,force)
            e=(spring*(x*x).sum()+(v*v).sum()/ACCELERATION)*0.5
            self.assertLess(abs(e-e0)/e0,0.003)

    @unittest.skipUnless(find_spec('rdkit'),"RDKit is an optional MD-runner dependency")
    def test_nve_stability_and_observer_stride_independence(self):
        mol,group=petox_cluster(2,"OH",6,2)
        kwargs=dict(seed=1234,equil_steps=0,production_steps=200,dt_fs=0.25,friction_ps=0.)
        freq=run_mmff(mol,stride=1,**kwargs); sparse=run_mmff(mol,stride=20,**kwargs)
        np.testing.assert_allclose(freq['positions_nm'][19::20],sparse['positions_nm'],rtol=0,atol=1e-11)
        e=freq['potential_kcal_mol']+freq['kinetic_kcal_mol']
        self.assertLess(np.ptp(e)/np.mean(freq['kinetic_kcal_mol']),0.01)
        self.assertEqual(len(group['carbonyl']),2)
        self.assertEqual(len(group['water_oxygen']),2)

    @unittest.skipUnless(find_spec('rdkit'),"RDKit is an optional MD-runner dependency")
    def test_dp4_chemistry_matches_expected_formulae(self):
        from rdkit.Chem import rdMolDescriptors
        for end,formula in (('OH','C21H40N4O5'),('N3','C21H39N7O4')):
            m,g=petox_cluster(4,end,109,0)
            self.assertEqual(rdMolDescriptors.CalcMolFormula(m),formula)
            atom=m.GetAtomWithIdx(int(g['omega'][0]))
            self.assertEqual(atom.GetSymbol(),'O' if end=='OH' else 'N')


G={"polymer":[0,1,2],"omega":[2],"carbonyl":[0],"water_oxygen":[3,4]}
MAP=dict(**G,sample="mapping unit test",force_field="analytic fixture",water_model="analytic fixture",
         temperature_K=298.15,source="unit test",dp=1,end="OH",replica="A",polymer_bonds=[[0,1],[1,2]],
         polymer_masses_Da=[16,12,16],carbonyl_names=["P:O1"],omega_names=["P:O2"],water_oxygen_names=["W:OW"])


def grofile(t,xyz):
    names=[(1,'P','O1'),(1,'P','C'),(1,'P','O2'),(2,'W','OW'),(3,'W','OW')]
    lines=[f'Analytic atom-coordinate fixture t= {t}',str(len(names))]
    for i,((res,r,n),pt) in enumerate(zip(names,xyz)):
        lines.append(f'{res:5d}{r:<5s}{n:>5s}{i+1:5d}'+''.join(f'{v:8.3f}' for v in pt))
    return '\n'.join(lines+[' 2.00000 2.00000 2.00000',''])


class Observables(unittest.TestCase):
    def setUp(self):
        self.xyz=np.array([[0.05,0,0],[0.2,0,0],[0.35,0,0],[1.9,0,0],[0.45,0,0]])
    def test_box_and_mass_radius(self):
        rg,co,omega=frame_observables(self.xyz,G,[16,12,16],0.23,[2,2,2])
        self.assertAlmostEqual(rg,0.15*math.sqrt(32/44))
        self.assertEqual((co,omega),(1.,1.))
        self.assertEqual(frame_observables(self.xyz,G,[16,12,16],0.23)[1],0.)
    def test_gro_mapping_names_timing_pbc_and_unwrapping(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fixture.gro'
            p.write_text(grofile(0,self.xyz)+grofile(1,self.xyz))
            rows=analyze_gro(p,MAP,0.23)
            self.assertEqual(len(rows),2)
            self.assertEqual(rows[1]['carbonyl_coordination'],1.)
            with self.assertRaisesRegex(ValueError,'mismatch'):analyze_gro(p,dict(MAP,omega_names=["P:WRONG"]))
            self.xyz[1,0]=1.9;p.write_text(grofile(0,self.xyz))
            with self.assertRaisesRegex(ValueError,'whole'):analyze_gro(p,MAP)
    def test_bundle_frame_and_replica_balance(self):
        groups,meta=md_summary('',4)
        self.assertEqual((len(groups['OH']),len(groups['N3'])),(3,3))
        self.assertEqual([len(a) for _,a in groups['OH']],[40]*3)
        out=calculate(dict(model="petox_md",pilot_dp=4))
        metrics={m['key']:m['value'] for m in out['metrics']}
        self.assertAlmostEqual(metrics['OH_omega'],1.4416667,6)
        self.assertAlmostEqual(metrics['N3_omega'],0.0583333,6)
    def test_mixed_simulations_and_duplicate_frames_rejected(self):
        from importlib.resources import files
        p=files('polymer_chain_lab').joinpath('data').joinpath('petox_cluster_md.csv')
        data=p.read_text()
        with self.assertRaisesRegex(ValueError,'matched'):
            md_summary(data.replace('MMFF94 water','incompatible water',1),2)
        with self.assertRaisesRegex(ValueError,'ordered'):
            md_summary(data+'\n'+data.splitlines()[1],2)
