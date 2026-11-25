"""Atom map and finite-cell checks used by the explicit-water MD runner."""
import unittest
import numpy as np
from rdkit import Chem
from openmm import app
from simulation import _integer_steps, _self_image, _solvent_count, _starting_geometry


class SimulationGeometryTest(unittest.TestCase):
    def test_periodic_copy_distance_excludes_own_molecule(self):
        box = np.diag([2.0, 3.0, 4.0])
        xyz = np.array([[0.25, 0.0, 0.0], [1.75, 0.0, 0.0]])
        self.assertAlmostEqual(_self_image(xyz, box), 0.5)
        self.assertAlmostEqual(_self_image(xyz + [4.0, -3.0, 8.0], box), 0.5)

    def test_whole_box_does_not_shrink_with_an_extended_chain(self):
        short = np.array([[0.0, 0, 0], [0.9, 0, 0]])
        extended = np.array([[0.0, 0, 0], [2.1, 0, 0]])
        self.assertEqual(_solvent_count(short, None), 500)
        self.assertGreater(_solvent_count(extended, None), 500)
        with self.assertRaises(ValueError):
            _solvent_count(extended * 3, None)

    def test_carbonyl_site_survives_reversed_force_field_atom_order(self):
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)N(C)CCO"))
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        topology = app.Topology()
        residue = topology.addResidue("LIG", topology.addChain())
        remap = {}
        for atom in reversed(list(mol.GetAtoms())):
            remap[atom.GetIdx()] = topology.addAtom(
                atom.GetSymbol(), app.Element.getBySymbol(atom.GetSymbol()), residue)
        for bond in mol.GetBonds():
            topology.addBond(remap[bond.GetBeginAtomIdx()], remap[bond.GetEndAtomIdx()])
        _, selected, _ = _starting_geometry(smiles, topology, 21, None)
        oxygen, _, _ = mol.GetSubstructMatch(Chem.MolFromSmarts("O=C-N"))
        self.assertEqual(selected, [remap[oxygen].index])

    def test_sample_intervals_are_integer_timestep_multiples(self):
        self.assertEqual(_integer_steps(0.2, 0.002), 100)
        with self.assertRaises(ValueError):
            _integer_steps(0.201, 0.002)


if __name__ == "__main__":
    unittest.main()
