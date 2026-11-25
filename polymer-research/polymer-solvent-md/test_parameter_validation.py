"""Model 08 native ligand mapping, force terms and solute–water mixing regressions."""

from collections import Counter
from pathlib import Path
import unittest

import numpy as np
import openmm as mm
from openmm import app, unit
from rdkit import Chem

from amber_bridge import create_solution_system, placement_forcefield


PARAMETERS = Path(__file__).resolve().parent / "pilot" / "parameters"
KJ = unit.kilojoules_per_mole
NM = unit.nanometer
E = unit.elementary_charge


def evaluate(system, xyz):
    """Compute each native term plus full gradient on OpenMM Reference."""
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(xyz * NM)
    energies = {
        f.getName(): context.getState(getEnergy=True, groups=1 << i)
                     .getPotentialEnergy().value_in_unit(KJ)
        for i, f in enumerate(system.getForces())
    }
    forces = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(KJ/NM)
    return energies, forces


def nb(system):
    return next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))


def particle_parameters(force):
    return np.array([[q.value_in_unit(E), s.value_in_unit(NM), eps.value_in_unit(KJ)]
                     for q, s, eps in (force.getParticleParameters(i)
                                       for i in range(force.getNumParticles()))])


class TestNativeGAFF(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prmtop = app.AmberPrmtopFile(str(PARAMETERS / "solute.prmtop"))
        coords = app.AmberInpcrdFile(str(PARAMETERS / "solute.inpcrd"))
        cls.xyz = np.asarray(coords.positions.value_in_unit(NM))
        cls.sdf = Chem.SDMolSupplier(str(PARAMETERS / "solute.sdf"), removeHs=False)[0]

    def test_atom_graph_and_signed_charge(self):
        atoms = list(self.prmtop.topology.atoms())
        self.assertEqual([a.element.symbol for a in atoms], [a.GetSymbol() for a in self.sdf.GetAtoms()])
        self.assertEqual(Counter(a.element.symbol for a in atoms), {"C": 15, "H": 29, "N": 3, "O": 3})
        self.assertEqual(len(set(a.name for a in atoms)), len(atoms))
        supplied = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in self.sdf.GetBonds()}
        imported = {tuple(sorted((a.index, b.index))) for a, b in self.prmtop.topology.bonds()}
        self.assertEqual(supplied, imported)
        system = self.prmtop.createSystem(constraints=None)
        parameters = particle_parameters(nb(system))
        self.assertLess(abs(parameters[:, 0].sum() - Chem.GetFormalCharge(self.sdf)), 1e-4)
        self.assertTrue(np.isfinite(parameters).all())
        self.assertTrue((parameters[:, 2] >= 0).all())
        self.assertTrue((parameters[:, 1] > 0).all())
        self.assertTrue((parameters[[a.index for a in atoms if a.element != app.element.hydrogen], 2] > 0).all())
        self.assertTrue(all(np.isfinite(self.xyz.ravel())))
        frcmod = (PARAMETERS / "solute.frcmod").read_text().upper()
        self.assertNotIn("ATTN", frcmod)
        self.assertNotIn("NEED REVISION", frcmod)

    def test_native_intramolecular_energy_and_forces(self):
        direct = self.prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
        bridge = create_solution_system(self.prmtop, self.prmtop.topology,
                                        nonbondedMethod=app.NoCutoff, constraints=None)
        nref, nbridge = nb(direct), nb(bridge)
        np.testing.assert_equal(particle_parameters(nbridge), particle_parameters(nref))
        self.assertGreater(nref.getNumExceptions(), 0)
        self.assertEqual(nref.getNumExceptions(), nbridge.getNumExceptions())
        for i in range(nref.getNumExceptions()):
            self.assertEqual(nref.getExceptionParameters(i), nbridge.getExceptionParameters(i))
        rng = np.random.default_rng(838)
        for xyz in (self.xyz, self.xyz + rng.normal(0, 0.002, self.xyz.shape)):
            reference_energy, reference_force = evaluate(direct, xyz)
            observed_energy, observed_force = evaluate(bridge, xyz)
            for term, value in reference_energy.items():
                self.assertAlmostEqual(observed_energy[term], value, places=8, msg=term)
            np.testing.assert_allclose(observed_force, reference_force, atol=1e-7, rtol=0)

    def test_tip3p_cross_energy_and_exclusions(self):
        oxygen = np.max(self.xyz, axis=0) + [0.27, 0.03, 0.02]
        # One equilibrium TIP3P water, nonperiodic single point.
        angle = 1.82421813418
        h1 = oxygen + [0.09572, 0, 0]
        h2 = oxygen + 0.09572*np.array([np.cos(angle), np.sin(angle), 0])
        water_xyz = np.array([oxygen, h1, h2])
        water_top = app.Topology()
        residue = water_top.addResidue("HOH", water_top.addChain())
        o, h1, h2 = [water_top.addAtom(name, elem, residue) for name, elem in
                     [("O", app.element.oxygen), ("H1", app.element.hydrogen), ("H2", app.element.hydrogen)]]
        water_top.addBond(o, h1)
        water_top.addBond(o, h2)
        model = app.Modeller(self.prmtop.topology, self.xyz*NM)
        model.add(water_top, water_xyz*NM)
        native = self.prmtop.createSystem(constraints=None)
        mixed = create_solution_system(self.prmtop, model.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        n = native.getNumParticles()
        m_nb = nb(mixed)
        self.assertEqual(m_nb.getNumExceptions(), nb(native).getNumExceptions()+3)
        water_distances = sorted(mixed.getConstraintParameters(i)[2].value_in_unit(NM)
                                 for i in range(mixed.getNumConstraints()))
        np.testing.assert_allclose(water_distances,
                [0.09572, 0.09572, 2 * 0.09572 * np.sin(angle/2)], atol=1e-9)
        actual_params = particle_parameters(m_nb)
        reference_params = particle_parameters(nb(native))
        np.testing.assert_allclose(actual_params[:n], reference_params, atol=0, rtol=0)
        np.testing.assert_allclose(actual_params[n:], [[-.834,.31507524065751241,.635968],
                                                          [.417,1.,0.],[.417,1.,0.]], atol=1e-12)
        distances = np.linalg.norm(self.xyz[:, None, :] - water_xyz[None, :, :], axis=-1)
        q1, s1, e1 = reference_params.T
        q2, s2, e2 = actual_params[n:].T
        ratio6 = (0.5*(s1[:, None] + s2[None, :])/distances)**6
        coulomb = 138.93545764438198 * q1[:, None]*q2[None, :]/distances
        lj = 4*np.sqrt(e1[:, None]*e2[None, :])*(ratio6**2-ratio6)
        native_e, _ = evaluate(native, self.xyz)
        mixed_e, _ = evaluate(mixed, np.vstack((self.xyz, water_xyz)))
        observed = mixed_e["NonbondedForce"] - native_e["NonbondedForce"]
        self.assertAlmostEqual(observed, np.sum(coulomb+lj), places=6)

    def test_periodic_placement_uses_native_sizes(self):
        model = app.Modeller(self.prmtop.topology, self.xyz*NM)
        model.addSolvent(placement_forcefield(self.prmtop), model="tip3p",
                         boxSize=mm.Vec3(2.4, 2.4, 2.4)*NM, neutralize=False)
        mixed = create_solution_system(self.prmtop, model.topology)
        natoms = self.prmtop.topology.getNumAtoms()
        nwater = (mixed.getNumParticles()-natoms)//3
        self.assertGreater(nwater, 100)
        self.assertEqual(mixed.getNumConstraints(), 29+3*nwater)
        self.assertEqual(nb(mixed).getNonbondedMethod(), mm.NonbondedForce.PME)


if __name__ == "__main__":
    unittest.main()
