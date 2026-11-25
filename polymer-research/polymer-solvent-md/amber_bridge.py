"""Use a native Amber ligand topology and OpenMM's stock three-site TIP3P water."""

import openmm as mm
from openmm import app, unit


def _force(system, kind):
    matches = [force for force in system.getForces() if isinstance(force, kind)]
    if len(matches) != 1:
        raise ValueError(f"Expected one {kind.__name__}; found {len(matches)}")
    return matches[0]


def _water_topology(prmtop, topology):
    """Verify ligand identity and extract only complete, neutral three-site waters."""
    reference = list(prmtop.topology.atoms())
    atoms = list(topology.atoms())
    n = len(reference)
    def identity(atom):
        return (atom.name, atom.element, atom.residue.index, atom.residue.name)

    if [identity(a) for a in atoms[:n]] != [identity(a) for a in reference]:
        raise ValueError("Solute atoms must be in the native Amber topology order")
    waters = app.Topology()
    chain = waters.addChain()
    water_atoms = {}
    for residue in topology.residues():
        members = list(residue.atoms())
        if members[0].index < n:
            if members[-1].index >= n:
                raise ValueError("A residue crosses the solute/solvent boundary")
            continue
        if sorted(a.element.symbol for a in members) != ["H", "H", "O"]:
            raise ValueError("Only explicit three-site water is supported after the solute")
        copied = waters.addResidue("HOH", chain)
        for a in members:
            water_atoms[a.index] = waters.addAtom(a.name, a.element, copied)
    solute_bonds = set()
    for a, b in topology.bonds():
        if a.index < n and b.index < n:
            solute_bonds.add(tuple(sorted((a.index, b.index))))
        elif a.index in water_atoms and b.index in water_atoms and a.residue == b.residue:
            waters.addBond(water_atoms[a.index], water_atoms[b.index])
        else:
            raise ValueError("Unexpected solute–water or inter-water covalent bond")
    expected = {tuple(sorted((a.index, b.index))) for a, b in prmtop.topology.bonds()}
    if solute_bonds != expected:
        raise ValueError("Solute bonds differ from the native Amber topology")
    return waters


def _append_water(system, water):
    """Append TIP3P parameters with index offsets to the intact GAFF2 System."""
    offset = system.getNumParticles()
    for i in range(water.getNumParticles()):
        if water.isVirtualSite(i):
            raise ValueError("Expected three-site TIP3P without virtual sites")
        system.addParticle(water.getParticleMass(i))
    for i in range(water.getNumConstraints()):
        a, b, distance = water.getConstraintParameters(i)
        system.addConstraint(offset+a, offset+b, distance)
    for source in water.getForces():
        if isinstance(source, mm.NonbondedForce):
            target = _force(system, mm.NonbondedForce)
            for i in range(source.getNumParticles()):
                target.addParticle(*source.getParticleParameters(i))
            for i in range(source.getNumExceptions()):
                a, b, charge, sigma, epsilon = source.getExceptionParameters(i)
                target.addException(offset+a, offset+b, charge, sigma, epsilon)
        elif isinstance(source, mm.HarmonicBondForce):
            if source.getNumBonds():
                target = _force(system, mm.HarmonicBondForce)
                for i in range(source.getNumBonds()):
                    a, b, distance, k = source.getBondParameters(i)
                    target.addBond(offset+a, offset+b, distance, k)
        elif isinstance(source, mm.HarmonicAngleForce):
            if source.getNumAngles():
                target = _force(system, mm.HarmonicAngleForce)
                for i in range(source.getNumAngles()):
                    a, b, c, angle, k = source.getAngleParameters(i)
                    target.addAngle(offset+a, offset+b, offset+c, angle, k)
        else:
            raise ValueError(f"Unsupported force in three-site water: {source.getName()}")


def create_solution_system(prmtop, solvated_topology, *, nonbondedMethod=app.PME,
                           nonbondedCutoff=0.9*unit.nanometer, constraints=app.HBonds,
                           rigidWater=True, ewaldErrorTolerance=0.0005):
    """Build a neutral GAFF2/TIP3P system; direct Amber import preserves all 1–4 exceptions."""
    methods = {app.NoCutoff: mm.NonbondedForce.NoCutoff, app.PME: mm.NonbondedForce.PME}
    if nonbondedMethod not in methods:
        raise ValueError("Use PME for dynamics or NoCutoff for verification")
    waters = _water_topology(prmtop, solvated_topology)
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=constraints,
                                rigidWater=False, removeCMMotion=True)
    supported = (mm.HarmonicBondForce, mm.HarmonicAngleForce, mm.PeriodicTorsionForce,
                 mm.NonbondedForce, mm.CMMotionRemover)
    if any(not isinstance(f, supported) for f in system.getForces()):
        raise ValueError("This bridge requires native Amber additive GAFF2 terms and standard LJ mixing")
    nonbonded = _force(system, mm.NonbondedForce)
    net_charge = sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                     for i in range(system.getNumParticles()))
    if abs(net_charge) > 1e-4:
        raise ValueError("This water model workflow currently requires a neutral solute")
    if waters.getNumAtoms():
        water = app.ForceField("tip3p.xml").createSystem(waters, nonbondedMethod=app.NoCutoff,
                      constraints=constraints, rigidWater=rigidWater, removeCMMotion=False)
        _append_water(system, water)
    box = solvated_topology.getPeriodicBoxVectors()
    if box is not None:
        system.setDefaultPeriodicBoxVectors(*box)
    elif nonbondedMethod == app.PME:
        raise ValueError("A periodic box is required for PME")
    nonbonded.setNonbondedMethod(methods[nonbondedMethod])
    nonbonded.setCutoffDistance(nonbondedCutoff)
    nonbonded.setEwaldErrorTolerance(ewaldErrorTolerance)
    nonbonded.setUseDispersionCorrection(True)
    return system


class _PlacementParameters:
    def __init__(self, prmtop):
        self.prmtop = prmtop

    def createSystem(self, topology, residueTemplates=None):
        if residueTemplates:
            raise ValueError("Explicit residue template overrides are unsupported")
        return create_solution_system(self.prmtop, topology, nonbondedMethod=app.NoCutoff,
                                      constraints=None)


def placement_forcefield(prmtop):
    """Provide Modeller.addSolvent() the same native ligand charges and LJ radii as dynamics."""
    return _PlacementParameters(prmtop)
