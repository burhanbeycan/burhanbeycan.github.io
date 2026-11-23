"""Classical, non-periodic all-atom trajectory pilot in Å, fs and kcal/mol.

Forces use the analytical MMFF94 gradients in RDKit. This module has no
dependency on an external MD executable and keeps molecular connectivity fixed.
"""
from __future__ import annotations

import math
import numpy as np

# 1 kcal/mol/Å per Da = 4.184e-4 Å/fs². Boltzmann constant in kcal/(mol K).
ACCELERATION = 4.184e-4
GAS_CONSTANT = 0.00198720425864083


def baoab_step(positions, velocities, force, masses, dt_fs, decay, noise_scale, rng, gradient_force):
    """B/A/O/A/B Langevin step. gradient_force returns kcal/(mol Å)."""
    inverse_mass = ACCELERATION / masses[:, None]
    velocities = velocities + (0.5 * dt_fs) * inverse_mass * force
    positions = positions + (0.5 * dt_fs) * velocities
    if decay < 1.0:
        velocities = decay * velocities + noise_scale * rng.standard_normal(velocities.shape)
    positions = positions + (0.5 * dt_fs) * velocities
    force = gradient_force(positions)
    velocities = velocities + (0.5 * dt_fs) * inverse_mass * force
    return positions, velocities, force


def run_mmff(mol, *, seed=1, temperature_K=298.15, equil_steps=500, production_steps=2000,
             stride=50, dt_fs=0.5, friction_ps=5.0):
    """Minimize a supplied explicit-H conformer and propagate an isolated cluster.

    Uses RDKit MMFF94 in vacuo (dielectric 1). Separate fragments interact.
    Trajectory time starts at the end of the initial thermalization interval.
    Rerunning with another velocity seed alone is not a new packing.
    """
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    ints = (seed, equil_steps, production_steps, stride)
    if any(isinstance(v,bool) or not isinstance(v,(int,np.integer)) for v in ints):
        raise ValueError("Seed, step counts and stride must be integers.")
    if not (0 <= seed < 2**32 and equil_steps >= 0 and production_steps >= 1 and
            1 <= stride <= production_steps and production_steps % stride == 0):
        raise ValueError("Select nonnegative seed/thermalization and positive steps divisible by stride.")
    therm = (temperature_K, dt_fs, friction_ps)
    if not all(math.isfinite(v) for v in therm) or not (0 < temperature_K <= 1000 and
            0 < dt_fs <= 0.5 and 0 <= friction_ps <= 1000):
        raise ValueError("Temperature, timestep or friction is outside the supported interval.")
    mol = Chem.Mol(mol)
    if not mol.GetNumAtoms() or mol.GetNumConformers()!=1:
        raise ValueError("A single initial 3D conformer is required.")
    if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in mol.GetAtoms()):
        raise ValueError("Represent all hydrogen atoms as explicit coordinate sites.")
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise ValueError("Initial coordinates must be finite.")
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError("Incomplete MMFF94 parameters for this graph.")
    properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    ff = AllChem.MMFFGetMoleculeForceField(mol, properties, nonBondedThresh=1e6,
                                           ignoreInterfragInteractions=False)
    initial_energy = ff.CalcEnergy()
    status = ff.Minimize(maxIts=2500, forceTol=1e-4, energyTol=1e-7)
    if status:
        raise ValueError("Initial minimization did not converge. Inspect packing or choose another seed.")
    positions = np.asarray(mol.GetConformer().GetPositions()).copy()
    minimum_energy = ff.CalcEnergy()
    if not math.isfinite(minimum_energy):
        raise ValueError("Nonfinite minimized energy.")
    masses = np.array([a.GetMass() for a in mol.GetAtoms()])
    rng = np.random.default_rng(seed)
    thermal_speed = np.sqrt(ACCELERATION * GAS_CONSTANT * temperature_K / masses[:, None])
    velocities = thermal_speed * rng.standard_normal(positions.shape)
    velocities -= np.average(velocities, axis=0, weights=masses)
    decay = math.exp(-friction_ps * dt_fs / 1000)
    noise_scale = math.sqrt(-math.expm1(-2 * friction_ps * dt_fs / 1000)) * thermal_speed

    def force(coordinates):
        point = tuple(coordinates.ravel())
        # RDKit energy refreshes the internal pair-distance cache before its
        # analytical gradient is evaluated at these external coordinates.
        ff.CalcEnergy(point)
        return -np.asarray(ff.CalcGrad(point)).reshape(-1, 3)

    forces = force(positions)
    times, snapshots, potential, kinetic = [], [], [], []
    for step in range(1, equil_steps + production_steps + 1):
        positions, velocities, forces = baoab_step(positions, velocities, forces, masses,
                                                  dt_fs, decay, noise_scale, rng, force)
        if step > equil_steps and (step-equil_steps) % stride == 0:
            if not np.isfinite(positions).all() or not np.isfinite(velocities).all():
                raise ValueError("Nonfinite trajectory. Inspect the system and timestep.")
            times.append((step - equil_steps) * dt_fs / 1000)
            snapshots.append(positions.copy() / 10)
            potential.append(ff.CalcEnergy(tuple(positions.ravel())))
            kinetic.append(float(np.sum(masses[:, None] * velocities**2) / (2 * ACCELERATION)))
    return dict(positions_nm=np.array(snapshots), time_ps=np.array(times), masses=masses,
                potential_kcal_mol=np.array(potential), kinetic_kcal_mol=np.array(kinetic),
                initial_energy_kcal_mol=initial_energy, minimized_energy_kcal_mol=minimum_energy,
                rdkit_version=rdBase.rdkitVersion, seed=int(seed), temperature_K=temperature_K,
                dt_fs=dt_fs, friction_ps=friction_ps, equil_steps=equil_steps,
                production_steps=production_steps, stride=stride)


def petox_cluster(dp, end, seed, waters=12):
    """Build α-N-methyl-PEtOx, with OH or organic azide terminal chemistry.

    The azide SMILES uses the neutral organic-azide resonance representation
    C–N=[N+]=[N−]. Water oxygen placement is independent for every packing seed.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from scipy.spatial.transform import Rotation
    if not isinstance(dp,int) or not 1 <= dp <= 20 or end not in {"OH","N3"}:
        raise ValueError("Choose DP 1–20 and OH or N3.")
    if not isinstance(waters,int) or not 0 <= waters <= 64:
        raise ValueError("Use 0–64 explicit waters.")
    if not isinstance(seed,int) or not 0 <= seed < 2**31-1:
        raise ValueError("Packing seed must be a nonnegative 31-bit integer.")
    smi = "C" + "N(C(=O)CC)CC" * dp + {"OH":"O","N3":"N=[N+]=[N-]"}[end]
    heavy = Chem.MolFromSmiles(smi)
    carbonyl = [a.GetIdx() for a in heavy.GetAtoms() if a.GetAtomicNum() == 8
                and any(b.GetBondTypeAsDouble()==2 for b in a.GetBonds())]
    terminal = heavy.GetNumAtoms()-1
    poly = Chem.AddHs(heavy)
    emb = AllChem.ETKDGv3()
    emb.randomSeed = seed
    if AllChem.EmbedMolecule(poly,emb) != 0:
        raise ValueError("Could not embed this oligomer; try another initialization seed.")
    AllChem.MMFFOptimizeMolecule(poly,maxIters=1200)
    poly_pos = poly.GetConformer().GetPositions()
    heavy_pos = poly_pos[:heavy.GetNumAtoms()]
    positions = list(poly_pos)
    water = Chem.AddHs(Chem.MolFromSmiles("O"))
    # Seed geometries only. MMFF94 freely minimizes and propagates all bonds.
    water_xyz = np.array([[0.,0.,0.],[0.9572,0,0],[-0.23999,0.92663,0.]])
    r = np.random.default_rng(seed)
    oxygen_pos = []
    mol = Chem.Mol(poly)
    for _ in range(waters):
        for _attempt in range(25000):
            origin = heavy_pos[r.integers(len(heavy_pos))]
            direction = r.normal(size=3)
            point = origin + direction/np.linalg.norm(direction)*r.uniform(2.9,4.1)
            if np.min(np.linalg.norm(heavy_pos-point,axis=1))<2.75:
                continue
            if oxygen_pos and np.min(np.linalg.norm(np.asarray(oxygen_pos)-point,axis=1))<2.8:
                continue
            placed = point + Rotation.random(random_state=r).apply(water_xyz)
            if np.min(np.linalg.norm(np.asarray(positions)[:,None,:]-placed[None,:,:],axis=2))<1.0:
                continue
            break
        else:
            raise ValueError("Could not pack waters without overlap.")
        oxygen_pos.append(point)
        positions.extend(placed)
        mol = Chem.CombineMols(mol,water)
    mol.RemoveAllConformers()
    conf=Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for index,point in enumerate(positions):
        conf.SetAtomPosition(index,point)
    mol.AddConformer(conf)
    Chem.SanitizeMol(mol)
    groups = dict(polymer=np.arange(poly.GetNumAtoms()), carbonyl=np.array(carbonyl,dtype=int),
                  omega=np.array([terminal]), water_oxygen=np.arange(poly.GetNumAtoms(),mol.GetNumAtoms(),3))
    return mol,groups


def frame_observables(coordinates_nm, groups, polymer_masses, cutoff_nm, box_nm=None):
    """One whole polymer, water oxygen distances, equal weight for sites within a class."""
    xyz = np.asarray(coordinates_nm, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.all(np.isfinite(xyz)):
        raise ValueError("Coordinates must be finite N × 3 positions in nm.")
    polymer, water = np.asarray(groups["polymer"]), np.asarray(groups["water_oxygen"])
    masses = np.asarray(polymer_masses,dtype=float)
    for name in ("polymer", "water_oxygen", "carbonyl", "omega"):
        selection = np.asarray(groups[name])
        if (selection.ndim != 1 or not len(selection) or not np.issubdtype(selection.dtype,np.integer)
                or len(np.unique(selection)) != len(selection) or np.any(selection<0) or np.any(selection>=len(xyz))):
            raise ValueError(f"Invalid or duplicate atom indices in {name}.")
    if not (set(groups["omega"]) | set(groups["carbonyl"])) <= set(polymer):
        raise ValueError("Chemical sites must belong to the selected polymer.")
    if set(groups["omega"]) & set(groups["carbonyl"]) or set(polymer) & set(water):
        raise ValueError("Site classes and polymer/water selections must be disjoint.")
    if masses.shape!=(len(polymer),) or not np.all(np.isfinite(masses)) or np.any(masses<=0):
        raise ValueError("A positive finite topology mass is required for each polymer atom.")
    if not math.isfinite(cutoff_nm) or not 0<cutoff_nm<=1:
        raise ValueError("Choose a contact distance in (0, 1] nm.")
    if box_nm is not None:
        box_nm=np.asarray(box_nm,dtype=float)
        if box_nm.shape!=(3,) or not np.all(np.isfinite(box_nm)) or min(box_nm)<=2*cutoff_nm:
            raise ValueError("Use an orthorhombic box larger than twice the contact cutoff.")
    chain_xyz=xyz[polymer]
    centre=np.average(chain_xyz,weights=masses,axis=0)
    rg=float(np.sqrt(np.average(np.sum((chain_xyz-centre)**2,axis=1),weights=masses)))

    def occupancy(sites):
        n=0
        for site in sites:
            diff=xyz[water]-xyz[site]
            if box_nm is not None:
                diff -= box_nm*np.rint(diff/box_nm)
            n+=np.count_nonzero(np.einsum("ij,ij->i",diff,diff)<=cutoff_nm**2)
        return n/len(sites)

    return rg, occupancy(groups["carbonyl"]), occupancy(groups["omega"])
