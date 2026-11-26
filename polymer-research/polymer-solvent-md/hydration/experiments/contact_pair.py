"""Simulate and analyze a deliberately contacting, finite PEtOx trimer pair in explicit water.

Run from the released project with --model path/to/polymer-solvent-md. Each run
uses two copies of the same native GAFF2 solute topology. TIP3P waters are
packed around the two chains and both chains are free throughout production.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
import openmm as mm
from openmm import app, unit


def checksum(path):
    with open(path, 'rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def save_json(path, value):
    pending = path.with_suffix('.pending')
    pending.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    os.replace(pending, path)


def nonbonded(system):
    return next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))


def duplicate_solute(native):
    """Copy all terms, constraints and exceptions, shifting every atom index."""
    s = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(native))
    n = s.getNumParticles()
    for i in range(n):
        if s.isVirtualSite(i):
            raise ValueError("Expected physical solute atoms")
        s.addParticle(s.getParticleMass(i))
    for j in range(s.getNumConstraints()):
        a,b,d = s.getConstraintParameters(j)
        s.addConstraint(n+a,n+b,d)
    for f in s.getForces():
        if isinstance(f, mm.NonbondedForce):
            for j in range(n):
                f.addParticle(*f.getParticleParameters(j))
            for j in range(f.getNumExceptions()):
                a,b,q,sigma,eps = f.getExceptionParameters(j)
                f.addException(n+a,n+b,q,sigma,eps)
        elif isinstance(f, mm.HarmonicBondForce):
            for j in range(f.getNumBonds()):
                a,b,d,k = f.getBondParameters(j)
                f.addBond(n+a,n+b,d,k)
        elif isinstance(f, mm.HarmonicAngleForce):
            for j in range(f.getNumAngles()):
                a,b,c,theta,k = f.getAngleParameters(j)
                f.addAngle(n+a,n+b,n+c,theta,k)
        elif isinstance(f, mm.PeriodicTorsionForce):
            for j in range(f.getNumTorsions()):
                a,b,c,d,m,phi,k = f.getTorsionParameters(j)
                f.addTorsion(n+a,n+b,n+c,n+d,m,phi,k)
        elif not isinstance(f, mm.CMMotionRemover):
            raise ValueError('Unrecognized solute force '+type(f).__name__)
    return s


def validate_pair(native,s):
    """Require equality to the native topology, per term and per copy."""
    n=native.getNumParticles()
    if s.getNumParticles()!=2*n or s.getNumConstraints()!=2*native.getNumConstraints():
        raise ValueError('Duplicate topology has the wrong size')
    def check(sequence, source, target, k, nindices):
        for i in range(k):
            old=list(source(i))
            for copy in (0,1):
                row=list(target(i+copy*k))
                if row[:nindices]!=[v+n*copy for v in old[:nindices]] or row[nindices:]!=old[nindices:]:
                    raise ValueError('Parameter mismatch for '+sequence)
    check('constraint',native.getConstraintParameters,s.getConstraintParameters,native.getNumConstraints(),2)
    signatures=[]
    for f,g in zip(native.getForces(),s.getForces(),strict=True):
        if isinstance(f, mm.CMMotionRemover):
            continue
        for label,nidx in [('Particle',0),('Exception',2),('Bond',2),('Angle',3),('Torsion',4)]:
            count=getattr(f,'getNum'+label+'s',None)
            if count and hasattr(f,'get'+label+'Parameters'):
                size=count()
                if getattr(g,'getNum'+label+'s')()!=2*size:
                    raise ValueError('Unexpected number of '+label)
                check(label,getattr(f,'get'+label+'Parameters'),getattr(g,'get'+label+'Parameters'),size,nidx)
                signatures.append([f.getName(),label,size])
    q=sum(nonbonded(s).getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(2*n))
    if abs(q)>2e-4:
        raise ValueError('Pair net charge inconsistent with neutral model')
    return dict(native_term_counts=signatures,net_pair_charge_e=q,native_atoms=n,constraints_per_chain=native.getNumConstraints())


def waters_from(topology, offset):
    wt=app.Topology()
    chain=wt.addChain()
    for r in topology.residues():
        atoms=list(r.atoms())
        if atoms[0].index<offset:
            continue
        if [x.element.symbol for x in atoms]!=['O','H','H']:
            raise ValueError('Invalid TIP3P atom or order')
        res=wt.addResidue('HOH',chain)
        o,h1,h2=[wt.addAtom(x.name,x.element,res) for x in atoms]
        wt.addBond(o,h1)
        wt.addBond(o,h2)
    return wt


class ContactSystem:
    def __init__(self,prmtop,append_water):
        self.single=prmtop.createSystem(nonbondedMethod=app.NoCutoff,constraints=app.HBonds,
                                        removeCMMotion=True,rigidWater=False)
        self.pair=duplicate_solute(self.single)
        self.validation=validate_pair(self.single,self.pair)
        self.water_ff=app.ForceField('tip3p.xml')
        self.append_water=append_water

    def createSystem(self,topology,residueTemplates=None):
        s=mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(self.pair))
        nw=waters_from(topology,s.getNumParticles())
        if nw.getNumAtoms():
            water=self.water_ff.createSystem(nw,nonbondedMethod=app.NoCutoff,constraints=app.HBonds,
                                             rigidWater=True,removeCMMotion=False)
            self.append_water(s,water)
        return s


def make_contact(xyz_a,xyz_b,heavy,seed):
    """Draw a relative orientation and translate along one axis to a prescribed nonoverlapping contact."""
    rng=np.random.default_rng(seed)
    x0=xyz_a-xyz_a[heavy].mean(axis=0)
    y0=xyz_b-xyz_b[heavy].mean(axis=0)
    for trial in range(300):
        x=Rotation.random(random_state=rng).apply(x0)
        y=Rotation.random(random_state=rng).apply(y0)
        yz=cdist(x[heavy,1:],y[heavy,1:])
        mask=yz<0.34
        if not mask.any():
            continue
        # Coming from +x, the first heavy contact is exactly 0.34 nm.
        dx=x[heavy,0,None]-y[None,heavy,0]
        translation=np.max((dx+np.sqrt(np.maximum(0,0.34**2-yz**2)))[mask])
        y[:,0]+=translation
        distances=cdist(x[heavy],y[heavy])
        if ((distances<=0.45).sum()>=5 and cdist(x,y).min()>=0.13
                and max(np.ptp(np.vstack([x,y]),axis=0))<1.85):
            return x,y,dict(generation_trial=trial,initial_heavy_contacts_045=int(np.sum(distances<=.45)),
                            target_minimum_heavy_distance_nm=.34,minimum_all_atom_distance_nm=float(cdist(x,y).min()))
    raise ValueError('Could not pack an eligible non-overlapping contacting orientation')


def mc(v,L):
    return v-L*np.rint(v/L)


def make_whole(xyz,indices,bonds,L):
    index={k:i for i,k in enumerate(indices)}
    graph=[[] for _ in indices]
    for a,b in bonds:
        if a in index and b in index:
            ai,bi=index[a],index[b]
            graph[ai].append(bi); graph[bi].append(ai)
    whole=np.empty((len(indices),3))
    whole[0]=xyz[indices[0]]
    visited={0}; stack=[0]
    while stack:
        j=stack.pop()
        for k in graph[j]:
            if k not in visited:
                whole[k]=whole[j]+mc(xyz[indices[k]]-xyz[indices[j]],L)
                visited.add(k); stack.append(k)
    if len(visited)!=len(indices):
        raise ValueError('Disconnected chain')
    return whole


def pair_dist(a,b,L):
    delta=mc(a[:,None,:]-b[None,:,:],L)
    return np.sqrt(np.einsum('ijk,ijk->ij',delta,delta))


def metrics(xyz,box,top):
    L=np.diag(box)
    if np.max(abs(box-np.diag(L)))>1e-6:
        raise ValueError('Expected orthorhombic box')
    A,B=top['chain_a_indices'],top['chain_b_indices']
    a=make_whole(xyz,A,top['bonds'],L); b=make_whole(xyz,B,top['bonds'],L)
    masses=np.array(top['mass_da'])[A]
    com=mc(np.average(a,axis=0,weights=masses)-np.average(b,axis=0,weights=masses),L)
    ha=np.array([i for i in A if top['elements'][i]!='H'])
    hb=np.array([i for i in B if top['elements'][i]!='H'])
    ab=pair_dist(xyz[ha],xyz[hb],L)
    wo=np.asarray(top['water_oxygens'])
    nearby_a=pair_dist(xyz[ha],xyz[wo],L).min(axis=0)<=.35
    nearby_b=pair_dist(xyz[hb],xyz[wo],L).min(axis=0)<=.35
    def hydration(oxygens):
        # All D–H···O acceptor triplets. Unlike O–O coordination, this requires angular alignment.
        acc=xyz[oxygens]
        da=mc(acc[:,None,:]-xyz[wo][None,:,:],L)
        close=np.sum(da**2,axis=2)<=.35**2
        hb=0
        for offset in (1,2):
            oh=mc(xyz[wo+offset]-xyz[wo],L)
            h_to_a=mc(da-oh[None,:,:],L)
            numerator=np.einsum('jk,ijk->ij',-oh,h_to_a)
            denominator=np.linalg.norm(oh,axis=1)[None,:]*np.linalg.norm(h_to_a,axis=2)
            hb+=np.count_nonzero(close&(numerator/np.maximum(denominator,1e-12)<=math.cos(math.radians(150))))
        return int(close.sum()), int(hb)
    ca,hba=hydration(top['carbonyl_oxygens_a']);cb,hbb=hydration(top['carbonyl_oxygens_b'])
    return dict(heavy_pair_contacts_045=int(np.sum(ab<=.45)),min_heavy_distance_nm=float(ab.min()),
                com_separation_nm=float(np.linalg.norm(com)),shared_water_035=int(np.sum(nearby_a&nearby_b)),
                unique_water_a_035=int(nearby_a.sum()),unique_water_b_035=int(nearby_b.sum()),
                carbonyl_water_pairs_a_035=ca,carbonyl_water_pairs_b_035=cb,
                water_donor_hbonds_a=hba,water_donor_hbonds_b=hbb)


def run(a):
    from simulation import _starting_geometry
    from amber_bridge import _append_water
    start=time.monotonic()
    p=a.model/'pilot/parameters'
    source=json.loads((a.model/'pilot/source_model.json').read_text())
    params=json.loads((p/'parameterization.json').read_text())
    prmtop=app.AmberPrmtopFile(str(p/'solute.prmtop'))
    # First use the same attributed conformation, second an independent ETKDG/MMFF warm start.
    smi=source['structure']['canonical_smiles']
    x,co,_=_starting_geometry(smi,prmtop.topology,a.seed,source)
    y,co2,geom=_starting_geometry(smi,prmtop.topology,a.seed+17,None)
    if co!=co2: raise ValueError('Carbonyl selections differ')
    elements=[atom.element.symbol for atom in prmtop.topology.atoms()]
    heavy=np.array([i for i,el in enumerate(elements) if el!='H'])
    left,right,packing=make_contact(x,y,heavy,a.seed)
    n=len(x)
    m=app.Modeller(prmtop.topology,left*unit.nanometer)
    m.add(prmtop.topology,right*unit.nanometer)
    factory=ContactSystem(prmtop,_append_water)
    m.addSolvent(factory,model='tip3p',numAdded=a.waters,neutralize=False)
    system=factory.createSystem(m.topology)
    system.setDefaultPeriodicBoxVectors(*m.topology.getPeriodicBoxVectors())
    f=nonbonded(system)
    f.setNonbondedMethod(mm.NonbondedForce.PME)
    f.setCutoffDistance(.8*unit.nanometer)
    f.setEwaldErrorTolerance(.0005)
    f.setUseDispersionCorrection(True)
    atoms=list(m.topology.atoms())
    nat=len(atoms)
    top=dict(elements=[i.element.symbol for i in atoms],atom_names=[i.name for i in atoms],
             bonds=[[i.index,j.index] for i,j in m.topology.bonds()],chain_a_indices=list(range(n)),
             chain_b_indices=list(range(n,2*n)),carbonyl_oxygens_a=co,carbonyl_oxygens_b=[i+n for i in co],
             water_oxygens=list(range(2*n,nat,3)),mass_da=[system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(nat)])
    if len(top['water_oxygens'])!=a.waters: raise ValueError('Incorrect water count')
    if a.out.exists() and any(a.out.iterdir()): raise ValueError('Use a new output directory')
    a.out.mkdir(parents=True,exist_ok=True)
    save_json(a.out/'topology.json',top)
    bar=mm.MonteCarloBarostat(1*unit.bar,298.15*unit.kelvin,25)
    bar.setRandomNumberSeed(a.seed+2);system.addForce(bar)
    # Allow water and box to relax without forcing an immediate, arbitrary pair displacement.
    restr=mm.CustomExternalForce('0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2')
    restr.addGlobalParameter('k',1000.)
    for v in ('x0','y0','z0'): restr.addPerParticleParameter(v)
    input_xyz=np.asarray(m.positions.value_in_unit(unit.nanometer))
    for i in np.r_[heavy,heavy+n]:
        restr.addParticle(int(i),input_xyz[i])
    system.addForce(restr)
    integrator=mm.LangevinMiddleIntegrator(298.15*unit.kelvin,1/unit.picosecond,.002*unit.picoseconds)
    integrator.setRandomNumberSeed(a.seed+1);integrator.setConstraintTolerance(1e-5)
    sim=app.Simulation(m.topology,system,integrator,mm.Platform.getPlatformByName('CPU'),{'Threads':'1','DeterministicForces':'true'})
    sim.context.setPositions(m.positions)
    sim.minimizeEnergy(tolerance=10*unit.kilojoules_per_mole/unit.nanometer,maxIterations=400)
    sim.context.setVelocitiesToTemperature(298.15*unit.kelvin,a.seed)
    sim.step(round(a.condition_ps/.002))
    sim.context.setParameter('k',0)
    manifest=dict(schema='contact-pair/1',seed=a.seed,temperature_K=298.15,pressure_bar=1,thermostat='LangevinMiddle 1/ps',
                  engine='OpenMM '+mm.version.version,ensemble='NPT',force_fields='GAFF2 + AM1-BCC / TIP3P',
                  parameter_id=params['parameter_id'],native_prmtop_sha256=checksum(p/'solute.prmtop'),input_artifact=source['artifact_id'],
                  water_molecules=a.waters,canonical_smiles=smi,initial_geometry_b=geom,packing=packing,
                  system_validation=factory.validation,condition_ps=a.condition_ps,condition_restraint_kj_mol_nm2=1000,
                  condition='Periodic heavy-atom Cartesian restraints during minimization and conditioning only',
                  production_ps=a.ps,frame_stride_ps=a.stride_ps,dt_ps=.002,real_space_cutoff_nm=.8,
                  production_frame_zero='After conditioning, positional restraint strength set to zero; before first free step.',
                  observables='All distances minimum-image. Cross-chain heavy pairs at d≤0.45 nm; shared waters are unique Ow within 0.35 nm of a heavy atom on each chain. Water donor H bonds: Ow···carbonyl O≤0.35 nm, Ow–H···O≥150°.',
                  scope='Two short chains in finite periodic water. Local contact reorganization, no bulk phase fraction or equilibrium solubility.')
    save_json(a.out/'protocol.json',manifest)
    print('Conditioned',a.seed,round(time.monotonic()-start,2),flush=True)
    pos=[];boxes=[];times=[];table=[]
    def save_track():
        tmp=a.out/'trajectory.part.npz'
        np.savez_compressed(tmp,positions_nm=np.asarray(pos,dtype=np.float32),box_nm=np.asarray(boxes,dtype=np.float32),
                            time_ps=np.asarray(times,dtype=float))
        os.replace(tmp,a.out/'trajectory.npz')
    for i in range(round(a.ps/a.stride_ps)+1):
        if i:
            sim.step(round(a.stride_ps/.002))
        state=sim.context.getState(positions=True,energy=True,enforcePeriodicBox=True)
        xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box=state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        energy=state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        if not np.isfinite(energy) or not np.isfinite(xyz).all() or min(np.diag(box))<=1.6: raise ValueError('Unstable frame')
        # Quantities are computed from exactly the coordinates written into the file.
        xyz=xyz.astype(np.float32);box=box.astype(np.float32)
        row=dict(time_ps=float(i*a.stride_ps),**metrics(xyz,box,top),density_g_ml=sum(top['mass_da'])/(602.214076*np.linalg.det(box.astype(float))),
                 potential_energy_kj_mol=float(energy))
        pos.append(xyz);boxes.append(box);times.append(row['time_ps']);table.append(row)
        with (a.out/'metrics.csv').open('a',newline='') as file:
            w=csv.DictWriter(file,fieldnames=list(row))
            if not i: w.writeheader()
            w.writerow(row)
        if i%10==0 or times[-1]==a.ps:
            save_track()
            print(a.seed,'free',times[-1],'ps, contacts',row['heavy_pair_contacts_045'],'shared waters',row['shared_water_035'],
                  'elapsed s',round(time.monotonic()-start),flush=True)
    manifest.update(complete=True,frames=len(times),wall_s=round(time.monotonic()-start,2),
                    files_sha256={s:checksum(a.out/s) for s in ('topology.json','trajectory.npz','metrics.csv')})
    save_json(a.out/'protocol.json',manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seed',type=int,required=True);p.add_argument('--waters',type=int,default=700)
    p.add_argument('--condition-ps',type=float,default=8.);p.add_argument('--ps',type=float,default=60.)
    p.add_argument('--stride-ps',type=float,default=1.)
    a=p.parse_args()
    if (not 1<=a.seed<2**31-20 or not 250<=a.waters<=10000 or not a.ps>0 or not a.condition_ps>0 or not a.stride_ps>0
            or abs(round(a.ps/a.stride_ps)*a.stride_ps-a.ps)>1e-8 or abs(round(a.stride_ps/.002)*.002-a.stride_ps)>1e-8):
        p.error('Incompatible duration, stride, size or seed')
    sys.path.insert(0,str(a.model.resolve()))
    run(a)

if __name__=='__main__':
    main()
