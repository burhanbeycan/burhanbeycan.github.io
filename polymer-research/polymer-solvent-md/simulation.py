"""Run one versioned, explicit-TIP3P, atomistic polymer solution replica."""
from contextlib import contextmanager
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import time

import numpy as np

SCHEMA = "polymer-solution-md/1"
CUTOFF_NM = 0.8
STEP_PS = 0.002
PILOT_SEEDS = (13031, 27011, 49033)
OFFSETS = np.asarray([(a,b,c) for a in (-1,0,1) for b in (-1,0,1) for c in (-1,0,1)
                      if (a,b,c) != (0,0,0)], dtype=float)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(obj):
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()


def file_hash(path):
    return digest(path.read_bytes())


def atomic_json(path, obj):
    raw = json_bytes(obj)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
        tmp.write(raw)
    os.replace(tmp.name, path)


@contextmanager
def _replica_lock(path):
    # The kernel releases the lock on exit, even after a process crash.
    f = path.open("a+b")
    try:
        if os.name == "nt":
            import msvcrt
            if f.tell() == 0:
                f.write(b"0"); f.flush()
            f.seek(0); msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        f.close()
        raise ValueError("Aynı koşu kimliği zaten hesaplanıyor.") from exc
    try:
        yield
    finally:
        f.close()


def _mol_graph(elements, bonds):
    from rdkit import Chem
    g = Chem.RWMol()
    for name in elements:
        atom=Chem.Atom(name); atom.SetNoImplicit(True); g.AddAtom(atom)
    for a,b in bonds:
        g.AddBond(int(a),int(b),Chem.BondType.SINGLE)
    return g.GetMol()


def _atom_map(source_elem, source_bonds, target_top):
    # Map only equivalent full connectivity graphs. Symmetry-related H atoms can permute.
    native_atoms=list(target_top.atoms())
    native_elem=[a.element.symbol for a in native_atoms]
    native_bonds=[(a.index,b.index) for a,b in target_top.bonds()]
    if len(source_elem)!=len(native_elem) or len(source_bonds)!=len(native_bonds):
        raise ValueError("Parametre paketi ile başlangıç geometrisinin atom veya bağ sayıları eşleşmiyor.")
    def edges(bs):
        return {tuple(sorted(pair)) for pair in bs}
    if source_elem == native_elem and edges(source_bonds)==edges(native_bonds):
        return list(range(len(source_elem)))
    source=_mol_graph(source_elem, source_bonds)
    order=source.GetSubstructMatch(_mol_graph(native_elem,native_bonds))
    if len(order)!=len(native_elem) or edges((order[i],order[j]) for i,j in native_bonds)!=edges(source_bonds):
        raise ValueError("Parametre paketi ile başlangıç geometrisinin bağ grafı uyuşmuyor.")
    return list(order)


def _starting_geometry(canonical, top, seed, record):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles(canonical))
    smi_elements=[a.GetSymbol() for a in mol.GetAtoms()]
    smi_bonds=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()]
    order=_atom_map(smi_elements,smi_bonds,top)
    source_to_native = {j:i for i,j in enumerate(order)}
    q=Chem.MolFromSmarts("[O;X1]=[C;X3]-[#7]")
    amide=sorted({source_to_native[a] for a,_,_ in mol.GetSubstructMatches(q)})
    if not amide:
        raise ValueError("Bu analiz için zincirde amid karbonil oksijeni bulunamadı.")
    if record is not None:
        source=record["structure"]["canonical_smiles"]
        if Chem.MolToSmiles(Chem.MolFromSmiles(source),isomericSmiles=True)!=canonical:
            raise ValueError("Başlangıç konformeri ve kuvvet alanı farklı yapı kimliklerine ait.")
        geom=record["energy"].get("lowest_geometry")
        if not geom or geom.get("unit")!="angstrom":
            raise ValueError("Kayıt, angström biriminde doğrulanmış başlangıç geometrisi içermiyor.")
        elements=[a["element"] for a in geom["atoms"]]
        bonds=[(b["a"],b["b"]) for b in geom["bonds"]]
        stored_order=_atom_map(elements,bonds,top)
        xyz=np.asarray([[a[k] for k in "xyz"] for a in geom["atoms"]],float)[stored_order]*0.1
        tag={"method":"previous module accepted isolated conformer", "artifact_id":record.get("artifact_id"),
             "model_id":record["energy"]["model_id"]}
    else:
        emb=AllChem.ETKDGv3();emb.randomSeed=seed;emb.enforceChirality=True;emb.numThreads=1
        if AllChem.EmbedMolecule(mol,emb)!=0 or not AllChem.MMFFHasAllMoleculeParams(mol):
            raise ValueError("Bu geometri tohumu için desteklenen başlangıç konformeri oluşturulamadı.")
        if AllChem.MMFFOptimizeMolecule(mol,mmffVariant="MMFF94s",maxIters=1200)!=0:
            raise ValueError("Başlangıç konformerinin MMFF minimizasyonu yakınsamadı.")
        xyz=mol.GetConformer().GetPositions()[order,:]*0.1
        tag={"method":"independent ETKDGv3 + MMFF94s starting geometry", "conformer_seed":seed}
    if not np.isfinite(xyz).all():
        raise ValueError("Sonlu koordinatlar gereklidir.")
    return xyz-xyz.mean(axis=0),amide,tag


def _self_image(whole,box):
    # For the compact registered oligomer all nearest other copies lie among 26 neighbors.
    shifts=OFFSETS @ box
    delta=whole[:,None,:]-whole[None,:,:]
    return math.sqrt(min(float(np.min(np.sum((delta-shift)**2,axis=2))) for shift in shifts))


def _integer_steps(length, step):
    n=round(length/step)
    if n<=0 or abs(n*step-length)>1e-8:
        raise ValueError("Zaman uzunluğu örnekleme adımının pozitif bir tam katı olmalıdır.")
    return n


def _solvent_count(reference_xyz,override):
    diam=float(np.max(np.linalg.norm(reference_xyz[:,None,:]-reference_xyz[None,:,:],axis=2)))
    nominal=max(500,math.ceil(33.4*(diam+CUTOFF_NM+0.32)**3))
    count=nominal if override is None else override
    if type(count) is not int or not 250<=count<=2200:
        raise ValueError(f"Önerilen kutu {count} su içeriyor; bu yerel protokol 250–2200 suyu destekler.")
    return count


def _protocol_checks(temperature,pressure,equil,prod,stride):
    if not all(math.isfinite(x) for x in (temperature,pressure,equil,prod,stride)):
        raise ValueError("Koşullar sonlu olmalıdır.")
    if not 273.15<=temperature<=350 or not 0.5<=pressure<=10:
        raise ValueError("Başlangıçtaki sıvı-su protokolü 273,15–350 K ve 0,5–10 bar aralığındadır.")
    if not 0<equil<=2000 or not 0<prod<=10000 or not 0.02<=stride<=5:
        raise ValueError("Örnekleme aralığı veya koşu süresi protokol dışında.")
    _integer_steps(equil,STEP_PS);_integer_steps(prod,stride);_integer_steps(stride,STEP_PS)


def run_replicate(packet_dir,out_dir,*,seed,temperature_K=298.15,pressure_bar=1.0,n_water=None,
                  equil_ps=12.,production_ps=24.,sample_ps=0.2,geometry_record=None):
    """Compute one independent NPT run. Repeated, identical requests return verified data."""
    from openmm import app, unit
    import openmm as mm
    from rdkit import Chem, rdBase
    from amber_bridge import placement_forcefield, create_solution_system
    import amber_bridge
    from importlib.metadata import version

    if type(seed) is not int or not 1<=seed<2**31-3:
        raise ValueError("Tohum pozitif bir 32-bit tam sayı olmalıdır.")
    _protocol_checks(temperature_K,pressure_bar,equil_ps,production_ps,sample_ps)
    packet=Path(packet_dir); out=Path(out_dir); out.mkdir(parents=True,exist_ok=True)
    md_top=packet/"solute.prmtop";md_crd=packet/"solute.inpcrd"
    sdf=packet/"solute.sdf"
    native_mol=Chem.MolFromMolFile(str(sdf),removeHs=True,sanitize=True)
    if native_mol is None:
        raise ValueError("Parametre paketindeki yapısal kimlik okunamadı.")
    canonical=Chem.MolToSmiles(native_mol,isomericSmiles=True)
    prmtop=app.AmberPrmtopFile(str(md_top))
    template_xyz=np.asarray(app.AmberInpcrdFile(str(md_crd)).positions.value_in_unit(unit.nanometer))
    count=_solvent_count(template_xyz,n_water)
    xyz,amide,init=_starting_geometry(canonical,prmtop.topology,seed,geometry_record)
    sig={"schema":SCHEMA, "smiles":canonical, "input_structure_sha256":file_hash(sdf),
         "topology_sha256":file_hash(md_top), "reference_coordinates_sha256":file_hash(md_crd),
         "runner_sha256":file_hash(Path(__file__)), "bridge_sha256":file_hash(Path(amber_bridge.__file__)),
         "water_parameters_sha256":file_hash(Path(app.__file__).parent/"data"/"tip3p.xml"),
         "initial_coordinates_sha256":digest(np.ascontiguousarray(xyz,dtype="<f8").tobytes()),
         "initial_geometry":init, "engine":"OpenMM", "engine_version":version("openmm"),
         "engine_build":mm.version.version, "rdkit_version":rdBase.rdkitVersion,
         "water_model":"TIP3P", "solute_model":"GAFF2 + AM1-BCC",
         "temperature_K":float(temperature_K),"pressure_bar":float(pressure_bar),"n_water":count,
         "equil_ps":float(equil_ps), "production_ps":float(production_ps), "sample_ps":float(sample_ps),
         "seed":seed, "velocity_seed":seed, "thermostat_seed":seed+1,"barostat_seed":seed+2,
         "timestep_ps":STEP_PS,"cutoff_nm":CUTOFF_NM, "pme_relative_tolerance":0.0005,
         "constraint_tolerance":1e-5,"solute_constraints":"HBonds","rigid_water":True,
         "platform":"CPU","threads":1, "barostat_interval_steps":25,"friction_per_ps":1.0,
         "minimization_steps":300,"minimization_tolerance_kj_mol_nm":10,
         "ensemble":"NPT","thermostat":"LangevinMiddleIntegrator","barostat":"MonteCarloBarostat"}
    request_id="sha256:"+digest(json_bytes(sig))
    npfile=out/f"seed_{seed}.npz";response_file=out/f"seed_{seed}.json";tf=out/"topology.json"
    with _replica_lock(out/f"seed_{seed}.lock"):
        if response_file.exists():
            old=json.loads(response_file.read_text("utf-8"))
            if old["run_id"]!=request_id:
                raise ValueError("Bu çıktı konumunda aynı tohumla farklı bir protokol var. Başka bir koşu dizini kullanın.")
            if npfile.is_file() and tf.is_file() and file_hash(npfile)==old["trajectory_sha256"] and file_hash(tf)==old["shared_topology_sha256"]:
                return old
            raise ValueError("Kayıtlı yörüngenin veya topolojinin dosya özeti uyuşmuyor.")
        start=time.monotonic()
        model=app.Modeller(prmtop.topology,xyz*unit.nanometer)
        model.addSolvent(placement_forcefield(prmtop),model="tip3p",numAdded=count,neutralize=False)
        system=create_solution_system(prmtop,model.topology,nonbondedMethod=app.PME,
                                      nonbondedCutoff=CUTOFF_NM*unit.nanometer, constraints=app.HBonds,
                                      rigidWater=True,ewaldErrorTolerance=0.0005)
        atoms=list(model.topology.atoms());ns=len(xyz); nat=len(atoms)
        waters=[r for r in model.topology.residues() if list(r.atoms())[0].index>=ns]
        if len(waters)!=count or nat!=ns+3*count:
            raise ValueError("Su sayısı veya çözücü atom sırası tutarsız.")
        wat_o=[]
        for r in waters:
            members=list(r.atoms())
            if [a.element.symbol for a in members]!=["O","H","H"]:
                raise ValueError("Beklenmeyen TIP3P su atom dizilimi.")
            wat_o.append(members[0].index)
        masses=[system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(nat)]
        topo={"schema":"solution-topology/1","canonical_smiles":canonical,"elements":[a.element.symbol for a in atoms],
              "atom_names":[a.name for a in atoms],"bonds":[[a.index,b.index] for a,b in model.topology.bonds()],
              "mass_da":masses,"solute_indices":list(range(ns)),"carbonyl_oxygens":amide,"water_oxygens":wat_o,
              "coordinate_unit":"nm","solvent":"TIP3P","n_water":count}
        with _replica_lock_wait(out/"topology.lock"):
            if tf.exists() and tf.read_bytes()!=json_bytes(topo):
                raise ValueError("Bu klasördeki topoloji, koşunun atom sırası veya bileşimi ile uyuşmuyor.")
            atomic_json(tf,topo)
        bar=mm.MonteCarloBarostat(pressure_bar*unit.bar,temperature_K*unit.kelvin,25);bar.setRandomNumberSeed(seed+2)
        system.addForce(bar)
        integrator=mm.LangevinMiddleIntegrator(temperature_K*unit.kelvin,1/unit.picosecond,STEP_PS*unit.picoseconds)
        integrator.setRandomNumberSeed(seed+1);integrator.setConstraintTolerance(1e-5)
        sim=app.Simulation(model.topology,system,integrator,mm.Platform.getPlatformByName("CPU"),
                           {"Threads":"1","DeterministicForces":"true"})
        sim.context.setPositions(model.positions)
        box0=np.array(model.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
        first_clearance=_self_image(xyz,box0)
        if first_clearance<CUTOFF_NM+0.05 or min(np.diag(box0))<2*CUTOFF_NM+0.15:
            raise ValueError("Başlangıçta periyodik görüntüye uzaklık yetersiz. Daha büyük su kutusu gerekir.")
        sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=300)
        sim.context.setVelocitiesToTemperature(temperature_K*unit.kelvin,seed)
        sim.step(_integer_steps(equil_ps,STEP_PS))
        print(f"{seed}: equilibration complete ({time.monotonic()-start:.1f} s)",flush=True)
        frames=_integer_steps(production_ps,sample_ps)
        coords=np.empty((frames,nat,3),dtype=np.float32); boxes=np.empty((frames,3,3),dtype=np.float32)
        etotal=np.empty(frames);densities=np.empty(frames);temperatures=np.empty(frames);distances=np.empty(frames)
        dt_steps=_integer_steps(sample_ps,STEP_PS);dof=3*nat-system.getNumConstraints()
        if any(isinstance(f,mm.CMMotionRemover) for f in system.getForces()):
            dof-=3
        for i in range(frames):
            sim.step(dt_steps)
            state=sim.context.getState(positions=True,energy=True,enforcePeriodicBox=True)
            pos=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            box=state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
            vol=np.linalg.det(box)
            if not np.isfinite(pos).all() or np.max(abs(box-np.diag(np.diag(box))))>1e-8 or min(np.diag(box))<=2*CUTOFF_NM:
                raise RuntimeError("Geçersiz periyodik kutu veya sonlu olmayan koordinat.")
            coords[i]=pos;boxes[i]=box
            etotal[i]=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            densities[i]=sum(masses)/(602.214076*vol)
            temperatures[i]=2*state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(0.00831446261815324*dof)
            distances[i]=_self_image(pos[:ns],box)
        if not all(np.isfinite(a).all() for a in (etotal,densities,temperatures,distances)):
            raise RuntimeError("Üretim sırasında sonlu olmayan bir fiziksel büyüklük saptandı.")
        times=np.arange(1,frames+1)*sample_ps
        with tempfile.NamedTemporaryFile(dir=out,delete=False) as tmp:
            np.savez_compressed(tmp,positions_nm=coords,box_nm=boxes,time_ps=times,
                                density_g_ml=densities,temperature_K=temperatures,cell_potential_kj_mol=etotal,
                                image_clearance_nm=distances)
        os.replace(tmp.name,npfile)
        concentrations=1/(np.linalg.det(boxes.astype(float))*0.602214076)
        result={"schema":SCHEMA,"run_id":request_id,"parameters":sig,"frames":frames,"solute_atoms":ns,"water_count":count,
                "trajectory":npfile.name,"trajectory_sha256":file_hash(npfile),"shared_topology_sha256":file_hash(tf),
                "initial_self_image_distance_nm":first_clearance, "minimum_self_image_distance_nm":float(distances.min()),
                "sample_frames_with_self_image_within_cutoff":int(np.sum(distances<CUTOFF_NM)),
                "mean_density_g_ml":float(densities.mean()),"mean_temperature_K":float(temperatures.mean()),
                "mean_chain_molarity_mol_l":float(concentrations.mean()),"solute_mass_fraction":sum(masses[:ns])/sum(masses),
                "wall_seconds":round(time.monotonic()-start,2),
                "limits":"Finite periodic box. Short production pilot. Sampling and size dependence require additional runs."}
        atomic_json(response_file,result)
        print(f"{seed}: {frames} production frames, {result['wall_seconds']} s",flush=True)
        return result


@contextmanager
def _replica_lock_wait(path):
    started=time.monotonic()
    while True:
        guard=_replica_lock(path)
        try:
            guard.__enter__();break
        except ValueError:
            if time.monotonic()-started>3:
                raise
            time.sleep(.05)
    try:
        yield
    finally:
        guard.__exit__(None,None,None)


def cli():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parameters",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--seed",type=int,required=True);p.add_argument("--temperature",type=float,default=298.15)
    p.add_argument("--pressure",type=float,default=1.);p.add_argument("--waters",type=int)
    p.add_argument("--equil-ps",type=float,default=12.);p.add_argument("--production-ps",type=float,default=24.)
    p.add_argument("--sample-ps",type=float,default=.2)
    p.add_argument("--initial-record",type=Path)
    a=p.parse_args()
    rec=json.loads(a.initial_record.read_text()) if a.initial_record else None
    run_replicate(a.parameters,a.out,seed=a.seed,temperature_K=a.temperature,pressure_bar=a.pressure,
                  n_water=a.waters,equil_ps=a.equil_ps,production_ps=a.production_ps,sample_ps=a.sample_ps,geometry_record=rec)


if __name__=="__main__":
    os.environ.setdefault("OPENMM_CPU_THREADS","1")
    cli()
