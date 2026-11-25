"""Extend a validated one-chain TIP3P trajectory from its final sampled coordinates.

The source contains positions and box but no velocities. A new Maxwell velocity
sample and a discarded NPT segment precede production; time_ps starts anew.
Each completed chunk also holds an engine checkpoint to support interrupted runs.
"""
import argparse
import gzip
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import numpy as np
import openmm as mm
from openmm import app, unit as u


def hash_file(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def save_bytes(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as file:
        file.write(value)
    os.replace(file.name, path)


def save_npz(path, **fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as file:
        np.savez_compressed(file, **fields)
    os.replace(file.name, path)


def load_source(source, params, seed):
    md = json.loads((source/f"seed_{seed}.json").read_bytes())
    traj = source/f"seed_{seed}.npz"
    top = source/"topology.json"
    checks = [(traj, md["trajectory_sha256"]), (top, md["shared_topology_sha256"]),
              (params/"solute.prmtop", md["parameters"]["topology_sha256"]),
              (params/"solute.inpcrd", md["parameters"]["reference_coordinates_sha256"])]
    for path, expected in checks:
        if hash_file(path) != expected:
            raise ValueError(f"Kaynak dosya özeti eşleşmiyor: {path.name}")
    if md["parameters"]["seed"] != seed:
        raise ValueError("Tohum ve kaynak kaydı uyuşmuyor.")
    topo = json.loads(top.read_bytes())
    with np.load(traj,allow_pickle=False) as saved:
        pos = saved["positions_nm"][-1].astype(float)
        box = saved["box_nm"][-1].astype(float)
        clock = float(saved["time_ps"][-1])
    if (pos.shape != (len(topo["elements"]),3) or box.shape != (3,3) or not np.isfinite(clock) or clock <= 0
            or not np.isfinite(pos).all() or not np.isfinite(box).all() or len(topo["elements"]) != len(pos)
            or not np.allclose(box,np.diag(np.diag(box))) or min(np.diag(box)) <= 2*md["parameters"]["cutoff_nm"]):
        raise ValueError("Geçersiz periyodik kutu veya koordinat kaydı.")
    return md, topo, pos, box, clock


def rebuild_topology(prmtop, topo, box):
    """Recreate named native atoms and stored TIP3P residues without placing new water."""
    full = prmtop.topology
    atoms = list(full.atoms())
    n = len(atoms)
    if topo["solute_indices"] != list(range(n)):
        raise ValueError("Çözünen atom sırası uyuşmuyor.")
    water_chain = full.addChain()
    for oxygen in topo["water_oxygens"]:
        if oxygen != len(atoms) or topo["elements"][oxygen:oxygen+3] != ["O","H","H"]:
            raise ValueError("Kaynak üç noktalı TIP3P su sırası bekleniyor.")
        res = full.addResidue("HOH", water_chain)
        for index in range(oxygen,oxygen+3):
            element = app.element.get_by_symbol(topo["elements"][index])
            atoms.append(full.addAtom(topo["atom_names"][index],element,res))
        full.addBond(atoms[oxygen],atoms[oxygen+1]); full.addBond(atoms[oxygen],atoms[oxygen+2])
    saved_edges = {tuple(sorted(pair)) for pair in topo["bonds"]}
    rebuilt_edges = {tuple(sorted((a.index,b.index))) for a,b in full.bonds()}
    if (saved_edges != rebuilt_edges or [a.name for a in atoms] != topo["atom_names"]
            or [a.element.symbol for a in atoms] != topo["elements"]):
        raise ValueError("Kaydedilmiş tüm-sistem atom kimliği veya bağları eşleşmiyor.")
    full.setPeriodicBoxVectors(tuple(mm.Vec3(*v) for v in box)*u.nanometer)
    return full


def steps(duration, dt):
    value = int(round(duration/dt))
    if duration <= 0 or value < 1 or not np.isclose(value*dt,duration,rtol=0,atol=1e-9):
        raise ValueError("Zaman uzunluğu, adımın pozitif tam katı olmalıdır.")
    return value


def seed_stream(seed):
    rng = np.random.SeedSequence([30971, seed]).generate_state(3)
    return [1+int(i)%(2**31-2) for i in rng]


def bind_system(base, params, topo, box, parent, seeds):
    sys.path.insert(0,str(base.resolve()))
    from amber_bridge import create_solution_system
    imported = app.AmberPrmtopFile(str(params/"solute.prmtop"))
    # Keep the ligand-only AmberPrmtopFile topology intact inside the native bridge.
    separate = app.AmberPrmtopFile(str(params/"solute.prmtop"))
    full = rebuild_topology(separate,topo,box)
    pp = parent["parameters"]
    system = create_solution_system(imported,full,nonbondedMethod=app.PME,nonbondedCutoff=pp["cutoff_nm"]*u.nanometer,
                                     constraints=app.HBonds,rigidWater=True,ewaldErrorTolerance=pp["pme_relative_tolerance"])
    masses = [system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())]
    if not np.allclose(masses,topo["mass_da"],rtol=0,atol=1e-7):
        raise ValueError("Kütleler önceki sistemle uyuşmuyor.")
    T=pp["temperature_K"]*u.kelvin
    pressure=pp["pressure_bar"]*u.bar
    bar=mm.MonteCarloBarostat(pressure,T,pp["barostat_interval_steps"])
    bar.setRandomNumberSeed(seeds[2])
    system.addForce(bar)
    integ=mm.LangevinMiddleIntegrator(T,pp["friction_per_ps"]/u.picosecond,pp["timestep_ps"]*u.picoseconds)
    integ.setConstraintTolerance(pp["constraint_tolerance"]); integ.setRandomNumberSeed(seeds[1])
    return full,system,integ


def ready_chunks(check_dir, identity):
    result = []
    for idx,f in enumerate(sorted(check_dir.glob("*.npz"))):
        if f.name != f"{idx:05d}.npz":
            raise ValueError("Koşu parçalarının sırası bozuldu.")
        with np.load(f,allow_pickle=False) as data:
            if str(data["run_id"].item()) != identity:
                raise ValueError("Önceki parçadaki protokol kimliği uyuşmuyor.")
        result.append(f)
    return result


def run(base,params,source,out,seed,production_ps=120.,discard_ps=12.,stride_ps=1.,chunk_frames=12):
    sys.path.insert(0,str(base.resolve()))
    from simulation import _replica_lock, _self_image

    out.mkdir(parents=True,exist_ok=True)
    parent,topo,pos,box,clock=load_source(source,params,seed)
    pp=parent["parameters"]; seeds=seed_stream(seed); T=pp["temperature_K"]*u.kelvin
    dt=pp["timestep_ps"];n_frames=steps(production_ps,stride_ps)
    stride_steps=steps(stride_ps,dt);drop_steps=steps(discard_ps,dt)
    if type(chunk_frames) is not int or not 1<=chunk_frames<=100:
        raise ValueError("Koşu parçası kare sayısı 1–100 olmalıdır.")
    runner=hash_file(Path(__file__)); bridge=hash_file(base/"amber_bridge.py")
    water_sha=hash_file(Path(app.__file__).parent/"data"/"tip3p.xml")
    if water_sha != pp["water_parameters_sha256"]:
        raise ValueError("Kurulu TIP3P parametreleri kaynak koşuyla uyuşmuyor.")
    if bridge != pp["bridge_sha256"]:
        raise ValueError("Kaynak koşu ile yeni kuvvet alanı bağlantı kodu uyuşmuyor.")
    sig={**pp,"schema":"polymer-hydration-extension/1","runner_sha256":runner,
         "engine_version":metadata.version("openmm"),"engine_build":mm.version.version,
         "velocity_seed":seeds[0],"thermostat_seed":seeds[1],"barostat_seed":seeds[2],
         "equil_ps":float(discard_ps),"production_ps":float(production_ps),"sample_ps":float(stride_ps),
         "minimization_steps":0,"minimization_tolerance_kj_mol_nm":None,
         "initial_coordinates_sha256":hashlib.sha256(np.ascontiguousarray(pos,dtype="<f8").tobytes()).hexdigest(),
         "initial_geometry":{"method":"parent final saved box and positions; fresh constrained Maxwell velocities",
            "parent_run_id":parent["run_id"],"parent_seed":seed,"source_production_time_ps":clock,
            "source_trajectory_sha256":parent["trajectory_sha256"],"source_metadata_sha256":hash_file(source/f"seed_{seed}.json")},
         "checkpoint_chunk_frames":chunk_frames}
    identity="sha256:"+hashlib.sha256(canonical(sig)).hexdigest()
    topfile=out/"topology.json"
    if topfile.exists() and hash_file(topfile)!=parent["shared_topology_sha256"]:
        raise ValueError("Çıktı dizini başka bileşime ait.")
    save_bytes(topfile,(source/"topology.json").read_bytes())
    final=out/f"seed_{seed}.json"; traj=out/f"seed_{seed}.npz"
    with _replica_lock(out/f"seed_{seed}.lock"):
        if final.is_file():
            old=json.loads(final.read_bytes())
            if old["run_id"]!=identity or hash_file(traj)!=old["trajectory_sha256"]:
                raise ValueError("Çıktı klasörü farklı ya da tutarsız bir hesap içeriyor.")
            return old
        start=time.monotonic()
        top,system,integrator=bind_system(base,params,topo,box,parent,seeds)
        sim=app.Simulation(top,system,integrator,mm.Platform.getPlatformByName("CPU"),{"Threads":"1","DeterministicForces":"true"})
        c=sim.context
        checkdir=out/"checkpoints"/f"seed_{seed}"
        checkdir.mkdir(parents=True,exist_ok=True)
        chunks=ready_chunks(checkdir,identity)
        if chunks:
            with np.load(chunks[-1],allow_pickle=False) as data:
                c.loadCheckpoint(data["checkpoint"].tobytes())
                completed=int(data["stop_frame"].item())
        else:
            c.setPeriodicBoxVectors(*(mm.Vec3(*v)*u.nanometer for v in box))
            c.setPositions(pos*u.nanometer)
            c.applyConstraints(pp["constraint_tolerance"])
            c.setVelocitiesToTemperature(T,seeds[0]);c.applyVelocityConstraints(pp["constraint_tolerance"])
            sim.step(drop_steps)
            completed=0
            print(f"seed {seed}: rethermalization {discard_ps:g} ps complete in {time.monotonic()-start:.1f} s",flush=True)
        mass=np.sum(topo["mass_da"])
        ns=len(topo["solute_indices"])
        dof=3*len(pos)-system.getNumConstraints()-sum(isinstance(f,mm.CMMotionRemover) for f in system.getForces())*3
        fields=("positions_nm","box_nm","time_ps","density_g_ml","temperature_K","cell_potential_kj_mol","image_clearance_nm")
        while completed<n_frames:
            end=min(completed+chunk_frames,n_frames)
            rows={k:[] for k in fields}
            for i in range(completed+1,end+1):
                sim.step(stride_steps)
                state=c.getState(positions=True,energy=True,enforcePeriodicBox=True)
                p=state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
                b=state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer)
                if not np.isfinite(p).all() or not np.allclose(b,np.diag(np.diag(b))) or min(np.diag(b))<=2*pp["cutoff_nm"]:
                    raise RuntimeError("Geçersiz kutu veya koordinatlar.")
                pe=state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                temp=2*state.getKineticEnergy().value_in_unit(u.kilojoule_per_mole)/(0.00831446261815324*dof)
                density=mass/(602.214076*np.linalg.det(b))
                dist=_self_image(p[:ns],b)
                if not np.isfinite([pe,temp,density,dist]).all():
                    raise RuntimeError("Sonlu olmayan enerji veya gözlem.")
                for k,v in zip(fields,(p.astype("f4"),b.astype("f4"),i*stride_ps,density,temp,pe,dist)):
                    rows[k].append(v)
            f=checkdir/f"{len(chunks):05d}.npz"
            save_npz(f,**rows,checkpoint=np.frombuffer(c.createCheckpoint(),np.uint8),run_id=identity,stop_frame=end)
            chunks.append(f);completed=end
            print(f"seed {seed}: {end}/{n_frames} frames in {time.monotonic()-start:.1f} s",flush=True)
        collected={k:[] for k in fields}
        for f in chunks:
            with np.load(f,allow_pickle=False) as data:
                for k in fields:collected[k].append(data[k])
        result_arrays={k:np.concatenate(v) for k,v in collected.items()}
        if not np.allclose(result_arrays["time_ps"],np.arange(1,n_frames+1)*stride_ps):
            raise ValueError("Çıktı zamanları protokolle uyuşmuyor.")
        save_npz(traj,**result_arrays)
        statefile=out/f"state_{seed}.xml.gz"
        state=mm.XmlSerializer.serialize(c.getState(positions=True,velocities=True,parameters=True))
        save_bytes(statefile,gzip.compress(state.encode(),mtime=0))
        vols=np.linalg.det(result_arrays["box_nm"].astype(float))
        result={"schema":sig["schema"],"run_id":identity,"parameters":sig,"frames":n_frames,
                "solute_atoms":ns,"water_count":len(topo["water_oxygens"]),"trajectory":traj.name,"trajectory_sha256":hash_file(traj),
                "shared_topology_sha256":hash_file(topfile),"initial_self_image_distance_nm":_self_image(pos[:ns],box),
                "minimum_self_image_distance_nm":float(result_arrays["image_clearance_nm"].min()),
                "sample_frames_with_self_image_within_cutoff":int((result_arrays["image_clearance_nm"]<pp["cutoff_nm"]).sum()),
                "mean_density_g_ml":float(result_arrays["density_g_ml"].mean()),"mean_temperature_K":float(result_arrays["temperature_K"].mean()),
                "mean_chain_molarity_mol_l":float(np.mean(1/(vols*.602214076))),"solute_mass_fraction":float(sum(topo["mass_da"][:ns])/mass),
                "final_portable_state":statefile.name,"final_portable_state_sha256":hash_file(statefile),
                "wall_seconds_this_invocation":round(time.monotonic()-start,2),
                "scope":"Warm positional restart. Fresh thermal velocities. Finite, explicit periodic cell with one chain."}
        save_bytes(final,canonical(result))
        return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for s in ("base","params","source","out"):p.add_argument("--"+s,type=Path,required=True)
    p.add_argument("--seed",type=int,required=True)
    p.add_argument("--production-ps",type=float,default=120);p.add_argument("--discard-ps",type=float,default=12)
    p.add_argument("--stride-ps",type=float,default=1);p.add_argument("--chunk-frames",type=int,default=12)
    args=p.parse_args();r=run(**vars(args));print(r["run_id"],flush=True)
