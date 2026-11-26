"""Reproduce the microsolvation pilot or analyze a matched explicit-water GROMACS run."""
import argparse
import csv
import hashlib
import json
import re
import io
import shutil
import subprocess
from collections import defaultdict
from importlib.resources import files
from pathlib import Path

import numpy as np

from .classical_md import frame_observables, petox_cluster, run_mmff

HEADER = ["sample","environment","dp","end","replica","frame","time_ps","rg_nm",
          "carbonyl_coordination","omega_coordination","cutoff_nm","water_count","force_field","water_model","temperature_K"]
PACKING_SEEDS=(109,211,353)


def save_csv(path, rows):
    with Path(path).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)


def run_pilot(directory, dp_values=(2,4), seeds=PACKING_SEEDS, *, steps=4000,
              thermalization=2000, stride=100, cutoff=0.35):
    from rdkit.Chem import rdMolDescriptors
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=False)
    rows=[]; records=[]
    for dp in dp_values:
        for end in ("OH","N3"):
            for seed in seeds:
                mol,group=petox_cluster(dp,end,seed)
                run=run_mmff(mol,seed=seed+100000,production_steps=steps,equil_steps=thermalization,stride=stride)
                mass=run["masses"][group["polymer"]]
                trial=f"petox_{dp}_{end}_{seed}"
                for i,(t,xyz) in enumerate(zip(run["time_ps"],run["positions_nm"])):
                    rg,n_co,n_end=frame_observables(xyz,group,mass,cutoff)
                    rows.append(dict(zip(HEADER,["isolated_12water_MMFF94","isolated microsolvated cluster",
                        dp,end,seed,i+1,round(t,6),round(rg,7),round(n_co,6),round(n_end,6),cutoff,12,"MMFF94","MMFF94 water",run["temperature_K"]])))
                z=run["positions_nm"][-1]
                with (directory/(trial+".xyz")).open("w",encoding="utf-8") as stream:
                    stream.write(f"{len(z)}\n{trial}, last production frame, angstrom, no periodic box\n")
                    for atom,point in zip(mol.GetAtoms(),z*10):
                        stream.write(f"{atom.GetSymbol()} {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")
                records.append(dict(id=trial,dp=dp,end=end,packing_seed=seed,velocity_seed=seed+100000,
                    formula_with_waters=rdMolDescriptors.CalcMolFormula(mol),n_sites=mol.GetNumAtoms(),
                    initial_energy_kcal_mol=run["initial_energy_kcal_mol"],
                    minimized_energy_kcal_mol=run["minimized_energy_kcal_mol"],
                    mean_total_kinetic_kcal_mol=float(np.mean(run["kinetic_kcal_mol"])),
                    mean_potential_kcal_mol=float(np.mean(run["potential_kcal_mol"]))))
    save_csv(directory/"petox_cluster_md.csv",rows)
    (directory/"pilot.json").write_text(json.dumps(dict(rdkit_version=run["rdkit_version"],
        numpy_version=np.__version__,force_field="RDKit MMFF94; all interfragment terms",
        integrator="BAOAB",temperature_K=run["temperature_K"],dt_fs=run["dt_fs"],
        friction_ps=run["friction_ps"],thermalization_steps=thermalization,
        production_steps=steps,stride=stride,cutoff_nm=cutoff,water_per_chain=12,
        results=records),indent=2)+"\n",encoding="utf-8")


def read_gro_frames(path):
    """Read fixed-width GROMACS multi-frame GRO coordinates in nm."""
    time_re=re.compile(r"\bt\s*=\s*([+\-0-9.eE]+)")
    with Path(path).open(encoding="utf-8") as stream:
        while title:=stream.readline():
            m=time_re.search(title)
            if not m:
                raise ValueError("GRO frame titles must include time t= in ps.")
            count=int(stream.readline())
            if not 1 <= count <= 2_000_000:
                raise ValueError("Invalid GRO atom count.")
            coords=[]; names=[]
            for _ in range(count):
                line=stream.readline()
                if len(line.rstrip("\n"))<44:
                    raise ValueError("Truncated GRO coordinate row.")
                # Standard GRO with three decimals. Refuse extended precision rather
                # than slicing a variable-width column at an incorrect location.
                if any(line[pos]!="." for pos in (24,32,40)):
                    raise ValueError("Use standard three-decimal GRO coordinates (nm).")
                names.append((line[5:10].strip(), line[10:15].strip()))
                coords.append([float(line[20:28]),float(line[28:36]),float(line[36:44])])
            box=np.fromstring(stream.readline(),sep=" ")
            if len(box)==9 and np.allclose(box[3:],0):
                box=box[:3]
            if len(box)!=3:
                raise ValueError("This analyzer supports orthorhombic boxes.")
            yield float(m.group(1)),names,np.asarray(coords),box


def analyze_gro(path, mapping, cutoff=0.35):
    required=("sample","source","force_field","water_model","temperature_K","dp","end","replica",
              "polymer","polymer_masses_Da","carbonyl","omega","water_oxygen","polymer_bonds",
              "carbonyl_names","omega_names","water_oxygen_names")
    if any(k not in mapping for k in required):
        raise ValueError("Mapping is incomplete; see the example and MD guide.")
    if any(not isinstance(mapping[k],str) or not mapping[k].strip()
           for k in ("sample","source","force_field","water_model","replica")):
        raise ValueError("Sample, provenance, force field, water model and independent replica are required.")
    if (type(mapping["dp"]) is not int or mapping["dp"]<1 or mapping["end"] not in ("OH","N3")
            or not np.isfinite(mapping["temperature_K"]) or mapping["temperature_K"]<=0):
        raise ValueError("Enter a positive DP, temperature and a supported terminal group.")
    poly=set(mapping["polymer"])
    if len(mapping["carbonyl"])!=mapping["dp"] or len(mapping["omega"])!=1:
        raise ValueError("Each repeat requires one carbonyl O; specify one distal omega site.")
    bonds=np.asarray(mapping["polymer_bonds"])
    if bonds.ndim!=2 or bonds.shape[1]!=2 or not np.issubdtype(bonds.dtype,np.integer) or not set(bonds.flat)<=poly:
        raise ValueError("Supply covalent polymer bonds as pairs of zero-based atom indices.")
    graph=defaultdict(set)
    for a,b in bonds:
        graph[a].add(b); graph[b].add(a)
    todo=[next(iter(poly))]; seen=set()
    while todo:
        now=todo.pop()
        if now not in seen:
            seen.add(now); todo.extend(graph[now]-seen)
    if seen!=poly:
        raise ValueError("Polymer atom selection must form one connected all-atom graph.")
    last_time=-np.inf; rows=[]; first_names=None
    for ix,(t,names,xyz,box) in enumerate(read_gro_frames(path)):
        if not np.isfinite(t) or t<=last_time:
            raise ValueError("Frame times must be finite and strictly increasing.")
        last_time=t
        if first_names is None:
            first_names=names
            for tag in ("omega","carbonyl","water_oxygen"):
                allowed=set(mapping[tag+"_names"])
                if not allowed or any(not isinstance(i,int) or not 0<=i<len(names) or
                                      ":".join(names[i]) not in allowed for i in mapping[tag]):
                    raise ValueError(f"Atom/residue name or index mismatch for {tag}.")
        elif names!=first_names:
            raise ValueError("Atom ordering or identity changed between frames.")
        obs=frame_observables(xyz,mapping,mapping["polymer_masses_Da"],cutoff,box)
        if np.max(np.linalg.norm(xyz[bonds[:,0]]-xyz[bonds[:,1]],axis=1))>0.25:
            raise ValueError("A polymer bond exceeds 0.25 nm. Make the molecule whole and verify the mapping.")
        rows.append(dict(zip(HEADER,[mapping["sample"],"periodic explicit water",
                mapping["dp"],mapping["end"],mapping["replica"],ix+1,t,*obs,cutoff,len(mapping["water_oxygen"]),
                mapping["force_field"],mapping["water_model"],mapping["temperature_K"]])))
    if len(rows)<2:
        raise ValueError("At least two trajectory frames are required.")
    return rows


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    choices=ap.add_subparsers(dest="action",required=True)
    pilot=choices.add_parser("pilot",help="Repeat all twelve finite-cluster trajectories")
    pilot.add_argument("directory",type=Path)
    for action in ("gro","gromacs"):
        cmd=choices.add_parser(action)
        cmd.add_argument("input",type=Path,help="GRO trajectory or a prepared production TPR")
        cmd.add_argument("mapping",type=Path)
        cmd.add_argument("directory",type=Path,help="New output directory")
        cmd.add_argument("--cutoff-nm",type=float,default=0.35)
        if action=="gromacs": cmd.add_argument("--discard-ps",type=float,required=True)
    args=ap.parse_args()
    if args.action=="pilot":
        run_pilot(args.directory)
        return
    input_file=args.input.resolve(strict=True)
    mapping=json.loads(args.mapping.read_text())
    args.directory.mkdir(parents=True,exist_ok=False)
    gro=input_file
    if args.action=="gromacs":
        engine=shutil.which("gmx")
        if engine is None: raise RuntimeError("Install GROMACS and provide gmx on PATH.")
        if not np.isfinite(args.discard_ps) or args.discard_ps<0: raise ValueError("Invalid first retained time.")
        subprocess.run([engine,"mdrun","-s",str(input_file),"-deffnm","production"],cwd=args.directory,check=True)
        candidates=[args.directory/"production.xtc",args.directory/"production.trr"]
        traj=next((f for f in candidates if f.is_file()),None)
        if not traj: raise RuntimeError("Enable coordinate trajectory output in the prepared TPR.")
        gro=(args.directory/"whole.gro").resolve()
        subprocess.run([engine,"trjconv","-s",str(input_file),"-f",str(traj.resolve()),"-o",str(gro),
                        "-pbc","whole","-b",str(args.discard_ps)],input="System\n",text=True,cwd=args.directory,check=True)
    rows=analyze_gro(gro,mapping,args.cutoff_nm)
    save_csv(args.directory/"observables.csv",rows)
    (args.directory/"run.json").write_text(json.dumps(dict(mapping=mapping,cutoff_nm=args.cutoff_nm,
        input_sha256=sha256(input_file),gro_sha256=sha256(gro),command=args.action,
        discard_ps=getattr(args,"discard_ps",None),frames=len(rows)),indent=2)+"\n")


def md_summary(text, target_dp):
    if not text.strip():
        text=files("polymer_chain_lab").joinpath("data").joinpath("petox_cluster_md.csv").read_text(encoding="utf-8")
    reader=csv.DictReader(io.StringIO(text))
    if reader.fieldnames!=HEADER:
        raise ValueError("Use an unmodified polymer-md CSV header; see the MD guide.")
    conditions=set(); runs=defaultdict(list); total=0; sample=""
    for row in reader:
        total += 1
        if total>20000:
            raise ValueError("Session supports a maximum of 20,000 trajectory summary rows.")
        if None in row or any(not isinstance(v,str) or not v.strip() for v in row.values()):
            raise ValueError("Missing value or unexpected column in trajectory CSV.")
        dp=int(row["dp"]); frame=int(row["frame"]); nwater=int(row["water_count"])
        if dp<1 or frame<1 or nwater<1 or row["end"] not in ("OH","N3"):
            raise ValueError("Invalid polymer group, frame number or water count.")
        if dp!=target_dp: continue
        t,rg,co,omega,cutoff,temp=map(float,[row[k] for k in
            ("time_ps","rg_nm","carbonyl_coordination","omega_coordination","cutoff_nm","temperature_K")])
        if (not np.all(np.isfinite([t,rg,co,omega,cutoff,temp])) or min(t,rg,co,omega)<0
                or not 0<cutoff<=1 or temp<=0 or max(co,omega)>nwater):
            raise ValueError("Invalid dimensions or contact counts.")
        context=row["environment"].strip()
        if context not in ("periodic explicit water","isolated microsolvated cluster"):
            raise ValueError("Unrecognized environment.")
        sample=row["sample"].strip()
        conditions.add((context,sample,nwater,cutoff,row["force_field"],row["water_model"],temp))
        runs[row["end"],row["replica"].strip()].append((frame,t,rg,co,omega))
    if not runs:
        raise ValueError("No observations for the selected DP.")
    if len(conditions)!=1:
        raise ValueError("Compare one matched sample series, water count/model, environment, temperature, force field and cutoff per analysis.")
    states={}
    for (end,identifier),vals in sorted(runs.items()):
        data=np.array(vals)
        if len(data)<2 or np.any(np.diff(data[:,0])<=0) or np.any(np.diff(data[:,1])<=0):
            raise ValueError("Replicas need at least two ordered, unique frames and times.")
        if not np.allclose(np.diff(data[:,1]),np.diff(data[:,1])[0],rtol=1e-4,atol=1e-6):
            raise ValueError("Use uniformly spaced saved times for time averages.")
        states.setdefault(end,[]).append((identifier,data))
    return states,conditions.pop()


def md_output(values):
    dp=values["pilot_dp"]
    states,condition=md_summary(values["trajectory_csv"],dp)
    env,sample,nwater,threshold,ff,solvent,temp=condition
    table=[]; metrics=[]; series=[]; curves={"OH":[],"N3":[]}
    for end,replicas in states.items():
        means=[]
        for identifier,data in replicas:
            rg,co,omega=np.mean(data[:,2:],axis=0)
            means.append([rg,co,omega]); curves[end].append(data)
            table.append([end,identifier,len(data),float(data[0,1]),float(data[-1,1]),float(rg),float(co),float(omega),float(np.sqrt(np.mean(data[:,2]**2)))])
        for k,(key,label,unit) in enumerate([
                ("rg","Polymer Rg","nm"),("co","Water per carbonyl O","H₂O/site"),("omega","Water per ω site","H₂O/site")]):
            arr=np.array(means)[:,k]
            metrics.append(dict(key=f"{end}_{key}",label=f"{end} · {label}",value=float(np.mean(arr)),unit=unit,
                note=f"Equal-weight mean of {len(arr)} trajectory means (range {min(arr):.3g}–{max(arr):.3g})."))
    for end,arr in curves.items():
        for (ident,_),data in zip(states.get(end,[]),arr):
            # Downsample plots only; table statistics above use every frame.
            if len(series)<12:
                idx=np.unique(np.linspace(0,len(data)-1,min(200,len(data)),dtype=int))
                series.append(dict(name=f"{end} · {ident}",x=data[idx,1].tolist(),y=data[idx,4].tolist()))
    warnings=["Independent starting packings are the units of replication. Consecutive frames are time-correlated.",
              "A water O within the cutoff counts as a neighbor. Donor geometry, lifetime and binding free energy require separate analyses."]
    if env.startswith("isolated"):
        warnings.append("Finite unconfined cluster in vacuum. This short pilot checks the calculation and site assignments; bulk hydration, long-chain scaling and equilibrium rankings remain untested.")
    return dict(metrics=metrics,plots=[dict(title="Water neighbors of the terminal site",x_label="Retained time (ps)",
        y_label="Water O atoms within cutoff",x_scale="linear",y_scale="linear",series=series)],
        tables=[dict(title="Trajectory means by independent starting configuration",
            columns=["ω end","Replica","Frames","First ps","Last ps","Mean Rg (nm)","Water / carbonyl","Water / ω site","RMS Rg (nm)"],rows=table),
            dict(title="Published periodic TIP3P reference · linear DP 30",columns=["Source","Force field","Waters","Duration (ns)","Published RMS Rg (nm)"],
                 rows=[["Galata & Kröger, Table 1","AMBER",1291,33,0.76],["Same source","DREIDING",1469,89,1.55],["Same source","UFF",7413,89,1.82]])],
        notes=[f"Series {sample}. {env}, {nwater} waters, {temp:g} K, force field {ff}, water model {solvent}, distance cutoff {threshold:g} nm.",
            f"DP {dp}, α-N-methyl. The ω site is the terminal hydroxyl O or the distal azide N. Carbonyl coordination is averaged across the {dp} repeat amide oxygens.",
            "The trajectory table reports arithmetic mean Rg and square root of mean Rg-squared separately. The published DP30 liquid simulations tabulate the latter statistic; the short DP2/4 cluster study has different environments and chain lengths. DOSY measures a hydrodynamic observable."],
        warnings=warnings,material_updates={},handoffs=[])


if __name__=="__main__":
    main()
