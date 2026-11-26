"""Write GROMACS/VMD-compatible GRO, PDB, XTC and index groups for one full replica.

Exports all atoms, persistent water labels, times and each sampled NPT box. The
original NPZ is the quantitative source; both GRO and XTC use rounded coordinates.
"""
from argparse import ArgumentParser
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import openmm as mm
from openmm import app, unit
from openmm.app.internal.xtc_utils import read_xtc


def file_digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()


def residue_groups(top):
    if "chain_a_indices" in top:
        groups=[("POA",top["chain_a_indices"]),("POB",top["chain_b_indices"])]
        indexes={"Chain_A":groups[0][1],"Chain_B":groups[1][1],"Amide_O_A":top["carbonyl_oxygens_a"],
                 "Amide_O_B":top["carbonyl_oxygens_b"]}
    else:
        groups=[("POL",top["solute_indices"])]
        indexes={"Polymer":top["solute_indices"],"Amide_O":top["carbonyl_oxygens"]}
    for oxygen in top["water_oxygens"]:
        if top["elements"][oxygen:oxygen+3] != ["O","H","H"]:
            raise ValueError("Water atom order does not match intact three-site water.")
        groups.append(("SOL",list(range(oxygen,oxygen+3))))
    if [i for _,idx in groups for i in idx] != list(range(len(top["elements"]))):
        raise ValueError("Atom IDs are not partitioned by ordered solute and waters.")
    indexes["Water_O"]=top["water_oxygens"]
    return groups,indexes


def make_graph(top, groups):
    graph=app.Topology()
    water_chain=None; atoms=[]
    for idx,(name,indices) in enumerate(groups,1):
        if name=="SOL":
            if water_chain is None: water_chain=graph.addChain("W")
            chain=water_chain;resname="HOH"
        else:
            chain=graph.addChain(chr(ord("A")+idx-1));resname=name
        res=graph.addResidue(resname,chain,str(idx))
        for i in indices:
            atoms.append(graph.addAtom(top["atom_names"][i],app.Element.getBySymbol(top["elements"][i]),res))
    for a,b in top["bonds"]:
        graph.addBond(atoms[a],atoms[b])
    return graph


def vec(box):
    return tuple(mm.Vec3(*axis) for axis in box)*unit.nanometer


def write_gro(output,top,groups,xyz,box,t):
    lines=[f"PEtOx–water | full periodic frame | t={t:.6f} ps",str(len(xyz))]
    for resid,(name,indices) in enumerate(groups,1):
        for i in indices:
            atom_name=top["atom_names"][i]
            lines.append(f"{resid:5d}{name:<5}{atom_name:>5}{i+1:5d}"+"".join(f"{v:8.3f}" for v in xyz[i]))
    lines.append(" ".join(f"{v:10.6f}" for v in np.diag(box)))
    output.write_text("\n".join(lines)+"\n",encoding="utf-8")


def export(raw,topology,manifest,out,name):
    rec=json.loads(manifest.read_bytes());top=json.loads(topology.read_bytes())
    if "parameters" in rec:
        expected=rec["trajectory_sha256"]
        top_digest=rec["shared_topology_sha256"];step_ps=rec["parameters"]["timestep_ps"]
    else:
        if not rec.get("complete"):
            raise ValueError("The contact-pair simulation must be complete before export.")
        expected=rec["files_sha256"][raw.name];top_digest=rec["files_sha256"][topology.name];step_ps=rec["dt_ps"]
    if file_digest(raw)!=expected or file_digest(topology)!=top_digest:
        raise ValueError("Source checksum mismatch.")
    with np.load(raw,allow_pickle=False) as dat:
        xyz,boxes,times=[dat[k] for k in ("positions_nm","box_nm","time_ps")]
    n=len(times);nat=len(top["elements"])
    if n<2 or xyz.shape!=(n,nat,3) or boxes.shape!=(n,3,3):
        raise ValueError("Invalid frame and topology dimensions.")
    delta=np.diff(times);time_step=int(np.rint(delta[0]/step_ps));first_step=int(np.rint(times[0]/step_ps))
    if (not np.allclose(delta,time_step*step_ps,atol=1e-7,rtol=0)
            or not np.isclose(times[0],first_step*step_ps,atol=1e-7,rtol=0) or time_step<1
            or not np.allclose(boxes, np.eye(3)[None]*np.diagonal(boxes,axis1=1,axis2=2)[:,None,:])):
        raise ValueError("Uniform time spacing and orthorhombic boxes are required.")
    groups,index=residue_groups(top);graph=make_graph(top,groups)
    graph.setPeriodicBoxVectors(vec(boxes[0]));out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        stage=Path(tmp); xfile=stage/(name+".xtc"); gfile=stage/(name+".gro")
        xtc=app.XTCFile(str(xfile),graph,dt=step_ps,firstStep=first_step,interval=time_step)
        for p,b in zip(xyz,boxes):xtc.writeModel(p*unit.nanometer,periodicBoxVectors=vec(b))
        write_gro(gfile,top,groups,xyz[0],boxes[0],times[0])
        with (stage/(name+".pdb")).open("w") as handle:
            app.PDBFile.writeFile(graph,xyz[0]*unit.nanometer,handle,keepIds=True)
        rows=[]
        for label,atoms in index.items():
            rows.append("[ "+label+" ]")
            for j in range(0,len(atoms),15):rows.append(" ".join(str(i+1) for i in atoms[j:j+15]))
        (stage/(name+".ndx")).write_text("\n".join(rows)+"\n")
        coords,cell,t,read_steps=read_xtc(os.fsencode(xfile))
        coords=coords.transpose(2,0,1);cell=cell.transpose(2,0,1)
        max_error=float(np.max(np.abs(coords-xyz)))
        gg=app.GromacsGroFile(str(gfile))
        gro_p=np.array(gg.positions.value_in_unit(unit.nanometer))
        if (max_error>0.00051 or not np.allclose(cell,boxes,atol=1e-7,rtol=0)
                or not np.allclose(t,times,atol=1.e-5,rtol=0) or not np.array_equal(read_steps,first_step+time_step*np.arange(n))
                or np.max(abs(gro_p-xyz[0]))>0.00051 or gg.atomNames!=top["atom_names"]):
            raise ValueError("GRO/XTC read-back failed the coordinate, indexing, timing or cell check.")
        pdb=app.PDBFile(str(stage/(name+".pdb")))
        read_atoms=list(pdb.topology.atoms())
        read_edges={tuple(sorted((a.index,b.index))) for a,b in pdb.topology.bonds()}
        if ([a.element.symbol for a in read_atoms] != top["elements"]
                or read_edges != {tuple(sorted(edge)) for edge in top["bonds"]}):
            raise ValueError("The PDB roundtrip changed elements or covalent topology.")
        record={"export_code_sha256":file_digest(Path(__file__)), "trajectory_sha256":expected,"topology_sha256":top_digest,"manifest_sha256":file_digest(manifest),
                "engine":mm.version.version,"frames":n,"atoms":nat,"waters":len(top["water_oxygens"]),
                "first_time_ps":float(times[0]),"last_time_ps":float(times[-1]),"stride_ps":float(delta[0]),
                "coordinate_precision_nm":0.001,"max_abs_component_error_nm":max_error,
                "coordinate_source":"Original periodic coordinates with per-frame box. No resampling or smoothing.",
                "files_sha256":{f.name:file_digest(f) for f in sorted(stage.iterdir())}}
        (stage/(name+".json")).write_text(json.dumps(record,indent=2)+"\n")
        for f in stage.iterdir():os.replace(f,out/f.name)
    return record


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    for arg in ("raw","topology","manifest","out"):parser.add_argument("--"+arg,type=Path,required=True)
    parser.add_argument("--name",required=True)
    print(json.dumps(export(**vars(parser.parse_args())),indent=2))
