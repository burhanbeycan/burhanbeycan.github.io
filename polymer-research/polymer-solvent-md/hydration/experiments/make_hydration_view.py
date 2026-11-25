"""Create a self-contained trajectory viewer from explicit-water two-chain MD."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

SHELL_NM = 0.35
CONTACT_NM = 0.45
HBOND_ANGLE_DEG = 150.0
HERE = Path(__file__).resolve().parent


def digest(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def minimum_image(vectors, lengths):
    return vectors - lengths * np.rint(vectors / lengths)


def whole_chain(xyz, indices, bonds, lengths):
    neighbors = {i: [] for i in indices}
    for a, b in bonds:
        if a in neighbors and b in neighbors:
            neighbors[a].append(b)
            neighbors[b].append(a)
    first = indices[0]
    pos = {first: xyz[first]}
    active = [first]
    while active:
        a = active.pop()
        for b in neighbors[a]:
            if b not in pos:
                pos[b] = pos[a] + minimum_image(xyz[b] - xyz[a], lengths)
                active.append(b)
    if len(pos) != len(indices):
        raise ValueError("Disconnected solute graph")
    return np.asarray([pos[a] for a in indices])


def as_pm(xyz):
    # Integer picometres keep the download small. All selections use full input precision.
    return np.rint(xyz * 1000).astype(int).reshape(-1).tolist()


def extract(topology, trajectory):
    elements = topology["elements"]
    a, b = topology["chain_a_indices"], topology["chain_b_indices"]
    atoms = a + b
    local = {ix: n for n, ix in enumerate(atoms)}
    local_bonds = [[local[x],local[y]] for x,y in topology["bonds"] if x in local and y in local]
    ref = topology["carbonyl_oxygens_a"] + topology["carbonyl_oxygens_b"]
    oi = np.array(topology["water_oxygens"])
    if not all(elements[i:i+3] == ["O", "H", "H"] for i in oi):
        raise ValueError("Expected explicit O,H,H water ordering")
    ai = np.array([i for i in a if elements[i] != "H"])
    bi = np.array([i for i in b if elements[i] != "H"])
    frames, extent = [], 0
    with np.load(trajectory, allow_pickle=False) as data:
        for t, xyz, cell in zip(data["time_ps"],data["positions_nm"],data["box_nm"],strict=True):
            xyz = xyz.astype(float)
            box = np.diag(cell)
            if not np.allclose(cell, np.diag(box)) or np.any(box<=0):
                raise ValueError("This projection requires an orthorhombic cell")
            pa = whole_chain(xyz,a,topology["bonds"],box)
            pb = whole_chain(xyz,b,topology["bonds"],box)
            # Translate an entire chain, preserving its covalent geometry.
            pb -= box*np.rint((pb.mean(0)-pa.mean(0))/box)
            polymer_xyz = np.vstack([pa,pb])
            center = (pa.mean(0) + pb.mean(0))/2
            inter = minimum_image(xyz[ai,None,:]-xyz[bi],box)
            n_contacts = int((np.sum(inter*inter,axis=-1)<=CONTACT_NM**2).sum())
            d_a = minimum_image(xyz[oi,None,:]-xyz[ai],box)
            d_b = minimum_image(xyz[oi,None,:]-xyz[bi],box)
            da2 = (d_a*d_a).sum(axis=-1)
            db2 = (d_b*d_b).sum(axis=-1)
            shared = (da2.min(1)<=SHELL_NM**2) & (db2.min(1)<=SHELL_NM**2)
            dr = minimum_image(xyz[oi,None,:]-xyz[ref],box)
            d2 = np.sum(dr*dr, axis=-1)
            shell = d2.min(1)<=SHELL_NM**2
            waters, hydrogen_bonds = [], []
            for n in np.flatnonzero(shell | shared):
                w = oi[n]
                if shell[n]:
                    r = int(np.argmin(d2[n]))
                    anchor = local[ref[r]]
                    oxygen = polymer_xyz[anchor] + dr[n,r]
                else:
                    i = int(np.argmin(da2[n]))
                    j = int(np.argmin(db2[n]))
                    if da2[n,i] < db2[n,j]:
                        oxygen = polymer_xyz[local[ai[i]]] + d_a[n,i]
                    else:
                        oxygen = polymer_xyz[local[bi[j]]] + d_b[n,j]
                oh = minimum_image(xyz[[w+1,w+2]]-xyz[w],box)
                sites = np.vstack([oxygen, oxygen+oh])
                waters.append([int(n)+1,int(shared[n]),*as_pm(sites-center)])
                extent=max(extent,float(np.max(np.linalg.norm(sites-center,axis=1))))
                for ri in np.flatnonzero(d2[n]<=SHELL_NM**2):
                    # At H, the two vectors point toward donor water O and acceptor amide O.
                    h_to_o = -oh
                    h_to_a = minimum_image(xyz[ref[ri]] - xyz[[w+1,w+2]],box)
                    cos = np.sum(h_to_o*h_to_a,axis=1)/(np.linalg.norm(h_to_o,axis=1)*np.linalg.norm(h_to_a,axis=1))
                    for h in np.flatnonzero(cos<=np.cos(np.deg2rad(HBOND_ANGLE_DEG))):
                        target = local[ref[ri]]
                        if abs(np.linalg.norm(polymer_xyz[target]-sites[1+h])-np.linalg.norm(h_to_a[h]))>0.003:
                            # Other periodic copy: draw exact image-local end point rather than a box-spanning line.
                            end = as_pm(sites[1+h]+h_to_a[h]-center)
                        else:
                            end = []
                        hydrogen_bonds.append([int(n)+1,int(h)+1,target,*end])
            extent=max(extent,float(np.max(np.linalg.norm(polymer_xyz-center,axis=1))))
            frames.append({"t":float(t),"p":as_pm(polymer_xyz-center),"w":waters,"hb":hydrogen_bonds,
                           "c":n_contacts,"s":int(shared.sum()),"near":int(shell.sum())})
    return {"atoms":[elements[i] for i in atoms],"chain_a_n":len(a),"bonds":local_bonds,"co":[local[i] for i in ref],
            "nwater":len(oi),"frames":frames,"radius_pm":round(1000*extent),"trajectory_sha256":digest(trajectory)}


def check_against_analysis(study, path):
    with path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(study["frames"]) != len(rows):
        raise ValueError("Coordinate frame and reference row counts disagree")
    for f,r in zip(study["frames"],rows):
        if (not np.isclose(f["t"],float(r["time_ps"])) or f["c"]!=int(r["heavy_pair_contacts_045"])
                or f["s"]!=int(r["shared_water_035"])
                or len(f["hb"]) != int(r["water_donor_hbonds_a"])+int(r["water_donor_hbonds_b"])):
            raise ValueError(f"Contact/water count mismatch at {f['t']} ps")


def add_tracked_water(case, topo, trajectory):
    """Track the first interposed water in every underlying frame, even outside the shell."""
    if not case["events"]:
        return
    event = case["events"][0]
    a,b=topo["chain_a_indices"],topo["chain_b_indices"]
    target_a=(a+b)[event["a"]]
    target_b=(a+b)[event["b"]]
    w=topo["water_oxygens"][event["w"]-1]
    with np.load(trajectory,allow_pickle=False) as raw:
        for f, xyz, cell in zip(case["frames"],raw["positions_nm"],raw["box_nm"],strict=True):
            xyz=xyz.astype(float)
            box=np.diag(cell)
            pa=whole_chain(xyz,a,topo["bonds"],box)
            pb=whole_chain(xyz,b,topo["bonds"],box)
            pb-=box*np.rint((pb.mean(0)-pa.mean(0))/box)
            center=(pa.mean(0)+pb.mean(0))/2
            oa=minimum_image(xyz[w]-xyz[target_a],box)
            ob=minimum_image(xyz[w]-xyz[target_b],box)
            o=pa[a.index(target_a)]+oa-center
            h=minimum_image(xyz[[w+1,w+2]]-xyz[w],box)+o
            f["tracked"]={"xyz":as_pm(np.vstack([o,h])),"da":round(float(np.linalg.norm(oa)),4),
                          "db":round(float(np.linalg.norm(ob)),4)}
    case["tracked_water"]=event["w"]


def build(input_root, out, temperature, pressure):
    runs = []
    for folder in sorted(input_root.glob("seed_*")):
        if not (folder/"trajectory.npz").is_file():
            continue
        topfile = folder/"topology.json"
        topo = json.loads(topfile.read_text())
        case = extract(topo,folder/"trajectory.npz")
        check_against_analysis(case,folder/"metrics.csv")
        protocol=json.loads((folder/"protocol.json").read_text())
        if not protocol.get("complete"):
            raise ValueError(f"Trajectory is still being written: {folder}")
        if protocol["temperature_K"]!=temperature or protocol["pressure_bar"]!=pressure:
            raise ValueError("Inconsistent ambient conditions")
        case.update(seed=int(folder.name.split("_")[-1]), topology_sha256=digest(topfile),
                    conditioning=protocol["condition_ps"],restraint=protocol["condition_restraint_kj_mol_nm2"],
                    stride=protocol["frame_stride_ps"],events=[])
        event_file = folder/"interposed_waters.csv"
        if event_file.exists():
            with event_file.open(encoding="utf-8") as stream:
                records = list(csv.DictReader(stream))
            atom_ids = topo["chain_a_indices"] + topo["chain_b_indices"]
            water_ids = list(topo["water_oxygens"])
            # All simultaneously qualifying triplets remain in the downloadable evidence CSV.
            for row in records:
                case["events"].append({"t":float(row["time_ps"]),"a":atom_ids.index(int(row["atom_a"])),
                        "b":atom_ids.index(int(row["atom_b"])),"w":water_ids.index(int(row["water_oxygen"]))+1,
                        "d0":float(row["initial_ab_nm"]),"d":float(row["current_ab_nm"]),"da":float(row["a_ow_nm"]),
                        "db":float(row["b_ow_nm"]),"angle":float(row["angle_a_ow_b_deg"])})

        add_tracked_water(case,topo,folder/"trajectory.npz")
        runs.append(case)
    if not runs:
        raise ValueError("No completed pair trajectories")
    data = {"runs":runs,"temperature":temperature,"pressure":pressure,"shell_nm":SHELL_NM,"contact_nm":CONTACT_NM,
            "min_angle":HBOND_ANGLE_DEG}
    template = (HERE/"hydration_view.template.html").read_text()
    out.write_text(template.replace('"__STUDY_JSON__"',json.dumps(data,allow_nan=False,separators=(",",":"))), encoding="utf-8")
    print(out, [(r["seed"],len(r["frames"]),r["nwater"]) for r in runs])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--temperature",type=float,default=298.15)
    p.add_argument("--pressure",type=float,default=1.0)
    args=p.parse_args()
    build(args.pair,args.out,args.temperature,args.pressure)
