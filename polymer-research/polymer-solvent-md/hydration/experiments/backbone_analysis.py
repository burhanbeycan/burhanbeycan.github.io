"""Measure two matched geometric observables on the nine shared PEtOx backbone heavy atoms."""
import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from statistics import mean, stdev
import numpy as np
from rdkit import Chem


SEEDS=(13031,27011,49033)
MATCH="[N;X3]-[C;X4]-[C;X4]-[N;X3]-[C;X4]-[C;X4]-[N;X3]-[C;X4]-[C;X4]"


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def backbone_atoms(sdf,top):
    mol=Chem.MolFromMolFile(str(sdf),removeHs=False)
    sol=top['solute_indices']
    if sol != list(range(mol.GetNumAtoms())):
        raise ValueError("Unexpected solute topology index range.")
    if [a.GetSymbol() for a in mol.GetAtoms()]!=[top['elements'][i] for i in sol]:
        raise ValueError("SDF and trajectory element order differ.")
    bonds={tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    if bonds!={tuple(sorted((a,b))) for a,b in top['bonds'] if a in sol and b in sol}:
        raise ValueError("SDF and trajectory graphs differ.")
    paths=mol.GetSubstructMatches(Chem.MolFromSmarts(MATCH))
    if len(paths)!=1:
        raise ValueError("Expected a unique N(alpha)–C–C–N–C–C–N–C–C(omega) backbone.")
    return paths[0]


def geometric_series(positions,boxes,index,mass):
    positions=np.asarray(positions,dtype=float)
    boxes=np.asarray(boxes,dtype=float)
    if boxes.shape!=(len(positions),3,3):
        raise ValueError("Expected one periodic box per trajectory frame.")
    lengths=np.diagonal(boxes,axis1=1,axis2=2)
    if not np.allclose(boxes,np.eye(3)[None,:,:]*lengths[:,None,:]) or np.any(lengths<=0):
        raise ValueError("This analysis is registered for orthorhombic boxes.")
    sub=positions[:,index,:]
    delta=np.diff(sub,axis=1)
    delta-=lengths[:,None,:]*np.rint(delta/lengths[:,None,:])
    bond=np.linalg.norm(delta,axis=2)
    if not np.isfinite(bond).all() or np.any(bond<.09) or np.any(bond>.20):
        raise ValueError("Unphysical covalent backbone bond after periodic correction.")
    # Anchor N(alpha) at the origin and reconstruct the entire chain in bond order.
    whole=np.concatenate((np.zeros((len(delta),1,3)),np.cumsum(delta,axis=1)),axis=1)
    mass=np.asarray(mass,dtype=float)
    center=np.average(whole,weights=mass,axis=1)
    rg=np.sqrt(np.average(np.sum((whole-center[:,None,:])**2,axis=2),axis=1,weights=mass))
    span=np.linalg.norm(whole[:,-1,:],axis=1)
    return rg,span


def summarize(values):
    return {"mean_nm":mean(values),"between_lineage_sd_nm":stdev(values)}


def run(cases,out):
    individual=[]
    records=[]
    metrics={}
    proven={}
    reference_time=None
    reference_mass=None
    for name,(folder,sdf) in cases.items():
        folder,sdf=Path(folder),Path(sdf)
        topfile=folder/'topology.json'
        top=json.loads(topfile.read_text())
        idx=backbone_atoms(sdf,top)
        mass=[top['mass_da'][i] for i in idx]
        if reference_mass is not None and mass!=reference_mass:
            raise ValueError("Backbone mass mapping differs between the two endpoints.")
        reference_mass=mass
        refs={"trajectory_sha256":{},"metadata_sha256":{},"sdf_sha256":sha(sdf),"topology_sha256":sha(topfile),
              "backbone_indices_alpha_to_omega":list(idx),"mass_da":mass}
        proven[name]=refs
        for seed in SEEDS:
            path=folder/f'seed_{seed}.npz'
            meta_path=folder/f'seed_{seed}.json'
            meta=json.loads(meta_path.read_text())
            refs["trajectory_sha256"][str(seed)]=sha(path)
            refs["metadata_sha256"][str(seed)]=sha(meta_path)
            if meta.get('trajectory_sha256')!=sha(path) or meta.get('shared_topology_sha256')!=sha(topfile):
                raise ValueError("Run data fail the native checksum.")
            if meta['parameters']['seed']!=seed:
                raise ValueError("Mismatching trajectory seed.")
            with np.load(path,allow_pickle=False) as data:
                time=np.array(data['time_ps'])
                if reference_time is not None and not np.array_equal(reference_time,time):
                    raise ValueError("Production time grids must match within and between conditions.")
                reference_time=time
                rg,span=geometric_series(data['positions_nm'],data['box_nm'],idx,mass)
            for t,r,s in zip(time,rg,span):
                individual.append((name,seed,float(t),float(r),float(s)))
            records.append({"condition":name,"seed":seed,"frames":len(time),
                            "mean_backbone_rg_nm":float(rg.mean()),"mean_endpoint_separation_nm":float(span.mean()),
                            "first_half_rg_nm":float(np.array_split(rg,2)[0].mean()),"second_half_rg_nm":float(np.array_split(rg,2)[1].mean()),
                            "first_half_span_nm":float(np.array_split(span,2)[0].mean()),"second_half_span_nm":float(np.array_split(span,2)[1].mean())})
        rows=[r for r in records if r['condition']==name]
        metrics[name]={"backbone_rg":summarize([r['mean_backbone_rg_nm'] for r in rows]),
                       "alpha_n_to_omega_backbone_c":summarize([r['mean_endpoint_separation_nm'] for r in rows])}
    summary={"schema":"petox-matched-backbone/1","metrics":metrics,"replicas":records,"provenance":proven,
             "analysis_sha256":sha(__file__),"mass_weighting":"GAFF2 masses of three amide N and six backbone C; all hydrogens, side chains and caps excluded.",
             "endpoint":"Norm of the reconstructed vector from the first amide N to the final backbone C. Both are common atoms; omit OH and the α-methyl carbon.",
             "uncertainty":"Arithmetic mean of three lineage means and sample SD across those means. SD is descriptive; no significance or converged conformational equilibrium inferred."}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'backbone.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n')
    with (out/'backbone_traces.csv').open('w',newline='') as h:
        w=csv.writer(h);w.writerow(('condition','seed','time_ps','backbone_rg_nm','n_alpha_to_backbone_c_omega_nm'));w.writerows(individual)
    with (out/'backbone_replicas.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=records[0]);w.writeheader();w.writerows(records)
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('hh-run','capped-run','hh-sdf','capped-sdf','out'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    s=run({'H/H':(a.hh_run,a.hh_sdf),'Me/OH':(a.capped_run,a.capped_sdf)},a.out)
    print(json.dumps(s["metrics"],ensure_ascii=False,indent=2))
