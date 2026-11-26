"""Report sampled O=C–N–C(backbone toward ω) signed angles in the six extended PEtOx runs.

Use --hh data/hh/run --caps data/caps/run --out results/amide_torsions.csv.
Basin changes skip intermediate samples. Counts depend on cadence and have no kinetic interpretation.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from compare_solutions import Model, SITES, file_sha


def select(top, chemistry):
    model=Model(top,chemistry)
    N=[model.o_to_n[o] for o in model.refs]
    parent={b:a for a,b in model.tree}
    path=[N[-1]]
    while path[-1]!=N[0]:
        path.append(parent[path[-1]])
    path.reverse()
    ccarbonyl,=model.graph[int(model.refs[-1])]
    c1,=[i for i in model.graph[N[-1]] if model.elem[i]=='C' and i not in (ccarbonyl,path[-2])]
    c2,=[i for i in model.graph[c1] if model.elem[i]=='C']
    path.extend([c1,c2])
    if len(path)!=9 or [path[0],path[3],path[6]]!=N or [model.elem[i] for i in path]!=list('NCCNCCNCC'):
        raise ValueError('Cannot recover α→ω repeat backbone')
    groups=[]
    for k,o in enumerate(model.refs):
        carbonyl,=model.graph[int(o)]
        groups.append([int(o),carbonyl,N[k],path[3*k+1]])
    return groups,path


def torsions(positions,boxes,atoms):
    # Unwrap each covalent edge; no global chain origin or wrapping convention needed.
    L=np.diagonal(boxes,axis1=1,axis2=2)
    if np.max(np.abs(boxes-np.eye(3)[None,:,:]*L[:,:,None]))>1e-6:
        raise ValueError('This implementation requires an orthorhombic cell')
    b=np.diff(positions[:,atoms,:].astype(float),axis=1)
    b-=L[:,None,:]*np.rint(b/L[:,None,:])
    axis=b[:,1,:]/np.linalg.norm(b[:,1,:],axis=-1)[:,None]
    start=-b[:,0,:]
    v=start-np.einsum('ij,ij->i',start,axis)[:,None]*axis
    w=b[:,2,:]-np.einsum('ij,ij->i',b[:,2,:],axis)[:,None]*axis
    if np.any(np.linalg.norm(v,axis=1)*np.linalg.norm(w,axis=1)<1e-8):
        raise ValueError('Collinear, undefined dihedral')
    a=np.degrees(np.arctan2(np.einsum('ij,ij->i',np.cross(axis,v),w),np.einsum('ij,ij->i',v,w)))
    return a,np.concatenate([np.zeros_like(b[:,:1]),np.cumsum(b,axis=1)],axis=1)


def verify(coords,angles):
    from rdkit import Chem
    from rdkit.Chem.rdMolTransforms import GetDihedralDeg
    # This independent coordinate-based oracle needs only the unwrapped four points.
    mol=Chem.MolFromSmiles('CCCC'); conf=Chem.Conformer(4);mol.AddConformer(conf)
    for frame in (0,len(angles)//2,len(angles)-1):
        for i,xyz in enumerate(coords[frame]):mol.GetConformer().SetAtomPosition(i,xyz)
        expected=GetDihedralDeg(mol.GetConformer(),0,1,2,3)
        delta=(expected-angles[frame]+180)%360-180
        if abs(delta)>1e-6:raise ValueError('Signed angle disagrees with RDKit')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('hh','caps','out'):parser.add_argument('--'+k,type=Path,required=True)
    a=parser.parse_args()
    rows=[]
    for folder,chem in ((a.hh,'HH'),(a.caps,'MeOH')):
        top=json.loads((folder/'topology.json').read_text())
        groups,path=select(top,chem)
        for seed in (13031,27011,49033):
            file=folder/f'seed_{seed}.npz'
            metadata=json.loads((folder/f'seed_{seed}.json').read_text())
            if file_sha(file)!=metadata['trajectory_sha256']:raise ValueError('Trajectory hash mismatch')
            with np.load(file,allow_pickle=False) as data:
                xyz=data['positions_nm'];box=data['box_nm'];times=data['time_ps']
            steps=np.diff(times)
            if len(steps)<1 or not np.allclose(steps,steps[0]) or steps[0]<=0:raise ValueError('Invalid time grid')
            for site,atoms in zip(SITES,groups,strict=True):
                phi,coordinates=torsions(xyz,box,atoms);verify(coordinates,phi)
                labels=np.where(abs(phi)<=60,1,np.where(abs(phi)>=120,2,0))
                qualified=labels[labels!=0]
                changes=int(np.count_nonzero(np.diff(qualified)))
                rows.append(dict(chemistry=chem,seed=seed,site=site,atom_o=atoms[0],atom_c=atoms[1],atom_n=atoms[2],
                                 atom_c_toward_omega=atoms[3],backbone_alpha_to_omega='-'.join(map(str,path)),frames=len(phi),
                                 first_sample_ps=times[0],last_sample_ps=times[-1],sample_step_ps=steps[0],
                                 frac_abs_phi_le_60=float(np.mean(labels==1)),frac_abs_phi_ge_120=float(np.mean(labels==2)),
                                 frac_intermediate=float(np.mean(labels==0)),sampled_basin_changes=changes,
                                 first_signed_deg=phi[0],last_signed_deg=phi[-1],trajectory_sha256=metadata['trajectory_sha256']))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':
    main()
