"""Short dry-cluster control for LDI diurethane and two PEtOx-derived dimers.

Run after installing Chemistry, Chain and Studio: python reproduce.py --output /tmp/result.json
Independent packing seeds are paired across partners. No periodic bulk or solvent.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from polymer_chain_lab.classical_md import run_mmff

# Maps identify the oxygen environment without relying on atom order.
LDI = "CCOC(=[O:31])[C@H](CCCC[NH:2]C(=[O:42])OCCCCCC(=[O:21])OC)[NH:1]C(=[O:41])OCC[O:11]CC[O:12]C"
PARTNERS = {
    "EtOx/EtOx": "CN(C(=[O:51])CC)CCN(C(=[O:52])CC)CC[OH:61]",
    "EtOx/EI": "CN(C(=[O:51])CC)CC[NH:71]CC[OH:61]",
}
GROUP = {11:"PEG ether",12:"PEG ether",21:"PCL ester carbonyl",31:"LDI ester carbonyl",
         41:"urethane carbonyl",42:"urethane carbonyl",51:"EtOx amide carbonyl",
         52:"EtOx amide carbonyl",61:"partner terminal OH",71:"neutral EI nitrogen"}
SEEDS = (11, 37, 71, 103)


def conformer(smi, seed):
    mol=Chem.AddHs(Chem.MolFromSmiles(smi))
    opt=AllChem.ETKDGv3()
    opt.randomSeed=seed
    if AllChem.EmbedMolecule(mol,opt)!=0:
        raise RuntimeError("Initial ETKDG embedding failed.")
    return mol


def prepare(partner,seed,initialization):
    host=conformer(LDI,seed)
    guest=conformer(partner,seed+200)
    a=host.GetConformer().GetPositions()
    b=guest.GetConformer().GetPositions()
    a-=a.mean(axis=0);b-=b.mean(axis=0)
    ha=[i.GetIdx() for i in host.GetAtoms() if i.GetAtomicNum()>1]
    hb=[i.GetIdx() for i in guest.GetAtoms() if i.GetAtomicNum()>1]
    direction=np.random.default_rng(seed).normal(size=3)
    direction/=np.linalg.norm(direction)
    for offset in np.arange(30.,0.,-.04):
        shifted=b+offset*direction
        if cdist(a[ha],shifted[hb]).min()<3.25:
            break
    if initialization=="carbonyl contact":
        donor=next(i for i in host.GetAtoms() if i.GetAtomMapNum()==2)
        h=next(i.GetIdx() for i in donor.GetNeighbors() if i.GetAtomicNum()==1)
        acceptor=next(i for i in guest.GetAtoms() if i.GetAtomMapNum()==51)
        carbon=acceptor.GetNeighbors()[0].GetIdx()
        # Place the conserved EtOx oxygen in front of the host N-H. Rotate
        # C=O to a 120° H...O=C approach and choose the least crowded azimuth.
        axis=a[h]-a[donor.GetIdx()];axis/=np.linalg.norm(axis)
        perp=np.cross(axis,np.eye(3)[np.argmin(np.abs(axis))]);perp/=np.linalg.norm(perp)
        oc=b[carbon]-b[acceptor.GetIdx()];oc/=np.linalg.norm(oc)
        base=.5*axis+np.sqrt(.75)*perp
        rot,_=Rotation.align_vectors([base],[oc])
        candidates=[]
        for angle in np.arange(0,360,15):
            twist=Rotation.from_rotvec(axis*np.deg2rad(angle))
            for torsion in [0,90,180,270]:
                spin=Rotation.from_rotvec(base*np.deg2rad(torsion))
                xyz=(twist*spin*rot).apply(b-b[acceptor.GetIdx()])+a[h]+axis*1.95
                ds=cdist(a[ha],xyz[hb])
                candidates.append((float(np.sum(np.maximum(3.05-ds,0)**2)),xyz))
        shifted=min(candidates,key=lambda t:t[0])[1]
    mol=Chem.CombineMols(host,guest)
    mol.RemoveAllConformers()
    conf=Chem.Conformer(mol.GetNumAtoms());conf.Set3D(True)
    for i,point in enumerate(np.vstack((a,shifted))):
        conf.SetAtomPosition(i,point)
    mol.AddConformer(conf);Chem.SanitizeMol(mol)
    return mol


def analyse(mol,pos):
    distance=Chem.GetDistanceMatrix(mol)
    donors=[a for a in mol.GetAtoms() if a.GetAtomMapNum() in (1,2)]
    by_group=defaultdict(list)
    for a in mol.GetAtoms():
        if a.GetAtomMapNum() in GROUP:by_group[GROUP[a.GetAtomMapNum()]].append(a.GetIdx())
    rows=[]
    for group,acceptors in by_group.items():
        counts=np.zeros(len(pos)); eligible=0
        for d in donors:
            h=next(n.GetIdx() for n in d.GetNeighbors() if n.GetAtomicNum()==1)
            for a in acceptors:
                if distance[d.GetIdx(),a]<=3:continue
                eligible+=1
                hd=pos[:,d.GetIdx()]-pos[:,h];ha=pos[:,a]-pos[:,h]
                cs=np.sum(hd*ha,axis=1)/(np.linalg.norm(hd,axis=1)*np.linalg.norm(ha,axis=1))
                da=np.linalg.norm(pos[:,d.GetIdx()]-pos[:,a],axis=1)
                counts+=(da<=.35)&(cs<=np.cos(np.deg2rad(150.)))
        rows.append(dict(group=group,available_pairs=eligible,
                         contact_count=int(counts.sum()),frames=len(counts),
                         contacts_per_donor=float(counts.mean()/len(donors)),
                         occupancy_per_eligible_pair=float(counts.mean()/eligible) if eligible else None))
    return rows


def run():
    rows=[]
    from itertools import product
    for initialization,(name,partner),seed in product(("undirected","carbonyl contact"),PARTNERS.items(),SEEDS):
        mol=prepare(partner,seed,initialization)
        traj=run_mmff(mol,seed=seed+1000,temperature_K=298.15,equil_steps=1000,
                      production_steps=6000,stride=50,dt_fs=.5,friction_ps=5.)
        rows.append(dict(partner=name,initialization=initialization,seed=seed,atoms=mol.GetNumAtoms(),
                         rdkit_version=traj['rdkit_version'],last_time_ps=float(traj['time_ps'][-1]),
                         contacts=analyse(mol,traj['positions_nm'])))
    return dict(system='isolated neutral LDI diurethane + dimer, dielectric 1, no water',
                host_smiles=LDI,partner_smiles=PARTNERS,n_h_donors_on_host=2,seeds=SEEDS,
                method='MMFF94; BAOAB Langevin, 298.15 K, 0.5 fs, 5/ps friction',
                initial_relaxation_ps=.5,production_ps=3,save_interval_ps=.025,
                criteria='D–A <= 0.35 nm and D–H–A >= 150 deg; exclude <=3 graph bonds',
                scope='Local geometry and initialization sensitivity only; no equilibrium affinity, film morphology or macroscopic ranking.',
                runs=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(run(),indent=2,ensure_ascii=False)+'\n',encoding='utf8')
