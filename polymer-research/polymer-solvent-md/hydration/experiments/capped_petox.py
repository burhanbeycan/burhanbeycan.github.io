"""Prepare an explicitly capped PEtOx trimer and generate the three pilot lineages.

Run prepare, then run the three pilot commands in parallel. extend_hydration.py
branches from each archived pilot endpoint using the same procedure as for H/H.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys

# Connectivity from the initiator-derived end toward the hydroxyl-terminated one.
EXPLICIT = "CN(C(=O)CC)CCN(C(=O)CC)CCN(C(=O)CC)CCO"
MODEL_ID = "petox-dp3-alpha-methyl-omega-hydroxyl"
CROP_DOI = "10.1039/D2PY00437B"
SEEDS = (13031, 27011, 49033)


def save_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n",
                    encoding="utf-8")


def identify_sites(mol):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    if rdMolDescriptors.CalcMolFormula(mol) != "C16H31N3O4":
        raise ValueError("Expected explicit methyl/hydroxyl terminated PEtOx DP3.")
    amides = mol.GetSubstructMatches(Chem.MolFromSmarts("[O;X1]=[C;X3]-[N;X3;H0]"))
    methyl = mol.GetSubstructMatches(Chem.MolFromSmarts("[CH3]-[N;X3]"))
    alcohol = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;X3]-[CH2]-[CH2]-[O;H1]"))
    if (len(amides), len(methyl), len(alcohol)) != (3, 1, 1):
        raise ValueError("Amide classes or end-group positions disagree with the intended chain.")
    alpha_methyl, alpha_n = methyl[0]
    omega_n, _, omega_end_c, hydroxyl = alcohol[0]
    ordered = sorted(amides, key=lambda site: len(Chem.GetShortestPath(mol,alpha_methyl,site[2])))
    if ordered[0][2] != alpha_n or ordered[-1][2] != omega_n:
        raise ValueError("Capping sites are misplaced along the chain.")
    return {
        "model_id":MODEL_ID,"end_groups":{"alpha":"N-methyl","omega":"2-hydroxyethyl backbone end"},
        "dp":3,"canonical_isomeric_smiles":Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True),
        "formula":rdMolDescriptors.CalcMolFormula(mol),"charge_e":Chem.GetFormalCharge(mol),
        "explicit_atoms":mol.GetNumAtoms(),"alpha_methyl_c":alpha_methyl, "omega_hydroxyl_o":hydroxyl,
        "omega_terminal_c":omega_end_c, "omega_hydroxyl_h":[a.GetIdx() for a in mol.GetAtomWithIdx(hydroxyl).GetNeighbors() if a.GetSymbol()=="H"],
        "ordered_sites_alpha_to_omega":[{"site":name,"amide_oxygen":o,"carbonyl_c":c,"amide_n":n,
                                         "class":"tertiary_amide"} for name,(o,c,n) in zip(("alpha","middle","omega"),ordered)],
        "index_base":0,"experimental_end_group_precedent_doi":CROP_DOI,
        "scope":"Neutral selected model of α-Me/ω-OH PEtOx; the polymerization reaction is not simulated."}


def prepare(out):
    from rdkit import Chem
    from upstream.polymer_identity import build_structure
    from upstream.energy_store import analyze
    from study import scientific_record
    from parameters import parameterize
    from openmm import app
    out.mkdir(parents=True, exist_ok=True)
    structure = build_structure({"mode":"finite", "smiles":EXPLICIT})
    energy = analyze(structure["canonical_smiles"],out/"conformers.sqlite",count=12,seed=2403)
    if not energy.get("lowest_geometry"):
        raise ValueError(energy.get("message","No geometry obtained"))
    source = scientific_record(structure,energy)
    save_json(out/"source_model.json",source)
    prm = parameterize(structure["canonical_smiles"],out/"cache",source_record=source["energy"])
    packet = out/"parameters"
    shutil.copytree(prm["directory"],packet,dirs_exist_ok=True)
    mol = Chem.MolFromMolFile(str(packet/"solute.sdf"),removeHs=False)
    native = app.AmberPrmtopFile(str(packet/"solute.prmtop")).topology
    if [a.element.symbol for a in native.atoms()] != [a.GetSymbol() for a in mol.GetAtoms()]:
        raise ValueError("Native topology and SDF atom order differ.")
    expected_edges = {tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    if {tuple(sorted((a.index,b.index))) for a,b in native.bonds()} != expected_edges:
        raise ValueError("Native topology and SDF graph differ.")
    sites = identify_sites(mol)
    if not len(sites["omega_hydroxyl_h"]) == 1:
        raise ValueError("Explicit omega O–H did not resolve.")
    sites.update(parameter_id=prm["parameter_id"],native_atom_order_checked=True,source_id=source["artifact_id"])
    save_json(out/"chemical_model.json",sites)
    return sites


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base",required=True,type=Path,help="Polymer–Solvent MD directory")
    p.add_argument("--out",required=True,type=Path,help="Case data directory")
    p.add_argument("--seed",type=int,choices=SEEDS,help="Run one preliminary lineage; omit to prepare")
    args = p.parse_args()
    sys.path.insert(0,str(args.base.resolve()))
    out = args.out.resolve()
    if args.seed is None:
        print(json.dumps(prepare(out), ensure_ascii=False, indent=2))
        return
    from simulation import run_replicate
    source = json.loads((out/"source_model.json").read_text())
    run_replicate(out/"parameters", out/"pilot", seed=args.seed,temperature_K=298.15,pressure_bar=1.0,
                  n_water=500,equil_ps=12.0,production_ps=24.0,sample_ps=.2,
                  geometry_record=source if args.seed==SEEDS[0] else None)


if __name__=="__main__":
    main()
