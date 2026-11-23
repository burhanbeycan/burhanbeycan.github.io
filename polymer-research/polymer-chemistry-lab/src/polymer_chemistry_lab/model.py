"""RDKit graph construction; no GUI dependency and no property calibration.

Templates describe the already polymerized repeat connectivity, not monomer
reaction substrates. Endpoint 2 of unit i joins endpoint 1 of unit i+1.
Terminal open valences are capped by one protium atom each. molzip handles
attachment-point substitution, including local stereochemical bookkeeping.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
import re

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors

MAX_UNITS = 64
MAX_ATOMS = 512
MAX_TEMPLATE_ATOMS = 64
H_MASS = Chem.GetPeriodicTable().GetAtomicWeight(1)
H_EXACT_MASS = Chem.GetPeriodicTable().GetMostCommonIsotopeMass(1)

TEMPLATES = {
    "pe": {"name": "Polyethylene", "smiles": "[*:1]CC[*:2]", "formula": "C2H4"},
    "peo": {"name": "Poly(ethylene oxide)", "smiles": "[*:1]CCO[*:2]", "formula": "C2H4O"},
    "pva": {"name": "Poly(vinyl alcohol)", "smiles": "[*:1]CC(O)[*:2]", "formula": "C2H4O"},
    "pmma": {"name": "Poly(methyl methacrylate)", "smiles": "[*:1]CC(C)(C(=O)OC)[*:2]", "formula": "C5H8O2"},
    "petox": {"name": "Poly(2-ethyl-2-oxazoline)", "smiles": "[*:1]CCN(C(=O)CC)[*:2]", "formula": "C5H9NO"},
    "pmetox": {"name": "Poly(2-methyl-2-oxazoline)", "smiles": "[*:1]CCN(C(=O)C)[*:2]", "formula": "C4H7NO"},
}


def parse_template(smiles: str) -> Chem.Mol:
    """Validate a connected, closed-valence repeat with two single-bond ports."""
    if not isinstance(smiles, str) or not 1 <= len(smiles) <= 2000:
        raise ValueError("Repeat SMILES must contain 1–2000 characters.")
    if re.search(r"\s", smiles):
        raise ValueError("Use a single SMILES string without a name, spaces or CXSMILES annotations.")
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    params.removeHs = False
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles, params)
    if molecule is None:
        raise ValueError("RDKit could not parse and sanitize the repeat SMILES.")
    if molecule.GetNumAtoms() > MAX_TEMPLATE_ATOMS:
        raise ValueError(f"A repeat template may contain at most {MAX_TEMPLATE_ATOMS} atoms, including ports.")
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("A repeat template must be a single connected fragment; counterions are not separate inputs here.")
    ports = [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(ports) != 2 or sorted(atom.GetAtomMapNum() for atom in ports) != [1, 2]:
        raise ValueError("Exactly two dummy endpoints [*:1] and [*:2] are required.")
    for atom in molecule.GetAtoms():
        if atom.GetNumRadicalElectrons():
            raise ValueError("Radical repeat templates are unsupported.")
        if atom.GetAtomicNum() and atom.GetAtomMapNum():
            raise ValueError("Only the two dummy endpoints may have atom-map numbers.")
    for port in ports:
        if port.GetDegree() != 1 or port.GetIsotope() or port.GetFormalCharge():
            raise ValueError("Each dummy endpoint must be neutral, nonisotopic and have exactly one neighbour.")
        bond = port.GetBonds()[0]
        if bond.GetBondType() != Chem.BondType.SINGLE or port.GetNeighbors()[0].GetAtomicNum() == 0:
            raise ValueError("Each endpoint must have one single bond to an ordinary atom.")
    return molecule


def validate_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise ValueError("Sequence must be an explicit A/B string.")
    # Spaces/newlines may delimit a long sequence; no random sequence is inferred.
    clean = "".join(sequence.split()).upper()
    if not 1 <= len(clean) <= MAX_UNITS or any(symbol not in "AB" for symbol in clean):
        raise ValueError(f"Enter 1–{MAX_UNITS} A/B symbols, for example AAAAAA or AAABBB.")
    return clean


def repeat_information(molecule: Chem.Mol) -> dict:
    """Mass/formula of the polymerized contribution, excluding dummy atoms."""
    # The exactly two massless dummy atoms retain the port valences. Removing
    # them first would add implicit terminal H and change the repeat formula.
    return {
        "smiles": Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
        "formula": rdMolDescriptors.CalcMolFormula(molecule, separateIsotopes=True).replace("*2", ""),
        "average_mass_g_mol": float(Descriptors.MolWt(molecule)),
        "exact_mass_Da": float(rdMolDescriptors.CalcExactMolWt(molecule)),
        "formal_charge": int(Chem.GetFormalCharge(molecule)),
    }


@dataclass(frozen=True)
class Oligomer:
    molecule: Chem.Mol
    sequence: str
    repeats: dict

    @property
    def summary(self) -> dict:
        mol = self.molecule
        assigned = Chem.FindMolChiralCenters(mol, includeUnassigned=True, includeCIP=True)
        return {
            "sequence": self.sequence,
            "degree_of_polymerization": len(self.sequence),
            "canonical_isomeric_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True),
            "average_molar_mass_g_mol": float(Descriptors.MolWt(mol)),
            "exact_mass_Da": float(rdMolDescriptors.CalcExactMolWt(mol)),
            "end_group_mass_g_mol": float(2 * H_MASS),
            "formal_charge": int(Chem.GetFormalCharge(mol)),
            "heavy_atoms": int(mol.GetNumHeavyAtoms()),
            "assigned_tetrahedral_centres": sum(tag != "?" for _, tag in assigned),
            "unassigned_tetrahedral_centres": sum(tag == "?" for _, tag in assigned),
            "rdkit_version": rdBase.rdkitVersion,
        }


def build_oligomer(template_a: str, sequence: str, template_b: str | None = None) -> Oligomer:
    sequence = validate_sequence(sequence)
    templates = {"A": parse_template(template_a)}
    if "B" in sequence:
        templates["B"] = parse_template(template_b)
    if sum(templates[key].GetNumAtoms() - 2 for key in sequence) > MAX_ATOMS:
        raise ValueError(f"The finite chain may contain at most {MAX_ATOMS} non-port atoms.")
    fragments = []
    for unit_index, key in enumerate(sequence):
        fragment = Chem.Mol(templates[key])
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() == 0:
                # head i+1, tail i+2; each interior map then appears twice.
                atom.SetAtomMapNum(unit_index + atom.GetAtomMapNum())
            else:
                atom.SetIntProp("unit_index", unit_index)
                atom.SetProp("repeat_symbol", key)
        fragments.append(fragment)
    # Do not infer initiator or termination chemistry from a polymer's name.
    # These are deliberately the two explicitly stated graph caps.
    caps = Chem.MolFromSmiles(f"[H][*:1].[H][*:{len(sequence) + 1}]", sanitize=False)
    fragments.append(caps)
    combined = fragments[0]
    for fragment in fragments[1:]:
        combined = Chem.CombineMols(combined, fragment)
    try:
        with rdBase.BlockLogs():
            mol = Chem.molzip(combined)
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("The supplied ports could not form a sanitized H-capped linear oligomer.") from exc
    if len(Chem.GetMolFrags(mol)) != 1 or any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ValueError("Construction did not resolve all attachment points into one molecule.")
    if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        raise ValueError("The constructed graph contains unsupported radicals.")
    repeats = {key: repeat_information(template) for key, template in templates.items()}
    expected = sum(repeats[key]["average_mass_g_mol"] for key in sequence) + 2 * H_MASS
    if not math.isclose(Descriptors.MolWt(mol), expected, abs_tol=1e-7, rel_tol=1e-10):
        raise ValueError("The terminal-capping mass balance failed; this template is outside the supported construction scope.")
    return Oligomer(mol, sequence, repeats)


def molecular_view(oligomer: Oligomer) -> dict:
    """Hydrogen-suppressed 2D connectivity, with wedge/dash and unit lineage."""
    mol = Chem.Mol(oligomer.molecule)
    rdDepictor.Compute2DCoords(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        point = conf.GetAtomPosition(atom.GetIdx())
        element = atom.GetSymbol()
        # Skeletal carbons use the empty label; isotope/charge/heteroatom
        # labels must remain visible. CIP is metadata, not inferred tacticity.
        label = ""
        if element != "C" or atom.GetIsotope() or atom.GetFormalCharge():
            isotope = str(atom.GetIsotope()) if atom.GetIsotope() else ""
            hydrogens = atom.GetTotalNumHs()
            hlabel = "" if not hydrogens else "H" + (str(hydrogens) if hydrogens > 1 else "")
            charge = atom.GetFormalCharge()
            clabel = "" if not charge else (str(abs(charge)) if abs(charge) > 1 else "") + ("+" if charge > 0 else "−")
            label = isotope + element + hlabel + clabel
        atoms.append({"index": atom.GetIdx(), "element": element, "x": float(point.x), "y": float(point.y),
                      "label": label, "charge": atom.GetFormalCharge(), "isotope": atom.GetIsotope(),
                      "chiral": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "",
                      "unit_index": atom.GetIntProp("unit_index") if atom.HasProp("unit_index") else -1})
    bonds = [{"begin": bond.GetBeginAtomIdx(), "end": bond.GetEndAtomIdx(),
              "order": float(bond.GetBondTypeAsDouble()), "stereo": str(bond.GetStereo()),
              "direction": str(bond.GetBondDir())} for bond in mol.GetBonds()]
    return {"title": "Finite oligomer · explicit connectivity", "caption": "Hydrogen-suppressed 2D connectivity diagram. Vertices and unlabelled ends denote carbon; implicit H complete valence. Drawing coordinates are arbitrary units, not an equilibrium conformation. Only specified stereochemistry is retained.",
            "smiles": oligomer.summary["canonical_isomeric_smiles"], "atoms": atoms, "bonds": bonds}


def sequence_statistics(sequence: str) -> dict:
    sequence = validate_sequence(sequence)
    counts = Counter(sequence)
    dyads = Counter(sequence[i:i + 2] for i in range(len(sequence) - 1))
    return {"counts": {key: counts[key] for key in "AB"},
            "fractions": {key: counts[key] / len(sequence) for key in "AB"},
            "dyad_counts": {key: dyads[key] for key in ("AA", "AB", "BA", "BB")},
            "dyad_fractions": {key: (dyads[key] / (len(sequence) - 1) if len(sequence) > 1 else None) for key in ("AA", "AB", "BA", "BB")}}
