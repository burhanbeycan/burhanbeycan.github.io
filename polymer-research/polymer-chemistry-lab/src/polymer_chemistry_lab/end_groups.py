"""Oriented PEtOx termini and prospective positive-ion MS targets."""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors

from .model import TEMPLATES, build_oligomer, molecular_view
from .sequences import integer

# This template runs C→N. Its leftmost backbone carbon is omega; the
# rightmost backbone nitrogen is alpha, opposite to the direction of CROP.
END_PAIRS = {
    "hh": ("H/H graph reference", "H", "H", None, None),
    "me_oh": ("α-methyl / ω-hydroxyl", "CH3", "OH", "C", "O"),
    "me_n3": ("α-methyl / ω-azide", "CH3", "N3", "C", "N=[N+]=[N-]"),
}
# Conventional LC-MS mass shifts (Da), Fiehn Lab MS Adduct Calculator.
PROTON_SHIFT = 1.007276
SODIUM_SHIFT = 22.989218


def petox_candidates(dp):
    dp = integer(dp, "PEtOx repeat count", 1, 64)
    base = build_oligomer(TEMPLATES["petox"]["smiles"], "A" * dp)
    ref = Chem.Mol(base.molecule)
    rdDepictor.Compute2DCoords(ref)
    # Each H-capped PEtOx graph has a single terminal backbone CH3 and a
    # single secondary amide NH. Side-chain methyls attach to carbonyl-adjacent
    # carbons instead of the backbone CH2–N.
    motifs = {"alpha": "[NH1](C=O)", "omega": "[CH3]-[CH2]-N"}
    sites = {}
    for role, smarts in motifs.items():
        matches = ref.GetSubstructMatches(Chem.MolFromSmarts(smarts))
        if len(matches) != 1:
            raise ValueError(f"Could not identify a unique {role} PEtOx attachment site.")
        sites[role] = matches[0][0]

    output = []
    for key, (title, alpha, omega, alpha_frag, omega_frag) in END_PAIRS.items():
        rw = Chem.RWMol(ref)
        for role, smiles in (("alpha", alpha_frag), ("omega", omega_frag)):
            rw.GetAtomWithIdx(sites[role]).SetProp("end_role", role)
            if smiles is None:
                continue
            fragment = Chem.MolFromSmiles(smiles)
            for atom in fragment.GetAtoms():
                atom.SetProp("end_role", role)
            first = rw.GetNumAtoms()
            rw.InsertMol(fragment)
            rw.AddBond(sites[role], first, Chem.BondType.SINGLE)
            for bond in rw.GetBonds():
                if max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) >= first:
                    bond.SetProp("end_role", role)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        if Chem.GetFormalCharge(mol) != 0:
            raise ValueError("A neutral covalent PEtOx graph was expected.")
        view = molecular_view(base, molecule=mol, alignment=ref)
        view.update(title=title, comparison_group="petox_end_groups",
                    caption="The N-terminal α site and carbon-terminal ω site share a common backbone orientation across all three drawings. Colours mark the sites, added atoms and their bonds.")
        avg = Descriptors.MolWt(mol)
        exact = rdMolDescriptors.CalcExactMolWt(mol)
        output.append({"key": key, "name": title, "alpha": alpha, "omega": omega,
                       "formula": rdMolDescriptors.CalcMolFormula(mol),
                       "smiles": Chem.MolToSmiles(mol),
                       "avg_g_mol": avg, "neutral_Da": exact,
                       "end_mass_g_mol": avg - dp * base.repeats["A"]["average_mass_g_mol"],
                       "mz_H": exact + PROTON_SHIFT, "mz_Na": exact + SODIUM_SHIFT,
                       "mz_2H": (exact + 2 * PROTON_SHIFT) / 2, "view": view})
    return output
