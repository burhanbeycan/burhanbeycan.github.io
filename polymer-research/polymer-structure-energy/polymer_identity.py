"""Finite, versioned chemical graphs for the polymer workbench (RDKit).

Two singly bonded ports encode oriented, already-polymerized repeat units.
molzip handles neighbor replacement and stereochemistry at the junctions.
"""
from __future__ import annotations

import hashlib
import json
import random
from collections import Counter

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

SCHEMA = "polymer-identity/1"
MAX_DP = 32
MAX_HEAVY_ATOMS = 150
CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
H_MASS = Chem.GetPeriodicTable().GetAtomicWeight(1)


def _id(label, payload):
    data = json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    return label + "-" + hashlib.sha256((SCHEMA + "\n" + data).encode()).hexdigest()


def _read(smiles):
    if not isinstance(smiles, str) or not 1 <= len(smiles) <= 1800 or any(ch.isspace() for ch in smiles):
        raise ValueError("Boşluk içermeyen tek bir SMILES girin (en çok 1800 karakter).")
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Kimyasal kimlik tek bir bağlı, geçerli molekül olmalıdır.")
    if mol.GetNumHeavyAtoms() > MAX_HEAVY_ATOMS:
        raise ValueError(f"Model en fazla {MAX_HEAVY_ATOMS} ağır atom içerebilir.")
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons():
            raise ValueError("Açık kabuklu yapılar bu modelin kapsamı dışında.")
        if atom.GetAtomicNum() and atom.GetAtomMapNum():
            raise ValueError("Atom eşleme numaraları bağlantı noktalarına ayrılmıştır.")
    return mol


def _ports(mol, *, core=False):
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    by_label = {a.GetAtomMapNum(): a for a in atoms}
    count = len(atoms)
    if (core and not 3 <= count <= 6) or (not core and count != 2):
        raise ValueError("Çekirdek 3–6, tekrar birimi 2 numaralı bağlantı noktası içermelidir.")
    if sorted(by_label) != list(range(1, count + 1)):
        raise ValueError("Bağlantı noktalarını 1’den başlayarak ardışık numaralandırın.")
    for atom in atoms:
        if (atom.GetDegree() != 1 or atom.GetFormalCharge() or atom.GetIsotope() or
                atom.GetNeighbors()[0].GetAtomicNum() <= 1 or
                atom.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE):
            raise ValueError("Her nötr bağlantı noktası bir ağır atoma tek bağla bağlanmalıdır.")
    return by_label


def _positive_dp(value):
    if type(value) is not int or not 1 <= value <= MAX_DP:
        raise ValueError(f"Tekrar sayısı 1–{MAX_DP} arasında bir tam sayı olmalıdır.")
    return value


def _tag_from_sign(order, target, sign):
    """Convert parity in a fixed reference frame to the RDKit-local neighbor ordering."""
    perm = [target.index(atom) for atom in order]
    inversions = sum(a > b for i, a in enumerate(perm) for b in perm[i + 1:])
    tag = (CW, CCW)[sign < 0]
    return ((CCW if tag == CW else CW) if inversions % 2 else tag)


def _vinyl_site(mol):
    ports = _ports(mol)
    path = Chem.GetShortestPath(mol, ports[1].GetIdx(), ports[2].GetIdx())
    msg = "Otomatik taktisite, tek yan gruplu yönlendirilmiş –CH₂–CH(R)– vinil tekrarları içindir."
    if len(path) != 4:
        raise ValueError(msg)
    atoms = [mol.GetAtomWithIdx(i) for i in path[1:3]]
    if (any(a.GetAtomicNum() != 6 or a.IsInRing() or a.GetFormalCharge() or
            any(b.GetBondType() != Chem.BondType.SINGLE for b in a.GetBonds()) for a in atoms)
            or sorted(a.GetTotalNumHs() for a in atoms) != [1, 2]):
        raise ValueError(msg)
    pos = 1 if atoms[0].GetTotalNumHs() == 1 else 2
    site = mol.GetAtomWithIdx(path[pos])
    side = [a.GetIdx() for a in site.GetNeighbors() if a.GetIdx() not in path]
    if len(side) != 1 or mol.GetAtomWithIdx(side[0]).GetAtomicNum() <= 1:
        raise ValueError(msg)
    order = [a.GetIdx() for a in site.GetNeighbors()]
    ref = [path[pos - 1], path[pos + 1], side[0]]
    # The single implicit H is last in both permutations.
    return site.GetIdx(), _tag_from_sign(order, ref, 1), tuple(ref)


def _double_bond_records(parts):
    """Remember relative E/Z reference substituents, including replaced ports."""
    swaps = {}
    records = []
    for piece_id, part in enumerate(parts):
        for atom in part.GetAtoms():
            if atom.GetAtomicNum():
                atom.SetIntProp("_prsAddr", piece_id * 2048 + atom.GetIdx())
        def reference(atom):
            if atom.GetAtomicNum():
                return atom.GetIntProp("_prsAddr")
            return (atom.GetAtomMapNum(), atom.GetNeighbors()[0].GetIntProp("_prsAddr"))
        for atom in part.GetAtoms():
            if not atom.GetAtomicNum():
                port, adjacent = reference(atom)
                swaps.setdefault(port, []).append(adjacent)
        for bond in part.GetBonds():
            sa = list(bond.GetStereoAtoms())
            if bond.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ) and len(sa) == 2:
                records.append((reference(bond.GetBeginAtom()), reference(bond.GetEndAtom()),
                                reference(part.GetAtomWithIdx(sa[0])), reference(part.GetAtomWithIdx(sa[1])),
                                bond.GetStereo() == Chem.BondStereo.STEREOE))
    def resolved(x):
        if isinstance(x, int):
            return x
        label, own = x
        others = [addr for addr in swaps[label] if addr != own]
        if len(others) != 1:
            raise ValueError("Çift bağ stereokimyası bağlantı noktasında izlenemedi.")
        return others[0]
    return [(lo, hi, resolved(a), resolved(b), trans) for lo, hi, a, b, trans in records]


def _restore_double_bonds(mol, records):
    # A directional single bond shared by adjacent polyenes cannot simply be
    # copied per repeat. Reconstruct E/Z from each pair of actual substituents.
    for bond in mol.GetBonds():
        if bond.GetBondDir() in (Chem.BondDir.ENDUPRIGHT, Chem.BondDir.ENDDOWNRIGHT):
            bond.SetBondDir(Chem.BondDir.NONE)
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)
    remove = Chem.RemoveHsParameters()
    remove.removeDefiningBondStereo = True
    mol = Chem.RemoveHs(mol, remove)
    address = {a.GetIntProp("_prsAddr"): a.GetIdx() for a in mol.GetAtoms()}
    potential = {d.centeredOn for d in Chem.FindPotentialStereo(mol)
                 if d.type == Chem.StereoType.Bond_Double}
    for u,v,ref_u,ref_v,trans in records:
        i,j = address[u], address[v]
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond.GetIdx() not in potential:  # e.g. terminal =CH2 after H capping
            continue
        refs = []
        for pos, other, target in ((i,j,ref_u), (j,i,ref_v)):
            attached = [a.GetIdx() for a in mol.GetAtomWithIdx(pos).GetNeighbors() if a.GetIdx() != other]
            if target in address and address[target] in attached:
                refs.append(address[target])
            elif len(attached) == 1:  # an ordinary H reference was removed
                refs.append(attached[0])
                trans = not trans
            else:
                raise ValueError("E/Z çift bağının referans atomları izlenemedi.")
        if i != bond.GetBeginAtomIdx():
            refs.reverse()
        bond.SetStereoAtoms(*refs)
        bond.SetStereo(Chem.BondStereo.STEREOE if trans else Chem.BondStereo.STEREOZ)
    Chem.SetDoubleBondNeighborDirections(mol)
    return mol


def _zip(parts, caps, expected_mass):
    caps_mol = Chem.MolFromSmiles(".".join(f"[H][*:{i}]" for i in caps), sanitize=False)
    records = _double_bond_records([caps_mol, *parts])
    mol = caps_mol
    for part in parts:
        mol = Chem.CombineMols(mol, part)
    try:
        mol = Chem.molzip(mol)
        Chem.SanitizeMol(mol)
        mol = _restore_double_bonds(mol, records) if records else Chem.RemoveHs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("Bağlantı noktaları geçerli bir H uçlu zincir oluşturmadı.") from exc
    if (len(Chem.GetMolFrags(mol)) != 1 or abs(Descriptors.MolWt(mol) - expected_mass) > 1e-5 or
            any(a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms())):
        raise ValueError("Sonlu grafın bağlantı veya atom kütlesi denetimi başarısız.")
    return mol


def _copy_unit(unit, left, right, unit_index):
    piece = Chem.Mol(unit)
    for atom in piece.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum({1: left, 2: right}[atom.GetAtomMapNum()])
        else:
            atom.SetIntProp("_prsUnit", unit_index)
            atom.SetIntProp("_prsLocal", atom.GetIdx())
    return piece


def _ordered_stereo(mol, center, ref, signs):
    """Check each surviving site's handedness on the *assembled* graph in head-to-tail order."""
    realized = {}
    for a in mol.GetAtoms():
        if not a.HasProp("_prsLocal") or a.GetIntProp("_prsLocal") != center or a.GetChiralTag() not in (CW, CCW):
            continue
        u = a.GetIntProp("_prsUnit")
        neighbors = list(a.GetNeighbors())
        if len(neighbors) != 3:
            raise ValueError("Sonlu zincirin stereomerkezi üç ağır komşu ile izlenemedi.")
        def role(n):
            v = n.GetIntProp("_prsUnit")
            if v < u:
                return 0
            if v > u:
                return 1
            return ref.index(n.GetIntProp("_prsLocal"))
        order = [role(n) for n in neighbors]
        positive = _tag_from_sign(order, [0, 1, 2], 1)
        actual = 1 if a.GetChiralTag() == positive else -1
        if signs[u] != actual:
            raise ValueError("Birleşme sonrasında taktisite paritesi korunmadı.")
        realized[u] = actual
    dyads = "".join("m" if realized[i] == realized[i + 1] else "r"
                    for i in sorted(realized) if i + 1 in realized)
    if not dyads:
        raise ValueError("H/H uçlarıyla en az bir stereojenik komşu ikilisi için model DP’sini artırın.")
    return {"requested_orientation": "".join("+" if v == 1 else "-" for v in signs),
            "retained_units": [i + 1 for i in sorted(realized)], "dyads": dyads,
            "meso_dyads": dyads.count("m"), "racemo_dyads": dyads.count("r")}


def _canonical(mol):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _repeat_signature(mol):
    # The two polymerization ports are equivalent for counting constitutional
    # repeat types. Unit orientation and stereochemical recipe remain separate.
    bare = Chem.Mol(mol)
    for a in bare.GetAtoms():
        a.SetAtomMapNum(0)
    Chem.RemoveStereochemistry(bare)
    return _canonical(bare)


def build_structure(spec):
    if not isinstance(spec, dict):
        raise ValueError("Yapı tanımı bir nesne olmalıdır.")
    mode = spec.get("mode", "finite")
    if mode not in ("repeat", "sequence", "star", "finite"):
        raise ValueError("Geçersiz yapı oluşturma yöntemi.")
    tacticity = spec.get("tacticity", "as_drawn")
    if tacticity not in ("as_drawn", "isotactic", "syndiotactic", "atactic"):
        raise ValueError("Taktisite seçimi geçersiz.")
    if tacticity != "as_drawn" and mode != "repeat":
        raise ValueError("Otomatik taktisite yalnız yönlendirilmiş bir tür vinil tekrarı için kullanılabilir.")
    with rdBase.BlockLogs():
        template = _read(spec.get("smiles"))
    recipe = {"mode": mode, "smiles": _canonical(template), "tacticity": tacticity}
    pattern = None
    sequence = None
    arms = None
    dp = None
    expected = None
    comps = None
    cap_count = None
    if mode == "finite":
        if any(a.GetAtomicNum() == 0 for a in template.GetAtoms()):
            raise ValueError("Tam yapıda açık bir bağlantı noktası kalamaz. Tekrar birimi modunu seçin.")
        mol = template
        architecture = "explicit_graph"
        composition = "unspecified"
    else:
        _ports(template)
        types = {"A": template}
        if mode == "sequence":
            with rdBase.BlockLogs():
                types["B"] = _read(spec.get("repeat_b"))
            _ports(types["B"])
            seq = spec.get("sequence")
            if not isinstance(seq, str) or not 2 <= len(seq) <= MAX_DP or set(seq) != {"A","B"}:
                raise ValueError(f"A ve B içeren açık bir dizi girin (2–{MAX_DP} tekrar).")
            sequence = seq
            dp = len(seq)
            recipe.update(repeat_b=_canonical(types["B"]), sequence=seq)
        else:
            dp = _positive_dp(spec.get("dp"))
            sequence = "A" * dp
        recipe["dp"] = dp
        pieces = []
        cap_labels = []
        repetitions = 1
        if mode == "star":
            with rdBase.BlockLogs():
                core = _read(spec.get("core_smiles"))
            ports = _ports(core, core=True)
            arms = repetitions = len(ports)
            recipe["core_smiles"] = _canonical(core)
            architecture = "star"
            pieces.append(Chem.Mol(core))
        else:
            architecture = "linear_backbone"
        num_heavy = (pieces[0].GetNumHeavyAtoms() if pieces else 0) + repetitions * sum(types[t].GetNumHeavyAtoms() for t in sequence)
        if num_heavy > MAX_HEAVY_ATOMS:
            raise ValueError(f"Bu zincir {num_heavy} ağır atom içeriyor. Sınır {MAX_HEAVY_ATOMS}.")
        site = None
        if tacticity != "as_drawn":
            center, positive, ref = _vinyl_site(template)
            site = template.GetAtomWithIdx(center)
            first = -1 if site.GetChiralTag() in (CW, CCW) and site.GetChiralTag() != positive else 1
            if tacticity == "atactic":
                seed = spec.get("stereo_seed", 2403)
                if type(seed) is not int or not 0 <= seed <= 2147483647:
                    raise ValueError("Rastgele dizi tohumu 0–2147483647 arasında bir tam sayı olmalıdır.")
                recipe["stereo_seed"] = seed
                rng = random.Random(seed)
                signs = [first] + [rng.choice((-1, 1)) for _ in range(dp - 1)]
            else:
                signs = [first * (-1 if tacticity == "syndiotactic" and i % 2 else 1) for i in range(dp)]
        next_label = repetitions + 1
        for arm in range(repetitions):
            left = arm + 1
            if mode != "star":
                cap_labels.append(left)
            for i, name in enumerate(sequence):
                right = next_label
                next_label += 1
                piece = _copy_unit(types[name], left, right, i)
                if site is not None:
                    tag = positive if signs[i] > 0 else (CW if positive == CCW else CCW)
                    piece.GetAtomWithIdx(center).SetChiralTag(tag)
                pieces.append(piece)
                left = right
            cap_labels.append(left)
        expected = sum(Descriptors.MolWt(m) for m in pieces) + len(cap_labels) * H_MASS
        cap_count = len(cap_labels)
        with rdBase.BlockLogs():
            mol = _zip(pieces, cap_labels, expected)
        if site is not None:
            pattern = _ordered_stereo(mol, center, ref, signs)
        comps = dict(Counter(sequence * repetitions))
        composition = "copolymer" if mode == "sequence" and _repeat_signature(types["A"]) != _repeat_signature(types["B"]) else "homopolymer_repeats"
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    unspec = [str(x.type) for x in Chem.FindPotentialStereo(mol) if x.specified != Chem.StereoSpecified.Specified]
    canonical = _canonical(mol)
    return {"schema": SCHEMA, "mode": mode, "recipe": recipe, "recipe_id": _id("recipe", recipe),
            "canonical_smiles": canonical, "structure_id": _id("structure", canonical),
            "formula": rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True),
            "molar_mass": round(Descriptors.MolWt(mol), 6), "net_charge": Chem.GetFormalCharge(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(), "architecture": architecture,
            "composition_class": composition, "repeat_counts": comps, "sequence": sequence,
            "dp": dp, "arm_count": arms, "h_caps": cap_count, "tacticity": tacticity,
            "stereo_pattern": pattern, "unspecified_stereo_count": len(unspec),
            "unspecified_stereo_types": sorted(unspec), "stereo_resolved": not unspec,
            "rdkit_version": rdBase.rdkitVersion}
