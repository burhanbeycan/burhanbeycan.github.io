"""Reproducible isolated-molecule conformer screening and exact-protocol cache."""
from contextlib import closing
import hashlib
import json
import math
import sqlite3

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors


SCHEMA_VERSION = "conformer-energy/1"


def _canonical_molecule(smiles):
    if not isinstance(smiles, str) or not smiles or len(smiles) > 4000 or any(char.isspace() for char in smiles):
        raise ValueError("Bağlı, sonlu bir molekülün SMILES dizisi gereklidir.")
    parser = Chem.SmilesParserParams()
    parser.parseName = False
    parser.allowCXSMILES = False
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles, parser)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Enerji hesabı tek bir bağlı ve geçerli molekül gerektirir.")
    if any(a.GetAtomicNum() == 0 or a.GetAtomMapNum() for a in mol.GetAtoms()):
        raise ValueError("Enerji hesabından önce bağlantı noktalarını gerçek bağlarla kapatın.")
    if mol.GetNumHeavyAtoms() > 150:
        raise ValueError("Bu yapı kaydı en fazla 150 ağır atom işler.")
    # Reparse the canonical graph so equivalent input traversals produce identical coordinates.
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    normalized = Chem.MolFromSmiles(canonical, parser)
    return canonical, normalized


def _serialize(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def _coverage(mol):
    if mol.GetNumHeavyAtoms() > 100:
        return "Yapı kaydı geçerli. Yerel konformer enerji hesabı en fazla 100 ağır atom için uygulanır."
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return "Açık kabuklu yapı bu MMFF94s protokolünün kapsamı dışındadır."
    if any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53) for a in mol.GetAtoms()):
        return "Yapıda bu organik MMFF94s protokolünün kapsamı dışında bir element var."
    if any(str(s.specified) != "Specified" for s in Chem.FindPotentialStereo(mol)):
        return "Enerji hesabı için olası stereomerkezlerin ve çift bağ geometrisinin kimliği tanımlanmalıdır."
    if not AllChem.MMFFHasAllMoleculeParams(Chem.AddHs(mol)):
        return "Bu sonlu yapıdaki atomların tümü için MMFF94s parametresi bulunamadı."
    return None


def _geometry(mol, cid):
    xyz = mol.GetConformer(cid).GetPositions()
    return {
        "atoms": [
            {"element": atom.GetSymbol(), "charge": atom.GetFormalCharge(),
             "x": float(x), "y": float(y), "z": float(z)}
            for atom, (x, y, z) in zip(mol.GetAtoms(), xyz)
        ],
        "bonds": [{"a": b.GetBeginAtomIdx(), "b": b.GetEndAtomIdx(), "order": b.GetBondTypeAsDouble()}
                  for b in mol.GetBonds()],
        "unit": "angstrom",
    }


def _screen(mol, record):
    reason = _coverage(mol)
    if reason:
        record["message"] = reason
        return record
    spec = record["protocol"]
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = spec["seed"]
    params.numThreads = spec["threads"]
    params.maxIterations = spec["embedding_iterations"]
    params.pruneRmsThresh = spec["prune_rms_thresh_angstrom"]
    params.enforceChirality = True
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=record["requested"], params=params))
    record["embedded"] = len(ids)
    if not ids:
        record["message"] = "ETKDGv3 bu protokol ile konformer oluşturamadı."
        return record
    # Keep the original normalized bond orders for serialization and stereo checks.
    graph = Chem.Mol(mol)
    fits = AllChem.MMFFOptimizeMoleculeConfs(
        mol, numThreads=spec["threads"], maxIters=spec["optimization_iterations"],
        mmffVariant=spec["force_field"], nonBondedThresh=spec["nonbonded_threshold"],
        ignoreInterfragInteractions=True,
    )
    stereo = bool(Chem.FindPotentialStereo(Chem.RemoveHs(graph)))
    rows = []
    for cid, (code, energy) in zip(ids, fits):
        if code != 0 or not math.isfinite(energy):
            continue
        conf = mol.GetConformer(cid)
        if not all(math.isfinite(value) for position in conf.GetPositions() for value in position):
            continue
        if stereo:
            candidate = Chem.Mol(graph)
            candidate.RemoveAllConformers()
            candidate.AddConformer(conf, assignId=True)
            Chem.AssignStereochemistryFrom3D(candidate, confId=0, replaceExistingTags=True)
            result_smiles = Chem.MolToSmiles(Chem.RemoveHs(candidate), isomericSmiles=True)
            if result_smiles != record["canonical_smiles"]:
                continue
        rg = rdMolDescriptors.CalcRadiusOfGyration(mol, confId=cid, useAtomicMasses=True, force=True)
        if math.isfinite(rg):
            rows.append({"id": cid, "energy_kcal_mol": float(energy), "rg_angstrom": float(rg)})
    rows.sort(key=lambda item: (item["energy_kcal_mol"], item["id"]))
    record["converged"] = len(rows)
    if not rows:
        record["message"] = "Geometri, stereokimya ve yakınsama kontrollerini geçen konformer bulunamadı."
        return record
    minimum = rows[0]["energy_kcal_mol"]
    for row in rows:
        row["delta_kcal_mol"] = row["energy_kcal_mol"] - minimum
    # MMFF may change aromaticity flags. Recover original graph with optimized coordinates.
    best = Chem.Mol(graph)
    best.RemoveAllConformers()
    best.AddConformer(mol.GetConformer(rows[0]["id"]), assignId=True)
    record["lowest_geometry"] = _geometry(best, 0)
    record["conformers"] = rows
    record["status"] = "ok" if len(rows) == record["requested"] else "partial"
    record["message"] = (f"{len(rows)}/{record['requested']} aday yakınsama ve yapı kontrolünü geçti. "
                         "Enerji farkları yalnız bu molekülün örneklenen konformerleri içindir.")
    return record


def _model_definition(smiles, count, seed):
    if type(count) is not int or not 1 <= count <= 36:
        raise ValueError("Konformer deneme sayısı 1–36 aralığında bir tam sayı olmalıdır.")
    if type(seed) is not int or not 1 <= seed <= 2**31 - 1:
        raise ValueError("Tohum 1–2147483647 aralığında bir tam sayı olmalıdır.")
    canonical, molecule = _canonical_molecule(smiles)
    protocol = {
        "rdkit_version": rdBase.rdkitVersion, "embedding": "ETKDGv3", "force_field": "MMFF94s",
        "environment": "isolated; default MMFF dielectric 1; no solvent", "seed": seed,
        "threads": 1, "count": count, "prune_rms_thresh_angstrom": -1.0,
        "optimization_iterations": 1500, "embedding_iterations": 300, "nonbonded_threshold": 100.0,
        "hydrogens": "explicit", "enforce_chirality": True, "energy_unit": "kcal/mol",
        "rg": "mass-weighted including explicit H; isolated coordinates; angstrom",
        "duplicates": "unpruned trials; minima may repeat; frequencies are not populations",
    }
    key = {"schema_version": SCHEMA_VERSION, "canonical_smiles": canonical, "protocol": protocol}
    digest = hashlib.sha256(_serialize(key).encode("utf-8")).hexdigest()
    model_id = "sha256:" + digest
    return molecule, key, model_id


def expected_model_id(smiles, *, count=12, seed=2403):
    """Identifier for exactly this installed method, canonical graph and sampling protocol."""
    return _model_definition(smiles, count, seed)[2]


def analyze(smiles, db_path, *, count=12, seed=2403):
    """Screen one fully specified finite graph. Invalid inputs raise ValueError.

    Results carry the full protocol. Completed and deterministically unsupported records
    can be reused. Two processes may compute simultaneously, but both consume one committed
    record and neither can expose partial writes.
    """
    molecule, key, model_id = _model_definition(smiles, count, seed)
    with closing(sqlite3.connect(db_path, timeout=30.0)) as db:
        db.execute("CREATE TABLE IF NOT EXISTS conformer_runs (model_id TEXT PRIMARY KEY, record TEXT NOT NULL)")
        cached = db.execute("SELECT record FROM conformer_runs WHERE model_id=?", (model_id,)).fetchone()
        if cached:
            return dict(json.loads(cached[0]), cache_hit=True)
        record = dict(key, model_id=model_id, status="unavailable", requested=count, embedded=0,
                      converged=0, conformers=[], lowest_geometry=None, message="")
        try:
            result = _screen(molecule, record)
        except (ValueError, RuntimeError, OverflowError) as exc:
            # Transient engine failures may be retried; do not persist failed attempts.
            record["message"] = "Konformer hesap motoru tamamlanamadı: " + str(exc)
            return dict(record, cache_hit=False)
        with db:
            inserted = db.execute("INSERT OR IGNORE INTO conformer_runs(model_id,record) VALUES(?,?)",
                                  (model_id, _serialize(result))).rowcount
            stored, = db.execute("SELECT record FROM conformer_runs WHERE model_id=?", (model_id,)).fetchone()
        return dict(json.loads(stored), cache_hit=not bool(inserted))
