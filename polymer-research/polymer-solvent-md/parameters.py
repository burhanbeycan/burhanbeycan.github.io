"""Generate and cache native Amber GAFF2/AM1-BCC topologies from an exact Stage 1 geometry."""
from collections import Counter
from decimal import Decimal
import errno
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from rdkit import Chem
from upstream.energy_store import _canonical_molecule, analyze, expected_model_id


SCHEMA = "polymer-studio/gaff2-parameterization/1"
BUNDLE = Path(__file__).resolve().parent / "pilot" / "parameters"
MANIFEST = "parameterization.json"
SQM = 'qm_theory="AM1", scfconv=1.d-10, grms_tol=0.0005, maxcyc=9999'
LEAP = """source leaprc.gaff2
loadamberparams solute.frcmod
POL = loadmol2 solute.mol2
check POL
saveamberparm POL solute.prmtop solute.inpcrd
quit
"""
PROTOCOL = {
    "force_field": "GAFF2", "charges": "Antechamber AM1-BCC, single MMFF94s selected conformer",
    "sqm_keywords": SQM, "multiplicity": 1, "formal_charge_e": 0,
    "charge_normalization": "equal offset on each explicit atom after Antechamber output; 10 decimal places",
    "formal_charge_policy": "all atoms formally neutral; single closed-shell C/H/N/O molecule; 2–100 heavy atoms; all potential stereo specified",
    "input_geometry": "Stage 1 lowest converged isolated conformer; canonical-order explicit H; MDL V2000, coordinates in angstrom to four decimals",
    "screening": "reject missing/marked parameters, score above 10, unconverged SQM, LEaP warnings/errors or graph/stereo/charge mismatch",
}
FILES = ("solute.sdf", "raw_charge.mol2", "solute.mol2", "solute.frcmod", "leap.in",
         "solute.prmtop", "solute.inpcrd", "sqm.in", "sqm.out", "antechamber.log", "tleap.log")
TOOLS = ("antechamber", "am1bcc", "bondtype", "atomtype", "parmchk2", "sqm", "tleap", "teLeap")
DATA = ("dat/leap/parm/gaff2.dat", "dat/leap/cmd/leaprc.gaff2", "dat/antechamber/PARMCHK.DAT",
        "dat/antechamber/ATOMTYPE_GFF2.DEF", "dat/antechamber/ATOMTYPE_BCC.DEF",
        "dat/antechamber/BCCPARM.DAT", "dat/antechamber/PARM_BLBA_GAFF2.DAT",
        "dat/antechamber/APS.DAT", "dat/antechamber/ATOM_EQU.TYPE")
ANTECHAMBER = ["antechamber", "-i", "solute.sdf", "-fi", "mdl", "-o", "raw_charge.mol2",
               "-fo", "mol2", "-c", "bcc", "-at", "gaff2", "-nc", "0", "-rn", "POL",
               "-s", "2", "-pf", "n", "-ek", SQM]
PARMCHK = ["parmchk2", "-i", "solute.mol2", "-f", "mol2", "-o", "solute.frcmod", "-s", "2"]
TLEAP = ["tleap", "-f", "leap.in"]


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _file_hash(path):
    return _digest(path.read_bytes())


def _canonical(smiles):
    canonical, mol = _canonical_molecule(smiles)
    if not 2 <= mol.GetNumHeavyAtoms() <= 100:
        raise ValueError("MD parametreleme 2–100 ağır atom içeren sonlu bir yapı gerektirir.")
    for a in mol.GetAtoms():
        if a.GetAtomicNum() not in (1, 6, 7, 8) or a.GetNumRadicalElectrons() or a.GetFormalCharge() or a.GetIsotope():
            raise ValueError("MD parametreleme yüksüz, açık kabuk ve izotop etiketi taşımayan C/H/N/O yapıları destekler.")
    if any(str(s.specified) != "Specified" for s in Chem.FindPotentialStereo(mol)):
        raise ValueError("Olası stereo merkezleri ve çift bağ geometrilerini tanımlayın.")
    return canonical


def _source_input(canonical, record):
    source = record.get("energy", record)
    if source.get("canonical_smiles") != canonical or not source.get("model_id") or source.get("status") not in ("ok", "partial"):
        raise ValueError("Stage 1 enerjisi ve kanonik kimyasal yapı eşleşmiyor.")
    g = source.get("lowest_geometry") or {}
    mol = Chem.AddHs(Chem.MolFromSmiles(canonical))
    atoms = g.get("atoms", [])
    if g.get("unit") != "angstrom" or len(atoms) != mol.GetNumAtoms():
        raise ValueError("Kaynak konformer atom sayısı veya uzunluk birimi uyuşmuyor.")
    expected = {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                 b.GetBondTypeAsDouble()) for b in mol.GetBonds()}
    actual = {(min(b["a"], b["b"]), max(b["a"], b["b"]), b["order"]) for b in g.get("bonds", [])}
    if actual != expected or len(g["bonds"]) != len(expected):
        raise ValueError("Kaynak konformer bağları uyuşmuyor.")
    conf = Chem.Conformer(len(atoms))
    for i, (a, row) in enumerate(zip(mol.GetAtoms(), atoms)):
        xyz = tuple(float(row[k]) for k in ("x", "y", "z"))
        if (a.GetSymbol(), a.GetFormalCharge()) != (row["element"], row["charge"]) or not all(map(math.isfinite, xyz)):
            raise ValueError("Kaynak konformer atom sırası, yükleri veya koordinatları uyuşmuyor.")
        conf.SetAtomPosition(i, xyz)
    mol.AddConformer(conf)
    _check_stereo(canonical, mol)
    mol.SetProp("_Name", "POL")
    return {
        "canonical_smiles": canonical, "source_model_id": source["model_id"],
        "source_geometry_sha256": _digest(_json(g).encode()), "protocol": PROTOCOL,
        "schema": SCHEMA,
    }, Chem.MolToMolBlock(mol)


def _check_stereo(canonical, molecule):
    if Chem.FindPotentialStereo(Chem.RemoveHs(molecule)):
        checked = Chem.Mol(molecule)
        Chem.AssignStereochemistryFrom3D(checked, replaceExistingTags=True)
        if Chem.MolToSmiles(Chem.RemoveHs(checked), isomericSmiles=True) != canonical:
            raise ValueError("3B koordinatlarda belirtilen stereokimya korunmadı.")


def _amber_home(explicit=None):
    if explicit:
        candidates = [Path(explicit)]
    else:
        candidates = []
        if os.getenv("AMBERHOME"):
            candidates.append(Path(os.environ["AMBERHOME"]))
        exe = shutil.which("antechamber")
        if exe:
            candidates.append(Path(exe).resolve().parent.parent)
        spec = importlib.util.find_spec("acpype")
        if spec and spec.origin:
            candidates.append(Path(spec.origin).resolve().parent / "amber_linux")
    for home in candidates:
        if all((home / "bin" / x).is_file() for x in ("antechamber", "sqm", "parmchk2", "tleap")) and (home / DATA[0]).is_file():
            return home.resolve()
    if explicit:
        raise ValueError("Belirtilen AMBERHOME kullanılabilir AmberTools program ve verilerini içermiyor.")
    return None


def _toolchain(home):
    tool_files = {}
    for name in TOOLS:
        for partial in ("bin/" + name, "bin/wrapped_progs/" + name):
            if (home / partial).is_file():
                tool_files[partial] = _file_hash(home / partial)
    # Amber distributions differ in how data and executables are packaged.
    data_files = {name: _file_hash(home / name) for name in DATA if (home / name).is_file()}
    with (home / DATA[0]).open() as handle:
        version = handle.readline().strip()
    return {"gaff2_version": version, "data_sha256": data_files, "executable_sha256": tool_files}


def _normalize(raw, dest):
    lines = raw.read_text().splitlines()
    low = lines.index("@<TRIPOS>ATOM") + 1
    high = lines.index("@<TRIPOS>BOND")
    charges = [Decimal(s.split()[8]) for s in lines[low:high]]
    net = sum(charges)
    offset = -net / len(charges)
    if abs(net) > Decimal("0.01") or abs(offset) > Decimal("0.0005"):
        raise ValueError("AM1-BCC toplam yük sapması yuvarlama sınırını aşıyor.")
    new_charges = [format(q + offset, ".10f") for q in charges]
    for i, q in enumerate(new_charges, low):
        lines[i] = lines[i].rsplit(maxsplit=1)[0] + "  " + q
    if abs(sum(map(Decimal, new_charges))) > Decimal("0.0000001"):
        raise ValueError("AM1-BCC yük normalizasyonu başarısız.")
    dest.write_text("\n".join(lines) + "\n")
    return {"raw_net_e": float(net), "uniform_offset_per_atom_e": float(offset),
            "normalized_mol2_net_e": float(sum(map(Decimal, new_charges))),
            "atom_count": len(charges)}


def _supplements(content):
    if re.search(r"ATTN|need(?:s)? revision|missing parameter|unable to find", content, re.I):
        raise ValueError("parmchk2 eksik veya düzeltilmesi gereken bir parametre raporladı.")
    scores = [float(v) for v in re.findall(r"penalty score\s*=\s*([\d.]+)", content, re.I)]
    if max(scores, default=0) > 10:
        raise ValueError("parmchk2 analog parametre puanı bu taramanın kabul sınırını aşıyor.")
    section = None
    counts = Counter()
    for line in content.splitlines()[1:]:
        row = line.strip()
        if row in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"):
            section = row
        elif row and section:
            counts[section] += 1
    return {"supplementary_rows": dict(counts), "max_reported_analogy_penalty": max(scores, default=0),
            "positive_penalty_rows": sum(s > 0 for s in scores),
            "policy_penalty_limit": 10, "score_is_error_bound": False}


def _native_checks(folder, canonical):
    from openmm import NonbondedForce, app, unit
    prm = app.AmberPrmtopFile(str(folder / "solute.prmtop"))
    inp = app.AmberInpcrdFile(str(folder / "solute.inpcrd"))
    system = prm.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
    nb, = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
    mol = Chem.AddHs(Chem.MolFromSmiles(canonical))
    if [a.element.symbol for a in prm.topology.atoms()] != [a.GetSymbol() for a in mol.GetAtoms()]:
        raise ValueError("Native Amber atom sırası kaynak molekülle uyuşmuyor.")
    native_bonds = {tuple(sorted((a.index, b.index))) for a, b in prm.topology.bonds()}
    expected_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    if native_bonds != expected_bonds:
        raise ValueError("Native Amber bağ bağlantıları kaynak molekülle uyuşmuyor.")
    qs = []
    for i in range(system.getNumParticles()):
        q, sigma, epsilon = nb.getParticleParameters(i)
        values = [q.value_in_unit(unit.elementary_charge), sigma.value_in_unit(unit.nanometer),
                  epsilon.value_in_unit(unit.kilojoule_per_mole)]
        if not all(map(math.isfinite, values)) or min(values[1:]) < 0:
            raise ValueError("Native Amber yük veya LJ parametresi geçersiz.")
        qs.append(values[0])
    if abs(sum(qs)) > 0.0001:
        raise ValueError("Native Amber sisteminin toplam kısmi yükü nötr değil.")
    conf = Chem.Conformer(len(qs))
    positions = inp.getPositions().value_in_unit(unit.angstrom)
    if len(positions) != len(qs):
        raise ValueError("Native Amber koordinat sayısı uyuşmuyor.")
    for i, row in enumerate(positions):
        if not all(map(math.isfinite, row)):
            raise ValueError("Native Amber koordinatları geçersiz.")
        conf.SetAtomPosition(i, tuple(row))
    mol.AddConformer(conf)
    _check_stereo(canonical, mol)
    return {"num_atoms": len(qs), "num_bonds": len(native_bonds), "native_charge_e": sum(qs),
            "nonbonded_exceptions": nb.getNumExceptions(), "force_classes": [f.__class__.__name__ for f in system.getForces()]}


def _key(description):
    return _digest(_json(description).encode())


def _manifest(folder, desc, chain, normalization):
    if not all((folder / f).is_file() for f in FILES):
        raise ValueError("Amber parametre çıktısı eksik.")
    sqm = (folder / "sqm.out").read_text()
    if "geometry converged" not in sqm or "Calculation Completed" not in sqm:
        raise ValueError("SQM AM1 geometrisi yakınsamadı.")
    leap = (folder / "tleap.log").read_text()
    if not re.search(r"Errors\s*=\s*0;\s+Warnings\s*=\s*0;", leap):
        raise ValueError("LEaP hata veya uyarı bildirdi.")
    prm_id = "sha256:" + _key({"input": desc, "toolchain": chain})
    record = dict(desc, parameter_id=prm_id, toolchain=chain, charge_balance=normalization,
                  native_validation=_native_checks(folder, desc["canonical_smiles"]),
                  parameter_analogies=_supplements((folder / "solute.frcmod").read_text()),
                  sqm_version=next((line.strip() for line in sqm.splitlines() if "AMBER SQM VERSION" in line), None),
                  files_sha256={name: _file_hash(folder / name) for name in FILES},
                  qualification="General transferable GAFF2 with single-conformer AM1-BCC. Parameter coverage and native topology checks do not establish property accuracy. No molecule-specific calibration or predictive error bound.")
    (folder / MANIFEST).write_text(json.dumps(record, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    return record


def _cached(folder, description, chain=None):
    if not (folder / MANIFEST).is_file():
        return None
    r = json.loads((folder / MANIFEST).read_text())
    if any(r.get(k) != val for k, val in description.items()) or (chain and chain != r.get("toolchain")):
        return None
    if r.get("parameter_id") != "sha256:" + _key({"input":description, "toolchain":r["toolchain"]}):
        raise ValueError("Parametre kaydının kimliği uyuşmuyor.")
    if set(r["files_sha256"]) != set(FILES) or not all(
            (folder / name).is_file() and _file_hash(folder / name) == r["files_sha256"][name] for name in FILES):
        raise ValueError("Parametre dosyası bütünlük denetiminden geçemedi.")
    return r


def _execute(folder, home, command, logname):
    env = dict(os.environ, AMBERHOME=str(home))
    env["PATH"] = str(home / "bin") + os.pathsep + env.get("PATH", "")
    with (folder / logname).open("w") as log:
        code = subprocess.run([str(home / "bin" / command[0]), *command[1:]], cwd=folder,
                              env=env, stdout=log, stderr=subprocess.STDOUT, timeout=900).returncode
    raw = (folder / logname).read_text(errors="replace").replace(str(home), "$AMBERHOME").replace(str(folder), "$WORKDIR")
    (folder / logname).write_text(raw)
    if code:
        raise ValueError(f"{command[0]} bilimsel parametreleme aşaması başarısız (çıkış {code}).")


def parameterize(canonical_smiles, workdir, *, source_record=None, amber_home=None, rebuild=False):
    """Return manifest fields, directory and cache_source; only new jobs need AmberTools."""
    canonical = _canonical(canonical_smiles)
    cache = Path(workdir).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    if source_record is None:
        source_path = BUNDLE.parent / "source_model.json"
        if source_path.is_file():
            bundled = json.loads(source_path.read_text())["energy"]
            if bundled["canonical_smiles"] == canonical and bundled["model_id"] == expected_model_id(canonical):
                source_record = bundled
        if source_record is None:
            source_record = analyze(canonical, cache / "source-conformers.sqlite3")
    description, sdf = _source_input(canonical, source_record)
    home = _amber_home(amber_home)
    chain = _toolchain(home) if home is not None else None
    base = cache / _key(description)
    if not rebuild:
        candidates = [("bundled", BUNDLE)]
        candidates.extend(("local", item) for item in sorted(base.glob("*")) if item.is_dir())
        for label, directory in candidates:
            entry = _cached(directory, description, chain)
            if entry is not None:
                return dict(entry, directory=str(directory), cache_source=label)
    if home is None:
        raise RuntimeError("Bu yeni yapı AmberTools gerektiriyor. AMBERHOME ayarlayın veya AmberTools içeren ACPYPE Linux paketini kurun.")
    base.mkdir(parents=True, exist_ok=True)
    destination = base / _key(chain)
    with tempfile.TemporaryDirectory(dir=cache, prefix="amber-work-") as tmp:
        folder = Path(tmp)
        (folder / "solute.sdf").write_text(sdf)
        _execute(folder, home, ANTECHAMBER, "antechamber.log")
        normalization = _normalize(folder / "raw_charge.mol2", folder / "solute.mol2")
        _execute(folder, home, PARMCHK, "parmchk2.log")
        _supplements((folder / "solute.frcmod").read_text())
        (folder / "leap.in").write_text(LEAP)
        _execute(folder, home, TLEAP, "tleap.log")
        record = _manifest(folder, description, chain, normalization)
        staged = folder / "ready"
        staged.mkdir()
        for f in (*FILES, MANIFEST):
            shutil.move(folder / f, staged / f)
        try:
            staged.rename(destination)
        except OSError as exc:
            if exc.errno not in (errno.EEXIST, errno.ENOTEMPTY):
                raise
            record = _cached(destination, description, chain)
            if record is None:
                raise ValueError("Eşzamanlı oluşturulan parametre kaydı uyuşmuyor.") from exc
    return dict(record, directory=str(destination), cache_source="computed")
