"""Identify a finite chain, run a registered water protocol and assemble its evidence."""
from hashlib import file_digest, sha256
import json
import math
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parent
PILOT = ROOT / "pilot"
SEEDS = (13031, 27011, 49033)
PROTOCOL = {"version": "polymer-water/1.0.0", "pressure_bar": 1.0, "equil_ps": 12.0,
            "production_ps": 24.0, "sample_ps": 0.2, "seeds": list(SEEDS),
            "water": "TIP3P", "parameters": "GAFF2/AM1-BCC", "water_count_rule": "automatic",
            "initial_states": "imported best conformer, two independent ETKDGv3/MMFF94s starts"}


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()


def digest(value):
    return sha256(encode(value)).hexdigest()


def file_hash(path):
    with Path(path).open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def write_record(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
        tmp.write(encode(value))
    os.replace(tmp.name, path)


def read(path):
    return json.loads(Path(path).read_bytes())


def scientific_record(structure, energy):
    energy = {k: v for k, v in energy.items() if k not in {"cache_hit", "cache_source"}}
    record = {"schema": "polymer-structure-energy/1", "structure": structure, "energy": energy}
    return {**record, "artifact_id": digest(record)}


def restore_export_numbers(value):
    """Restore decimal fields in the v1 scientific schema after JS stringification.

    Numeric values stay identical. The original full SHA-256 is verified afterwards.
    """
    def decimals(row, names):
        for field in names:
            if type(row.get(field)) is int:
                row[field] = float(row[field])

    from copy import deepcopy
    result = deepcopy(value)
    decimals(result["structure"], ("molar_mass",))
    energy = result["energy"]
    decimals(energy.get("protocol", {}), ("prune_rms_thresh_angstrom", "nonbonded_threshold"))
    for row in energy.get("conformers", []):
        decimals(row, ("delta_kcal_mol", "energy_kcal_mol", "rg_angstrom"))
    geometry = energy.get("lowest_geometry") or {}
    for atom in geometry.get("atoms", []):
        decimals(atom, ("x", "y", "z"))
    for bond in geometry.get("bonds", []):
        decimals(bond, ("order",))
    return scientific_record(result["structure"], energy)


def check_export(source):
    from upstream.polymer_identity import build_structure
    if not isinstance(source, dict) or source.get("schema") != "polymer-structure-energy/1":
        raise ValueError("İlk uygulamadan dışa aktarılan yapı ve enerji kaydını seçin.")
    structure, energy = source.get("structure"), source.get("energy")
    if not isinstance(structure, dict) or not isinstance(energy, dict):
        raise ValueError("Yapı veya konformer kaydı eksik.")
    record = scientific_record(structure, energy)
    if source.get("artifact_id") != record["artifact_id"]:
        try:
            record = restore_export_numbers(record)
        except (AttributeError, TypeError, OverflowError) as error:
            raise ValueError("Konformer dosyasının sayısal alanları geçersiz.") from error
    if source.get("artifact_id") != record["artifact_id"]:
        raise ValueError("Dışa aktarılan dosyanın içerik özeti uyuşmuyor.")
    rebuilt = build_structure(structure.get("recipe"))
    for key in ("structure_id", "canonical_smiles", "recipe_id", "formula", "dp"):
        if rebuilt[key] != structure.get(key):
            raise ValueError("Yapı tarifi ve taşınan kimyasal kimlik uyuşmuyor.")
    energy = record["energy"]
    if (energy.get("canonical_smiles") != rebuilt["canonical_smiles"] or energy.get("status") not in ("ok", "partial")
            or not energy.get("lowest_geometry") or not energy.get("model_id")):
        raise ValueError("Bu yapıya ait yakınsamış bir başlangıç konformeri gerekli.")
    return record


def check_chemistry(structure):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(structure["canonical_smiles"])
    if (structure["heavy_atoms"] > 100 or not structure["stereo_resolved"]
            or any(a.GetAtomicNum() not in (1, 6, 7, 8) or a.GetFormalCharge() or a.GetNumRadicalElectrons()
                   or a.GetIsotope() for a in mol.GetAtoms())):
        raise ValueError("Bu su protokolü stereokimyası belirli, yüksüz C/H/N/O oligomerlerini işler (en çok 100 ağır atom).")
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[O;X1]=[C;X3]-[#7]")):
        raise ValueError("Su–karbonil çalışması için yapıda en az bir amid karbonili gerekli.")


def inputs(payload):
    from upstream.polymer_identity import build_structure
    if not isinstance(payload, dict) or payload.get("solvent", "water") != "water":
        raise ValueError("Bu sürümde çözücü sudur.")
    value = payload.get("temperature_K", 298.15)
    if type(value) not in (float, int) or not math.isfinite(value) or not 273.15 <= value <= 350:
        raise ValueError("Su sıcaklığı 273,15–350 K arasında olmalıdır.")
    if payload.get("material") is not None:
        source = check_export(payload["material"])
        structure = source["structure"]
    else:
        source = None
        structure = build_structure({"mode": "repeat", "smiles": payload.get("repeat_smiles"), "dp": payload.get("dp")})
    check_chemistry(structure)
    return structure, source, float(value)


def request_key(structure, source, temperature):
    from upstream.energy_store import expected_model_id
    # Imported coordinate ensembles with different seeds/rounding remain distinct.
    model = source["energy"]["model_id"] if source else expected_model_id(structure["canonical_smiles"])
    return digest({"protocol": PROTOCOL, "recipe_id": structure["recipe_id"], "structure_id": structure["structure_id"],
                   "temperature_K": temperature, "conformer_model": model,
                   "geometry": digest(source["energy"]["lowest_geometry"]) if source else "generate"})


def runtime_fingerprint():
    from importlib.metadata import version
    modules = ("study.py", "simulation.py", "parameters.py", "amber_bridge.py", "trajectory_analysis.py",
               "upstream/polymer_identity.py", "upstream/energy_store.py")
    return digest({"files": {s: file_hash(ROOT/s) for s in modules},
                   "versions": {p: version(p) for p in ("numpy", "openmm", "rdkit")}})


def get_source(structure, cache):
    from upstream.energy_store import analyze, expected_model_id
    pilot = read(PILOT / "source_model.json")
    if (structure["recipe_id"] == pilot["structure"]["recipe_id"]
            and expected_model_id(structure["canonical_smiles"]) == pilot["energy"]["model_id"]):
        return pilot
    energy = analyze(structure["canonical_smiles"], Path(cache) / "conformers.sqlite")
    if energy.get("status") not in ("ok", "partial"):
        raise ValueError(energy.get("message") or "Konformer taraması tamamlanamadı.")
    return scientific_record(structure, energy)


def assemble(folder, source, temperature, params, seed_runs):
    folder = Path(folder)
    analysis = read(folder/"analysis/observables.json")
    s = source["structure"]
    run0 = seed_runs[0]["parameters"]
    if ([m["parameters"]["seed"] for m in seed_runs] != list(SEEDS)
            or any(m["parameters"]["temperature_K"] != temperature for m in seed_runs)
            or any(m["parameters"]["smiles"] != s["canonical_smiles"] for m in seed_runs)):
        raise ValueError("Koşuların tohumları, kimyası veya sıcaklıkları protokolle uyuşmuyor.")
    hashes = analysis["provenance"]["trajectories_sha256"]
    if any(hashes.get(m["trajectory"]) != m["trajectory_sha256"] for m in seed_runs):
        raise ValueError("Analiz ve simülasyon yörünge özetleri uyuşmuyor.")
    files = ["analysis/rg_trace.csv", "analysis/rdf.csv", "analysis/rg_distribution.csv", "analysis/observables.json",
             "analysis/snapshot.json", "analysis/site_hydration.csv", "source_model.json", "run/topology.json", "parameters/parameterization.json",
             "parameters/solute.prmtop", "parameters/solute.inpcrd"]
    for seed in SEEDS:
        files.extend([f"run/seed_{seed}.json", f"run/seed_{seed}.npz"])
    file_manifest = {name: file_hash(folder/name) for name in files}
    identity = digest({"source": source["artifact_id"], "parameter_id": params["parameter_id"],
                       "protocol": PROTOCOL, "run_ids": [m["run_id"] for m in seed_runs], "files": file_manifest})
    pilot_id = read(PILOT / "source_model.json")["structure"]["structure_id"]
    label = "PEtOx · H/H · DP 3" if s["structure_id"] == pilot_id else s["formula"]
    return {"schema": "polymer-water-study/1", "status": "completed", "study_id": identity,
            "request_key": request_key(s, source, temperature), "system": {
                "label": label, "structure_id": s["structure_id"], "canonical_smiles": s["canonical_smiles"],
                "temperature_K": temperature, "pressure_bar": run0["pressure_bar"], "dp": s["dp"],
                "formula": s["formula"], "water_count": seed_runs[0]["water_count"], "solvent": "water",
                "source_artifact_id": source["artifact_id"], "solute_mass_fraction": seed_runs[0]["solute_mass_fraction"]},
            "method": {**PROTOCOL, "engine": "OpenMM "+run0["engine_version"], "force_field": "GAFF2 + AM1-BCC",
                       "water_model": "TIP3P", "replicas": len(SEEDS), "parameter_id": params["parameter_id"],
                       "ensemble": "NPT", "cutoff_nm": run0["cutoff_nm"], "timestep_ps": run0["timestep_ps"]},
            "analysis": analysis, "viewer": read(folder/"analysis/snapshot.json"),
            "files_sha256": file_manifest,
            "downloads": [{"label": path.rsplit("/", 1)[-1], "url": f"/api/data/{identity}/{path}"} for path in files],
            "quality": {"image_contacts_at_cutoff": sum(m["sample_frames_with_self_image_within_cutoff"] for m in seed_runs),
                        "minimum_image_clearance_nm": min(m["minimum_self_image_distance_nm"] for m in seed_runs),
                        "replicas": [{k: m[k] for k in ("run_id", "mean_density_g_ml", "mean_temperature_K",
                                                         "mean_chain_molarity_mol_l")} for m in seed_runs],
                        "scope": "Kısa, sonlu kutuda model örneklemesi. Denge ve deneysel doğrulama ayrı sınanır."}}


def run_study(folder, source, temperature, cache, update):
    import shutil
    from parameters import parameterize
    from simulation import run_replicate
    from trajectory_analysis import analyze_replicas
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    write_record(folder/"source_model.json", source)
    update("GAFF2 atom tipleri ve AM1-BCC yükleri hazırlanıyor.")
    params = parameterize(source["structure"]["canonical_smiles"], Path(cache)/"parameter-cache", source_record=source["energy"])
    packet = Path(params["directory"])
    shutil.copytree(packet, folder/"parameters", dirs_exist_ok=True)
    runs = []
    for i, seed in enumerate(SEEDS):
        update(f"Bağımsız MD koşusu {i+1}/{len(SEEDS)} çalışıyor.")
        runs.append(run_replicate(folder/"parameters", folder/"run", seed=seed, temperature_K=temperature,
                                  pressure_bar=PROTOCOL["pressure_bar"], equil_ps=PROTOCOL["equil_ps"],
                                  production_ps=PROTOCOL["production_ps"], sample_ps=PROTOCOL["sample_ps"],
                                  geometry_record=source if i == 0 else None))
    update("Yörüngelerden zincir boyutu ve hidratasyon hesaplanıyor.")
    analyze_replicas(folder/"run/topology.json", [folder/"run"/m["trajectory"] for m in runs], folder/"analysis")
    result = assemble(folder, source, temperature, params, runs)
    write_record(folder/"study.json", result)
    return result
