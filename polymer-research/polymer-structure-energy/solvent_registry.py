"""Versioned, locally calculated molecular starting geometries for common solvents.

Read the bundled registry without running a simulation. Regenerate it with
python solvent_registry.py --write solvent_components.json [--db models.sqlite3].
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

REGISTRY_PATH = Path(__file__).with_name("solvent_components.json")
SCHEMA_VERSION = "solvent-components/1"
COUNT = 12
SEED = 2403

# Identifiers curated in this project. All formulae, masses, coordinates and
# energies are produced by the same model code used for an entered molecule.
MOLECULES = (
    ("water", "Su", "Water", "O"),
    ("ethanol", "Etanol", "Ethanol", "CCO"),
    ("acetonitrile", "Asetonitril", "Acetonitrile", "CC#N"),
    ("acetone", "Aseton", "Acetone", "CC(=O)C"),
    ("thf", "Tetrahidrofuran (THF)", "Tetrahydrofuran", "C1CCOC1"),
    ("dmf", "N,N-Dimetilformamid (DMF)", "N,N-Dimethylformamide", "CN(C)C=O"),
    ("dmso", "Dimetil sülfoksit (DMSO)", "Dimethyl sulfoxide", "CS(C)=O"),
)


def load_registry(path: str | Path = REGISTRY_PATH) -> dict:
    """Load stored model outputs. Coordinates and energies require no recomputation."""
    with Path(path).open(encoding="utf-8") as stream:
        document = json.load(stream)
    if document.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Desteklenmeyen çözücü kayıt şeması.")
    return document


def get_component(component_id: str, registry: dict | None = None) -> dict:
    """Resolve a selected solvent, including its geometry and calculation provenance."""
    registry = load_registry() if registry is None else registry
    for component in registry["components"]:
        if component["id"] == component_id:
            return component
    raise KeyError(component_id)


def lookup_energy(smiles: str, *, count: int = COUNT, seed: int = SEED, registry: dict | None = None) -> dict | None:
    """Reuse a stored result only when the entire installed calculation identity matches."""
    from energy_store import expected_model_id

    key = expected_model_id(smiles, count=count, seed=seed)
    registry = load_registry() if registry is None else registry
    for component in registry["components"]:
        record = component["energy_record"]
        if key == component["model_id"] == record["model_id"]:
            result = deepcopy(record)
            result.update(cache_hit=True, cache_source="bundled_solvent_registry")
            return result
    return None


def analyze_with_registry(smiles: str, db_path: str | Path, *, count: int = COUNT, seed: int = SEED) -> dict:
    """Check shipped results, then use the usual SQLite-backed analysis."""
    from energy_store import analyze

    stored = lookup_energy(smiles, count=count, seed=seed)
    if stored is not None:
        return stored
    result = analyze(smiles, db_path, count=count, seed=seed)
    return {**result, "cache_source": "sqlite" if result["cache_hit"] else "computed"}


def build_registry(db_path: str | Path) -> dict:
    """Generate numeric records through the regular identity and energy pipeline."""
    from energy_store import analyze
    from polymer_identity import build_structure

    components = []
    for component_id, name, name_en, smiles in MOLECULES:
        identity = build_structure({"mode": "finite", "smiles": smiles})
        result = analyze(identity["canonical_smiles"], db_path, count=COUNT, seed=SEED)
        if result["status"] not in ("ok", "partial") or not result["conformers"]:
            raise RuntimeError(f"{name} için yakınsamış konformer bulunamadı: {result['message']}")
        if not result["lowest_geometry"] or result["canonical_smiles"] != identity["canonical_smiles"]:
            raise RuntimeError(f"{name} geometri veya kimyasal kimlik doğrulaması başarısız.")
        lowest = min(record["energy_kcal_mol"] for record in result["conformers"])
        components.append(
            {
                "id": component_id,
                "name": name,
                "name_en": name_en,
                "smiles": smiles,
                "canonical_smiles": identity["canonical_smiles"],
                "structure_id": identity["structure_id"],
                "formula": identity["formula"],
                "molar_mass_g_mol": identity["molar_mass"],
                "net_charge": identity["net_charge"],
                "model_id": result["model_id"],
                "energy_record": {key: value for key, value in result.items() if key not in ("cache_hit", "cache_source")},
                "energy_schema": result["schema_version"],
                "status": result["status"],
                "message": result["message"],
                "requested": result["requested"],
                "embedded": result["embedded"],
                "converged": result["converged"],
                "energy_unit": "kcal/mol",
                "energy_reference": "Bu tek molekülde aynı protokolle örneklenen konformerler.",
                "lowest_energy_kcal_mol": lowest,
                "geometry": result["lowest_geometry"],
                "conformers": result["conformers"],
                "protocol": result["protocol"],
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "scope": "Tekil çözücü molekülünün kimyasal kimliği ve izole konformer hesapları.",
        "source": {
            "smiles": "Proje içinde tanımlanmış moleküler yapılar.",
            "numerical_results": "polymer_identity.build_structure + energy_store.analyze",
            "conformer_count_requested": COUNT,
            "seed": SEED,
        },
        "protocol_use": (
            "Kimyasal kimlik ve model geometrisi yeniden kullanılabilir. "
            "MMFF94s enerjileri aynı molekülün aynı protokoldeki konformerlerini karşılaştırır. "
            "Sıvı solvent kutusu ve polimer–çözücü karışımı kendi kuvvet alanı, "
            "yük modeli, bileşimi, sıcaklığı, basıncı ve dengelemesiyle hesaplanır."
        ),
        "water_model_note": (
            "Su kartı tek bir H₂O molekülü içerir. Su içeren MD için polimer kuvvet alanıyla "
            "uyumlu su topolojisi ve o su modelinin geometri ve yükleri kullanılır."
        ),
        "coordinate_unit": "Å",
        "components": components,
    }


def write_registry(document: dict, target: Path) -> None:
    """Write the complete export atomically after successful model generation."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(mode="w", encoding="utf-8", dir=target.parent, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(document, stream, ensure_ascii=False, allow_nan=False, indent=2)
            stream.write("\n")
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", required=True, type=Path, metavar="JSON")
    parser.add_argument("--db", type=Path, help="İsteğe bağlı kalıcı model önbelleği")
    args = parser.parse_args()
    if args.db and args.db.resolve() == args.write.resolve():
        parser.error("JSON dosyası ve SQLite veritabanı için farklı yollar seçin.")
    if args.db:
        document = build_registry(args.db)
    else:
        with TemporaryDirectory() as temp:
            document = build_registry(Path(temp) / "models.sqlite3")
    write_registry(document, args.write)
    print(f"{len(document['components'])} çözücü kaydı yazıldı.")


if __name__ == "__main__":
    main()
