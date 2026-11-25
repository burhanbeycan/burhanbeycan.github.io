"""Serve the local Polymer Structure & Energy application or compute one record."""
import argparse
from contextlib import closing
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import re
import sqlite3
from urllib.parse import urlsplit

from rdkit import rdBase
from polymer_identity import build_structure
from solvent_registry import analyze_with_registry, get_component, load_registry

ROOT = Path(__file__).resolve().parent
VERSION = "1.0.0"
DEFAULT_DB = Path.home() / ".polymer-research-studio" / "structure-energy.sqlite3"
ARTIFACT_SCHEMA = "polymer-structure-energy/1"
ASSETS = {"editor.js": "text/javascript", "editor.css": "text/css", "ocl-vendor.js": "text/javascript"}


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False, ensure_ascii=False)


def connect(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(path, timeout=30)
    db.execute("CREATE TABLE IF NOT EXISTS calculation_records (id TEXT PRIMARY KEY, content TEXT NOT NULL)")
    return db


def calculate(spec, database, *, count=12, seed=2403):
    structure = build_structure(spec)
    Path(database).parent.mkdir(parents=True, exist_ok=True)
    energy = analyze_with_registry(structure["canonical_smiles"], database, count=count, seed=seed)
    # Acquisition (new run, local cache, shipped calculation) is request-specific.
    # Keep the portable scientific record identical on a cache hit or a fresh run.
    result = {k: v for k, v in energy.items() if k not in {"cache_hit", "cache_source"}}
    record = {"schema": ARTIFACT_SCHEMA, "structure": structure, "energy": result}
    artifact_id = hashlib.sha256(encode(record).encode()).hexdigest()
    record["artifact_id"] = artifact_id
    with closing(connect(database)) as db, db:
        db.execute("INSERT OR IGNORE INTO calculation_records(id,content) VALUES (?,?)", (artifact_id, encode(record)))
    return dict(record, energy=energy)


def find_record(database, artifact_id):
    with closing(connect(database)) as db:
        row = db.execute("SELECT content FROM calculation_records WHERE id=?", (artifact_id,)).fetchone()
    return json.loads(row[0]) if row else None


class AppHandler(BaseHTTPRequestHandler):
    server_version = "PolymerStructureEnergy/" + VERSION

    def emit(self, status, data, content_type="application/json"):
        if content_type == "application/json":
            data = encode(data).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type + "; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        p = urlsplit(self.path).path
        if p in {"/", "/index.html"}:
            return self.emit(200, (ROOT / "ui.html").read_bytes(), "text/html")
        if p.startswith("/assets/") and p[8:] in ASSETS:
            asset = p[8:]
            return self.emit(200, (ROOT / "assets" / asset).read_bytes(), ASSETS[asset])
        if p == "/api/status":
            return self.emit(200, {"status": "ready", "version": VERSION, "rdkit_version": rdBase.rdkitVersion,
                                   "default_protocol": {"count": 12, "seed": 2403, "method": "ETKDGv3 + MMFF94s"}})
        if p == "/api/liquid-reference":
            reference = ROOT / "pure_liquid" / "water_npt.json"
            if reference.is_file():
                return self.emit(200, json.loads(reference.read_text("utf-8")))
            return self.emit(404, {"error": "Saf sıvı için tamamlanmış bir pilot kayıt bulunamadı."})
        if p == "/api/solvents":
            return self.emit(200, load_registry())
        if p.startswith("/api/solvents/"):
            try:
                return self.emit(200, get_component(p.split("/")[-1]))
            except KeyError:
                pass
        if re.fullmatch(r"/api/records/[0-9a-f]{64}", p):
            item = find_record(self.server.database, p.rsplit("/", 1)[1])
            if item:
                return self.emit(200, item)
        return self.emit(404, {"error": "Kayıt veya adres bulunamadı."})

    def do_POST(self):
        p = urlsplit(self.path).path
        if p not in {"/api/analyze", "/api/structure"}:
            return self.emit(404, {"error": "İşlem bulunamadı."})
        try:
            length = int(self.headers.get("Content-Length") or 0)
            if not 0 < length <= 48000:
                raise ValueError("Geçerli bir yapı tanımı girin.")
            spec = json.loads(self.rfile.read(length))
            if not isinstance(spec, dict):
                raise ValueError("İstek bir yapı kaydı içermelidir.")
            if p == "/api/structure":
                result = build_structure(spec)
            else:
                result = calculate(spec, self.server.database)
        except (ValueError, TypeError, OverflowError) as error:
            return self.emit(400, {"error": str(error)})
        except (RuntimeError, sqlite3.Error) as error:
            return self.emit(500, {"error": "Hesap veya kayıt tamamlanamadı.", "detail": str(error)})
        return self.emit(200, result)


def main():
    parser = argparse.ArgumentParser(description="Polimer yapısı ve konformer enerjisi")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--database", default=DEFAULT_DB, type=Path)
    parser.add_argument("--spec", type=Path, help="JSON yapı tanımından bir sonuç üret")
    args = parser.parse_args()
    if args.spec:
        print(encode(calculate(json.loads(args.spec.read_text("utf-8")), args.database)))
        return
    with ThreadingHTTPServer(("127.0.0.1", args.port), AppHandler) as server:
        server.database = args.database
        print(f"Polymer Structure & Energy · http://127.0.0.1:{server.server_port}", flush=True)
        server.serve_forever()


if __name__ == "__main__":
    main()
