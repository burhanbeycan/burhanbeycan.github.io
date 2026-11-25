"""Local dashboard for explicit-water polymer MD. Start with python app.py."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import json
from pathlib import Path
import re
from threading import Lock
from urllib.parse import urlsplit

from study import PILOT, ROOT, encode, file_hash, get_source, inputs, read, request_key
from study import run_study, runtime_fingerprint, write_record

DEFAULT_STORE = Path.home() / ".polymer-research-studio" / "polymer-water"
HEX = r"[a-f0-9]{64}"


def archived(folder):
    data = read(folder/"study.json")
    if not all((folder/p).is_file() and file_hash(folder/p) == h for p,h in data["files_sha256"].items()):
        raise ValueError("Kaydedilmiş çalışmanın dosyaları değişmiş veya eksik.")
    return data


class Studies:
    def __init__(self, store):
        self.store = Path(store)
        self.store.mkdir(parents=True, exist_ok=True)
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.lock = Lock()
        self.jobs = {}

    def get_folder(self, study_id):
        if study_id == read(PILOT/"study.json")["study_id"]:
            return PILOT
        idx = self.store/"index"/study_id
        if idx.is_file():
            key = read(idx)
            if re.fullmatch(HEX,key):
                return self.store/"studies"/key
        raise KeyError(study_id)

    def submit(self, payload):
        structure, source, temperature = inputs(payload)
        if source is None:
            # Exact shipped conformer can be selected before an expensive new workflow starts.
            p = read(PILOT/"source_model.json")
            from upstream.energy_store import expected_model_id
            if (structure["recipe_id"] == p["structure"]["recipe_id"]
                    and expected_model_id(structure["canonical_smiles"]) == p["energy"]["model_id"]):
                source = p
        req = request_key(structure, source, temperature)
        pilot = read(PILOT/"study.json")
        if req == pilot["request_key"]:
            return dict(archived(PILOT), acquisition="bundled")
        from study import digest
        fingerprint = runtime_fingerprint()
        key = digest([req, fingerprint])
        alias = self.store/"aliases"/key
        if alias.is_file():
            return dict(archived(self.get_folder(read(alias))), acquisition="cache")
        folder = self.store/"studies"/key
        with self.lock:
            if (folder/"study.json").is_file():
                item = archived(folder)
                write_record(self.store/"index"/item["study_id"], key)
                return dict(item, acquisition="cache")
            active = self.jobs.get(key)
            if active and active["status"] in ("running", "queued"):
                return active.copy()
            initial = {"job_id": key, "status": "queued", "message": "Koşu kuyruğa alındı."}
            self.jobs[key] = initial
            self.executor.submit(self.execute, key, structure, source, temperature, folder, fingerprint)
            return initial.copy()

    def execute(self, key, structure, source, temperature, folder, fingerprint):
        def update(message):
            with self.lock:
                self.jobs[key] = {"job_id": key, "status": "running", "message": message}
        try:
            from simulation import _replica_lock
            folder.mkdir(parents=True, exist_ok=True)
            with _replica_lock(folder/"job.lock"):
                if (folder/"study.json").is_file():
                    item = archived(folder)
                else:
                    update("Sonlu zincirin başlangıç konformeri hazırlanıyor.")
                    source = source or get_source(structure, self.store)
                    item = run_study(folder, source, temperature, self.store, update)
                write_record(self.store/"index"/item["study_id"], key)
                from study import digest
                write_record(self.store/"aliases"/digest([item["request_key"], fingerprint]), item["study_id"])
            response = dict(item, acquisition="calculated")
        except Exception as error:
            response = {"status": "failed", "error": str(error) or type(error).__name__, "job_id": key}
        with self.lock:
            self.jobs[key] = response

    def job(self, key):
        with self.lock:
            value = self.jobs.get(key)
        if value is not None:
            return value
        folder = self.store/"studies"/key
        if (folder/"study.json").is_file():
            return dict(archived(folder), acquisition="cache")
        raise KeyError(key)


class Handler(BaseHTTPRequestHandler):
    def send(self, code, data, mime="application/json", filename=None):
        if mime == "application/json" and not isinstance(data, bytes):
            data = encode(data)
        self.send_response(code)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        if filename:
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        p = urlsplit(self.path).path
        try:
            if p in ("/", "/index.html"):
                return self.send(200,(ROOT/"ui.html").read_bytes(),"text/html; charset=utf-8")
            research_pages = {
                "/hydration/contact-viewer.html": "contact-viewer.html",
                "/hydration/comparison.html": "comparison.html",
            }
            if p in research_pages:
                return self.send(200,(ROOT/"hydration"/research_pages[p]).read_bytes(),"text/html; charset=utf-8")
            if match := re.fullmatch(f"/api/trajectory/({HEX})/([0-9]{{1,10}})",p):
                from hydration.player import trajectory_frames
                folder = self.server.studies.get_folder(match[1])
                return self.send(200,trajectory_frames(folder,read(folder/"study.json"),int(match[2])))
            if p == "/api/status":
                return self.send(200, {"version": "1.1.0", "engines": {k: bool(importlib.util.find_spec(k))
                                     for k in ("rdkit","numpy","openmm")}})
            if p == "/api/pilot":
                return self.send(200, dict(archived(PILOT), acquisition="bundled"))
            if match := re.fullmatch(f"/api/jobs/({HEX})",p):
                return self.send(200,self.server.studies.job(match[1]))
            if match := re.fullmatch(f"/api/data/({HEX})/([a-zA-Z0-9_./-]+)",p):
                folder = self.server.studies.get_folder(match[1])
                path = match[2]
                record = read(folder/"study.json")
                if path in record["files_sha256"]:
                    raw = (folder/path).read_bytes()
                    from hashlib import sha256
                    if sha256(raw).hexdigest() != record["files_sha256"][path]:
                        raise ValueError("Dosya özeti uyuşmuyor.")
                    mime = {".csv":"text/csv", ".json":"application/json"}.get(Path(path).suffix, "application/octet-stream")
                    return self.send(200,raw,mime,Path(path).name)
        except KeyError:
            pass
        except (OSError, ValueError) as error:
            return self.send(409,{"error":str(error)})
        return self.send(404,{"error":"Çalışma veya dosya bulunamadı."})

    def do_POST(self):
        if self.path != "/api/studies":
            return self.send(404,{"error":"İşlem bulunamadı."})
        try:
            origin = self.headers.get("Origin")
            if origin and urlsplit(origin).netloc != self.headers.get("Host"):
                return self.send(403,{"error":"Bu istek yerel arayüzden gönderilmelidir."})
            if self.headers.get_content_type() != "application/json":
                return self.send(415, {"error":"JSON veri gerekli."})
            length = int(self.headers.get("Content-Length") or 0)
            if not 0 < length <= 5*1024*1024:
                raise ValueError("İstek boyutu 5 MB sınırını aşıyor veya veri eksik.")
            reply = self.server.studies.submit(json.loads(self.rfile.read(length)))
        except (ValueError, ArithmeticError, TypeError, KeyError) as error:
            return self.send(400, {"error":str(error)})
        except (ImportError, OSError) as error:
            return self.send(503, {"error":"Yerel hesap bağımlılığı veya dosyası kullanılamadı. "+str(error)})
        return self.send(202 if reply["status"] in ("running","queued") else 200, reply)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--port", type=int, default=8766)
    ap.add_argument("--store", type=Path, default=DEFAULT_STORE)
    args = ap.parse_args()
    service = Studies(args.store)
    with ThreadingHTTPServer(("127.0.0.1",args.port),Handler) as server:
        server.studies = service
        print(f"Polymer–Solvent MD · http://127.0.0.1:{server.server_port}", flush=True)
        try:
            server.serve_forever()
        finally:
            service.executor.shutdown(wait=True)


if __name__ == "__main__":
    main()
