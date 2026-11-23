"""Loopback-only research GUI and API. No external services or credentials."""
from __future__ import annotations

import argparse
import importlib
import importlib.util
from pathlib import Path
from typing import Any
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from starlette.concurrency import run_in_threadpool
from starlette.middleware.trustedhost import TrustedHostMiddleware

from . import __version__
from .exchange import MAX_JSON_BYTES, default_material, ensure_json, handoff, make_result, validate_exchange, validate_material

APPS = {
    "chemistry": "polymer_chemistry_lab.plugin",
    "synthesis": "polymer_synthesis_lab.plugin",
    "chain": "polymer_chain_lab.plugin",
    "transport": "polymer_transport_lab.plugin",
    "network": "polymer_network_lab.plugin",
    "electrolyte": "polymer_electrolyte_lab.plugin",
}


def discover(selected: str | None = None) -> dict:
    found = {}
    for app_id, module in APPS.items():
        if selected and app_id != selected:
            continue
        try:
            available = importlib.util.find_spec(module.split(".")[0])
        except ModuleNotFoundError:
            available = None
        if available:
            found[app_id] = importlib.import_module(module)
    if selected and selected not in found:
        raise ValueError(f"{selected} is not installed. Install its project in this environment.")
    return found


class Calculation(BaseModel):
    inputs: dict[str, Any] = Field(default_factory=dict)
    material: dict[str, Any] | None = None


class Transfer(BaseModel):
    document: dict[str, Any]
    target: str


def create_app(selected: str | None = None) -> FastAPI:
    plugins = discover(selected)
    app = FastAPI(title="Polymer Research Studio", version=__version__)
    app.add_middleware(TrustedHostMiddleware, allowed_hosts=["127.0.0.1", "localhost", "[::1]", "testserver"])

    @app.middleware("http")
    async def boundary(request: Request, call_next):
        # A local app must not silently become a cross-origin computation service.
        origin = request.headers.get("origin")
        if origin and origin != str(request.base_url).rstrip("/"):
            return JSONResponse({"detail": "Cross-origin requests are disabled."}, status_code=403)
        if request.method == "POST":
            body = await request.body()
            if len(body) > MAX_JSON_BYTES:
                return JSONResponse({"detail": "Request exceeds 4 MB."}, status_code=413)
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["Content-Security-Policy"] = "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; object-src 'none'; frame-ancestors 'none'; base-uri 'self'"
        response.headers["Cache-Control"] = "no-store"
        return response

    @app.exception_handler(ValueError)
    async def invalid(request, exc):
        return JSONResponse({"detail": str(exc)}, status_code=422)

    @app.get("/api/catalog")
    def catalog():
        return {"version": __version__, "mode": "standalone" if selected else "integrated", "selected": selected, "apps": [x.APP for x in plugins.values()], "material": default_material()}

    @app.get("/api/health")
    def health():
        return {"status": "ok", "version": __version__, "installed": list(plugins)}

    @app.post("/api/calculate/{app_id}")
    async def calculate(app_id: str, request: Calculation):
        if app_id not in plugins:
            raise HTTPException(404, "Application is not installed in this environment.")
        plugin = plugins[app_id]
        ensure_json(request.inputs)
        allowed = {f["key"] for f in plugin.APP["fields"]}
        unknown = set(request.inputs) - allowed
        if unknown:
            raise ValueError("Unknown model inputs: " + ", ".join(sorted(unknown)))
        inputs = {f["key"]: f.get("default") for f in plugin.APP["fields"]}
        inputs.update(request.inputs)
        material = validate_material(request.material)
        output = await run_in_threadpool(plugin.calculate, inputs, material)
        return make_result(plugin.APP, inputs, output, material)

    @app.post("/api/exchange/validate")
    def validate(document: dict):
        return validate_exchange(document)

    @app.post("/api/exchange/transfer")
    def transfer(request: Transfer):
        if request.target not in plugins:
            raise HTTPException(404, "Target application is not installed.")
        return handoff(request.document, request.target, plugins[request.target].APP["fields"])

    static = Path(__file__).parent / "static"
    app.mount("/static", StaticFiles(directory=static), name="static")

    @app.get("/")
    def index():
        return FileResponse(static / "index.html")

    return app


def main(selected: str | None = None):
    parser = argparse.ArgumentParser(description="Local scientific polymer applications")
    parser.add_argument("--port", type=int, default=8765 if selected is None else {"chain": 8766, "transport": 8767, "network": 8768, "electrolyte": 8769, "synthesis": 8770, "chemistry": 8771}[selected])
    parser.add_argument("--no-browser", action="store_true", help="Do not open the local GUI automatically.")
    args = parser.parse_args()
    if not 1024 <= args.port <= 65535:
        parser.error("Choose a port between 1024 and 65535.")
    application = create_app(selected)
    if not args.no_browser:
        import threading
        import webbrowser
        threading.Timer(1.0, lambda: webbrowser.open(f"http://127.0.0.1:{args.port}")).start()
    import uvicorn
    uvicorn.run(application, host="127.0.0.1", port=args.port, log_level="info")

