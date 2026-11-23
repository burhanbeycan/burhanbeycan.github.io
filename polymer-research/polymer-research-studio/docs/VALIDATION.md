# Verification record · 0.2.0

**Executed on 2026-11-23: 26 unittest tests passed, with all six laboratories available.** JavaScript syntax was checked with `node --check`.

From the installed repository:

```bash
python -m unittest discover -s tests -v
node --check src/polymer_workbench/static/app.js
```

The local development run made the seven `src` directories importable using `PYTHONPATH`; no remote CI success is implied.

## Checks performed

- 12 exchange-contract tests verify identity preservation, checksums, versioned records, explicit handoffs, population semantics and malformed input rejection.
- 4 molecular-view tests verify graph round trips, unique atom indices, bounded finite coordinates, valid bond endpoints/orders, repeat memberships, malformed signed imports and added chemistry metadata.
- 5 local HTTP tests verify calculation records, error responses, origin/host boundaries, request size and static assets.
- 5 integration tests exercise all six applications' defaults and reference scenarios, each independent launcher, Chemistry→Chain repeat/end-mass transfer, Network→Transport geometry and curve agreement, population round trips and the absence of an implicit electrolyte-tracer→penetrant-D mapping.

These verify software contracts and selected equations. They do not establish experimental prediction accuracy for a named polymer. Molecular formulas and masses are calculations on the supplied, valence-checked finite chemical graph; physical demonstration coefficients remain illustrative.

## Executed environment

Python 3.12.14, Linux.

| Dependency | Version |
|---|---|
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |
| FastAPI | 0.141.1 |
| Uvicorn | 0.52.1 |
| Pydantic | 2.13.4 |
| HTTPX | 0.28.1 |
| RDKit | 2026.9.2 |

Python 3.10 is the declared minimum; the local checks above executed on Python 3.12.14. GitHub Actions supplies additional configuration but its result must be inspected on the actual uploaded commit.

## Interface boundary

Static assets, graph payload validation and JavaScript syntax were checked. The SVG molecular renderer builds elements from validated numerical atom/bond data, offers carbon labels, atom indices, repeat colours, zoom, structure SVG and graph CSV exports. Drawing coordinates have arbitrary depiction units and carry no conformational or energy interpretation. Live desktop/mobile rendering and user interactions should be reviewed in a browser before a public research demonstration; these Python tests are not browser visual tests.

Each laboratory owns its model-specific evidence in its own `docs/VALIDATION.md`, scientific documentation and tests. Calibration and independent observations remain necessary for material-specific physical predictions.
