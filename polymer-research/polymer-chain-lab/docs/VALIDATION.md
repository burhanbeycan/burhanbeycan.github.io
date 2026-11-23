# Verification record · 0.2.0

**Executed on 2026-11-23: 51 unittest tests passed.**

Command, from this repository after installing its dependencies:

```bash
python -m unittest discover -s tests -v
```

## What was checked

End-corrected chemical mass balance, terminal-mass rejection, population-moment invariance, common terminal correction of weighted dimensions, explicit unsupported-architecture rejection; population normalization and moments; independent scattering weights; WLC rigid/long-chain limits and numerical integral; exact 3^4 RIS enumeration; chain geometry identities; seeded sampling; malformed inputs and finite result contracts.

The tests ran against the delivered code. They verify stated equations,
selected analytical/numerical limits, input handling and software contracts.
They do **not** establish experimental prediction accuracy for a named polymer.
Demonstration inputs remain synthetic and explicitly labelled.

## Observed environment

Python 3.12.14 on Linux. These are the versions actually used in this
run, not a claim of exhaustive compatibility across all permitted versions.

| Dependency | Tested version |
|---|---|
| numpy | 2.3.5 |
| scipy | 1.17.0 |
| fastapi | 0.141.1 |
| uvicorn | 0.52.1 |
| pydantic | 2.13.4 |
| httpx | 0.28.1 |

Declared compatibility starts at Python 3.10; only Python 3.12.14 was
executed locally here. GitHub Actions is configured for additional checks but
remote CI status must be inspected after upload. No remote CI run is implied
by this record.

## Interface verification boundary

The shared runtime's local HTTP endpoints, static asset responses, result
contracts, default calculations and every synthetic reference scenario passed.
JavaScript syntax was checked with `node --check`. Browser rendering and visual
interaction review could not be completed in the build environment: the local
browser process was restricted and the cloud browser could not reach the local
server. Responsive CSS and the full GUI implementation are included; desktop,
mobile, keyboard and export interactions still require a live browser review.

## Next material-specific evidence

Use traceable experimental inputs with units, conditions, sample distributions,
preparation history and model-compatible reference states. Compare predictions
against independent observations. Read SCIENCE.md and REFERENCES.md in the
owning laboratory before assigning material-specific meaning to a calculation.
