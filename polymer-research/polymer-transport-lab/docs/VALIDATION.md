# Current release verification · 0.3.0

All 33 tests in this project passed in the suite run on 24 November 2026. The [release report](../../VALIDATION_TR.md) lists the installed-wheel and interface checks. The following notes also retain the earlier release record where identified.

# Verification record · 0.2.0

**Executed on 2026-11-23: 24 unittest tests passed.**

Command, from this repository after installing its dependencies:

```bash
python -m unittest discover -s tests -v
```

## What was checked

Full-thickness slab limits and similarity; independent finite-volume convergence and integrated flux balance; fixed-geometry D recovery and uninformative-data diagnostics; Henry-law units; invalid inputs and finite exports.

The tests ran against the delivered code. They verify stated equations,
selected analytical/numerical limits, input handling and software contracts.
They do **not** establish experimental prediction accuracy for a named polymer.
Demonstration inputs remain synthetic and explicitly labelled.

## Current numerical review

Re-executed on 2026-11-23 with Python 3.12.14, NumPy 2.3.5 and SciPy
1.17.0 on Linux. All laboratory tests listed above passed. The scientific
kernels are unchanged by the 0.2.0 compatibility update. This review checked
the full-thickness, dry-reference and polymer-fixed conventions where relevant.

GUI and integrated-workflow verification belongs to the shared runtime's
current verification record. Remote CI status must be inspected after upload.
No remote run is implied here.

## Next material-specific evidence

Use traceable experimental inputs with units, conditions, sample distributions,
preparation history and model-compatible reference states. Compare predictions
against independent observations. Read SCIENCE.md and REFERENCES.md in the
owning laboratory before assigning material-specific meaning to a calculation.
