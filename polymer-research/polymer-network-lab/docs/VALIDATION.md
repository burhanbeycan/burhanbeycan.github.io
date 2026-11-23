# Current release verification · 0.4.0

42 tests passed in the local suite run on 24 November 2026. See the [release verification report](../../VALIDATION_TR.md) and the [worked scientific examples](../INCELEME.md).

# Verification record · 0.2.0

**Executed on 2026-11-23: 21 unittest tests passed.**

Command, from this repository after installing its dependencies:

```bash
python -m unittest discover -s tests -v
```

## What was checked

Consistent free-energy derivatives, curvature and stable root selection; dry/dilute limits and dimensionless collapse; the published-model numerical equilibrium example; independent image-series release comparison; handoff units and provenance.

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

## WBPU formulation checks

The WBPU tests check the twentyfold difference in ideal PEG 20k and PCL 1k diol OH equivalents per gram, independent assay conversions, amine/water accounting, operational gel/sol limits, frequency-specific DMA grouping and suppression of incompatible handoffs. Tests confirm that a held-out replicate's response cannot change its own predictions. Reference measurements, source table mapping, extrapolated folds and the interpretation of the three experimental errors are documented in [WBPU](WBPU.md).
