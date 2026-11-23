# Current release verification · 0.3.0

All 34 tests in this project passed in the suite run on 24 November 2026. The [release report](../../VALIDATION_TR.md) lists the installed-wheel and interface checks. The following notes also retain the earlier release record where identified.

# Verification and measured benchmarks

Run `python -m unittest discover -s tests -v` after installation or run `PYTHONPATH=src python -m unittest discover -s tests -v` from this project. All 34 tests passed locally for this revision.

## Material-specific checks

| Check | Result | Evidence |
|---|---|---|
| Pesko 5 kg/mol mixture concentration | 1.58765 mol/L from the measured mixture density and repeat-based composition; rounds to published 1.59 | Table I |
| Pesko measured conductivity versus tracer baseline | 2 mS/cm (rounded EIS) vs ≈1.84 mS/cm (calculated); estimated figure-reading uncertainty propagated | Text with Fig. 2a and Fig. 4 |
| Zheng 100 kg/mol PEO | Computed tracer fraction 0.192 rounds to published 0.19 | Independent Table 3 |
| Zheng volumetric comparison | Suppressed without density or molarity for this specimen | Cross-study specimen guard |
| Weighed illustrative film | 100 µm, 1 cm² and 100 Ω give 0.10 mS/cm | Unit conversion only |

The tests check analytical units, recovered source quantities, sample and temperature matching, missing density, phase metadata, optional end-group correction and uncertainty propagation. These checks do not establish a chemistry-to-conductivity predictor. See the [worked comparison](PEO_COMPARISON.md) for provenance and replication limits.

## Numerical checks

The generic Nernst–Einstein and Onsager kernels are checked against independently assembled species currents and number-density expressions, positive-semidefinite boundary conditions, singular matrices, negative transference, coefficient scaling and finite output contracts. The GUI uses the same calculation entry points as the tests.
