# Current release verification · 0.3.0

All 24 tests in this project passed in the suite run on 24 November 2026. The [release report](../../VALIDATION_TR.md) lists the installed-wheel and interface checks. The following notes also retain the earlier release record where identified.

# Numerical verification

Version 0.2.0: 14 unittest cases passed on 23 November 2026. Tests cover balanced Carothers behavior, the published 99/100 full-conversion case Xn=199, A/B relabeling, functional-group end count, mass conservation, p=0 and the singular balanced p=1 limit. Initiated-chain checks cover finite chain count, residue/end mass conservation, zero-addition chains and invalid efficiency or noninteger counts. Terminal-model checks cover Bernoulli/alternating limits, ordered-dyad normalization, stationarity, A/B relabeling and invalid feeds/ratios. Presets are checked for finite JSON exports; all five also passed the runtime output-contract validator.

These are analytical, conservation and software checks. They are not wet-lab validation or prediction-error estimates. The original balance-model release used analytical checks only. GitHub Actions configuration is supplied; remote CI success must be checked in the actual repository after publication.

The PEtOx extension adds separate experimental cohorts and tests for terminal mass, Arrhenius exposure, zero-time behavior, grouped predictions, leakage, method matching, missing conversions and numerical agreement with the reported cohort scores. See [the PEtOx results](PETOX.md).
