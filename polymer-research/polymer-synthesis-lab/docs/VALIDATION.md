# Numerical verification

Version 0.2.0: 14 unittest cases passed on 23 November 2026. Tests cover balanced Carothers behavior, the published 99/100 full-conversion case Xn=199, A/B relabeling, functional-group end count, mass conservation, p=0 and the singular balanced p=1 limit. Initiated-chain checks cover finite chain count, residue/end mass conservation, zero-addition chains and invalid efficiency or noninteger counts. Terminal-model checks cover Bernoulli/alternating limits, ordered-dyad normalization, stationarity, A/B relabeling and invalid feeds/ratios. Presets are checked for finite JSON exports; all five also passed the runtime output-contract validator.

These are analytical, conservation and software checks. They are not wet-lab validation or prediction-error estimates. No material-specific kinetic training dataset is bundled. GitHub Actions configuration is supplied; remote CI success must be checked in the actual repository after publication.
