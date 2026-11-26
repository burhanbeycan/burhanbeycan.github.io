# Contributing

Before adding a physical model, specify its equations, variable definitions, units, reference state, domain and primary sources. Add independent analytical, conservation or convergence checks. Preserve the separation of numerical verification and experimental validation.

Do not silently reinterpret material metadata as calibrated parameters. Do not compare coefficients from different reference frames without an explicit transformation. Changes to the exchange schema need a versioned compatibility decision.

The GUI and numerical model must remain independent. Keep calculations deterministic where appropriate and record random seeds for stochastic algorithms. Do not commit credentials, proprietary data or undocumented material constants.
