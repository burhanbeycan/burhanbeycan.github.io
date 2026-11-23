# Scientific model and scope

## Implemented system

Isothermal, homogeneous, electroneutral binary **1:1 monovalent salt** in a neutral polymer matrix. The velocity reference is the **polymer-fixed frame** throughout. User-supplied transport coefficients must refer to the same composition, temperature, normalization and frame. A solvent-fixed frame is equivalent only when that solvent has the same mean velocity as the polymer. Independent solvent motion is not modelled here. The application does not convert between barycentric, volume-fixed, solvent-fixed or laboratory frames.

Architecture, molar mass, tacticity, sequence, crystallinity, water content and end groups are scientifically important metadata. This release preserves them but **does not infer physical coefficients from them**. No force field, molecular simulation or learned material-parameter database is included.

## Ideal baseline

For c+ = c− = c, ion charges z+ = +1 and z− = −1:

\[\sigma_{NE}=\frac{F^2c}{RT}(D_++D_-),\qquad t_{+,ideal}=\frac{D_+}{D_++D_-}.\]

The input salt concentration is mol/L; internal c is mol/m³ (multiply by 1000). D is tracer/self diffusion in m²/s. This approximation omits ion-ion correlations. It does not supply a salt chemical diffusivity. The concentration plot holds D fixed and is a conditional sensitivity plot, not a concentration-dependent material prediction.

The exact SI defining constants are N_A = 6.02214076e23 mol⁻¹, e = 1.602176634e-19 C and k_B = 1.380649e-23 J/K. F = N_A e and R = N_A k_B are computed from them. 1 S/m = 10 mS/cm.

## Collective Onsager model

Define molar fluxes and molar electrochemical potentials by

\[J_i=-\sum_j L_{ij}\nabla\widetilde\mu_j.\]

J is mol/(m² s), mu is J/mol, and **L is mol²/(J m s)**. The input 2×2 matrix is symmetric. Nonnegative entropy production requires it to be positive semidefinite: L++ ≥ 0, L−− ≥ 0 and |L+−| ≤ sqrt(L++ L−−). Invalid matrices are rejected.

At the PSD boundary the square-root calculation can differ by a few floating-point ULPs. The implementation therefore permits a relative rounding tolerance of `8 × machine epsilon`, with **no dimensional absolute tolerance**. This is a numerical tolerance, not uncertainty in the physical coefficients. If either diagonal is exactly zero, the cross coefficient must be exactly zero. A reported boundary eigenvalue can differ slightly from zero due to floating-point rounding.

Under a purely electrical force, with no chemical-potential gradient:

\[\sigma=F^2(L_{++}+L_{--}-2L_{+-}),\]
\[t_+=\frac{L_{++}-L_{+-}}{L_{++}+L_{--}-2L_{+-}}.\]

The denominator is nonnegative for admissible matrices. A zero charge-conduction mode makes t+ undefined; the app returns null, not zero. The conditioning decision uses the matrix's own scale: a charge response at most `10⁻¹⁴` of the sum of its diagonal coefficients is treated as numerically unresolved for t+. Uniformly scaling a well-conditioned matrix does not make t+ undefined. A current fraction outside [0,1] is not automatically an error. The current fraction reported here is not automatically equivalent to an experimentally inferred Bruce–Vincent value.

**Collective diagonal coefficients include same-species distinct correlations.** Their values cannot generally be reconstructed from tracer D alone. The NE limit has L++ = cD+/(RT), L−− = cD−/(RT), L+− = 0. In collective mode users supply all three independently. A sweep of L+−/sqrt(L++L−−) explores mathematically admissible matrices with the diagonals fixed; it does not assert a material mechanism or infer cluster populations.

If a diagonal coefficient is zero, the cross-correlation ratio is undefined and its sensitivity graph is omitted. The permitted scalar conductivity and current-fraction results remain available. A singular matrix is not necessarily nonconducting: only a zero response in the charge direction makes conductivity vanish.

## Demonstrations and numerical verification

Every bundled preset is synthetic. The negative-transference preset uses a positive-definite matrix to demonstrate a mathematically admissible outcome, not a measurement on PEO or any other named material.

Tests compare NE to an independent number-density/elementary-charge expression, recover the ideal limit from Onsager coefficients, and compare t+ with independently assembled species fluxes and currents. They also check positive quadratic dissipation, rank-one PSD boundaries over different scales, genuine PSD violations at small coefficient magnitudes, zero-diagonal behavior, scale invariance of t+, singularity handling, SI conversion, strict finite payloads and invalid inputs. These are numerical/model checks, not experimental validation.

## Computational ranges

Backend validation enforces the same numeric ranges as the form, including fields that are currently hidden: temperature 1–5000 K, salt concentration `10⁻⁹–100 mol/L`, each self-diffusivity `10⁻³⁰–10⁻⁴ m²/s`, diagonal L values 0–1 and cross L values −1–1 in the stated Onsager units. These are computational bounds, not a claim that a polymer electrolyte physically exists or that the model is valid at every point in these ranges. The concentration sensitivity graph remains inside the supported concentration range. NaN, infinity, booleans, unrepresentably large numbers, unknown input keys and malformed input containers are rejected with readable errors before calculation.

## Interoperability

Material identity and lineage may travel to any laboratory. This app intentionally provides no coefficient handoff to the Fickian Transport Lab: self-diffusion, collective mobility and chemical diffusion are different quantities. No predicted cell performance, battery lifetime, safety, dendrite behavior or electrochemical stability is claimed.
