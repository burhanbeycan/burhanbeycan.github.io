# Scientific scope and numerical verification

## What version 0.2.0 calculates

Polymer Chain Lab is a reproducible calculation tool for three explicitly separated models. Its numerical kernels are independent of the workbench and do not call RDKit, RadonPy, PolyChemMD, an atomistic force field, or an external simulation engine.

| Mode | Implemented observables | Required physical information |
| --- | --- | --- |
| Gaussian population | Molar-mass moments and tails; number-, mass-, and intensity-weighted chain dimensions; normalized Debye mixture | Molar-mass distribution, common repeat-unit molar mass, contour increment, Kuhn length |
| Worm-like population | The same distribution analysis; exact continuous WLC mean-square dimensions across the rod-to-coil crossover | Molar-mass distribution, contour mapping, persistence length |
| Finite-chain RIS | Exact torsional-state and neighboring-pair marginals; independent 3D conformational samples; mean-square dimensions and sampling standard errors | Bond count, fixed geometry, three-state energies, opposite-gauche pair energy, temperature |

There is no chemistry-specific fitted parameter database. All shipped presets are **synthetic demonstrations**. Naming a polymer, sequence, tacticity, architecture, end group, or solvent in a material record does not supply a validated Hamiltonian. WLC scattering, excluded volume, branching, chain dynamics, rheology, phase equilibria, parameter fitting, and prediction of transport coefficients are outside this release.

## 1. Molar-mass distributions and experimental weighting

### Accepted input

The two-column CSV contains molar mass in g/mol and a nonnegative amount. An optional `molar_mass_g_mol,amount` header, blank lines, and `#` comments are accepted; tab-separated tables also work. At most 512 rows and 200,000 characters are accepted. Discrete duplicate masses are combined; zero-weight components are discarded. A distribution must retain at least one positive amount.

Four amount conventions are available:

1. **Discrete number amounts:** counts or number fractions; normalization gives `n_i`.
2. **Discrete mass amounts:** masses or mass fractions; `n_i ∝ amount_i/M_i`.
3. **Number density per log10 M:** samples of a density with respect to `x=log10(M)`.
4. **Mass density per log10 M:** samples of a mass density with respect to the same coordinate.

For a log-density curve, sorted masses must be distinct. Trapezoidal node widths are `h_0=(x_1−x_0)/2`, `h_last=(x_last−x_previous)/2`, and `h_i=(x_(i+1)−x_(i−1))/2` internally. The effective amount is `density_i h_i`, followed by the selected number/mass conversion. These are quadrature nodes, not inferred physically discrete species. The supplied interval is the whole integration range; no missing tail is reconstructed. Refine the grid and extend the measured range to assess convergence, especially for Mz and scattering.

Raw SEC retention times or uncorrected detector signals are not molar masses or universally valid mass densities. Calibration, concentration-response corrections, and composition-dependent detector sensitivity must be established before import. Amounts per unit M and amounts per unit log M are different; the app does not silently convert between these density conventions.

An imported material MWD can be selected explicitly. It must contain `basis: "number_fraction"`, equal-length `molar_mass_g_mol` and `number_fraction` arrays. Missing or incompatible records raise an error. Imported records are not considered experimentally verified.

### Moments

For normalized number fractions `Σ_i n_i=1`, the implemented moments are

\[
M_n=\sum_i n_i M_i,\qquad
M_w=\frac{\sum_i n_i M_i^2}{\sum_i n_i M_i},\qquad
M_z=\frac{\sum_i n_i M_i^3}{\sum_i n_i M_i^2},\qquad
Đ=\frac{M_w}{M_n}.
\]

`Mw` here denotes mass-average molar mass, also written `Mm` in IUPAC terminology. These definitions and the use of *dispersity* follow the [IUPAC number-average molar mass entry](https://goldbook.iupac.org/terms/view/12216) and [molar-mass dispersity entry](https://goldbook.iupac.org/terms/view/12224). Dimensionless relative molecular mass and dimensional molar mass must not be conflated.

Mass and zero-angle coherent intensity weights are

\[
w_i=\frac{n_iM_i}{\sum_j n_jM_j},\qquad
s_i=\frac{n_iM_i^2}{\sum_j n_jM_j^2}.
\]

The latter follows because single-chain amplitude at zero angle is proportional to mass when contrast per unit mass is common, and the independent-chain intensity is proportional to squared amplitude. It assumes a common composition/contrast and no interchain structure-factor contribution. It is not a universal weight for arbitrary copolymers, aggregation, mixtures with different scattering-length densities, or concentrated interacting solutions.

Every reported RMS radius averages the **squared** observable before taking its square root: for example `Rg,number = sqrt(Σ n_i <Rg_i²>)`. It is not `Σ n_i sqrt(<Rg_i²>)`. The intensity-weighted RMS radius is the Guinier radius of the independent-chain mixture under the stated contrast assumptions.

Fractions above the user threshold use the strict criterion `M > threshold`. The graph shows component/quadrature-bin fractions, not a density per unit mass. Lines joining bins are guides.

### Why the two default populations matter

| Synthetic population | M values (g/mol) | Number fractions | Mn | Mw | Mz | Đ |
| --- | --- | --- | --- | --- | --- | --- |
| A | 20,000; 80,000 | 0.5; 0.5 | 50,000 | 68,000 | 76,470.588235 | 1.36 |
| B | 10,000; 50,000; 250,000 | 0.09375; 0.8875; 0.01875 | 50,000 | 68,000 | 118,823.529412 | 1.36 |

These are exact finite populations. They demonstrate that matching Mn and dispersity does not identify the high-mass tail or intensity-weighted dimensions. They do not represent measured polymer samples.

## 2. Mapping molar mass to contour length

For the population models,

\[
DP_i=(M_i-M_{\mathrm{end}})/M_0,\qquad L_i=DP_i a,
\]

where `M0` is the molar-mass contribution of a constitutional repeat, `Mend` is the common **net terminal mass correction**, and `a` is contour increment per repeat. This is the rearrangement of the chemical mass balance `M=DP M0+Mend`; it does not infer end-group identity. The [IUPAC degree-of-polymerization definition](https://goldbook.iupac.org/terms/view/D01569) counts monomeric units; terminal groups therefore do not become extra repeating units. The chosen monomeric-unit versus constitutional-repeat convention must be stated, especially for condensation polymers.

`Mend` counts only atoms not already included in `DP` copies of the chosen repeat contribution. For example an H-terminated polymethylene representation built from `−CH2−CH2−` contributions uses a net `H2` correction; the molar mass of ethane is not the repeat mass. This module accepts a nonnegative common correction, defaults explicitly to zero, and requires every positive-weight component to have `M>Mend`. Choosing zero neglects terminal mass; it is not a claim that real chains have no termini. Terminal length is not added to L, so end-dominated oligomers still require discrete chemical geometry.

Fractional **effective** DP is retained for continuum distribution bins, not rounded to a falsely precise integer structure. The number-average effective DP is reported, as is each component's DP. Molar-mass moments, mass fractions and intensity fractions continue to use total measured chain molar masses; subtracting end groups is only the contour mapping. A common M0, Mend and contour increment require a homopolymer or another justified common effective mapping. A single composition-average repeat mass is insufficient for composition-dependent or sequence-dependent populations. No such mapping is inferred from metadata. Explicitly declared cyclic, branched, star, brush or crosslinked architectures are rejected because the available models describe open linear chains.

All displayed geometry is in nm and q in nm⁻¹; squared dimensions use nm². Molar masses use g/mol. Temperature is K. The app does not derive `b(T)` or `lp(T)` for Gaussian/WLC models: the supplied stiffness must apply to the entered temperature and relevant material state.

## 3. Gaussian chain and form factor

For contour length L and Kuhn length b,

\[
\langle R_e^2\rangle=bL,\qquad
\langle R_g^2\rangle=bL/6.
\]

This is the long, ideal, linear-chain model. The app flags components with fewer than ten Kuhn segments as a practical diagnostic, not a sharp universal validity criterion. Nonlocal excluded-volume expansion and electrostatic stiffening are not computed. User-specified effective stiffness cannot by itself reproduce every solvent effect.

The normalized Gaussian continuum form factor is the Debye function,

\[
P_i(q)=\frac{2[e^{-x_i}-1+x_i]}{x_i^2},\qquad
x_i=q^2\langle R_{g,i}^2\rangle.
\]

The [official SasView model documentation](https://www.sasview.org/docs/user/models/mono_gauss_coil.html) gives the formula and its ideal-chain conditions, citing Debye (1947). Polymer Chain Lab implements this analytic formula independently; it does not import SasView code.

For the supplied arbitrary discrete population,

\[
I(q)/I(0)=\sum_i s_iP_i(q).
\]

No Schulz–Zimm distribution is assumed, and no absolute intensity is returned. Small x is evaluated by the Taylor series `1−x/3+x²/12−x³/60+x⁴/360`; the remaining range uses `expm1` to reduce cancellation. At zero q, P=1 exactly in the implementation. The app flags `q_max b ≥ 0.3`; this is a conservative diagnostic for approaching microscopic scales, not a universal cutoff.

When the terminal correction is zero, `L ∝ M` and b is common, the three weighted Rg² values are proportional to Mn, Mw, and Mz, respectively. With a common nonzero terminal correction, each relationship has the same subtraction in molar-mass units: the weighted value is proportional to `(Mn−Mend)`, `(Mw−Mend)` or `(Mz−Mend)`, respectively. This supplies an independently useful invariant for checking a population analysis.

## 4. Continuous worm-like chain

The model has uniform persistence length lp, contour length L, and exponential tangent correlation along contour. For `x=L/lp`,

\[
\langle R_e^2\rangle=2l_p^2[x-1+e^{-x}],
\]

\[
\langle R_g^2\rangle=l_p^2\left[\frac{x}{3}-1+\frac{2}{x}-\frac{2(1-e^{-x})}{x^2}\right].
\]

The gyration expression is the Benoit–Doty result, reproduced with its limiting behavior in [AFRC's official WLC documentation](https://afrc.readthedocs.io/en/stable/polymer_models/worm_like_chain_obrien.html). The original study is [Benoit and Doty (1953)](https://doi.org/10.1021/j150510a025). Only the displayed moment formulas are implemented here; no end-to-end-distribution approximation or sampling code from AFRC is used.

For `x<0.1`, direct subtraction becomes ill-conditioned, especially for Rg². The implementation uses the convergent expansions

\[
\langle R_e^2\rangle=L^2\sum_{k=0}^{10}\frac{2(-x)^k}{(k+2)!},\qquad
\langle R_g^2\rangle=L^2\sum_{k=0}^{10}\frac{2(-x)^k}{(k+4)!}.
\]

The rod limit is `Re²→L²`, `Rg²→L²/12`. The long-chain limit is `Re²→2lpL`, `Rg²→lpL/3`; the corresponding Kuhn length is `b=2lp`. Finite contour length is retained, but nonlocal excluded volume, variable stiffness along a sequence, and finite chain thickness are omitted.

**This mode returns no form-factor curve.** The Gaussian Debye expression is not assigned to a finite-stiffness WLC across the rod-to-coil crossover. The zero-angle population weighting of Rg remains meaningful under the common-contrast assumptions even though a full WLC P(q) is absent.

## 5. Finite-chain rotational-isomeric-state model

### Hamiltonian and geometry

RIS replaces continuous torsion angles with a finite set of local conformations; see the [IUPAC definition](https://goldbook.iupac.org/terms/view/12187). The classic chemical literature includes [Abe, Jernigan and Flory (1966)](https://doi.org/10.1021/ja00956a001). Those references motivate the model class; **their numerical chemistry parameters have not been imported or claimed to validate this implementation**.

An n-bond chain has n+1 equal-mass sites and m=n−2 torsions. Bond lengths and internal angles are fixed. States are trans=π, gauche+=π/3, and gauche−=−π/3. The trans state has energy zero and both gauche states have energy ΔEg. Consecutive unlike-gauche pairs receive an additional user-specified energy J; all other pair energies are zero:

\[
U(s_1,\ldots,s_m)=\sum_{i=1}^m E(s_i)+\sum_{i=1}^{m-1}J(s_i,s_{i+1}),
\qquad p(s)=e^{-U(s)/(RT)}/Z.
\]

Input energies are kJ/mol and `R=0.00831446261815324 kJ mol⁻¹ K⁻¹`. Negative supplied energies are allowed within bounded ranges, but their chemical interpretation remains the user's responsibility. The default length, angle, and energies define a synthetic model with no assigned polymer identity.

The model supports chemically motivated local energy changes, not an automatic chemical structure interpreter. Tacticity, stereocenters, distinct skeletal bond types, side groups, and sequence-dependent transfer matrices are not resolved in v0.1. Different named polymers with identical numerical inputs produce identical predictions.

### Exact finite-chain sampling

Backward partition messages are

\[
B_m(a)=1,\qquad
B_i(a)=\sum_b e^{-[E(b)+J(a,b)]/(RT)}B_{i+1}(b).
\]

The first state is sampled with probability proportional to `exp(−E(a)/(RT)) B1(a)`. Each subsequent state is sampled conditional on the previous state, with weights `exp(−[E(b)+J(a,b)]/(RT)) B_(i+1)(b)`. These finite-chain messages include the remaining chain and the open boundaries. Simply normalizing each row of the bare transfer matrix would generally give the wrong finite-chain distribution.

Forward messages provide exact state and adjacent-pair marginals. Both recurrences use logarithms and `logsumexp`, so allowed low-temperature/high-energy inputs remain finite. The implementation is an original direct recurrence of the stated partition sum. Independent random chains are drawn from this exact discrete torsion distribution; there is no Metropolis chain, burn-in, or kinetic clock.

Positions are reconstructed from the specified bond length, internal bond angle, and torsion. Each chain starts in the same fixed orientation; reported distances are rotationally invariant. The plot projects selected chains onto XY and is not an isotropically oriented ensemble visualization or a proof that nonlocal overlaps are absent.

### Reported estimates and uncertainty

For each reconstructed chain with N=n+1 equal-mass sites,

\[
R_g^2=N^{-1}\sum_i|\mathbf r_i-\mathbf r_{cm}|^2,
\qquad R_e^2=|\mathbf r_N-\mathbf r_1|^2.
\]

The app returns sample means, their RMS dimensions, and standard errors of the **means of squared dimensions**, `sample_std/sqrt(number of chains)`. These are Monte Carlo sampling errors. They exclude parameter uncertainty, model discrepancy, experimental uncertainty, and systematic errors from neglecting nonlocal interactions. A reported small sampling error is not evidence of chemical accuracy.

The finite-chain characteristic ratio is `Cn=<Re²>/(n l²)`. It is not an extrapolated C∞. Exact mean state fractions are compared with sampled fractions; for m=1, the adjacent-pair fraction is undefined and returned as null. RIS exports no invented molar mass for its skeletal sites and does not consume the population CSV.

## 6. Meaningful research questions supported now

- How much can the high-mass tail or zero-angle scattering contribution differ for populations with matched Mn and dispersity?
- Under a stated contour mapping, which chain fractions lie outside a plausible long-Gaussian-chain regime, and how large is the finite-WLC correction?
- How do local gauche energy, neighboring-torsion coupling, and temperature affect finite-chain dimensions under the stated RIS Hamiltonian?
- Is an apparent numerical change larger than Monte Carlo sampling uncertainty, and does it persist as the RIS sample count or MWD quadrature resolution is increased?

Each question requires recording the complete inputs, units, assumptions, and software version. The workbench supplies run export; the scientific API can be called independently. These questions do not establish that one chemistry has better drug delivery, conductivity, detection, or energy-device performance.

## 7. Verification, limits, and release gates

The delivered unittest suite checks:

- Exact mono-/bidisperse moments; number/mass conversion; squared-mass scattering weights; unequal log-grid quadrature; equal-Mn/equal-Đ/different-Mz populations.
- Gaussian dimensions, P(0)=1, P(q²Rg²=1)=2/e, positive decreasing Debye curves, and the mixture Guinier coefficient.
- WLC rod and long-chain limits, continuity across the stable-evaluation branch, contour bounds, and Rg² against an independently evaluated pair-distance integral.
- Exact finite-chain RIS partition sums and state/pair marginals against enumeration of all 3⁴ sequences; uncoupled Boltzmann states; reconstructed bond lengths/angles; planar all-trans end-to-end dimensions; and the all-pair identity for Rg².
- Reproducible sampling, sampled-versus-exact local probabilities, finite extreme-parameter recurrences, invalid inputs, compute limits, material import behavior, and finite JSON output for all presets.

Run from the project root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The delivered suite passed **51 tests**. Passing these tests establishes agreement with the encoded model and selected independent numerical/analytical references. It does not establish experimental validation, clinical usefulness, or predictive accuracy for a named polymer.

Bounded inputs include ≤512 distribution nodes, ≤512 scattering points, 3–500 RIS bonds, 2–2,000 independent chains, and ≤200,000 total reconstructed sites. Exported plots remain below the workbench's 1,200 points per series limit. Numeric fields must be finite and within documented ranges; integer fields reject fractional values. Finite numerical ranges are resource and stability bounds, not statements that every allowed parameter combination is chemically plausible.

Before a chemistry-specific claim, additional gates are needed: justify the structure-to-parameter mapping; calibrate to relevant independent evidence; test transferability across chain length, temperature, and material state; propagate parameter uncertainty; compare held-out observables; and document model failure regimes. Sequence/tacticity-sensitive RIS extensions require explicit chemically justified state/transition definitions and their own enumeration and conformational benchmarks. These gates remain **open** for every named chemistry.
