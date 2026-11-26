# Scientific scope and numerical verification

Version 0.2.0 implements a transparent, independently usable mathematical core for **constant-diffusivity transport through a planar polymer film**. It is a starting point for research software, not a validated material database. All bundled examples use synthetic parameters. The model, spatial discretization, boundary-flux accounting, fitting residuals and workbench adapter are implemented in this project. NumPy and SciPy supply general numerical operations, integration and optimization; no polymer application is a runtime dependency.

## 1. Geometry, concentration and physical assumptions

The film occupies `0 ≤ x ≤ L`. **L is the full thickness**, measured between its two exposed faces. The UI accepts `thickness_um` in micrometres and converts it to metres before calculation. This distinction matters: replacing full thickness by half thickness changes diffusion times by a factor of four.

The governing equation is

\[
\frac{\partial c}{\partial t}=D\frac{\partial^2c}{\partial x^2}.
\]

`D` is a positive, scalar, effective diffusivity in m²/s. It is constant in position, time and concentration. The model assumes one-dimensional transport, uniform initial concentration, fixed geometry, imposed surface concentrations and no external transfer resistance. It excludes binding, degradation, erosion, reactions, finite-bath depletion, convection, concentration-dependent diffusivity, moving boundaries and dynamic swelling. It cannot predict a release curve from repeat-unit chemistry, molar mass or mesh size alone.

**Release:** the film is initially uniformly loaded at `c₀`; both faces are placed in perfect sinks, so `c(0,t)=c(L,t)=0` for `t>0`. The displayed fraction is released mass divided by initial mass.

**Uptake:** the initial concentration is zero; both faces are held at a common concentration `c∞`. The displayed fraction is mass gained divided by equilibrium mass. The normalized progress is the same function as for release, while the concentration profiles are complements. This normalization does not determine an absolute mass or solubility.

The exact initial profile and the imposed boundary conditions are incompatible at the instant of switching. The plot displays the uniform initial state at `t=0`, followed by the imposed face concentrations for `t>0`. A finite spatial grid cannot resolve an arbitrarily thin initial boundary layer.

Temperature is recorded as context for the supplied D. No temperature dependence, Arrhenius prefactor or activation energy is inferred.

## 2. Analytical full-slab solution

Define the Fourier number with the **full thickness**:

\[
\mathrm{Fo}=Dt/L^2.
\]

For either normalized uptake or release, the transferred fraction is

\[
f(\mathrm{Fo})=1-\frac{8}{\pi^2}\sum_{n=0}^{\infty}
\frac{\exp[-(2n+1)^2\pi^2\mathrm{Fo}]}{(2n+1)^2}.
\]

The official [PolyKin slab documentation](https://hugomvale.github.io/polykin/reference/hmt/diffusion/uptake_constc_sheet/) states the same solution with half-thickness `a=L/2` and cites Crank (1975), page 48. PolyKin is a reference, not an imported dependency.

Direct truncation followed by `1 − remaining` is inaccurate very near zero time. This implementation uses an equivalent integrated image series for `0<Fo≤0.05`:

\[
f=4\sqrt{\frac{\mathrm{Fo}}{\pi}}
+\sum_{k=1}^{\infty}(-1)^k\left[
8\sqrt{\frac{\mathrm{Fo}}{\pi}}e^{-k^2/(4\mathrm{Fo})}
-4k\,\operatorname{erfc}\left(\frac{k}{2\sqrt{\mathrm{Fo}}}\right)
\right].
\]

Eight image terms suffice over this branch; omitted exponential factors start below `exp(-405)`. The long-time branch uses 40 odd Fourier terms. At exactly zero time, the initial condition gives `f=0` explicitly. Unit tests compare the image expression with a separately evaluated long Fourier sum near their overlap. The short-time limit is `f≈4 sqrt(Fo/π)`, representing two exposed faces. Tests also verify full-thickness scaling and long-time convergence to unity.

The first-mode relaxation time is `τ₁=L²/(π²D)`. It is an exponential time constant, not the completion time. Times to 50%, 90% and 99% transfer are obtained by scalar root finding on the analytical fraction. Complete transfer is asymptotic.

## 3. Independent finite-volume verification

The film is partitioned into N equal-volume cells with centres half a cell from each face. The solver evolves a normalized remaining concentration `r`; for uptake this is a concentration deficit. Initially every cell has `r=1`, and both face values are zero.

In dimensionless time the interior equations are

\[
\dot r_i=N^2(r_{i-1}-2r_i+r_{i+1}).
\]

The left boundary equation is `ṙ₀=N²(r₁−3r₀)`, with the symmetric expression at the right boundary. The factor three follows from the half-cell boundary distance. It is not the stencil for a boundary node.

A separate state integrates transfer through both faces:

\[
\dot q=2N(r_0+r_{N-1}),\qquad q(0)=0.
\]

The mass-balance invariant is `mean(r)+q=1`. Cumulative flux is integrated as a state, **not calculated by subtracting mean concentration from one**. Consequently the balance checks the assembled flux accounting as well as the time solve. The raw balance is evaluated before tiny concentration overshoots from numerical roundoff are clipped for display. A major bound violation raises an error.

SciPy's generic BDF initial-value integrator solves the sparse linear ODE system with its own adaptive time steps (`rtol=2e-9`, `atol=2e-11`). Output points are uniformly spaced in physical time; they are not integration time steps. The analytical series is not used to advance the finite-volume state.

For requests beyond `Fo=8`, the computed terminal numerical state is held to avoid unbounded integration over extreme scales. At this point the exact continuum remaining fraction is at most `exp(-8π²)≈5.2×10⁻³⁵`; numerical error is much larger than that physical tail. This approximation is disclosed in results. It does not overwrite the raw mass-balance calculation.

The graph compares analytical fraction, fraction from the spatial average and independently integrated boundary transfer. Their agreement is a numerical verification of the stated idealized problem. It is **not experimental validation**. The test suite checks initial and long-time limits, profile symmetry and bounds, monotonicity, conservation, grid convergence and agreement with the analytical solution. Early-time discrepancy must decrease with refinement; increasing output points alone does not refine the spatial grid.

## 4. Fitting an effective diffusivity

The optional fit accepts 5–1000 `time_s,fraction` observations, with at least three distinct positive times. Fractions must be in `[0,1]` and must mean **transferred fraction**, not fraction remaining. The optional CSV header is `time_s,fraction`. Repeated observation times are permitted; input is sorted without discarding replicates.

The objective is an equal-weight sum of squared fraction residuals. D is optimized in `log10(D)` within explicit user-supplied positive bounds. A coarse objective scan initializes bounded least squares and is returned for inspection. Thickness, normalization, surface conditions and time zero remain fixed. No burst release, lag time, baseline or instrument response is silently fitted.

**Identifiability:** the curve depends on `D/L²`. If L is unknown, fitting the curve cannot separately identify D and L. With a supplied L, the result is conditional on that thickness. A small relative thickness error produces approximately twice that relative error in D. Model mismatch, external mass-transfer resistance, inaccurate equilibrium normalization, swelling or an uncertain time origin can also bias D.

The app reports RMSE, the RMS local derivative with respect to `log10(D)`, bound warnings and a conditional local standard error when it is informative. The latter uses residual variance and the single-parameter Jacobian, assumes independent equal-variance fraction errors, and excludes thickness and model uncertainty. It is omitted at active bounds, at low sensitivity or with essentially zero residuals. It is not a robust uncertainty estimate for correlated or heteroscedastic measurements, and the objective scan is not a confidence interval. These diagnostics do not guarantee identifiability.

The bundled **Synthetic D recovery** dataset is generated from this same analytical model without noise. Recovery is a software check only. A second test fits independently computed finite-volume data. Neither substitutes for external experimental validation.

## 5. Separate steady permeation calculation

An optional calculation evaluates

\[
c=Sp,\qquad P=DS,\qquad J_{\mathrm{ss}}=P\frac{\Delta p}{L}.
\]

S must be supplied on the explicit molar Henry basis, in `mol m⁻³ Pa⁻¹`. Therefore P has units `mol m⁻¹ s⁻¹ Pa⁻¹`, and flux has units `mol m⁻² s⁻¹`. A positive pressure difference gives positive upstream-to-downstream flux. The official [Sandia reference](https://www.sandia.gov/app/uploads/sites/158/2021/12/8100TechRef_polymers.pdf) gives `P=DS` and `c=Sf` with fugacity f; this app uses pressure only where the ideal-gas approximation `f≈p` is appropriate. No tabulated material parameters are copied from that reference.

This is a **different steady experiment** from transient two-face release or uptake. It assumes constant D and S, local interfacial equilibrium, a homogeneous defect-free film and negligible external resistance. Nonlinear sorption, dual-mode sorption, competitive sorption, pressure-dependent D, gas nonideality and plasticization are outside scope. An uptake or release fit by itself does not identify S.

## 6. Interchange and provenance

The mathematical core has no workbench dependency. The adapter accepts material metadata without using chemical names or molecular-weight distributions to manufacture coefficients. Molar-mass populations are not reweighted or transformed. Imported metadata is not marked experimentally verified.

Each exported parameter records its quantity, unit, origin and model scope. A fitted D records whether its input data were labelled synthetic or user supplied and unverified. A network handoff can supply a **fully equilibrated wet full thickness** plus an explicit solute diffusivity. The subsequent calculation remains a fixed-geometry slab; it is not a swelling–diffusion coupling model.

Model identifier: `constant-D-two-face-full-slab/1.0`. Application version: `0.2.0`. Reproducible comparisons should retain the application version, input values, full-thickness convention, boundary conditions, data provenance and numerical grid. Default constants must not be relabelled as measured properties.

## 7. Computational limits and release tests

Supported input ranges bound computation and floating-point behavior; they do not establish physical validity over those ranges. Full thickness is 1 nm–1 m, D is `10⁻²⁴–10⁻⁴ m²/s`, duration is `10⁻⁹–10¹² s`, output size is 20–1000 points, and the verification grid is 8–160 cells. A continuum slab model at molecular dimensions requires independent physical justification. All inputs are validated by the backend; NaN, infinity, nonintegral grid sizes and invalid choices raise readable errors.

Run numerical tests from the repository root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The scope of these tests is mathematical and numerical behavior. A publication would additionally need a stated research contribution, independently obtained material data or a justified benchmark, reproducibility evidence and appropriately quantified uncertainty. Implementing a known slab equation alone is not a new transport theory.
