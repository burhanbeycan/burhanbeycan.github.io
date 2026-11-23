# Scientific scope and verification

## What this application calculates

Polymer Network Lab implements homogeneous equilibrium swelling of a **neutral,
permanently connected Gaussian network in a pure solvent reservoir**. The network
is free to swell isotropically and has no applied mechanical load. Its mixing and
elastic energies are evaluated consistently; equilibrium states are roots of the
derivative of that same energy. An optional, independent diffusion calculation
uses the resulting **fully equilibrated slab geometry**.

The constitutive convention follows the free-energy forms in Hong, Zhao, Zhou and
Suo (2008), equations (20)–(24), and the free-swelling equation in Hong, Liu and Suo
(2009), equation (19). Primary-source links are in [REFERENCES.md](REFERENCES.md).
The implementation does not call a molecular dynamics engine, PolyChemMD, another
application's numerical model, or a commercial simulator.

This is a research calculation tool with numerical checks. No experimental
dataset has been fitted or used to establish material-specific predictive
accuracy. The demonstration inputs are synthetic scenarios.

## Definitions, units and reference state

| Quantity | Definition | Units |
|---|---|---|
| `nu_e_mol_m3`, νe | Effective coefficient of the selected network elastic energy, per dry polymer volume | mol/m³ dry polymer |
| `solvent_molar_volume_m3_mol`, Vs | Incompressible solvent molar volume | m³/mol |
| `temperature_K`, T | Absolute temperature | K |
| `chi`, χ | Constant Flory–Huggins interaction parameter for one calculation | dimensionless |
| φ | Polymer volume fraction in the swollen network | dimensionless |
| Q = J | Swollen volume divided by dry polymer volume | dimensionless |
| λ | Isotropic stretch relative to the dry reference | dimensionless |
| W | Free-energy change per dry polymer volume | J/m³ dry polymer |
| Gdry, Gswollen | Ideal incremental shear moduli in the stated convention | Pa |
| `dry_thickness_um` | **Full** slab thickness in the dry reference | µm in the interface |
| `diffusivity_m2_s`, D | Independently supplied effective diffusivity of a dilute solute in the equilibrated gel | m²/s |

The gas constant is R = 8.31446261815324 J/(mol K). Molar quantities replace the
per-molecule quantities in the cited continuum model: `νe Vs` corresponds to its
dimensionless `Nv`. Constituent incompressibility and retention of the connected
polymer imply Q = 1/φ. Isotropic swelling gives λ = Q^(1/3).

Q is **not** a mass swelling ratio. Converting to solvent mass uptake additionally
requires appropriate solvent and dry-polymer densities and a clearly defined
mass basis. This application does not invent those inputs.

### Effective density and network defects

νe is an explicit, effective elastic parameter. It is not a crosslinker feed
concentration, a count of chemical junctions, a nominal crosslink percentage, or
an automatically inferred strand count. Its interpretation must be compatible
with the selected free energy and reference state. Loops, dangling ends, soluble
fractions, entanglements and heterogeneous junction environments can complicate
that interpretation. The application records parameter and topology provenance
without calculating a correction from unsupported information.

No extra phantom-network functionality factor is applied. Applying such a factor
to an already effective elastic density could double-count a correction. A
network prepared in solvent may require a preparation-state theory; merely
reporting a dry thickness does not establish that this constitutive reference
state is appropriate. Chain-length distributions, synthesis history, tacticity
and chemistry metadata are useful sample information but are not sufficient to
determine χ, νe or D here.

## One consistent thermodynamic model

Define a = νe Vs. In dimensionless form, relative to the dry state,

```text
w(Q) = W / (RT/Vs)
     = (Q − 1) ln(1 − 1/Q) + χ(1 − 1/Q)
       + (a/2) [3Q^(2/3) − 3 − 2 ln Q],       Q ≥ 1.
```

The first line of contributions is mixing and the second is elastic. The dry
mixing limit is defined by continuity and W(1) = 0. The mixing expression differs
from the displayed Hong form by an additive constant; derivatives and selected
equilibria are unchanged. The connected network has no independent polymer-chain
translational mixing term.

The reservoir chemical potential is zero relative to pure solvent. Thus the
relevant homogeneous potential equals this free-energy change. With φ = 1/Q,

```text
g(φ) = dw/dQ
     = ln(1 − φ) + φ + χφ² + a(φ^(1/3) − φ).

Equilibrium: g(φ) = 0.

g'(φ) = φ[2χ − 1/(1 − φ)] + a[(1/3)φ^(−2/3) − 1].
d²w/dQ² = −φ² g'(φ).
```

The elastic energy contains **−2 ln Q**. Consequently the elastic term in the
equilibrium equation contains **−φ**. A frequently used historical
Flory–Rehner equation contains **−φ/2**, corresponding to a different elastic
free-energy convention. Substituting one equation's fitted parameters into the
other without reconsidering that convention is not justified. This application
implements only the stated Hong/Suo convention.

At fixed χ, νe and Vs, temperature multiplies the full potential by RT. It changes
energies and ideal moduli but leaves the ideal equilibrium Q unchanged. A
temperature-responsive material requires a justified χ(T,φ), another constitutive
change, or an additional mechanism; changing the temperature field alone is not
a model of an LCST or UCST transition.

### Root search and homogeneous stability

The solver searches 10⁻¹⁴ ≤ φ ≤ 1 − 10⁻¹⁴. It combines 1,025 logit-spaced
points with the extrema of g. The extrema are generated from the degree-eight
polynomial obtained after setting y = φ^(1/3):

```text
3 y² (1 − y³) g'(φ)
    = a − 3a y² − a y³ + (3a + 6χ − 3)y⁵ − 6χ y⁸.
```

Brent bracketing in logit coordinates resolves sign-changing roots; candidate
extrema are also checked for tangent roots. A Taylor representation avoids
subtracting nearly equal first-order terms in the dilute mixing residual.
The free-energy implementation uses a series for the very small elastic energy
close to Q = 1.

Each reported root includes its residual, energy, curvature and local stability.
Positive curvature means a local minimum with respect to homogeneous volume
change. The lowest-energy locally stable homogeneous root is selected. The
reported relative residual is `|g|` divided by the sum of the absolute terms in
the equilibrium equation. It is a numerical diagnostic, not experimental error.
Extremely small elastic terms in a poor solvent can be below floating-point
resolution after the mixing terms cancel; the report therefore retains the
absolute residual as well.

This procedure does not calculate spatial phase coexistence, gradients,
nucleation barriers, surface instabilities, mechanical buckling or kinetic
hysteresis. A parameter scan independently equilibrates every point and must
not be interpreted as a time trajectory.

## Shear modulus in the same elastic convention

For a small isochoric shear superposed on an isotropically swollen state,

```text
Gdry = νe RT
Gswollen = Gdry λ²/J = Gdry/λ.
```

The factor λ² comes from the elastic strain energy, while division by J converts
dry-volume energy to current swollen-volume energy. Mixing does not change under
this isochoric perturbation. This is an ideal incremental **shear** modulus at
fixed solvent content; it is not a drained tensile/Young's modulus. No glassy,
crystalline, entanglement, viscoelastic or finite-extensibility contribution has
been added.

## Solute release after full swelling equilibration

The optional preview is a second, explicitly limited boundary-value problem.
The full dry thickness Ldry is converted to the full equilibrium thickness
L = Ldry λ. At time zero the **already swollen** slab contains a uniform dilute
solute. The two faces are exposed to perfect sinks, and the user supplies an
effective, constant solute diffusivity D:

```text
∂c/∂t = D ∂²c/∂x²,   0 < x < L
c(0,t) = c(L,t) = 0,   t > 0
c(x,0) = c0,          0 < x < L.
```

Separation of variables and integration over the full slab give

```text
Mt/M∞ = 1 − (8/π²) Σ[n=0..∞]
             exp[−(2n+1)² π² D t/L²] / (2n+1)².
```

`L` is **full thickness**, not half-thickness. For the Fourier number
Fo = Dt/L² ≤ 0.005, the implementation uses the short-time expression
4 sqrt(Fo/π); the first face-overlap correction is below 10⁻²⁰ in this interval.
For larger Fo it evaluates 128 odd modes. The value at t = 0 is exactly zero.
Target release times are solved independently of the displayed observation
window. Tests compare the Fourier branch against the independently transformed
image/erfc representation, as well as the short-time limit and diffusion
similarity laws.

This preview does not calculate dynamic network swelling, diffusion–relaxation
coupling, drug–polymer binding, concentration-dependent D, dissolution, erosion,
chemical degradation, finite bath depletion or external mass-transfer resistance.
It assumes the solute does not change the equilibrium solvent/network state.
Changing νe affects the estimated geometry; it does not automatically establish
how solute diffusivity changes. D is never inferred from a mesh size or νe.

### Transport handoff

The application can emit a handoff with these exact keys:

```text
target: transport
mode: release
thickness_um: FULL equilibrium slab thickness
diffusivity_m2_s: explicit user input
temperature_K: explicit input temperature
duration_s: observation duration
n_points: number of sampled times
```

The handoff includes the same two-face perfect-sink and fully equilibrated,
fixed-geometry assumptions. No handoff is generated if the swollen thickness is
outside Transport Lab's supported range. Each repository retains an independent
numerical implementation.

## Supported input domain and resource limits

The scientific API accepts −5 ≤ χ ≤ 5, 10⁻⁶ ≤ νe ≤ 10⁶ mol/m³ dry,
10⁻⁸ ≤ Vs ≤ 10⁻³ m³/mol and 1 ≤ T ≤ 3000 K. The interface restricts T to
2000 K for compatibility with the transport handoff. These are bounded numerical
domains, not evidence that the physical model is valid at every allowed state.

Scans contain 2–160 increasing values. Output series contain at most 1,200 points;
the release interface allows 20–600 times. Interface diffusivity is restricted
to 10⁻²⁴–10⁻⁴ m²/s, observation duration to 10⁻⁹–10¹² s, and full dry thickness
to 0.001–10⁶ µm. Every backend input is checked independently of the interface.
Non-finite values, booleans used as numbers, unsupported selections, inconsistent
scan limits and excessive sizes raise readable `ValueError` exceptions.

## Provenance, interpretation and next validation work

The material exchange records each exported quantity, unit, origin and model
scope. Model outputs additionally carry the source inputs' demonstration/user
origin. Imported sample chemistry, molar-mass distributions and conditions are
not silently converted into model parameters. The calculation does not mutate
the imported record.

The test suite checks energy derivatives, stable homogeneous minima, root
residuals, dry and dilute limits, dimensional collapse, temperature scaling,
monotonic parameter trends in specified regimes, bounded release, analytical
limits, an independent image-series comparison, invalid inputs, finite JSON and
handoff units. The published-model check reproduces the rounded stretch
λ ≈ 3.215 for νeVs = 0.001 and χ = 0.2 from Hong et al. (2008). That is a
**numerical model benchmark**, not reproduction of experimental measurements.

For material research, the next step is a traceable sample-specific dataset:
equilibrium volumes at defined solvent activity and temperature, preparation and
thermal history, independently constrained elastic response, and solute release
under documented geometry/bath conditions. Fitting χ and νe to a single swelling
measurement generally does not identify both parameters independently. A scan
shows parameter sensitivity; it does not provide statistical uncertainty or
prove transfer to another composition or sample.

Ionic/pH-sensitive, constrained, degrading and dynamically swelling networks
require additional validated equations and inputs. They are outside version
0.2.0 and are not represented by synthetic correction factors.

## Run numerical checks

From this repository:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Minimal independent Python use:

```python
from polymer_network_lab import NetworkParameters, solve_equilibrium

parameters = NetworkParameters(
    chi=0.2,
    nu_e_mol_m3=55.55555555555556,
    solvent_molar_volume_m3_mol=18e-6,
    temperature_K=298.15,
)
equilibrium = solve_equilibrium(parameters)
print(equilibrium.selected.isotropic_stretch)
# Approximately 3.21502, for this model demonstration only.
```
