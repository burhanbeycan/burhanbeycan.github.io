# Water-vapour transport in a fixed PVA film

The experiment starts with a film at equilibrium at initial vapour activity. A humidity step at fixed temperature changes the interfacial water loading according to an independently measured equilibrium branch. The observable is the mean mass of water per unit mass of dry polymer over time. The material identity links the PVA hydrolysis grade, molecular mass, crystallinity or thermal history, and crosslink treatment across runs. Those descriptors constrain which records can be pooled. They do not supply a diffusivity.

## Equilibrium and kinetics

[Equilibrium CSV template](../src/polymer_transport_lab/data/pva_eq_demo.csv) · [Kinetics CSV template](../src/polymer_transport_lab/data/pva_kinetics_demo.csv)

| Field | Meaning |
| --- | --- |
| record_id / curve_id | Provenance ID for an independent equilibrium record or one physical film humidity step. Include an article/table/figure locator or local experiment key. The two ID sets must be disjoint. |
| material_id, temperature_K, branch | Use one material grade and processing history at one temperature, for either sorption or desorption. Encode batch/process and conditioning history consistently. |
| aw, w_g_g | Vapour water activity (RH as a fraction) and equilibrium grams water per gram dry polymer. |
| thickness_um, faces | Full nominal film thickness in µm, plus one exposed face with impermeable backing or two equally exposed faces. Group replicate films by their nominal thickness. |
| aw_start, aw_end, time_s, w_g_g | Initial and final activity, seconds since the humidity step, and observed water per dry polymer at this time. At least three positive times per curve. |

The parser requires increasing, unique activity values after sorting and monotonically nondecreasing equilibrium loading. Linear interpolation is limited to the input range. Sort and control a desorption branch separately when hysteresis is present. A dry-basis mass ratio can exceed one. To convert a measured wet-basis water mass fraction x, use w = x/(1−x).

## Models and geometry

Let w denote g water/g dry polymer. In the fixed geometry and uniform dry mass density approximation, solve ∂w/∂t = ∂[D(w)∂w/∂x]/∂x. The initial loading is spatially uniform. A surface water activity step sets the new loading at the exposed face. Zero flux applies at the backing for a supported film and at the midplane for a symmetrically exposed film. The corresponding one-way diffusion distances are L and L/2. At equal L and D(w), a one-face film takes four times as long to reach the same water loading as its two-face counterpart.

The baseline is the constant positive D slab series. The variable model uses D(w) = Dref exp[B(w−wmid)/(wmax−wmin)], where wmid is the midpoint of the **independently supplied** isotherm loading range. B is the natural log of the ratio of D at its two endpoints. The conservative volume scheme differences the exact integral of D across neighbouring loadings. BDF integrates the resulting cell balances. Equal dry mass per cell gives the specimen mean. Fitting bounds and residual units appear in the result.

## Thickness validation

For an analyst-selected nominal thickness, retain all its specimens, humidity steps, replicates and time points in the test set. Freeze the independent equilibrium isotherm, fit the two transport models on the remaining kinetic curves, then compute held-out RMSE in g water/g dry polymer. The objective gives equal total weight to each curve. Inside the training set, include at least two different humidity steps and sample the transient region. Curve provenance alone cannot establish independent specimens or guarantee a reliable parameter estimate. Use actual thickness grouping consistently. The chosen split evaluates thickness transfer for the supplied material, temperature and branch.

The synthetic example uses Dref = 1.5 × 10⁻¹² m²/s and B = 2. It was generated with 96 cells, then fitted with 32 cells. Its two 50 µm backed training curves and two 90 µm free-film holdout curves cover two activity steps each. This fixture verifies numerical implementation under an imposed transport law. A measured thickness-holdout error is pending for a time-resolved PVA dataset with independent equilibrium and specimen thickness labels.

## Published material check

El Shafee and Naguib report vapour-equilibrium sorption and dry-film sorption kinetics at 25 °C for HDI-crosslinked PVA. For a nominal 5 mol% HDI feed per vinyl alcohol repeat, Table 2 reports GAB parameters hm = 0.062 g/g, f = 0.83 and c = 1.70. Equation 2 gives w = hm·c·f·aw/[(1−f·aw)(1+(c−1)·f·aw)]. The GUI evaluates that fitted relation. Table 3 reports apparent early-time D values from gravimetry at aw = 0.38, 0.61 and 0.74. They correspond to 9.4 × 10⁻¹², 1.4 × 10⁻¹¹ and 2.4 × 10⁻¹¹ m²/s. These are whole-step apparent coefficients. The local exponential D(w) requires calibration to full kinetic records. [Article](https://doi.org/10.1016/S0032-3861(02)00865-0)

Jeck and co-workers measured equilibrium uptake independently and assessed a concentration-dependent diffusion law in **physically** crosslinked PVA. Their material differs from the HDI-crosslinked films above. No parameters are transferred between these studies. A change in crystallinity, pronounced swelling, a finite gas-side transfer coefficient or slow polymer relaxation can affect the observed kinetics. The fixed-film workflow addresses vapour steps in an existing film. A solution-cast wet layer undergoing evaporation requires moving boundaries, solvent-air mass transfer and a mass/density balance. [Article](https://doi.org/10.1016/j.cep.2010.09.004)

### A one-activity interpolation check

The published-paper preset fits only the apparent D at aw = 0.38 and 0.74, then compares each model's prediction at the reserved aw = 0.61 against Table 3. An equal-weight, squared-error fit to the two endpoint coefficients gives 1.67 × 10⁻¹¹ m²/s for a constant D. For the exponential local diffusivity the numerical comparator solves the **semi-infinite early-time problem**, with the GAB loading as its surface concentration in each activity step. For a dry initial film, D(w)=Ddry exp(kw); set y=w/w_surface, η=x/(2√(D_surface t)), and j=(D/D_surface)dy/dη. The similarity system is dy/dη=j D_surface/D and dj/dη=−2η dy/dη, with y(0)=1 and y(∞)=0. Mass balance converts its surface flux into the apparent early-time coefficient D_app = π D_surface j(0)²/4. This recovers the constant-D gravimetric slope when k=0. The exponential model's two parameters match two endpoint values.

At the held-out middle humidity the published apparent D is 1.40 × 10⁻¹¹ m²/s. The exponential model gives 1.509 × 10⁻¹¹ m²/s (7.82% absolute relative error) and constant D gives 1.670 × 10⁻¹¹ m²/s (19.29%). The authors estimated slopes over an initial finite window. The computed slope is a short-time limit. There is one intermediate point, from the same paper, and no reported estimate of its uncertainty in that table. These results describe a within-study interpolation of apparent coefficients.

### Why the measured thickness-holdout is still open

Horkay et al. studied nominal 0.1, 0.2 and 0.5 µm attached PVA films. The authors discuss long relaxation times, the chamber vapour path and substrate constraints. Their Figure 1 is a single activity step across thicknesses, whereas identifying the concentration dependence in this workflow requires multiple training steps. Hasimi et al. also varied thickness and reported relaxation and differential swelling stresses. Their liquid-water uptake changes the interfacial conditions and can strongly change geometry. These studies identify separate validation regimes and the limits of the fixed-film, imposed-vapour-activity approximation. [Horkay et al.](https://doi.org/10.1021/bm049332c) · [Hasimi et al.](https://doi.org/10.1016/j.eurpolymj.2008.09.011)
