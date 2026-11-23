# Primary references and implementation provenance

The equations were implemented independently from their mathematical
definitions. These references establish the selected theory and its historical
context. They do not validate this application against experiments, endorse its
software, or supply a measured polymer-property database.

## Implemented free-energy convention

**Hong, W.; Zhao, X.; Zhou, J.; Suo, Z. (2008).** A theory of coupled diffusion and
large deformation in polymeric gels. *Journal of the Mechanics and Physics of
Solids* **56**, 1779–1793.

- DOI: [10.1016/j.jmps.2007.11.010](https://doi.org/10.1016/j.jmps.2007.11.010)
- [Author-hosted full paper](https://suo.seas.harvard.edu/sites/g/files/omnuum4271/files/suo/files/200.pdf)
- Equations (20)–(24) specify additive mixing and Gaussian elastic energies.
  The elastic logarithm has coefficient −2. The paper's freely swollen initial
  state for `Nv=10^-3`, `χ=0.2` gives the rounded stretch 3.215 used in one
  numerical test. This application implements homogeneous free equilibrium,
  not the paper's full coupled deformation–diffusion field problem.

**Hong, W.; Liu, Z.; Suo, Z. (2009).** Inhomogeneous swelling of a gel in
equilibrium with a solvent and mechanical load. *International Journal of Solids
and Structures* **46**, 3282–3289.

- DOI: [10.1016/j.ijsolstr.2009.04.022](https://doi.org/10.1016/j.ijsolstr.2009.04.022)
- [Author-hosted full paper](https://suo.seas.harvard.edu/sites/g/files/omnuum4271/files/suo/files/204.pdf)
- Equation (19) directly supplies the unconstrained isotropic equilibrium
  relation. This application sets the reservoir chemical potential to zero.
  Its small-shear expression is obtained by differentiating the same elastic
  energy about the homogeneous swollen state. The constrained and
  inhomogeneous examples in the source are outside the implemented scope.

## Historical theory

**Flory, P. J.; Rehner, J. Jr. (1943).** Statistical Mechanics of Cross-Linked
Polymer Networks II. Swelling. *The Journal of Chemical Physics* **11**, 521–526.

- DOI: [10.1063/1.1723792](https://doi.org/10.1063/1.1723792)
- Historical basis for balancing mixing and network elasticity. The current
  software states the Hong/Suo constitutive convention explicitly; it does
  not mix its elastic `−φ` term with a historical `−φ/2` relation.

## Diffusion law; independently derived slab solution

**Fick, A. (1855).** On liquid diffusion. *The London, Edinburgh, and Dublin
Philosophical Magazine and Journal of Science*, Series 4, **10**, 30–39.

- DOI: [10.1080/14786445508641925](https://doi.org/10.1080/14786445508641925)
- [Publisher record](https://www.tandfonline.com/doi/abs/10.1080/14786445508641925)
- Establishes the diffusion-law provenance. The particular uniformly loaded,
  two-face, perfect-sink slab series in this repository follows by separation
  of variables for the stated boundary-value problem; it is not presented as
  a new transport theory or an experimental result from this historical paper.

## Validation status

Unit tests provide numerical and analytical verification, including comparison
of two mathematical representations of slab release. The checks do not establish
prediction accuracy for any named hydrogel, drug carrier, membrane, sensor,
electrolyte, biological response or clinical application. Any material-specific
calibration should record the sample, preparation state, solvent, temperature,
measurement protocol, uncertainty and held-out validation conditions.

## WBPU experiments and mass calculations

[WBPU methods](WBPU.md) identifies the study, source tables, individual measurements, mass denominators and matched test conditions for the independent WBPU modes.

- Beycan et al., LDI/PEG/PCL nanofiber preparation, alternative extenders and aqueous L-lysine stage. [Original study](https://doi.org/10.1021/acs.iecr.4c00759)
- Hao et al., PEG/PCL waterborne polyurethanes with L-lysine ethyl ester diisocyanate and PDO. [Original study](https://doi.org/10.1039/C6BM00588H)
