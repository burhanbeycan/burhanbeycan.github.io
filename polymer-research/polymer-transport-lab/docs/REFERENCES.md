# References and source boundaries

The following official sources were checked for the implementation. External polymer applications are references only; none is a required runtime dependency and no source code was copied from them.

1. **PolyKin official documentation — fractional mass uptake in a plane sheet.**
   https://hugomvale.github.io/polykin/reference/hmt/diffusion/uptake_constc_sheet/

   States the two-face, constant-surface-concentration slab series and its half-thickness convention. It cites J. Crank, *The Mathematics of Diffusion*, Oxford University Press (1975), p. 48. This project uses full thickness `L=2a`, so its Fourier number is `Dt/L²`. The official documentation was used to verify that convention and equation; the package is not imported. The book publisher landing page was access-limited, so it is not presented as a directly inspected source.

2. **SciPy official documentation — `scipy.integrate.solve_ivp`.**
   https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

   Describes BDF integration, sparse Jacobians, tolerance semantics and requested output times. The finite-volume spatial operator and independently integrated boundary-flux equation are implemented in this repository.

3. **SciPy official documentation — `scipy.optimize.least_squares`.**
   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html

   Defines bounded nonlinear least squares and its residual/Jacobian interface. The log-D parameterization, slab objective, fit-data validation and interpretation are implemented here.

4. **C. San Marchi, Sandia National Laboratories — Technical Reference on Hydrogen Compatibility of Materials: Nonmetals, Polymers (code 8100).**
   https://www.sandia.gov/app/uploads/sites/158/2021/12/8100TechRef_polymers.pdf

   Section 2, equations 1–3, states the permeability–diffusivity–solubility relation, molar concentration/fugacity convention and dimensions. This app uses Pa rather than the MPa convention in the reference. Pressure may replace fugacity only under an appropriate ideal-gas approximation. The app does not use the reference's material tables or extrapolate hydrogen properties to other species.

## Validation terminology

Analytical–numerical comparison, mass conservation, grid convergence and synthetic parameter recovery are **verification** of this implementation under a stated ideal model. They do not establish that constant-D Fickian transport describes a particular polymer, penetrant, solvent, temperature or experiment. No clinical, membrane-performance or experimental validation claim is made.
