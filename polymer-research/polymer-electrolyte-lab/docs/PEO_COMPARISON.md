# Relating a PEO/LiTFSI formulation to EIS and PFG-NMR

Open **PEO/LiTFSI · Pesko 2017** to reproduce the 5 kg/mol calculation, or **PEO/LiTFSI · Zheng 2018** to inspect a separate 100 kg/mol specimen. **PEO/LiTFSI · weighed film** is an editable synthetic example for a sample characterised by thickness, electrode contact area and EIS bulk resistance.

## Measurements at EO:Li = 12.5 and 90 °C

| Source and PEO | Salt concentration | Bulk impedance | Lithium tracer D | TFSI tracer D | Computed tracer fraction |
|---|---|---|---|---|---|
| Pesko et al., 5 kg/mol | 1.59 mol/L, Table I, calculated from measured mixture density | 2 mS/cm, Results accompanying Fig. 2a | ≈7.5 × 10⁻¹² m²/s | ≈3.0 × 10⁻¹¹ m²/s | ≈0.20 |
| Zheng et al., 100 kg/mol | Unavailable for this record | 1.5 mS/cm, Table 3 | 6.4 × 10⁻¹² m²/s, Table 3 | 2.7 × 10⁻¹¹ m²/s, Table 3 | 0.192; published rounded value 0.19 |

The 1 g PEO calculation basis in each preset scales the recipe. It does not represent a reported experimental weighing. State IDs refer to matched formulation and temperature across the methods; the instruments used separate specimens.

### Reproducing the figure reading

For the 5 kg/mol material, read the markers at r = Li/EO = 0.08 in Fig. 4, p. E3573. Use the figure's logarithmic diffusion axis in cm²/s. The lithium and fluorine markers lie approximately −0.125 and +0.48 decades from 10⁻⁷ cm²/s. Converting cm²/s to m²/s multiplies both readings by 10⁻⁴. A reading allowance of ±0.035 decades (approximately ±9%) represents about three pixels for a rendering with ~85 pixels per decade. This sensitivity allowance is a reader's estimate. It excludes experimental uncertainties, including uncertainty in conductivity and the paper's estimated 2% density uncertainty. The original numerical tracer measurements are unavailable in tabulated form in this article.

Using these approximate central inputs gives **σNE ≈ 1.84 mS/cm** and **σEIS/σNE ≈ 1.09**. The reading-only interval reaches approximately 1.00–1.20. No inference about a significant excess above unity follows from this comparison. Replicate-level uncertainty would require the underlying observations. [Pesko et al., Journal of The Electrochemical Society (2017), Table I and Figs. 2a and 4](https://doi.org/10.1149/2.0581711jes)

For the 100 kg/mol control, the app reproduces the published tracer fraction from the experiment columns of Table 3. The same table reports a steady-state current fraction of 0.08. Imperfect correspondence between these measurements is expected in a concentrated electrolyte. An independently measured final volume, mixture density or concentration for this state is required to calculate its σNE. [Zheng et al., Macromolecules (2018), Table 3](https://doi.org/10.1021/acs.macromol.7b02706)

## Entering your own formulation

1. Give the material state an ID. Record the PEO molar mass and whether it is Mn, Mw or a nominal vendor grade, the phase, temperature and thermal history. A measured crystalline fraction accepts a method and a mass basis, either PEO or whole electrolyte. Enter `null` if crystallinity is unknown.
2. Enter dry PEO mass and either weighed LiTFSI mass or target ether oxygens per Li. We assume one Li per anhydrous LiTFSI formula unit and calculate EO from 44.05 g/mol. For appreciable end-group mass, set `eo_mass_fraction_in_peo` to the mass fraction actually due to EO repeats.
3. Enter volume evidence at the measurement temperature. Supported alternatives are `density_g_mL` of the dry **mixture**, measured `final_volume_mL` of the entire weighed mixture, or `salt_mol_L`. Omit evidence to retain nonvolumetric results. Component-volume additivity is not assumed. The dry binary model excludes liquid plasticizer, retained solvent, filler and absorbed water.
4. Enter ⁷Li and ¹⁹F self-diffusivities in m²/s. For EIS enter `sigma_mS_cm` or finished film thickness in µm, contacted area in cm² and fitted bulk resistance in Ω. Exclude electrode polarization and contact resistance from the supplied bulk resistance. PFG-NMR self-diffusion and restricted chemical diffusion use different observables.
5. Attach each observation's source and use its actual state ID and temperature. The program rejects unequal IDs or temperatures. Matching labels alone do not establish experimental equivalence; inspect specimen preparation and phase.
6. Run and export the CSV tables or calculation JSON. Each result retains the sources, schema inputs and reference conditions for independent review.

In Python the same GUI presets are available under `polymer_electrolyte_lab.records`. Run a preset through the public plugin entry point:

```python
from polymer_electrolyte_lab.plugin import calculate
from polymer_electrolyte_lab.records import PESKO, gui_inputs

result = calculate(gui_inputs(PESKO))
for metric in result["metrics"]:
    print(metric["label"], metric["value"], metric["unit"])
```
