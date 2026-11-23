# Intrinsic viscosity of linear PEtOx

How well does SEC–MALS mass-average molar mass explain the measured dilute-solution intrinsic viscosity within one homologous series? This module fits each temperature independently. It reports an interpolation for a comparable linear PEtOx sample and checks the effect of withholding an entire nominal mass class.

## Measurements and protocol

The seven samples P1–P7 come from Gubarev and colleagues. The [source data in this project](../src/polymer_chain_lab/data/petox_pbs.csv) preserve the published units. Each viscosity column is in cm³/g, equivalent to mL/g. Blank cells denote unreported measurements. The seven nominal synthesis targets define mass classes. The regression uses measured SEC–MALS **Mw**, converted from kDa to g/mol. The nominal target, Mn and chromatographic peak molar mass are excluded from the predictor. The series spans Mw of 11.0–330.5 kDa and reported dispersities of 1.04–1.18. [Study](https://doi.org/10.1039/C8PY00255J) · [SI, Tables S1 and S2 and experimental methods](https://www.rsc.org/suppdata/c8/py/c8py00255j/c8py00255j1.pdf)

| Measurement | Medium and conditions |
|---|---|
| SEC–MALS Mw | N,N-dimethylacetamide containing 50 mM LiCl, 50 °C columns, flow 0.500 mL/min |
| Intrinsic viscosity [η] | Aqueous PBS prepared with one Sigma-Aldrich tablet per 200 mL. The authors report pH 7.60 ± 0.04 for the prepared buffer. Rolling-ball dilution viscometry with Huggins and Kraemer intercepts averaged |
| Temperature | Fit one measured viscometry temperature per analysis. The temperature of the SEC column stays at 50 °C |

## Mass-class cross-validation

For each class, remove every sample in that class, fit ln[η] = ln K + a ln Mw on the remaining classes and predict the excluded samples. The table in the application reports absolute percentage error for every omitted sample. “Outside” marks a test sample beyond that fold's training range. Median and 90th percentile summarize sample-wise absolute percentage errors. The fitted K assumes M in g/mol and [η] in mL/g.

| Viscometry temperature | Measured samples / mass classes | All-data exponent a | Held-out median absolute relative error | Held-out 90th percentile | Held-out MAE (mL/g) |
|---:|---:|---:|---:|---:|---:|
| 15 °C | 4 / 4 | 0.6351 | 6.73% | 9.57% | 5.27 |
| 25 °C | 7 / 7 | 0.6273 | 1.72% | 7.47% | 3.28 |
| 37 °C | 7 / 7 | 0.5921 | 2.33% | 7.12% | 2.57 |
| 45 °C | 4 / 4 | 0.5472 | 9.32% | 13.34% | 4.22 |
| 55 °C | 4 / 4 | 0.4793 | 10.91% | 15.53% | 3.62 |

At 37 °C the full-data fit predicts 34.24 mL/g at Mw = 50,000 g/mol. The seven-fold test comes from one published series. It gives no interlaboratory error estimate and does not establish transfer to other buffers, temperatures, end groups, architectures or broad distributions. Measurement errors in Mw and [η] are not propagated. The fitted slope also depends on the mass estimator. The original article's exponent of 0.62 at 37 °C was based on sedimentation–diffusion masses and the authors' scaling analysis. The fit above pairs SEC–MALS Mw with the same viscosity observations. [Study, Figure 1](https://doi.org/10.1039/C8PY00255J)

## Supply a new series

Paste a long CSV with the header shown below. All rows must represent linear PEtOx homopolymers with comparable end chemistry and dispersity, SEC–MALS Mw and dilute-solution [η]. The fit selects the requested viscometry temperature and requires consistent solutions and SEC measurement conditions across the selected rows. Use stable nominal mass-class labels across batches. Combine technical replicates before entry, keep independent preparations as distinct samples and record a retrievable primary source or lab record for each value. A minimal example uses four of the published rows.

```csv
sample_id,mass_class,Mw_g_mol,eta_mL_g,medium,viscosity_temp_C,sec_eluent,sec_temp_C,source
P1,10 kDa,11000,13.8,PBS,37,DMAc + 50 mM LiCl,50,10.1039/C8PY00255J S1/S2
P2,20 kDa,18800,19.0,PBS,37,DMAc + 50 mM LiCl,50,10.1039/C8PY00255J S1/S2
P3,40 kDa,43400,31.4,PBS,37,DMAc + 50 mM LiCl,50,10.1039/C8PY00255J S1/S2
P4,60 kDa,62500,39.6,PBS,37,DMAc + 50 mM LiCl,50,10.1039/C8PY00255J S1/S2
```

## Apply a calibrated component law to a distribution

For chemically equivalent independent chains at infinite dilution, dilute-solution viscosity contributions are additive by mass. Given a component law [η]ᵢ = K Mᵢᵃ, [η]mix = Σwᵢ[η]ᵢ and Mᵥ = (Σwᵢ Mᵢᵃ)^(1/a). K and a must describe that polymer architecture, solvent and temperature over the complete populated mass range. The GUI accepts number or mass fractions, distributions per log₁₀ M, and imported material number fractions. It lists each mass-weighted contribution. An empirical fit against whole-sample Mw does not establish the corresponding single-component law. Supply that calibration independently and avoid assigning the synthetic default to PEtOx. The viscosity-average molar mass definition uses mass fractions. [IUPAC definition](https://goldbook.iupac.org/terms/view/12220)

For example, two equally numerous components of 10 and 40 kDa have mass fractions 0.2 and 0.8. With illustrative K = 0.1 in (mL/g)/(g/mol)^0.5 and a = 0.5, each component has [η] of 10 and 20 mL/g, respectively. The mixture has [η] = 18 mL/g and Mᵥ = 32.4 kDa.
