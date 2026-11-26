# Waterborne polyurethane characterization and validation

## Functional-group ledger

The GUI accepts one row for each component and addition stage. Use the mass of dry active material and a matching basis for every assay. With nominal molar mass, the calculation uses the number-average molar mass of a polymeric feed and its number-average reactive functionality. Unit conversions below give mmol of reactive groups for mass in grams.

| Basis | Required quantity | Calculation |
|---|---|---|
| `molar_mass` | g/mol and reactive functionality $f$ | $1000mf/M_n$ |
| `OH_number` | mg KOH/g of polyol | $m I_{OH}/56.1$ |
| `NCO_wt_pct` | percent NCO on the stated mass basis | $1000 m w_{NCO}/42.02$ (mass fraction) |
| `equivalent_weight` | g per mol of reactive groups | $1000m/EW$ |
| `net_urea` | independently quantified mass of water reacted to CO₂ and urea | $2000m/18.01528$ consumed NCO equivalents |

One primary or secondary amine **site** contributes one NCO equivalent in ordinary urea formation. Tertiary amine neutralizer contributes no such site. The net-water basis applies only to the amount that underwent both reactions in the net urea pathway. Dispersion water generally dominates the liquid mass and is excluded. Enter multifunctional compounds carrying several different reactive groups as separate group contributions and label the basis for each. The output sums group equivalents; the displayed component masses should not be summed across duplicate entries.

The recipe preset uses one gram each of nominal PEG 20,000 and PCL 1,000 diols. Their contributions are 0.10 and 2.00 mmol OH. BDO contributes another 1.00 mmol with the selected example mass, giving 3.10 mmol total OH. The 0.35 g example IPDI charge supplies 3.1490 mmol NCO. Index is 101.58. These illustrative masses support an equivalent-balance check. Conversion, allophanate formation and elastically active strand density require independent evidence. Distinguish prepolymer assays from total feed balance. [Formulation calculations](https://www.dollmar.com/wp-content/uploads/2022/09/Brochure-Polymeg.pdf)

## Operational gel, sol and swelling

Use the same unfilled specimen throughout the extraction. Record the initial dry mass $m_0$, blotted swollen residue $m_w$ and mass redried to constant weight $m_d$. The display reports four mass ratios derived from these weighings.

| Quantity | Expression (%) |
|---|---|
| Insoluble fraction under the protocol | $100m_d/m_0$ |
| Extractable fraction | $100(1-m_d/m_0)$ |
| Net swelling on initial dry mass | $100(m_w-m_0)/m_0$ |
| Solvent on retained dry gel | $100(m_w-m_d)/m_d$, defined for $m_d>0$ |

A 100 mg specimen blotted at 220 mg and redried to 80 mg gives 80% insoluble fraction, 20% extractables, 120% initial-mass swelling and 175% solvent on retained gel. Unknown solvent quality, incomplete extraction and persistent crystalline or hydrogen-bonded domains affect interpretation. Establish matrix solubility using a chemically comparable linear control. Record ionic form and neutralization, fillers, thermal history, extent of cure and evidence for hard-domain or crystallite melting. The free-energy solver's dry reference and `νe` require independent justification for each material.

A DMA sweep describes tensile storage modulus over temperature at the recorded frequency. Record strain amplitude, linear viscoelastic range, conditioning and heating rate alongside the measurements. Thermoplastic polyurethanes can exhibit an apparent rubbery plateau supported by physical junctions. A finite-frequency plateau on its own cannot assign covalent topology. The module retains DMA results in their measured units. [Segmented TPU study](https://doi.org/10.3390/polym11121910)

## Experimental holdout

Saalah et al. report freestanding jatropha-based anionic WBPU films. The net swelling target comes from a seven-day toluene immersion at 23 °C. The authors averaged three measurements. The protocol defines net swelling relative to initial dry mass. Table 4 supplies the study's hard-segment contents and Table 6 supplies tensile storage modulus at 25 °C and 1 Hz, independently of the extraction mass. Tables 3, 4 and 6 use JPU 188; Table 1 prints JPU 181. The present benchmark uses matching labels from the former tables. [Primary data and methods](https://doi.org/10.3390/polym13050795)

| Sample | Hard segment (wt%) | E′ (MPa) | Measured net swelling (%) | Observed insoluble fraction (%) |
|---|---:|---:|---:|---:|
| JPU 161 | 38 | 92.8 | 231.9 ± 7.6 | 51.2 ± 5.1 |
| JPU 188 | 41 | 247.7 | 167.8 ± 17.9 | 55.4 ± 8.3 |
| JPU 217 | 45 | 373.8 | 128.2 ± 14.9 | 74.2 ± 1.1 |

Every fold withholds one entire formulation. The remaining formulation means have equal weight in training. The fold error uses percentage points of initial dry mass. Comparing the training-set mean, linear hard-segment model, and linear model against `ln(E′/1 MPa)` gives MAE of **55.93, 33.31, and 17.36 percentage points**, respectively. Both endpoints require extrapolation. These are within-paper exploratory results from three films, with two observations to fit each single-variable fold. Modulus is an empirical covariate that brings an additional measurement cost. The comparison cannot establish that modulus generalizes to new WBPU families. The same experiment's gel fraction or swelling-derived `νe` would leak the outcome if used as predictors here.

For a new series, specify specimen family and processing route, a common definition of hard-segment wt%, solvent, temperature, exposure time and full conditioning protocol. Replace the dataset source and DMA protocol. Multiple rows per formulation are averaged across specimens before the fold assignment. With at least four formulations and full-rank training designs, the tool also evaluates `S=a+b w_hard+c ln(E′/1 MPa)`. Relative humidity, ionic form, casting substrate and sample thickness should be specified in the experimental protocol when relevant. Data from separate protocols require separate runs. A paper-level external test is a distinct validation exercise.

The measured preset concerns jatropha WBPU cast films. For PEG/PCL WBPU, enter the original sample records and keep film and fiber-mat preparation routes distinct. [Beycan et al., electrospun WBPU research](https://doi.org/10.1021/acs.iecr.4c00759)
