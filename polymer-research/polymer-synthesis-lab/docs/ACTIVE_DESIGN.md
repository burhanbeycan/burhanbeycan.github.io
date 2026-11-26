# Selecting the next EtOx/MeOTs/MeCN measurement

At fixed reactor, assay, terminal chemistry, temperature and hold time, choose from a user-supplied pool of feasible concentration/feed combinations. Existing observations come from the same CSV schema as the [PEtOx kinetic audit](PETOX.md). Assign repeated samples to their original preparation, study and series. Use actual initial solution molarity. Define feed as the initial EtOx:MeOTs molar ratio. Relative SEC masses retain their eluent, detector and calibration basis.

## Design criterion

Define a dimensionless feature vector $\phi=(1,\log_2([M]_0/(3\,\mathrm{M})),\log_2(\mathrm{feed}/60))$. The fixed reference levels specify a coordinate convention. Each represented reaction series receives total weight one within the selected condition. A series with $n_s$ measurements contributes $1/n_s$ for each feature vector. Construct $A=I+\sum_i w_i\phi_i\phi_i^T$. Rank each unmeasured pool member by $g=\tfrac12\ln(1+\phi^TA^{-1}\phi)$. By the matrix determinant lemma this is the half log-determinant gain for one new independent measurement with unit standardized precision. The identity prior and equal measurement precision are explicit design choices. The score quantifies information under this local log-linear design assumption. It carries no units of prediction error. The ranking depends on candidate bounds and the prior. Check process feasibility before admitting a candidate. [Optimal experimental design](https://doi.org/10.1017/S0962492924000023)

## Conditional residual model

Fit the cohort-specific quantity $\ln(M_{n,\mathrm{SEC}}/M_{n,\mathrm{kin}})=\phi^T\beta$ by series-weighted least squares. Coefficients describe an assay-specific empirical intercept and the incremental change in log ratio for a doubling of concentration or feed. Neither observed conversion nor candidate Mn enters acquisition. Observed conversion stays available as an independent diagnostic of the kinetic assumptions.

Fit eligibility requires a minimum of four independent preparation series, four paired recipe coordinates, variation in concentration and in feed, a full-rank design and weighted design condition number at most 1000. These are release-specific safeguards against unidentified coefficients and unstable matrices. They do not establish scientific validity. Prospective adjusted SEC estimates are shown only inside the convex hull of the matched observed designs in the two log coordinates. Predictions outside that hull are displayed under the fixed kinetic baseline. No predictive interval is claimed from either bundled cohort. If different series share sample material, enter a common series ID before validation.

## Independent predictions

The audit withholds entire publications and reports a separate stringent transfer test that withholds both concentration and feed, plus any preparation series represented in the test group. Each residual surface uses only the remaining measurements at the selected condition. Holdouts can test extrapolation. Results show eligible training and scored test counts, median and 90th percentile absolute relative error, and a kinetic baseline on exactly the same scored rows as the residual. An empty fold leaves the residual result blank. The data from studies supplying the kinetic constants are excluded by the CSV parser.

## Current evidence

At 160 °C and 5 min the four published flow observations all used 4 M EtOx. In the supplied pool, 2 M, feed 120 ranks first with gain 0.3939 nat. The corresponding model input [MeOTs]₀ is 0.01667 M. The nominal kinetic mass of 8452 g/mol uses an Arrhenius extrapolation beyond the kinetic source's 60–140 °C measurement range. No published kinetic measurement or SEC measurement is assigned to this candidate. Single-paper and joint holdouts score four kinetic baseline rows (median error 49.01%) and zero residual estimates. At the already measured 4 M the earlier one-scale leave-feed-out comparison gave median error 28.54%. A second slope cannot be learned for concentration until concentration varies independently. [Baeten et al., Table S1](https://doi.org/10.1039/C5CC04319K) · [Kinetic coefficients, Table 2](https://doi.org/10.1021/acs.macromol.2c01930)

After an independent measurement, append its traceable SEC Mn to the measured CSV and rerun. Optional conversion, the full chromatogram and authentic end-group evidence help discriminate propagation, branching and measurement bias. The apparent-mass residual alone cannot distinguish these mechanisms. [Arráez et al., SEC evidence for macropropagation](https://doi.org/10.1021/acs.macromol.9b00544)


## Visible recovery control

Open **PEtOx · synthetic residual recovery**. Two explicitly simulated studies use the same nine paired recipes at 100 °C and 10 min. Their synthetic log residual has known coefficients (0.10, 0.20, −0.30). Inspect the recovered coefficient table and group coverage. An unmeasured candidate inside the observed log-coordinate hull receives a residual estimate. The candidate at 5 M, feed 90 receives only a theoretical kinetic value. This exact recovery exercise tests computation and provenance handling; all displayed masses in this preset are synthetic.
