# Recipe-only PEtOx molar-mass prediction

In the Synthesis Lab select **PEtOx · CROP from feed, temperature and time**, or start with either measured-cohort preset. Concentrations refer to total initial solution volume. The model assumes EtOx homopolymer, MeOTs in acetonitrile, constant mixture volume, an isothermal hold, immediate quantitative initiation and a constant number of propagating chains.

## Kinetics and dry-chain mass

For time in seconds and temperature in kelvin, use $k_p=A\exp[-E_a/(RT)]$, $X=1-\exp(-k_p[I]_0t)$, $DP_n=X[M]_0/[I]_0$ and $M_n=99.133DP_n+M_{ends}$ in g/mol. The combined terminal masses are 32.042 g/mol for methyl/hydroxyl and 57.056 g/mol for methyl/azide under the uncapped C5H9NO repeat convention. The full-incorporation comparison sets $X=1$ irrespective of time. The benchmark scores each calculation against the particular reported Mn method in its own cohort.

The MeCN Arrhenius fit uses $E_a=71.6$ kJ/mol and $A=1.06\times10^8$ L/(mol s), from Table 2. It yields 0.01006 L/(mol s) at 100 °C, while Table 1 lists 0.0131 L/(mol s) for the individual experiment. The fit uses five kinetic temperatures and its central prediction need not coincide with each measured point. The MeCN series in this paper is attributed to the original 2015 work. [Arrhenius fit](https://doi.org/10.1021/acs.macromol.2c01930) · [Original kinetics](https://doi.org/10.1021/acsmacrolett.5b00392)

## Independent observations

| Cohort | Conditions and measurement basis | Median / P90 absolute relative Mn error |
|---|---|---|
| Batch, four observations from one time course | 82 °C, 1–6 h, feed 113, azide termination. Crude PEtOx measured by RI SEC in DMAc with PMMA calibration. | Fixed kinetics **10.30% / 16.75%**. Full incorporation **54.06% / 176.18%**. |
| Microflow, four feed ratios in one study | 4 M EtOx, 160 °C, 5 min residence, feed 15/30/60/120. SEC-RI in DMAc–LiCl at 50 °C, PMMA standards. | Fixed kinetics **49.01% / 56.50%**. Leave-one-feed-out assay scaling **28.54% / 46.66%**. |

The batch recipe provides 24.6 mg MeOTs, 1.3 mL EtOx, 1.9 mL MeCN and nominal feed 113. Using 186.23 g/mol for MeOTs and 3.2 mL nominal volume gives $[I]_0\approx0.04128$ M. The score uses $[M]_0=113[I]_0$. No final solution volume or uncertainty was reported. Concentrations are reconstructions. Figure S2D gives the four rounded SEC masses; pair them with the listed sampling times in increasing conversion order. [Batch experiment, SI pp. S2 and S5](https://www.rsc.org/suppdata/c9/cc/c9cc00407f/c9cc00407f1.pdf)

The flow study reports four concentration ratios and their corresponding SEC masses directly in Table S1. We use a nominal Me/OH mass term after the reported KOH/MeOH quench. End groups were not independently assigned for these four samples. The calculation equates residence time to a uniform isothermal exposure. Reactor dispersion and SEC calibration can contribute to discrepancies. One cannot assign the error to either effect uniquely. [Microflow study, SI pp. S2–S5](https://www.rsc.org/suppdata/c5/cc/c5cc04319k/c5cc04319k1.pdf)

## Add and audit a cohort

Copy either [batch CSV](../src/polymer_synthesis_lab/data/batch_azide.csv) or [flow CSV](../src/polymer_synthesis_lab/data/flow_feed.csv) into **Experimental PEtOx observations**, preserving the header. Record the originating paper, reaction series, sample, reactor, Mn assay (solvent, calibration and detector), end-group convention, initial molarities, temperature in °C, hold in minutes and observed Mn in g/mol. Optional measured conversion is a fraction in [0,1] used for display only. Each CSV represents the same chemistry, reactor mode and measurement basis. Assign each reprinted observation its original source and deduplicate it before import.

The full-conversion and fixed-kinetics models use no test or training outcomes. For the third model, each eligible training series supplies its median ratio of measured to kinetic Mn. The median across these series gives one assay-scale factor. Apply that factor to the held-out kinetic Mn. Every model, including apparent-Mn scaling, is scored against the same measurements. The factor is an empirical correction tied to the cohort. It cannot separate initiation efficiency, chain transfer or instrumental bias. Report the number of scored rows when some folds lack training data.

Withhold by publication, reaction series, concentration or feed. The joint option excludes every training row sharing either the held-out monomer concentration or the held-out feed ratio. The bundled time course and flow series each cover one initial monomer concentration, one assay and one publication. Neither dataset supports within-cohort testing on unseen concentrations. Differences in reactor, quench and SEC conditions preclude forming a matched cross-paper benchmark by pooling these files. Four samples per cohort do not support an accurate population-tail estimate. The observed errors and sample counts are narrow diagnostics.

The software omits free fits for initiation efficiency, transfer and branching. Relative SEC Mn and conversion do not uniquely identify these effects. Arráez et al. used full SEC traces to constrain macropropagation. Additional matched independent data and identifiable parameters are needed for those mechanisms and for the proposed sub-10% cross-condition target. [SEC constraints on macropropagation](https://doi.org/10.1021/acs.macromol.9b00544)

## Source inventory outside the two scored Mn cohorts

| Source | Relevant evidence | Use in this release |
|---|---|---|
| [Vergaelen et al., 2015](https://doi.org/10.1021/acsmacrolett.5b00392) and [2023](https://doi.org/10.1021/acs.macromol.2c01930) | Same underlying MeCN kinetic experiments, reported again in the 2023 comparison. The authors excluded the earlier MeCN time-series measurements from the Mn-versus-conversion panels in Fig. 1. | One coefficient provenance chain. These papers are excluded from independent audit cohorts. |
| [Arráez et al.](https://doi.org/10.1021/acs.macromol.9b00544) | Microwave CROP across temperatures at target DP 100; SEC distributions used to constrain branching. | Mechanistic reference. Conversion-only digitizations cannot be counted as Mn measurements. |
| [Weber et al.](https://doi.org/10.1163/156855509X412090) | A-series hydroxy oligomers at diverse feed ratios, PS-calibrated SEC. A 1–5 min hold range is specified rather than an unambiguous hold for each Table 1 sample. | Useful end-group comparison. Excluded from time-resolved benchmark. |
| [Wendler et al.](https://doi.org/10.1039/C6PY00033A) | Azide-terminated PEtOx at 4 M measured with PS-SEC and MALDI. Work-up and mass basis differ from the crude DMAc-RI batch cohort. | Excluded. Converted and maleimide-functionalized materials must not be relabelled as the precursor PEtOx-N3. |
