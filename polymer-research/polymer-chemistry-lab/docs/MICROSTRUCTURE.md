# NMR-constrained PVA and PEtOx/EI microstates

## Experimental PVA comparison

McCrary et al. measured acetylated PVA prepared by homogeneous solution acetylation and heterogeneous acetylation of semicrystalline films in supercritical CO₂. We include the six newly prepared grades for which Table 1 reports a numeric block character η. Labels such as rAcPVOH16 retain the nominal acetyl percentage encoded in the sample name. η and run lengths are taken from Table 1. The commercial Poval comparison comes from a different preparation route. Block character was unavailable for bAcPVOH5 and bAcPVOH25. [Primary experimental paper, Table 1 and Eqs. 1–3](https://doi.org/10.1021/acsapm.6c00317)

Each fold withholds one whole grade. The predictor uses the preparation route and the mean η of the other two grades from that route. The Bernoulli reference fixes η to one for every grade. Neither predictor uses held-out η, dyad intensities or held-out run lengths. Mean absolute errors are reported in η units across all six folds. This small within-study comparison does not establish generalization to a new study or route.

A separate fixed-composition comparison uses η = 0.95 from rAcPVOH16 and η = 0.30 from bAcPVOH15. At a common model fraction of acetate x, the mixed dyad expectation is 2ηx(1–x). The stationary first-order Markov ensemble then samples a chain population for each measured η. At x = 0.15 the boundary probabilities are 0.24225 and 0.07650. These are computed scenarios. The actual materials have different nominal compositions, and input-constrained dyads do not provide an independent experimental endpoint. Tacticity, water coordination, water transport and crystallinity require their own records or models.

## PEtOx hydrolysis and ion states

The chain is oriented from α-methyl to ω-hydroxyl. Type one symbol per backbone repeat:

| Symbol | Chemical state | H-bond site |
|---|---|---|
| E | N-propionyl ethyleneimine, intact EtOx-derived repeat | Side-chain carbonyl oxygen accepts H bonds |
| I | Neutral linear EI (secondary amine) | Secondary amine donates and accepts |
| P | Singly protonated linear EI | Secondary ammonium donates and carries +1 |

The covalent model gives hydrolysis fraction (I + P)/DP and protonation fraction P/(I + P). The latter is undefined for an intact chain. Each P needs one monovalent counterion (chloride, bromide or acetate). The counterion is shown once with its multiplicity. Counterion masses are excluded from the covalent-chain mass, and the salt-entry SMILES includes all counterions. Graphs for intact, hydrolysed neutral and charged structures share the same DP and termini. Repeated identical scenarios appear only once. Different patterns with identical overall composition can have different covalent connectivity.

This is an explicit microstate calculation, with protonation supplied site by site. A measured pH alone does not fix the microscopic pattern in a polyelectrolyte. The user should record sample-specific pH, ionic strength, titration, chain length and termination evidence. Morgan et al. used NMR and tandem MS to characterize partially hydrolysed oligomers and reported end-group-dependent behavior under their conditions. Their measured sequences are not assigned to a WBPU sample here. [Primary oligomer study](https://doi.org/10.1039/D2PY00437B)

For a WBPU blend simulation, retain the input pattern, ion counts, solvent, mass basis and force-field version together. Equivalent polymer mass fractions can imply different salt concentrations and repeat counts. Compare prepared replicas and local hydration under defined experimental conditions. The application provides molecular identity and stoichiometric input for such comparisons. It does not generate validated polymer/ion force-field parameters or a trajectory. [Beycan, MACRO 2024 abstract, p. 176](https://static1.squarespace.com/static/62c326b96a09b762cf8302ac/t/667e8cd4578728784d3601c9/1719569633844/AbstractBook-MACRO2024Final.pdf)
