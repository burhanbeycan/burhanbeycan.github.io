# Scientific specification

## AA+BB acyclic step growth

For nAA and nBB moles of exactly bifunctional monomers, define A0=2nAA, B0=2nBB and r=min(A0,B0)/max(A0,B0). The input p is the fraction of the *initial limiting groups* reacted. The number of intermolecular bond equivalents is b=p min(A0,B0). Every bond reduces the number of retained molecules by one, so n=nAA+nBB−b and Xn=(nAA+nBB)/n=(1+r)/(1+r−2rp).

Remaining groups are A=A0−b, B=B0−b. Their sum is 2n, independently checking the open-chain end count. Each original AA or BB monomer contributes one unit to Xn. The chemical paired-repeat count would use a different convention and is not silently substituted.

Input masses MA and MB describe the original monomer molecules. With loss L g/mol of bonds, initial mass is nAA MA+nBB MB, condensate mass bL, retained population mass m=initial−bL and Mn=m/n. This includes unreacted monomers and all retained oligomers. Removing residual monomers or isolating only a polymer fraction changes the population and its Mn. L=0 describes polyaddition; a chemically appropriate nonzero loss can describe condensation accounting. Validate 2L<MA+MB so the paired residue has positive mass.

At exact stoichiometric balance and p=1, n=0 in the ideal continuous population limit and Xn/Mn diverge. The implementation returns null with an infinite-chain-limit explanation instead of a floating-point infinity or a finite fabricated mass. A finite acyclic system cannot consume every end group; finite-size/cyclization corrections are outside this model.

For balanced ideal equal-reactivity linear step growth, the Flory chain-*length* distribution has Xw/Xn=1+p. The app reports this as a reference length dispersity only at exact balance. It is not generally equal to Mw/Mn with alternating unequal monomer masses, condensation mass corrections and end-group contributions. No dispersity is inferred for imbalanced feeds.

No branching, monofunctional stopper, cyclization, unequal reactivity, diffusion limitation, gelation or conversion-time model is implemented. The basic Carothers expression is literature theory, not a claim of a novel polymerization model.

## Ideal initiated-chain accounting

Given nM initial monomer moles, nI initiator moles, supplied conversion p, initiation efficiency η and s *separate chains per initiator molecule*, nchains=η s nI. Reacted monomer is p nM; DPn=p nM/nchains. Every incorporated monomer contributes the supplied residue molar mass M0. With a chemically derived combined terminal/core correction Mend per separate chain, Mn=DPn M0+Mend and total chain-population mass=p nM M0+nchains Mend.

All initiated chains, including zero-addition chains, are counted. This assumption matters at low conversion/DP. A star with multiple growing arms joined to one initiator core is one molecule, and cannot be represented by counting each attached arm as a separate molecular chain for Mn. Initiator and quencher atoms cannot be inferred from their concentrations: derive Mend using the explicit reaction/end-group convention.

The calculation supplies no kinetics, no conversion from time or temperature, no initiation-efficiency fit, no termination/transfer/coupling, no chain-length distribution and no GPC response curve. It cannot establish <10% prediction error for PEtOx or another polymer. Errors in supplied conversion and chain count pass directly into Mn.

## Terminal binary copolymerization

At one fixed free-monomer feed fA+fB=1, rA=kAA/kAB and rB=kBB/kBA. The instantaneous incorporated fraction is FA=(rA fA²+fA fB)/(rA fA²+2fA fB+rB fB²), FB=1−FA. These ratios are supplied at the intended mechanism and conditions; no values are assigned automatically to chemical names.

Transition probabilities are PAA=rA fA/(rA fA+fB), PAB=fB/(rA fA+fB), PBA=fA/(fA+rB fB), PBB=rB fB/(fA+rB fB). Stationary ordered dyads are FAA=FA PAA, FAB=FA PAB, FBA=FB PBA, FBB=FB PBB. FAB=FBA and their sum with like dyads is one. Mean run lengths are 1/PAB and 1/PBA. These are stationary infinite-chain statistics at the fixed feed, not exact finite-oligomer sequence measurements.

Strict interior feeds avoid absent-monomer undefined transitions. Both r=0 is the alternating stationary limit; both r=1 is a Bernoulli chain with FA=fA. The model assumes terminal reactivity, stationary propagating-end statistics and irreversible additions. Penultimate effects, slow-initiation transients, depropagation, feed depletion and cumulative batch composition require different models. In particular, applying the instantaneous equation directly to high-conversion bulk composition can bias estimated reactivity ratios.

## EtOx recipe kinetics and SEC benchmarks

The EtOx-specific calculation uses temperature, time and initial MeOTs concentration to model conversion. Its source coefficients, terminal masses, cohort definitions and validation design are detailed in [the PEtOx study](PETOX.md).
