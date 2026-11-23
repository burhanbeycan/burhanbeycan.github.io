# Chemical model and interpretation

## Scope and molecular identity

The general H/H constructor accepts one connected polymerized repeat graph for A and, when used, a second for B. Two degree-one neutral nonisotopic dummies have maps 1 and 2 and single bonds to real atoms. Internal side groups and rings are allowed; the *sequence of repeat units* is open and linear. Ports are the only atom maps. Radical templates, unresolved ports, disconnected counterions, invalid valence and non-single port bonds are rejected.

Each input repeat is copied without a chemical reaction transformation. Tail of unit i joins head of unit i+1 through RDKit molzip. Two separate H/port fragments cap the global head and tail. Sanitization verifies the resulting graph; an independent additive mass balance must agree. This verifies bookkeeping, not synthesis feasibility. The custom template may carry charge, but the output does not model counterions, solvation, acid–base equilibria or ion association.

## Chemical examples

PEtOx uses the sequence of backbone atoms –CH₂–CH₂–N(C(=O)CH₂CH₃)–. Its repeat formula is C₅H₉NO. H capping introduces a carbon-terminal H and N-terminal H, so the terminal N is a secondary amide in this deliberately simplified graph. The interior PEtOx nitrogens remain tertiary amides. The PEtOx comparison mode explicitly substitutes the α and ω sites as described below. PMeOx replaces the propionyl substituent by acetyl.

PEO and PVA both contribute C₂H₄O per chosen repeat; their backbone connectivity differs. H-capped PEO contains one terminal OH, while each PVA repeat contributes a pendant OH. PMMA uses a carbon backbone with pendant methyl and methyl ester groups; stereocentres are unspecified in the preset. These examples are structural templates rather than measured material records.

## Formula and mass conventions

The two dummy atoms carry zero mass and reserve valence. Mass is computed while those ports are present, so the repeat contribution does not gain terminal implicit H. Removing ports before calculating mass would instead create the capped molecule and overestimate the repeating contribution by 2H.

For unit sequence s₁…sₙ:

\[
M_{\mathrm{avg}}=\sum_{i=1}^n M_{s_i,\mathrm{avg}}+2M_{H,\mathrm{avg}}.
\]

RDKit's average atomic-weight convention is used unless an isotope is explicitly labelled. With the tested version, the terminal correction is 2.016 g/mol. Exact molecular mass is reported separately in Da and uses explicit isotopes or most-abundant isotope masses. Formula is generated with isotope separation. The software records the RDKit version because atomic-weight tables and algorithms are dependencies.

For a homopolymer, the corresponding graph-based mass-to-repeat relation is:

\[
DP=(M-2M_H)/M_0.
\]

That equation is conditional on these H/H caps. For a real sample, actual end groups and molecular mass basis must be supplied. A sample with one oligomer structure per row could be used to construct a distribution, but the current application does not claim to infer one from a single SMILES. No Mn, Mw or dispersity is estimated here.

## Stereochemistry and drawing

Input is isomeric SMILES. molzip substitutes attachment points and performs local stereochemical bookkeeping. The product is sanitized; stereochemistry is reassigned and exported in canonical isomeric SMILES. If capping removes a distinction between two substituents, a former template stereocentre may cease to be stereogenic. That is a molecular identity change, not loss of a valid stereocentre. Tests compare a fluorinated chain against an independently written isomeric reference.

The constructor does not generate isotactic, syndiotactic or atactic sequences. In particular, unspecified stereochemistry is not equivalent to an atactic ensemble. Unit ordering and stereochemical ordering are separate questions. Explicit repeat templates can encode stereochemistry, but no stereoselectivity mechanism is inferred.

RDKit generates 2D coordinates and wedge/dash cues. The molecular view exports atoms with element, label, isotope, charge, optional CIP annotation and originating unit index; bonds include order, stereo and direction. These are connectivity drawing coordinates in arbitrary units. They must not be used as atomistic initial coordinates, Rg, end-to-end distance, a force-field topology, optimized geometry or evidence of tacticity.

## Sequence statistics

The explicit sequence contains n units and n−1 directed adjacent pairs. Unit fractions are N_A/n and N_B/n. Each dyad XY has fraction N_XY/(n−1). For n=1 the dyad fraction is undefined and exported as null. A block chain and an alternating chain with identical composition can have equal mass but different dyad counts. The explicit-graph calculation gives dyads of one chain. The PVA analysis below uses a separate population model.

## Interoperability

A homopolymer may transfer average repeat mass and H/H end correction to Chain Lab. No mass distribution, stiffness or contour increment is transferred. For a copolymer, a common repeat mass is not inferred. Material metadata preserve the explicit sequence, mapped repeats, finite graph SMILES and stated cap assumption. Open linear repeat sequence does not imply that pendant rings are absent.

## PEtOx end-group graphs and mass spectra

The mapped PEtOx repeat runs from backbone carbon to nitrogen, opposite to the propagation direction. The methyl initiator fragment forms a covalent N–C bond at the α nitrogen. The terminal backbone carbon bears OH or an organic azide, represented as C–N=[N⁺]=[N⁻]. Each extra heavy-atom bond replaces one terminal H. The graph additions relative to the uncapped repeat sum are H₂, CH₄O and CH₃N₃ respectively. These graphs have the same DP, backbone molecular identity and formal net charge zero. RDKit's sanitization and explicit atom correspondence preserve a common backbone orientation in the drawings. The α/ω annotations refer to propagation, independently of the drawing direction.

Positive-ion targets use the calculated neutral monoisotopic mass. [M+H]⁺ adds 1.007276 Da, [M+Na]⁺ adds 22.989218 Da (conventional LC–MS adduct shifts), and [M+2H]²⁺ divides the doubly protonated mass by two. Ion abundances, ionization sites, isotope envelopes, reaction conversion and stability are outside this calculation. Precursor mass cannot distinguish isomeric end groups with the same elemental composition.

## NMR-constrained PVA sequence population

A is –CH₂–CH(OH)– and B is –CH₂–CH(OCOCH₃)–. For normalized quantitative dyads O = OH–OH, X = OH–OAc in either orientation and V = OAc–OAc, the bulk alcohol fraction is h = O + X/2. Equal forward and reverse interface flux in a stationary binary chain gives transition probabilities P(B|A) = X/(2h) and P(A|B) = X/[2(1−h)]. Long-chain mean run lengths are the reciprocals of those transition probabilities. For 0<h<1 the index η = X/[2h(1−h)] equals one for a Bernoulli sequence, is smaller for block-rich dyads, and is larger for alternating enrichment. Hydrolysis fraction alone supplies a Bernoulli *reference* when no pair measurements exist.

The NMR preset transcribes the primary dyad fractions of Poval 217. The ensemble samples a stationary first-order Markov model consistent with the measured dyads in expectation. It cannot reconstruct higher-order sequence correlations, actual chains, experimental DP, polydispersity or stereochemistry. The seed, number and length of chains are recorded. Changes to the default dyad constraints require a sample-specific citation or an explicit hypothetical-data label. Quantitative integration, signal assignment, normalization and chain-end contributions determine whether an NMR input is appropriate.

In a pure population η is undefined. With no mixed dyads a mixture of two homopolymer populations is an admissible stationary limit. The simulation draws each chain's initial state from the bulk composition and no switches occur. For the pure/zero-interface limits, the run-length table describes the limiting model separately from finite simulated sequence lengths.
