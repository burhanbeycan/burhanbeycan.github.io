# Chemical model and interpretation

## Scope and molecular identity

The constructor accepts one connected polymerized repeat graph for A and, when used, a second for B. Two degree-one neutral nonisotopic dummies have maps 1 and 2 and single bonds to real atoms. Internal side groups and rings are allowed; the *sequence of repeat units* is open and linear. Ports are the only atom maps. Radical templates, unresolved ports, disconnected counterions, invalid valence and non-single port bonds are rejected.

Each input repeat is copied without a chemical reaction transformation. Tail of unit i joins head of unit i+1 through RDKit molzip. Two separate H/port fragments cap the global head and tail. Sanitization verifies the resulting graph; an independent additive mass balance must agree. This verifies bookkeeping, not synthesis feasibility. The custom template may carry charge, but the output does not model counterions, solvation, acid–base equilibria or ion association.

## Chemical examples

PEtOx uses the sequence of backbone atoms –CH₂–CH₂–N(C(=O)CH₂CH₃)–. Its repeat formula is C₅H₉NO. H capping introduces a carbon-terminal H and N-terminal H, so the terminal N is a secondary amide in this deliberately simplified graph. The interior PEtOx nitrogens remain tertiary amides. Neither end is inferred from methyl tosylate initiation or a termination reagent. PMeOx replaces the propionyl substituent by acetyl.

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

The explicit sequence contains n units and n−1 directed adjacent pairs. Unit fractions are N_A/n and N_B/n. Each dyad XY has fraction N_XY/(n−1). For n=1 the dyad fraction is undefined and exported as null. A block chain and an alternating chain with identical composition can have equal mass but different dyad counts. This module does not estimate reactivity ratios, conversion-dependent composition drift or a distribution of chains.

## Interoperability

A homopolymer may transfer average repeat mass and H/H end correction to Chain Lab. No mass distribution, stiffness or contour increment is transferred. For a copolymer, a common repeat mass is not inferred. Material metadata preserve the explicit sequence, mapped repeats, finite graph SMILES and stated cap assumption. Open linear repeat sequence does not imply that pendant rings are absent.
