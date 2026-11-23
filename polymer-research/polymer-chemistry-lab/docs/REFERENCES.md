# References and implementation provenance

Official documentation was checked on 23 November 2026. Reference links identify the external algorithms; the finite-chain protocol, input limits, sequence ledger, handoff and chemistry scope are implemented in this repository.

1. [RDKit: Getting Started in Python](https://www.rdkit.org/docs/GettingStartedInPython.html). Molecular graphs, sanitization, isomeric SMILES and graph editing.
2. [RDKit rdmolops](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html). `molzip`, attachment-point maps, H handling and stereochemical assignment.
3. [RDKit rdMolDescriptors](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) and [Descriptors](https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html). Formula, exact molecular mass and average molar mass conventions.
4. [The RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html). Representation and interpretation of stereochemistry.
5. Chiu, T. T.; Thill, B. P.; Fairchok, W. J. *Poly(2-ethyl-2-oxazoline): A New Water- and Organic-Soluble Adhesive.* In *Water-Soluble Polymers: Beauty with Performance*; Advances in Chemistry 213; American Chemical Society, 1986; Chapter 23. [DOI:10.1021/ba-1986-0213.ch023](https://doi.org/10.1021/ba-1986-0213.ch023). Chemical identity of the N-propionyl-substituted ethyleneimine repeat. Publisher retrieval was restricted during this build; metadata and the publisher-indexed abstract were available to the scientific review.

No RDKit algorithm is claimed as newly invented here. The software makes the polymer-specific distinction between polymerized repeat contributions, finite capped molecules, explicit sequences and sample-level distributions operational and reviewable.
