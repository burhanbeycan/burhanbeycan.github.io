# Current release verification · 0.3.0

All 26 tests in this project passed in the suite run on 24 November 2026. The [release report](../../VALIDATION_TR.md) lists the installed-wheel and interface checks. The following notes also retain the earlier release record where identified.

# Verification record

## Chemistry and sequence additions

The 0.3.0 suite passed 26 tests locally with RDKit 2026.09.2, including all legacy cases. Independent DP 1 reference SMILES check the orientation of both termini. DP 6 elemental formulas, the azide motif and formal charge, expected proton/sodium/double proton ion shifts, aligned backbone coordinates, assigned-terminus transfers, and DP limits are tested.

The Poval 217 example recomputes the published ¹³C-NMR hydrolysis fraction, blockiness index and bulk OH/OAc mean run lengths from measured dyads. The comparison checks arithmetic against the original example. It does not validate a unique sequence or identify the physical chain-length distribution. Tests separately check composition-only Bernoulli variance, Markov dyads, reproducibility, ensemble boundaries and an acetate substructure count in the drawn model chain. PEtOx examples are calculated prospective masses. The repository contains no sample-specific experimental mass spectrum.

## Earlier graph constructor

Verified on 23 November 2026 using Python 3.12.14 and RDKit 2026.09.2 in the build environment. Command:

```bash
PYTHONPATH=../polymer-chemistry-lab/src:../polymer-research-studio/src python -m unittest discover -s tests -v
```

**14 tests passed, zero failures, zero skips** in the recorded local run. The checks cover:

- All six preset polymerized repeat formulas.
- Independent PE alkane formula, atom and bond counts across several DP values.
- Average and exact mass balances with 2H terminal correction.
- PEO/PVA isomeric connectivity and hydroxyl-count distinction.
- PEtOx N-propionyl motifs and one secondary terminal amide under H/H capping.
- Equal-composition block/alternating chains: equal mass, distinct graph and dyads.
- Explicit internal tetrahedral stereo against a separately written reference; correct terminal loss of a now-achiral site.
- No stereochemistry invented for unspecified PVA centres.
- Isotope-labelled mass balance.
- Invalid ports, radicals, disconnected fragments, input/atom limits and malformed sequences.
- Finite 2D graph coordinates, repeat lineage and cross-unit bond count.
- Homopolymer-only mass handoff, all GUI presets and runtime output contract.

The initial tests verify chemical bookkeeping and software behavior. They do not establish an experimental property model, synthetic accessibility, force-field validity or broad coverage of all stereochemical edge cases. The 2D renderer is shared with Workbench; its visual review is documented in that repository. The included GitHub Actions workflow is a configuration, not evidence that hosted CI has already run.
