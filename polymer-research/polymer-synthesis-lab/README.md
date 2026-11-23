# Polymer Synthesis Lab

**Functional groups, chain counts, end groups and copolymer sequence statistics.**

A local scientific application by **Burhan Beycan** for explicit polymer-chemistry bookkeeping, independently launchable and connected to Polymer Research Studio through its shared workbench. Version **0.3.0**.

## Implemented chemistry

| Calculation | What is calculated | Necessary boundary |
|---|---|---|
| PEtOx CROP in MeCN | Recipe → conversion and conditional Mn from published MeOTs kinetics. Compare two experimental SEC cohorts and audit custom CSV by entire study, series, concentration, feed or both | Relative SEC masses, end-group identity, reactor and concentration basis must be tracked |
| AA+BB step growth | Limiting/excess functional-group ratio, Carothers Xn, remaining A/B groups, condensate and whole-population Mn | Exactly bifunctional, acyclic, A–B reactions; conversion supplied |
| Ideal chain growth | Initiated-chain count, incorporated residues per chain, end/core contribution and conditional theoretical Mn | Supplied conversion/efficiency, fixed chain count, no transfer or termination |
| Terminal binary copolymer | Instantaneous composition, stationary ordered dyads, stationary mean run lengths | Supplied reactivity ratios at one fixed feed; no cumulative conversion drift |

The AA+BB **Xn counts original AA and BB molecules separately**. It is not the number of paired AA+BB chemical repeat units. The full mass balance therefore uses individual molecule masses and subtracts the stated byproduct mass per bond. Unreacted monomers remain included in this model population.

The balanced Flory **length** dispersity is reported only at exact AA/BB balance and is explicitly distinguished from molar-mass dispersity when monomer masses or terminal masses differ. A multiarm star joined to one core is one molecule: the “chains per initiator” field counts separate chains and must not be used to turn arms into molecules.

## Install and launch

Python 3.10+ and sibling source folders:

```bash
git clone --single-branch --branch polymer-research-review-v0.2.0 https://github.com/burhanbeycan/burhanbeycan.github.io.git
cd burhanbeycan.github.io/polymer-research
python -m pip install ./polymer-research-studio ./polymer-synthesis-lab
polymer-synthesis-lab
```

Open `http://127.0.0.1:8770`. Use `--no-browser` to suppress automatic opening or `--port 8780` to choose another loopback port. Run `polymer-workbench` to see all laboratories installed in the same environment. No API key is required.

## Chemical examples and evidence

The nylon-6,6 example uses approximate adipic-acid and hexamethylenediamine molar masses and water loss per formed bond. Its conversion of 0.99 is a stipulated ideal condition, not measured experimental data or a prediction of a reaction recipe. Other presets illustrate analytical limits and use synthetic coefficients.

The PEtOx module includes two small, separately scored experimental datasets. Read [the conditions, data sources and numerical results](docs/PETOX.md). The recipe-only model exceeded the target of 10% median absolute relative error on both cohorts. End-group masses must be derived for the actual finite-chain convention. Chemistry Lab provides explicit finite oligomer formulae; a single Mn does not determine a full molar-mass distribution.

## Verification

```bash
cd polymer-synthesis-lab
python -m unittest discover -s tests -v
```

See [SCIENCE](docs/SCIENCE.md), [REFERENCES](docs/REFERENCES.md), and [VALIDATION](docs/VALIDATION.md). Numerical kernels are independent of the GUI. General numerical plotting grids use NumPy; reaction balances and terminal-model equations are implemented here. MIT license; developed with AI assistance and intended for scientific review.
