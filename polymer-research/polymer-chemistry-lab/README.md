# Polymer Chemistry Lab

**Explicit repeat connectivity, finite chains, sequence chemistry, and end-group-aware mass bookkeeping.**

A local scientific application by **Burhan Beycan** for a polymer-chemistry research workflow. It builds atom/bond graphs from mapped polymerized repeat units and connects to [Polymer Research Studio](../polymer-research-studio). No PolyChemMD installation is required.

## Compare PEtOx termini at fixed DP

Select **PEtOx · same DP, three end groups**. The workbench shows the same backbone with H/H model caps and the two literature end pairs expected from methyl tosylate initiation with hydroxyl or azide termination. Inspect the N-bound α-methyl and the carbon-bound ω substituent, then compare neutral mass with the selected ion composition. The structural assignments are grounded in the end-group and tandem-MS study by [Morgan et al.](https://doi.org/10.1039/D2PY00437B).

![Three aligned PEtOx trimer structures with their α and ω sites highlighted](docs/assets/petox-ends.svg)

| DP 6 candidate | Formula | Neutral exact mass (Da) | [M+H]⁺ (m/z) |
|---|---|---:|---:|
| H/H reference | C₃₀H₅₆N₆O₆ | 596.426134 | 597.433410 |
| α-Me / ω-OH | C₃₁H₅₈N₆O₇ | 626.436698 | 627.443974 |
| α-Me / ω-N₃ | C₃₁H₅₇N₉O₆ | 651.443181 | 652.450457 |

Calculated ion targets do not identify a measured chain. Record diagnostic end-retaining fragments and the spectrum identifier alongside any proposed assignment. **Compare only** leaves the sample record intact. Selecting one end pair enables transfer of its specific molar-mass correction to Chain Lab.

## Distinguish PVA samples with the same degree of hydrolysis

Open **PVA · Poval 217 (¹³C NMR dyads)**. The published fractions are 0.837 OH–OH, 0.088 mixed and 0.075 OAc–OAc. They give 88.1 mol% OH, η = 0.420 and bulk mean run lengths of 20.02 OH and 2.70 OAc units. Changing to the **equal-hydrolysis Bernoulli reference** preserves the target OH fraction and changes the mixed dyad fraction to 0.2097. The GUI displays bulk target statistics, empirical dyads from a seeded finite ensemble, the distribution of chain compositions and the covalent graph of one simulated chain. The primary NMR measurements and run-length definitions are given in [Budhlall et al., equations 1–8](https://bpb-us-w2.wpmucdn.com/sites.uml.edu/dist/4/136/files/2018/09/Budhlall_2000_Macromol.Symposia_PVOHviaNMRandHPLC-22fn4h7.pdf).

## Install and run

Use Python 3.10 or newer. The reviewed sources for all seven projects share one GitHub branch.

```bash
git clone --single-branch --branch polymer-research-review-v0.2.0 https://github.com/burhanbeycan/burhanbeycan.github.io.git
cd burhanbeycan.github.io/polymer-research
python -m venv .venv
```

Activate using `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` on Windows PowerShell, then:

```bash
python -m pip install ./polymer-research-studio ./polymer-chemistry-lab
polymer-chemistry-lab
```

The local GUI opens at **http://127.0.0.1:8771**. Use `--no-browser` to suppress browser opening. For integrated use, install the other laboratories in the same environment and run `polymer-studio`. RDKit is an external chemistry dependency installed from PyPI; the software is not itself published to PyPI.

## Chemically explicit inputs

Enter the **already polymerized repeat**, with exactly two neutral dummy ports: `[*:1]` and `[*:2]`. Endpoint 2 of each unit joins endpoint 1 of the next. Every unit has the same declared orientation; the constructor does not decide regioselectivity. The general-purpose graph builder uses explicit H/H model caps. Use **PEtOx end groups & ion masses** for chemically assigned α and ω candidates.

| Preset | Polymerized mapped repeat | Repeat contribution |
|---|---|---|
| PE | `[*:1]CC[*:2]` | C₂H₄ |
| PEO | `[*:1]CCO[*:2]` | C₂H₄O |
| PVA | `[*:1]CC(O)[*:2]` | C₂H₄O |
| PMMA | `[*:1]CC(C)(C(=O)OC)[*:2]` | C₅H₈O₂ |
| PEtOx | `[*:1]CCN(C(=O)CC)[*:2]` | C₅H₉NO |
| PMeOx | `[*:1]CCN(C(=O)C)[*:2]` | C₄H₇NO |

PEtOx contains an N-propionyl-substituted ethyleneimine backbone; its cyclic monomer is not used as a polymer repeat. PEO and PVA illustrate why identical elemental formulas do not establish identical polymer chemistry. The H-capped PEO model has one terminal hydroxyl; PVA carries a hydroxyl on each repeat.

The sequence `AAAAAA` constructs six A units. `AAABBB` and `ABABAB` contain the same numbers of A and B but have different directed dyads and connectivity. These are individual explicit chains, not random-copolymer ensembles. Bounds are 64 units and 512 non-port atoms. Stereochemistry is retained where explicitly encoded; leaving it unspecified does not produce an atactic population.

![Finite four-repeat molecular graphs](docs/assets/finite-chain-gallery.png)

The gallery shows H/H-capped finite graphs. Unspecified stereocentres are not an atactic ensemble. Reproduce it with `python scripts/generate_gallery.py`.

## Mass bookkeeping and connected calculations

In the general H/H builder the repeat contribution excludes both terminal atoms.

\[
M_{\mathrm{chain}}=\sum_i M_{\mathrm{repeat},i}+2M_H.
\]

Formula, average molar mass and exact molecular mass are reported separately. A finite molecule's mass is not a sample's measured Mn or Mw. The homopolymer handoff transfers the average repeat mass and the end correction for the selected model. It does not invent a measured MWD, contour increment, persistence length or Kuhn length. Mixed-repeat chains have no automatic one-repeat mapping.

## Verify and inspect

```bash
cd polymer-chemistry-lab
python -m unittest discover -s tests -v
```

The numerical module runs without the workbench. Read [SCIENCE](docs/SCIENCE.md), [REFERENCES](docs/REFERENCES.md) and [VALIDATION](docs/VALIDATION.md) for definitions, verification evidence and limits. Source code and tests distinguish our attachment/sequence protocol and chemistry ledger from RDKit's molecular parsing, graph substitution, mass and depiction algorithms. Developed with AI assistance; MIT license.
