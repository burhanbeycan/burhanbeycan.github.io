# Polymer Research Studio

**Chemical structure, polymer populations and physical response in one connected research environment.**

Polymer Research Studio is developed for Burhan Beycan's polymer-chemistry research workflow. Six independent projects share a local academic interface and a versioned sample record. The GitHub repository is `polymer-research-studio`; the Python distribution remains `polymer-research-workbench` for compatibility. Start the user interface with `polymer-studio`. PolyChemMD is a separate project and is not a dependency.

## Independent laboratories

| Project | Implemented work | Independent command | Local port |
|---|---|---|---|
| [Polymer Chemistry Lab](https://github.com/burhanbeycan/polymer-chemistry-lab) | Valence-checked finite oligomers from mapped repeat units, explicit A/B sequences, H caps, molecular formula and molar mass, atom/bond diagrams | `polymer-chemistry-lab` | 8771 |
| [Polymer Synthesis Lab](https://github.com/burhanbeycan/polymer-synthesis-lab) | AA+BB Carothers balances, ideal initiated-chain DPn/Mn, terminal Mayo–Lewis composition, stationary dyads and run lengths | `polymer-synthesis-lab` | 8770 |
| [Polymer Chain Lab](https://github.com/burhanbeycan/polymer-chain-lab) | Number/mass MWD, Mn/Mw/Mz/Đ, end-group-aware DP mapping, Gaussian/WLC dimensions, Gaussian scattering and finite-chain RIS | `polymer-chain-lab` | 8766 |
| [Polymer Transport Lab](https://github.com/burhanbeycan/polymer-transport-lab) | Slab diffusion, finite-volume verification, fixed-geometry D fitting and Henry-law permeation | `polymer-transport-lab` | 8767 |
| [Polymer Network Lab](https://github.com/burhanbeycan/polymer-network-lab) | Neutral-network equilibrium swelling, consistent mixing/elastic free energy, homogeneous stability and equilibrated-film release | `polymer-network-lab` | 8768 |
| [Polymer Electrolyte Lab](https://github.com/burhanbeycan/polymer-electrolyte-lab) | Binary Nernst–Einstein baseline, supplied collective Onsager coefficients, PSD checks and polymer-fixed transference | `polymer-electrolyte-lab` | 8769 |

Version **0.2.0** distinguishes a finite chemical graph from a polymer sample and from a statistical chain model. Molecular Chemistry performs explicit structural calculations on the supplied templates. Physics-model demonstration coefficients are illustrative and do not establish predictive accuracy for a named polymer. Each laboratory documents equations, sources, tests and model boundaries.

## Install all six

Use Python 3.10 or newer and clone the repositories into sibling folders:

```bash
git clone https://github.com/burhanbeycan/polymer-research-studio.git
git clone https://github.com/burhanbeycan/polymer-chemistry-lab.git
git clone https://github.com/burhanbeycan/polymer-synthesis-lab.git
git clone https://github.com/burhanbeycan/polymer-chain-lab.git
git clone https://github.com/burhanbeycan/polymer-transport-lab.git
git clone https://github.com/burhanbeycan/polymer-network-lab.git
git clone https://github.com/burhanbeycan/polymer-electrolyte-lab.git
python -m venv .venv
```

Activate with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell, then run:

```bash
python -m pip install ./polymer-research-studio ./polymer-chemistry-lab ./polymer-synthesis-lab ./polymer-chain-lab ./polymer-transport-lab ./polymer-network-lab ./polymer-electrolyte-lab
polymer-workbench
```

The Studio opens at **http://127.0.0.1:8765**. `polymer-studio` is an equivalent launcher. Use `--no-browser --port 8775` for a different local port. Installation is from source; these instructions do not imply a PyPI release. A source archive uses the same install command without cloning.

To use one laboratory independently, install this shared runtime and the chosen laboratory, then use its own command. Its launcher exposes only that laboratory. Separate environments provide fully isolated installations. Scientific kernels do not depend on the GUI or on another laboratory.

## A chemist's workflow

1. **Define the chemistry.** Supply polymerized repeat units with explicit head/tail attachment points. Chemistry Lab connects actual atoms and covalent bonds, validates valence and retains specified stereochemistry. H capping is an explicit model of end groups.
2. **Resolve the population.** A finite oligomer structure is not a molar-mass distribution. Chain Lab accepts a distribution with a stated number/mass basis and distinguishes Mn, Mw, Mz and scattering averages. End-group mass is explicit in the DP mapping.
3. **Choose the physical model.** Select conformation, synthesis, network, transport or electrolyte assumptions for the intended question. Structure metadata does not silently supply unknown Kuhn lengths, interaction parameters or diffusion coefficients.
4. **Record evidence.** Preserve sample identity, mapped SMILES, architecture, branch topology, tacticity evidence, sequence, functional groups, route, end groups and conditions with each calculation.
5. **Transfer defined quantities.** The Workflow view offers explicit input mappings and preserves calculation lineage. Independently installed apps exchange the same result JSON.

The molecular view has repeat-unit colours, carbon labels, atom indices, zoom and SVG export. These are **2D connectivity drawings**, not relaxed conformations or molecular-dynamics trajectories. The RIS visualization uses equal-mass statistical sites; those sites are not assigned chemical atom identities.

## Scientific interoperability

A homopolymer Chemistry → Chain handoff can carry the calculated average repeat-unit molar mass and the explicitly chosen cap mass. It does not determine contour increments or stiffness. Copolymer chemistry requires explicit composition/sequence treatment rather than an unqualified homopolymer mapping.

Network → Transport transfers the **full equilibrium swollen thickness** and the independently supplied effective solute diffusivity. The receiver solves diffusion in a fixed, already swollen geometry. Electrolyte tracer D is never silently treated as a Fickian penetrant diffusivity.

Each result contains inputs, units, material metadata, model sources, UTC timestamp, run ID and integrity checksum. The checksum detects accidental edits; it is not a signature or scientific validation. See [EXCHANGE](docs/EXCHANGE.md).

## Reproducibility

```bash
python -m pip install './polymer-research-studio[test]'
cd polymer-research-studio
python -m unittest discover -s tests -v
```

Run each laboratory's tests from its own repository. Full-suite integration tests skip when any laboratory is absent. [VALIDATION](docs/VALIDATION.md) states the checks performed and their limits. GitHub Actions configuration is supplied; remote CI success must be checked on the actual commit.

The interface and all calculations run locally with no external CDN, API key or commercial simulator. RDKit provides chemical graph handling and depiction coordinates; NumPy/SciPy supply numerical primitives. Original code joins these tools through explicit polymer sample semantics and independently implemented scientific models. Literature models and third-party tools retain their attribution. MIT licence; see [CONTRIBUTING](CONTRIBUTING.md).

## Türkçe

Studio, polimer kimyasını tekrar birimi, uç grup, kopolimer dizisi ve gerçek atom–bağ bağlantılarıyla ele alır; numunenin molar kütle dağılımını ve fiziksel model varsayımlarını ayrı tutar. Altı uygulama ayrı ayrı veya birlikte çalışır. Yapı çizimi termodinamik denge geometrisi sayılmaz; fiziksel tahminler için kullanılan katsayıların kaynağı ve geçerlilik aralığı açık tutulur.
