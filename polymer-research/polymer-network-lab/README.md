# Polymer Network Lab

**Polyurethane functional-group balances, film extraction, DMA and condition-matched swelling analysis.**

An independent, local scientific application developed for a polymer-chemistry research workflow. It also connects to [Polymer Research Workbench](https://github.com/burhanbeycan/burhanbeycan.github.io/tree/polymer-research-review-v0.2.0/polymer-research/polymer-research-studio) through a versioned material and result protocol. There is no PolyChemMD dependency.

## Scientific status

Choose a WBPU preset in the GUI to review the component ledger, an extraction/DMA example, or a measured formulation series. The legacy equilibrium mode calculates free isotropic swelling of an ideal neutral network with supplied parameters.

| WBPU analysis | Input | Readout |
|---|---|---|
| Functional groups | Component masses, lot-specific OH number or NCO content, equivalent weights or nominal molar masses and functionality | OH, NCO, amine and independently reacted-water equivalents by component and stage, NCO/OH, overall reactive index |
| Extraction and DMA | Initial, blotted and redried masses, solvent, exposure and temperature, optional dry-film tensile DMA | Operational gel/sol, net swelling on initial dry mass, solvent uptake on retained gel and E′ traces by frequency |
| Withheld formulation | Matched recipe and swelling observations, optional independent DMA | Fold predictions and MAE in percentage points for a training mean, recipe fraction, or independent modulus |

The experimental preset contains three jatropha WBPU films measured in toluene. The recipe example illustrates equal masses of PEG 20k and PCL 1k diols. Each gram supplies 0.10 and 2.00 mmol of OH, respectively. Keep sample-level results for your PEG/PCL films and electrospun mats linked to their original records. [WBPU methods and worked example](docs/WBPU.md)

## Install and run independently

Use Python 3.10 or newer. Check out the review branch containing both source folders:

```bash
git clone --single-branch --branch polymer-research-review-v0.2.0 https://github.com/burhanbeycan/burhanbeycan.github.io.git
cd burhanbeycan.github.io/polymer-research
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell. Then:

```bash
python -m pip install ./polymer-research-studio ./polymer-network-lab
polymer-network-lab
```

The GUI opens at `http://127.0.0.1:8768`. To suppress automatic browser opening, use `polymer-network-lab --no-browser`. The server binds to loopback only. No API key or commercial simulation package is required. This release is installed from the repositories; no PyPI publication is implied.

## Use with other laboratories

Install any other laboratory into the same environment, then run `polymer-workbench`. The workbench discovers installed laboratories. Separate virtual environments provide independent installations.

The shared material record includes architecture, tacticity, sequence, end groups, conditions and provenance. These descriptors are **not automatically converted into unknown model parameters**. Parameter transfer is explicit and dimensionally defined. Export a result JSON to connect separately installed applications, or use the integrated Workflow view.

## GUI capabilities

- Parameter forms with units and model-specific scope.
- Publication-friendly SVG plots and CSV data export.
- Complete calculation JSON with inputs, references and checksum.
- Material record and calculation lineage.
- Model equations, assumptions and sources alongside results.
- Responsive, local academic interface; no external CDNs.

## Verify

```bash
cd polymer-network-lab
python -m unittest discover -s tests -v
```

The numerical module in `src/polymer_network_lab/model.py` is independent of the GUI. Read the tests to see the specific verified properties. A test count is not a claim of experimental validation.

## Türkçe kullanım özeti

Bu uygulama tek başına çalışır; ortak çalışma ortamında diğer laboratuvarlarla birlikte de kullanılabilir. Model, birimler, parametre kaynakları ve geçerlilik sınırları sonuçların yanında gösterilir. Varsayılan değerler gösterim içindir. Gerçek bir polimer için tahmin iddiası, uygun parametreler ve bağımsız deneysel doğrulama gerektirir.

## Development and attribution

Scientific algorithms are implemented in this repository; NumPy/SciPy provide general numerical primitives. Literature models are attributed in the scientific documentation. The implementation was developed with AI assistance and must be reviewed like other research software. See [CONTRIBUTING](CONTRIBUTING.md). MIT license.

## LDI polyurethane, molecular contacts and dry-film comparison

Open the three **LDI / PEG / PCL** presets. Inspect the effect of replacing BDO by glycerol, measure urethane N–H contact preference for ether and ester oxygens in typed coordinates, then compare formulation-held-out FTIR predictions with and without the contact descriptor. For GROMACS output the module reads multi-frame GRO files directly. [Methods, worked examples and input guide](docs/LDI.md)

A local replica adapter runs prepared systems through optional Polyply coordinate generation, GROMACS preprocessing and MD. The module records source files, commands and outputs per packing. [Two-replica configuration](examples/ldi-replica-plan.json) · [Operation and input requirements](docs/LDI.md#run-independent-gromacs-replicas-locally)
