# Polymer Electrolyte Lab

**Nernst–Einstein estimates and reference-frame-explicit collective Onsager transport.**

An independent, local scientific application developed for a polymer-chemistry research workflow. It also connects to [Polymer Research Workbench](https://github.com/burhanbeycan/polymer-research-studio) through a versioned material and result protocol. There is no PolyChemMD dependency.

## Scientific status

Version **0.2.0** retains the independently verified scientific kernels and connects to the expanded Polymer Research Studio. It implements explicitly scoped models, numerical checks, synthetic reference scenarios and reproducible exports. Initial example coefficients are **illustrative, not measured material constants**. Numerical verification does not establish accuracy for a particular polymer. See [SCIENCE](docs/SCIENCE.md), [REFERENCES](docs/REFERENCES.md) and [VALIDATION](docs/VALIDATION.md).

## Install and run independently

Use Python 3.10 or newer. Clone the common runtime and this application as sibling folders:

```bash
git clone https://github.com/burhanbeycan/polymer-research-studio.git
git clone https://github.com/burhanbeycan/polymer-electrolyte-lab.git
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell. Then:

```bash
python -m pip install ./polymer-research-studio ./polymer-electrolyte-lab
polymer-electrolyte-lab
```

The GUI opens at `http://127.0.0.1:8769`. To suppress automatic browser opening, use `polymer-electrolyte-lab --no-browser`. The server binds to loopback only. No API key or commercial simulation package is required. This release is installed from the repositories; no PyPI publication is implied.

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
cd polymer-electrolyte-lab
python -m unittest discover -s tests -v
```

The numerical module in `src/polymer_electrolyte_lab/model.py` is independent of the GUI. Read the tests to see the specific verified properties. A test count is not a claim of experimental validation.

## Türkçe kullanım özeti

Bu uygulama tek başına çalışır; ortak çalışma ortamında diğer laboratuvarlarla birlikte de kullanılabilir. Model, birimler, parametre kaynakları ve geçerlilik sınırları sonuçların yanında gösterilir. Varsayılan değerler gösterim içindir. Gerçek bir polimer için tahmin iddiası, uygun parametreler ve bağımsız deneysel doğrulama gerektirir.

## Development and attribution

Scientific algorithms are implemented in this repository; NumPy/SciPy provide general numerical primitives. Literature models are attributed in the scientific documentation. The implementation was developed with AI assistance and must be reviewed like other research software. See [CONTRIBUTING](CONTRIBUTING.md). MIT license.
