# Polymer Transport Lab

**Slab sorption and release, independently discretized diffusion, parameter fitting and Henry permeation.**

An independent, local scientific application developed for a polymer-chemistry research workflow. It also connects to [Polymer Research Workbench](https://github.com/burhanbeycan/burhanbeycan.github.io/tree/polymer-research-review-v0.2.0/polymer-research/polymer-research-studio) through a versioned material and result protocol. There is no PolyChemMD dependency.

## Scientific status

Version **0.4.0** retains the independently verified scientific kernels and connects to the expanded Polymer Research Studio. It implements explicitly scoped models, numerical checks, synthetic reference scenarios and reproducible exports. Initial example coefficients are **illustrative, not measured material constants**. Numerical verification does not establish accuracy for a particular polymer. See [SCIENCE](docs/SCIENCE.md), [REFERENCES](docs/REFERENCES.md) and [VALIDATION](docs/VALIDATION.md).

## Install and run independently

Use Python 3.10 or newer. Check out the review branch containing both source folders:

```bash
git clone --single-branch --branch polymer-research-review-v0.2.0 https://github.com/burhanbeycan/burhanbeycan.github.io.git
cd burhanbeycan.github.io/polymer-research
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell. Then:

```bash
python -m pip install ./polymer-research-studio ./polymer-transport-lab
polymer-transport-lab
```

The GUI opens at `http://127.0.0.1:8767`. To suppress automatic browser opening, use `polymer-transport-lab --no-browser`. The server binds to loopback only. No API key or commercial simulation package is required. This release is installed from the repositories; no PyPI publication is implied.

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
cd polymer-transport-lab
python -m unittest discover -s tests -v
```

The numerical module in `src/polymer_transport_lab/model.py` is independent of the GUI. Read the tests to see the specific verified properties. A test count is not a claim of experimental validation.

## Türkçe kullanım özeti

Bu uygulama tek başına çalışır; ortak çalışma ortamında diğer laboratuvarlarla birlikte de kullanılabilir. Model, birimler, parametre kaynakları ve geçerlilik sınırları sonuçların yanında gösterilir. Varsayılan değerler gösterim içindir. Gerçek bir polimer için tahmin iddiası, uygun parametreler ve bağımsız deneysel doğrulama gerektirir.

## Development and attribution

Scientific algorithms are implemented in this repository; NumPy/SciPy provide general numerical primitives. Literature models are attributed in the scientific documentation. The implementation was developed with AI assistance and must be reviewed like other research software. See [CONTRIBUTING](CONTRIBUTING.md). MIT license.

## PVA films under a humidity step

Select **PVA film · independent isotherm + humidity step** to follow the mean water content of an equilibrated film after changing vapour activity. Enter equilibrium and kinetic measurements on a dry-polymer mass basis and select a supported film (one exposed face) or suspended film (two equally exposed faces). Compare the analytical constant-D baseline with a conservative numerical model in which D varies with local water loading.

Select **PVA film · compare models at a held-out thickness** to reserve every time series at one nominal thickness before fitting either model. Each curve carries its thickness, exposed faces, grade/history identifier, temperature and branch. Separate plots and RMSEs show the effect of model choice on the reserved thickness. The bundled 50/90 µm example uses simulated data.

Select **PVA · inspect published equilibrium and apparent D** for the 25 °C HDI-crosslinked PVA results of El Shafee and Naguib. The table distinguishes the measured-curve apparent coefficients from equilibrium loadings calculated using the authors’ GAB fit. Numerical methods, CSV templates and the interpretation of each result are in [PVA water transport](docs/PVA_WATER.md).

## Sequence-aware PVA vapour sorption

Open **PVA · evidence for sequence and uptake** to see which primary data actually match. **PVA · withhold complete grades (synthetic)** runs three nested models on identical observation sets. The form accepts separate specimen and adsorption CSVs, with selection of grade, publication or joint holdout. **PVA · oxygen coordination across a periodic boundary** accepts selected oxygen coordinates from an external MD trajectory and calculates water contacts per site. See [methods, CSV units and benchmark](docs/PVA_SEQUENCE.md).
