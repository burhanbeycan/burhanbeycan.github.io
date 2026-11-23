# Polymer Electrolyte Lab

**Nernst–Einstein estimates and reference-frame-explicit collective Onsager transport.**

An independent, local scientific application developed for a polymer-chemistry research workflow. It also connects to [Polymer Research Workbench](https://github.com/burhanbeycan/burhanbeycan.github.io/tree/polymer-research-review-v0.2.0/polymer-research/polymer-research-studio) through a versioned material and result protocol. There is no PolyChemMD dependency.

## Scientific status

The PEO/LiTFSI experiment view calculates dry salt loading from EO:Li or weighed components and compares measured impedance with the Nernst–Einstein baseline from PFG-NMR. The two primary-source presets retain different PEO molar masses, histories and independent transport measurements. The numerical Onsager and ideal-ion demonstration presets remain available alongside them.

Start with **PEO/LiTFSI · Pesko 2017**. At 90 °C and EO:Li = 12.5, measured impedance gives 2 mS/cm. Approximate readings from the paper's tracer plot give σNE ≈ 1.84 mS/cm. The GUI reports the sensitivity to plot reading and records all source locations. **PEO/LiTFSI · Zheng 2018** checks the tracer fraction of a separate 100 kg/mol material using tabulated values. [Worked comparison and input guide](docs/PEO_COMPARISON.md) · [Scientific model](docs/SCIENCE.md) · [Validation](docs/VALIDATION.md).

## Install and run independently

Use Python 3.10 or newer. Check out the review branch containing both source folders:

```bash
git clone --single-branch --branch polymer-research-review-v0.2.0 https://github.com/burhanbeycan/burhanbeycan.github.io.git
cd burhanbeycan.github.io/polymer-research
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

**PEO/LiTFSI · Pesko 2017** örneğini açın. EO:Li bileşimi, 90 °C’de karışım yoğunluğu ve PFG-NMR özdifüzyon katsayılarından hesaplanan iletkenliği, aynı bileşimde EIS ile ölçülen değerle karşılaştırın. **PEO/LiTFSI · weighed film** örneğinde kendi filminizin kalınlığını, elektrot alanını ve yığın direncini girebilirsiniz. Kristallik ölçümünüz yoksa alanı `null` bırakın. İletkenlik oranından serbest iyon yüzdesi elde edilemez.

## Development and attribution

Scientific algorithms are implemented in this repository; NumPy/SciPy provide general numerical primitives. Literature models are attributed in the scientific documentation. The implementation was developed with AI assistance and must be reviewed like other research software. See [CONTRIBUTING](CONTRIBUTING.md). MIT license.

## Review oxygen coordination and charge transport

Choose **PEO shell · original versus scaled force field** to compare published oxygen counts at EO:Li = 10 or 18. Two further modes calculate continuous contact persistence and collective/tracer conductivity from your own independently prepared trajectory records. The contact and transport presets contain synthetic checks. The methods page includes Python adapters for orthorhombic coordinates and unwrapped ion displacements. [Trajectory analysis and worked examples](docs/TRAJECTORY_ANALYSIS.md)
