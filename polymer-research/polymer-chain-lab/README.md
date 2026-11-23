# Polymer Chain Lab

## Automated end-group MD

Compare [twelve executed PEtOx microsolvation trajectories](docs/PEtOx_MD.md), three independent packings per DP and terminus. Source files and a CLI reproduce the pilot from covalent structures. A GROMACS adapter reads condition-matched periodic trajectories supplied with validated end-group topologies. The graphical readout distinguishes finite clusters from bulk-water runs.


**PEtOx dilute-solution viscometry, molar-mass distributions, Gaussian/WLC dimensions and finite-chain RIS conformations.**

An independent, local scientific application developed for a polymer-chemistry research workflow. It also connects to [Polymer Research Studio](../polymer-research-studio) through a versioned material and result protocol. There is no PolyChemMD dependency.

## Scientific status

The PEtOx viscometry benchmark uses 26 published experimental values from seven sample classes. The Gaussian, WLC, RIS and component-viscosity presets use illustrative parameters. Source tables, model assumptions and the per-temperature validation results are in [PEtOx viscometry](docs/PEtOx_VISCOMETRY.md). Other models are documented in [SCIENCE](docs/SCIENCE.md) and [VALIDATION](docs/VALIDATION.md).

## Examine the PEtOx series

1. Select **Measured PEtOx · PBS at 37 °C**. The chart plots measured intrinsic viscosity against SEC–MALS mass-average molar mass. Each sample row reports the prediction made while its entire nominal mass class was withheld.
2. Change the viscometry temperature to 15, 25, 45 or 55 °C. Each fit uses only measurements at that temperature. Set the target Mw inside the mass range for an interpolation.
3. Select **Intrinsic viscosity from a complete MWD** to study a broad distribution. Enter a component-level calibration, its condition and mass limits, and the mass distribution. The result lists each bin's mass fraction and contribution to the mixture's intrinsic viscosity.

## Install and run independently

Use Python 3.10 or newer. From the directory containing both the Studio and the Chain Lab source folders, create a virtual environment:

```bash
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell. Then:

```bash
python -m pip install ./polymer-research-studio ./polymer-chain-lab
polymer-chain-lab
```

The GUI opens at `http://127.0.0.1:8766`. To suppress automatic browser opening, use `polymer-chain-lab --no-browser`. The server binds to loopback only. No API key or commercial simulation package is required. This release is installed from the repositories; no PyPI publication is implied.

## Use with other laboratories

Install any other laboratory into the same environment, then run `polymer-workbench`. The workbench discovers installed laboratories. Separate virtual environments provide independent installations.

The shared material record includes architecture, tacticity, sequence, end groups, conditions and provenance. These descriptors are **not automatically converted into unknown model parameters**. Parameter transfer is explicit and dimensionally defined. Export a result JSON to connect separately installed applications, or use the integrated Workflow view.

## Chemical mass mapping

Supply the repeat contribution and a common net end-group mass correction:
`DP=(M−Mend)/Mrepeat`, then `L=DP×contour increment`. End groups are not repeating units. The app reports each component's effective DP and its number average while keeping total chain mass in Mn, Mw, Mz and scattering weights. Zero end correction explicitly neglects it. A capped monomer's mass is not a repeat-unit mass. Terminal geometry, composition-resolved copolymer mapping and unsupported architectures require other models; declared cyclic/branched/network architectures are rejected.

RIS conformations are equal-mass skeletal sites. They are not assigned atom identities, substituents, stereochemistry or force-field topology; use Polymer Chemistry Lab for an explicit atom-bond oligomer graph.

## GUI capabilities

- Parameter forms with units and model-specific scope.
- Publication-friendly SVG plots and CSV data export.
- Complete calculation JSON with inputs, references and checksum.
- Material record and calculation lineage.
- Model equations, assumptions and sources alongside results.
- Responsive, local academic interface; no external CDNs.

## Verify

```bash
cd polymer-chain-lab
python -m unittest discover -s tests -v
```

The numerical module in `src/polymer_chain_lab/model.py` is independent of the GUI. Read the tests to see the specific verified properties. A test count is not a claim of experimental validation.

## Türkçe kullanım özeti

Bu uygulama tek başına çalışır; ortak çalışma ortamında diğer laboratuvarlarla birlikte de kullanılabilir. Model, birimler, parametre kaynakları ve geçerlilik sınırları sonuçların yanında gösterilir. PEtOx viskozimetri ekranı Gubarev ve arkadaşlarının deneysel verilerini kullanır. Her kütle sınıfı sırayla dışarıda bırakılır. Karışım ekranında ise kullanıcı kendi polimerine, çözücüsüne ve sıcaklığına ait bileşen katsayılarını belirtir.

## Development and attribution

Scientific algorithms are implemented in this repository; NumPy/SciPy provide general numerical primitives. Literature models are attributed in the scientific documentation. The implementation was developed with AI assistance and must be reviewed like other research software. See [CONTRIBUTING](CONTRIBUTING.md). MIT license.
