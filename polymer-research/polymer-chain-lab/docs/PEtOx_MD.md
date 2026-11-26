# PEtOx end chemistry in explicit water

Open **Chain Lab → Executed MD · PEtOx termini with 12 waters**. Select DP 2 or 4. The graph plots actual neighbor counts from independently generated starting configurations. The table gives trajectory-mean values for each packing and separate arithmetic and RMS radii of gyration. The program weights replicas equally when calculating the displayed grand mean.

## What ran

Each neutral chain has a covalent α-N-methyl group and a carbon-bound ω-OH or C–N=[N⁺]=[N⁻] group. The terminal comparison counts water oxygens within 0.35 nm of the hydroxyl oxygen or the **distal** azide nitrogen. Backbone coordination averages over every pendant amide carbonyl O. One water can participate in the environments of several different sites. This geometric occupancy has no angular hydrogen-bond criterion.

Each independent system contains one oligomer and 12 explicit waters in vacuum. RDKit ETKDGv3 embeds each polymer; a seeded procedure places and orients the waters around its heavy atoms. All 12 systems passed the MMFF94 atom-typing and minimization checks, with interfragment interactions enabled. The production trajectories use an unconstrained BAOAB Langevin integrator at 298.15 K, 0.5 fs timestep and friction 5 ps⁻¹. The initial interval lasts 1 ps and the observed interval 2 ps. Forty equally spaced frames are saved per run. Each condition has three independent starting conformers, solvent packings and initial velocity draws. There is no periodicity, water chemical potential, bulk density, applied pressure or confining potential.

| DP | ω site | Mean Rg (nm) | Water O per amide O | Water O per ω site | Range of the three replica means (ω) |
|---:|---|---:|---:|---:|---|
| 2 | OH oxygen | 0.326 | 0.962 | 1.975 | 1.625–2.225 |
| 2 | Distal azide nitrogen | 0.352 | 1.183 | 0.275 | 0.175–0.450 |
| 4 | OH oxygen | 0.438 | 1.071 | 1.442 | 1.000–1.975 |
| 4 | Distal azide nitrogen | 0.463 | 0.929 | 0.058 | 0–0.100 |

The oxygen and nitrogen at the two ω sites probe different chemical environments. The early-time occupancies depend on initialization, force field and cutoff. Three short replicates do not establish equilibrium ordering, bulk hydration, chain-length scaling or the accuracy of the group charges. Missing water around a distal N does not imply zero hydration at the other azide nitrogens. The saved positions are model results; the panel has no experimental or hydrodynamic radius label.

## Reproduce the 12 trajectories

From the suite directory, install the Studio runtime and the chain package with its MD extra. Then run:

```sh
python -m pip install ./polymer-research-studio './polymer-chain-lab[md]'
polymer-petox-md pilot petox-md-rerun
```

Use a new directory. The command writes the 480 retained frame summaries, one final XYZ per trajectory and a run manifest including software versions, seeds, energies and formulas. Paste the CSV into **Analyzed trajectory CSV** to view the rerun. The release contains the [frame summaries](../src/polymer_chain_lab/data/petox_cluster_md.csv), [manifest](../examples/md-pilot/pilot.json), and the seed 109 final structures at DP4 for [OH](../examples/md-pilot/dp4_OH_final.xyz) and [N₃](../examples/md-pilot/dp4_N3_final.xyz). [RDKit ForceField](https://www.rdkit.org/docs/source/rdkit.ForceField.rdForceField.html)

Forces are negative analytical MMFF94 gradients in kcal mol⁻¹ Å⁻¹. Multiplying force divided by atom mass in Da by 0.0004184 gives acceleration in Å/fs². A B/A/O/A/B step performs two half kicks and half drifts around an exact Ornstein–Uhlenbeck velocity increment. Energies, coordinate-frame independence from the output stride and an analytical harmonic oscillator are tested. [Langevin integration, Leimkuhler & Matthews](https://doi.org/10.1063/1.4802990)

## Compare to a periodic system

Galata & Kröger Table 1 supplies a reference for linear DP30 PEtOx in TIP3P water. Reported RMS Rg values at DP30 are 0.76 nm (AMBER), 1.55 nm (DREIDING) and 1.82 nm (UFF); the reported water counts and run durations differ across these setups. Use these as protocol-specific simulation comparison targets. A physical replication needs the matching atom types, charges, end chemistry, architecture, solvent and sampling protocol. The cluster runs do not reproduce this reference system. A DOSY diffusion coefficient, hydrodynamic radius and mass-weighted Rg are different observables. [Primary report, Table 1](https://doi.org/10.1021/acs.macromol.4c00174)

The GROMACS adapter executes a prepared production TPR, then converts the trajectory to a whole-molecule GRO and calculates Rg and site-specific water contacts. The TPR should contain an equilibrated, orthorhombic system with validated PEtOx and end-group parameters, and store the whole system in XTC or TRR. Use the protocol's actual start/retained simulation times; `--discard-ps` is an absolute trajectory time. The supplied zero-based mapping must identify one complete chain and all water oxygens with the same force field and temperature across compared replicas. For a previously completed run you can also analyze a multiframe whole-molecule GRO directly.

```sh
polymer-petox-md gromacs production.tpr sites.json periodic-run --discard-ps 10000
polymer-petox-md gro whole.gro sites.json analyzed-run
```

The conversion uses `gmx trjconv -pbc whole` and the fixed-width default GRO coordinate precision of 0.001 nm. The reader checks known residue/atom names, connectedness of selected polymer atoms, bond lengths and frame ordering. Distances use minimum images for water contacts and whole coordinates for the chain. The result records SHA-256 input hashes. Local tests cover the GRO parser and analysis using analytical coordinates. The external GROMACS executable was unavailable for an end-to-end periodic production rerun in the release environment. [GROMACS trjconv options](https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html)

### Site mapping

Use an explicit index and mass list from your topology. `polymer` contains every atom in one chain; `polymer_masses_Da` follows the same order. `polymer_bonds` connects this chain and `carbonyl` holds one amide oxygen per repeat. `omega` holds a single end-group oxygen or distal azide nitrogen. `water_oxygen` enumerates all solvent oxygen atoms in that system. Supply the allowed `RES:ATOM` names from the GRO in `carbonyl_names`, `omega_names` and `water_oxygen_names`. The required scalar fields are `sample` (matched series label), `source` (topology reference), `force_field`, `water_model`, `temperature_K`, `dp`, `end` (OH or N3) and `replica` (a distinct packing/seed identifier). Compare identical water counts and thermodynamic conditions within each series.

For a no-download format and parser example, [test_gro_mapping_names_timing_pbc_and_unwrapping](../tests/test_md.py) generates a small analytical GRO fixture. Replace both its atom assignment and conditions before processing scientific data.
