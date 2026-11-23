# LDI–PEG/PCL polyurethane chemistry and molecular contacts

## What the original experiments establish

Beycan et al. prepared WBPU from PEG and PCL diols with lysine diisocyanate, then emulsified the material with aqueous L-lysine. The alternative extenders included 1,3-propanediol (PDO), 1,4-butanediol (BDO) and glycerol. Supporting Information Figure S1 presents reaction schemes. The study describes electrospun specimens. In a separate study, Hao et al. explicitly used L-lysine ethyl ester diisocyanate, PDO and aqueous L-lysine with PEG/PCL. The reported ranges of tensile strength, 12–18 MPa, and elongation at break, 1400–2500%, belong to that study's series. Neither range is a specimen-resolved ATR-FTIR band area. For this release, a matched set of LDI dry-film FTIR areas and replicated bulk trajectories could not be verified. [Beycan et al.](https://doi.org/10.1021/acs.iecr.4c00759) · [Hao et al.](https://doi.org/10.1039/C6BM00588H)

## Compare an illustrative feed

In **LDI / PEG / PCL · feed and branching**, quantities are mmol of dry molecules. At nominal diol functionality, 0.05 mmol PEG and 1 mmol PCL supply 0.10 and 2.00 mmol OH. Adding 0.50 mmol BDO gives 3.10 mmol OH. A 1.65 mmol LDI feed supplies 3.30 mmol NCO and yields a nominal NCO/OH index of 106.45%. With identical molecules of glycerol the result is 3.60 mmol OH and an index of 91.67%. To compare at *fixed* nominal index and OH count, replace 0.50 mmol BDO by 1/3 mmol glycerol. Check the hard-segment mass fraction separately, because equifunctional feeds need not have the same component mass.

The preset is a worked stoichiometric example. Verify the LDI ester form for the actual batch before building an atomistic topology. Protected termini, off-stoichiometry and incomplete conversion matter. Presence of glycerol permits junctions; formation of a covalent network requires additional evidence. In the dispersion stage the upper bound for lysine–isocyanate additions, in mmol, is min(measured pre-dispersion residual NCO, 2 × added L-lysine). A blank assay leaves the bound unspecified. Amine protonation, accessibility and competition with water can lower incorporation.

## Count urethane N–H···O contacts

Use the **NH contact geometry** mode for a short trajectory, or the streaming GRO input below for larger trajectories. Atom IDs begin at zero. All coordinates and box lengths are in nm, time in ps. Frames must be equally spaced. Supply the complete fixed molecular topology as an `atom_chain` array, a list of bonded urethane `donor_NH` atom pairs, and one labeled oxygen entry per `acceptors` site. Use the same atom order in every frame. Labels are `PEG_ether`, `PCL_ester`, `LDI_ester`, `urethane_carbonyl` and `urea_carbonyl`. For the ester groups, include the carbonyl oxygen only. Count urea carbonyls only where the modeled topology actually contains a urea. Assign each *molecule* a unique chain ID. With covalently percolated networks, a molecule may contain the entire component. This interchain analysis cannot recover distinct covalent strands without a new, explicitly defined exclusion policy.

At every donor, the algorithm considers acceptors on other chains, evaluates minimum-image N···O distances and minimum-image H→N and H→O vectors, then filters using ∠N–H···O. The default cutoffs are distance ≤ 0.35 nm, angle ≥ 150°. These explicit geometric choices should be sensitivity-tested. Orthorhombic, frame-dependent boxes are supported. For group k, the eligible inventory counts all cross-chain donor–acceptor pairs of its class. The enrichment is its observed share of contact observations divided by its share of the eligible inventory. Enrichment 1 denotes proportionality to that inventory. This denominator corrects stoichiometric abundance, not solvent-accessible area or phase segregation. A donor may have several simultaneous contacts.

The constructed four-frame example has three ether and two PCL-ester triplet observations. There are two eligible ether pairs, one PCL-ester pair and one urethane-carbonyl pair. Resulting enrichments are 1.20, 1.60 and 0.00, respectively. The observation count and enrichment rankings therefore differ. The seventh atom in the example is a same-chain oxygen deliberately excluded from the denominator. Replace these coordinates and their source when interpreting polymer data.

The persistence table reports the fraction of saved-frame pair observations also present one frame later and the span between first and last consecutive contact snapshots. The latter is zero for a one-frame event. Episodes reaching either trajectory boundary are omitted from the complete-episode mean. Listed maximum spans may be censored. Reformation between saved frames and continuous-time lifetimes cannot be resolved. For real data, examine equilibrium, analyze independently packed systems separately and use between-system uncertainty. Hydration, box size, force-field and saved-frame interval must accompany the result.

### Importing a GROMACS trajectory

Export a multi-frame GRO with whole-system atom order using GROMACS (e.g. `gmx trjconv -f production.xtc -s production.tpr -o trajectory.gro`, selecting the whole system). Use an orthorhombic simulation box and save the typed atom map from the same topology as `sites.json`. A minimal atom map is:

```json
{"atom_chain":["A","A","B"],"donor_NH":[[0,1]],"acceptors":[[2,"PCL_ester"]]}
```

For the real system replace the atom map and pass the trajectory file to:

```bash
python -m polymer_network_lab.trajectory sites.json trajectory.gro > contacts.json
```

Add `--spacing-ps 1` only for regularly spaced GRO frames lacking time headers and use the actual output interval. `--cutoff-nm` and `--angle-deg` change the geometric thresholds. GRO coordinates have finite output precision. Use the Python `contact_analysis(meta, frames, ...)` function with original nm coordinates for higher precision. Input is streamed frame by frame; the output retains per-frame counts. The analyzer performs no parameter assignment or MD integration.

## Test the extra value of one MD descriptor

In **recipe vs MD / FTIR**, supply a single matched protocol and one row per specimen. The target is a *measured and consistently fitted* bonded carbonyl area divided by the total area under the chosen carbonyl model. Fix preprocessing, component definitions, fit window, history and specimen type. Record peak overlap from PCL/LDI ester and urethane/urea, and check crystallinity independently. Do not import a mat tensile measurement as a dry-film spectrum. Each source and trajectory has a specimen/formulation assignment, and rows from one formulation always remain together. Treat correlated technical repeats from one physical preparation as one specimen. In a real study predefine donor and acceptor definitions and do not select descriptors using held-out spectra.

The baseline inputs are PEG mass fraction among soft segments and nominal dry-feed NCO/OH. The augmented model adds the interchain NH→ether enrichment. A ridge coefficient of 1 is fixed before splitting; the intercept is unpenalized and all scaling uses the training formulations only. The split withholds each entire formulation after averaging its specimen rows. Predictions outside the physical fraction range are reported as such. Fit the same training folds and compare held-out MAE in percentage points. This is a within-protocol test; prediction across studies needs an external study with matching definitions. In eight **constructed** formulations, MAE is 5.59 and 1.45 percentage points. This numerical example checks the analysis and does not support a claim of predictive improvement in measured films.

## Run independent GROMACS replicas locally

The installed `polymer_network_lab.bulk` module can run a supplied, parameterized system. A topology includes the actual ester and ends, chosen urethane/urea linkages, atom types, charges, bonded and nonbonded terms and matching molecule counts. Supply MDP files appropriate for this representation and state. Validate the force field against the chemistry and observables of interest. Explicit urethane hydrogens and chemical atom mapping are required for the contact analysis above. A Martini topology prepared with PolySMart can run with its own validated CG inputs; atomistic N–H contact criteria are inapplicable to its beads.

Start from the [two-replica configuration example](../examples/ldi-replica-plan.json). Its paths and 6 nm example box are placeholders. Choose a box size consistent with the actual system, and put each replica's complete inputs—including local force-field directories, relative include files and stage MDP files—in the declared dedicated input folder. Use independent prepacked coordinates or independent Polyply packing calls. Record velocity generation and any seeds in each replica's MDP/evidence. A different folder name alone provides no evidence of independence. Specify an optional per-replica `index` file if using custom groups. For a restrained stage an optional `restraint_reference` is the coordinate file used by `grompp -r`. A replica may override the top-level stages.

```bash
python -m polymer_network_lab.bulk replicas.json results/run-01 --dry-run
python -m polymer_network_lab.bulk replicas.json results/run-01
```

With prebuilt coordinates, only GROMACS is needed. The optional packing route runs `polyply gen_coords -p system.top -o packed.gro -name REPLICA -box X Y Z`. Velocity seeds and the stochastic packing algorithm can vary across runs. The exact generated or imported starting coordinates are retained for replay. Each ordered stage runs `gmx grompp` then `gmx mdrun`. The previous final GRO starts the next stage. Setting `continue_from_previous_checkpoint` also supplies the previous checkpoint to `grompp -t`; configure matching continuation and velocity-generation behavior in that stage's MDP. Use `--gmx` or `--polyply` for alternate binaries.

A dry run validates input paths and prints the planned commands, working directories and SHA-256 hashes. Execution requires a fresh output directory and preserves individual inputs, generated coordinates, processed topology/MDP, TPR, trajectory outputs selected by the MDP, and step logs. The manifest records execution and file hashes. For exact reproduction archive the GROMACS installation and external force-field include paths if the input topologies refer to them. When a command fails, completed files and the error are retained; use a fresh run directory for a changed configuration. Tool completion establishes that preprocessing and integration finished. Equilibration, finite-size effects and replicate agreement still need evaluation from the trajectories. The adapter was tested with a command-line test double; GROMACS and Polyply were unavailable for executing a WBPU bulk cell during this release. [Polyply coordinate generation](https://github.com/marrink-lab/polyply_1.0/blob/master/docs/tutorials/martini-polymers.md) · [GROMACS preprocessing](https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html) · [GROMACS dynamics](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html)
