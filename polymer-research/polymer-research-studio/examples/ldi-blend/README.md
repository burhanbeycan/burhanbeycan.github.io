# LDI diurethane and PEtOx dimer contact pilot

In the local Studio, open **Research pilot** and change starting pose, partner and acceptor site. The run record also includes the reproducible mapped SMILES, atom counts, packing seeds and RDKit version. Open [the run results](results.json). The displayed host is a finite, calculated structure with truncated PEG-like and PCL-like arms.

The selected lysine scaffold is ethyl (2S)-2,6-diisocyanatohexanoate transformed at both isocyanates to urethanes. The PEG arm contains two ether oxygens and the PCL arm contains a methyl-capped caprolactone unit. This explicit model uses ethyl-LDI chemistry with precedent in the PEG/PCL literature. The accessible Beycan abstract identifies lysine diisocyanate; the ester form of that specific reagent remains unconfirmed. The molecular graphs of the partners contain two originally polymerized EtOx units with α-Me/ω-OH ends. The second partner has one neutral EI site following formal deacylation. These two fixed protonation states in vacuo do not represent an aqueous pH or the user’s experimental film. [Beycan et al.](https://doi.org/10.1021/acs.iecr.4c00759) · [Hao et al.](https://doi.org/10.1039/C6BM00588H) · [Ethyl-LDI molecule](https://pubchem.ncbi.nlm.nih.gov/compound/71440643)

Each starting geometry was embedded using ETKDGv3. We used four paired seeds for each of two partners and two placement rules. Undirected placements bring the heavy atoms of independently embedded molecules to approximately 3.25 Å separation along a random approach. Contact-initialized placements direct an oxygen of the conserved EtOx carbonyl towards the host urethane N–H and select a low-overlap orientation. Starting with such a contact intentionally biases sampling. Comparing the two preparations tests dependence on initialization.

RDKit supplies MMFF94 analytical gradients with interacting fragments; Chain Lab integrates BAOAB Langevin dynamics. After minimization, each trajectory used 0.5 ps initial thermalization, 3 ps recorded dynamics at 298.15 K, and 0.5 fs steps. Statistics use 120 saved production frames per run. Named acceptors and the two urethane donors are fixed by atom maps. Geometric contacts require N–A ≤ 0.35 nm and N–H–A ≥ 150°, excluding atom pairs separated by three or fewer covalent bonds. Nonperiodic, short trajectories and a small solvent-free cluster limit interpretation to local geometry. Pair counts are normalized per donor or per eligible donor–acceptor pair. Independent packing results should be examined individually.

For EtOx/EtOx, the four undirected trajectories have zero carbonyl contacts under this definition. The contact-seeded trajectories give 0, 0.383, 0.279 and 0.396 mean contacts per host donor. For EtOx/EI the same comparison gives four zero values and 0.404, 0.258, 0.421 and 0.375. The influence of placement precludes assigning either dimer a thermodynamic binding preference. For comparison with films, match validated bulk force-field parameters, sample preparation, experimental composition, independent simulation cells and an experimental observable.

Run from this folder in an environment with the released Studio, Chemistry and Chain Lab packages:

```bash
python reproduce.py --output rerun.json
```
