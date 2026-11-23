# Which oxygen occupies the Li⁺ shell in PEO/LiTFSI?

## Published force-field benchmark

Choose **PEO shell · original versus scaled force field** and select the EO:Li ratio. At 300 K the authors integrated Li–O radial distribution functions to 0.27 nm. The values below were reported in the discussion of Figure 3. Each sulfonyl oxygen that satisfies the cutoff contributes to the count, including two oxygens belonging to one TFSI anion. [Primary article and Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11483941/)

| EO:Li | Model | PEO ether O/Li | TFSI O/Li | Sum |
|---:|---|---:|---:|---:|
| 10 | Original OPLS | 3.03 | 2.37 | 5.40 |
| 10 | Ionic charge scaling 0.79 + modified Li–O LJ | 4.20 | 1.13 | 5.33 |
| 18 | Original OPLS | 4.53 | 1.21 | 5.74 |
| 18 | Ionic charge scaling 0.79 + modified Li–O LJ | 5.32 | 0.10 | 5.42 |

The redistribution between polymer and anion is larger than the change in the total shell. The second model changes electrostatic and short-range interaction parameters. These aggregate MD results do not contain the timestamps and pair identities needed for residence estimates. A different salt ratio, sample molar mass, thermal history or force field needs a separate calculation.

## Count contacts from a trajectory

The contact-series GUI accepts uniformly spaced, consecutive frames with stable zero-based IDs for all Li ions and each oxygen population. Empty lists explicitly record a frame with zero contacts. Use positions from an equilibrated bulk interval. For a neutral PEO/LiTFSI box, supply every ether O site and all four sulfonyl O atoms per TFSI. State ID, sampling interval, temperature, population and cutoff are shared across the runs of a comparison. Each replica supplies its force-field source, trajectory origin and unique identifier.

For an orthorhombic box, the included adapter computes minimum-image contacts from coordinates in nm. For a triclinic box, compute pairs with a triclinic-aware trajectory engine, then use the same JSON input. Terminal hydroxyl oxygens have a different chemical environment. Select ether oxygens explicitly and state the convention for your EO:Li ratio. The GUI ratio for imported events is the count of supplied ether oxygen sites per Li.

```python
from polymer_electrolyte_lab.shell import orthorhombic_contacts, summarize_shell

# Each atom array has shape (frames, selected atoms, 3).
# lengths_nm is (3,) or (frames, 3), never a full cell matrix.
frames = orthorhombic_contacts(li_nm, ether_o_nm, tfsi_o_nm, lengths_nm, 0.27)
record = dict(state_id="PEO-film-01", temperature_K=300, n_li=li_nm.shape[1],
              n_ether_o=ether_o_nm.shape[1], dt_ps=output_interval_ps, cutoff_nm=0.27,
              max_lag=10, runs=[dict(state_id="PEO-film-01", model=model_label,
              replica=independent_run_id, force_field_source=parameter_file_and_version,
              trajectory_source=trajectory_and_equilibrated_interval, frames=frames)])
report = summarize_shell(record)
```

For a chosen lag, continuous survival divides the number of atom pairs present at every intervening frame by the number of pairs occupied at eligible time origins. Re-entry after departure starts a new contact. Trapezoidal integration ends at the selected lag. It reports a finite observation-window persistence area; loss and re-formation of contacts between saved frames cannot be resolved. Compare frame spacing and cutoff sensitivity. Mean coordination and areas are aggregated with equal weights per independent replica. For two or more independent estimates the sample SD is shown. Null marks an unestimable between-run spread or absent time-origin contact.

## Resolve tracer versus collective motion

Use an equilibrated, homogeneous, fixed-volume cell, physical ionic charges ±e in the current estimator, and a shared time-origin grid. Molecular TFSI positions must be unwrapped centers of mass of whole anions; unwrap components before evaluating the polymer center of mass. All times are ps and coordinates nm. The double sum implicit in the signed charge moment accounts for distinct-ion correlations. Common translation cancels from the total charge displacement in a neutral box. Subtracting the neutral polymer's translation is important for the individual tracer terms.

```python
from polymer_electrolyte_lab.displacements import ion_moments, displacement_result
moments = ion_moments(unwrapped_ions_nm, unwrapped_polymer_com_nm, output_interval_ps,
                      n_li, lag_steps=[10, 20, 30, 40, 50, 60, 70, 80])
run = dict(state_id=state_id, model=model_label, replica=independent_run_id,
           force_field_source=parameter_file_and_version,
           trajectory_source=trajectory_and_equilibrated_interval, moments=moments)
record = dict(state_id=state_id, temperature_K=T_K, volume_nm3=V_nm3, n_li=n_li,
              frame="polymer_com", transport_charge="formal_1e",
              fit_ps=[10*output_interval_ps, 80*output_interval_ps], runs=[run])
report = displacement_result(record)
```

Each row contains the lag, Li per-ion MSD, TFSI per-ion MSD, and the **extensive** charge-displacement moment ⟨|Σᵢ zᵢ Δrᵢ|²⟩ in nm². In 3D, each tracer diffusion coefficient is one sixth of its MSD slope. Collective conductivity equals the slope of the extensive moment times e²/(6VkBT). The NE baseline uses N e²(DLi + DTFSI)/(VkBT) for N LiTFSI formula units. For example, an analytic linear-moment fixture gives DLi=1×10⁻¹¹ m²/s and DTFSI=2×10⁻¹¹ m²/s in a 20 nm³ box containing four formula units at 300 K. Its constructed charge correlation gives 0.186 mS/cm, compared with a 0.372 mS/cm NE baseline. These fixture values are a unit check.

The user chooses the fit window before comparing outcomes. At least six lag points are required. The screening rule requires a log–log exponent of 0.8–1.2, two positive half-window slopes and a difference of at most 30% of the full-window slope. Analyses suppress individual estimates failing this heuristic. Overlapping origins are correlated. Review long-time plots, run-length and size sensitivity, alternate fit windows, and independent trajectory estimates before comparison with matched PFG-NMR and EIS measurements. The printed 300 K reference shell is a separate cohort from the 90 °C EIS/tracer examples.
