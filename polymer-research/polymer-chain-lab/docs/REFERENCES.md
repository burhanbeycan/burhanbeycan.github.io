# References and provenance

Sources were checked on 23 November 2026. Literature and official documentation support the equations, terminology, and scope. Source availability is not a benchmark result; no external experimental dataset has been fitted or validated in this release.

| ID | Source | Use in this project |
| --- | --- | --- |
| `iupac_dp` | [IUPAC Gold Book: degree of polymerization](https://goldbook.iupac.org/terms/view/D01569), DOI 10.1351/goldbook.D01569 | Counts monomeric units, not end groups. The implemented end-corrected mapping is the direct mass-balance rearrangement M=DP M0+Mend; common repeat and terminal chemistry must be supplied. |
| `iupac_mn` | [IUPAC Gold Book: number-average molar mass](https://goldbook.iupac.org/terms/view/12216), DOI 10.1351/goldbook.12216, 2025 online edition | Molar-mass terminology, units, number/mass weighting |
| `iupac_dispersity` | [IUPAC Gold Book: molar-mass dispersity](https://goldbook.iupac.org/terms/view/12224), DOI 10.1351/goldbook.12224, 2025 online edition | Mass-average/number-average ratio; preferred term dispersity |
| `sasview_debye` | [SasView official documentation: mono_gauss_coil](https://www.sasview.org/docs/user/models/mono_gauss_coil.html) | Debye formula and ideal Gaussian-chain assumptions; documentation cites P. Debye, J. Phys. Colloid Chem. 51 (1947), 18 |
| `wlc_dimensions` | [AFRC official documentation: Worm-like chain (O'Brien)](https://afrc.readthedocs.io/en/stable/polymer_models/worm_like_chain_obrien.html) | Explicit Benoit–Doty mean-square gyration formula and rod/coil limits; no AFRC code or end-to-end-distribution approximation is imported |
| `benoit_doty_1953` | [H. Benoit and P. Doty, “Light Scattering from Non-Gaussian Chains”](https://doi.org/10.1021/j150510a025), J. Phys. Chem. (1953) | Primary source for the non-Gaussian-chain/WLC context; bibliographic record verified; this release implements dimensions, not the paper's complete scattering treatment |
| `iupac_ris` | [IUPAC Gold Book: rotational isomeric state](https://goldbook.iupac.org/terms/view/12187), DOI 10.1351/goldbook.12187, 2025 online edition | Definition of a discrete rotational state; does not prescribe chemistry parameters |
| `abe_jernigan_flory_1966` | [A. Abe, R. L. Jernigan, and P. J. Flory, “Conformational Energies of n-Alkanes and the Random Configuration of Higher Homologs Including Polymethylene”](https://doi.org/10.1021/ja00956a001), J. Am. Chem. Soc. 88 (1966), 631–639 | Primary RIS chemical context; bibliographic record verified; no fitted energy constants copied into this app |

The finite-chain forward/backward recurrences follow directly from the explicitly stated finite partition sum in SCIENCE.md. The implementation and its enumeration tests are original project code. Neither a citation nor an imported material record is treated as evidence that the synthetic presets represent a real polymer.
