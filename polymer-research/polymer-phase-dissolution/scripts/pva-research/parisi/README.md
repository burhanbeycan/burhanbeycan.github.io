# PVA solution steady-shear records from Parisi et al.

The tables contain original numerical readings of identifiable plotted symbols in Figures 4 and 5 of D. Parisi, C. D. Ditillo, A. Han, S. Lindberg, M. W. Hamersky and R. H. Colby, *Journal of Rheology* 66 (2022), 1141–1150. Physical measurements belong to the cited study. Extracted values are approximate. [Primary article](https://doi.org/10.1122/8.0000435) · [University copy](https://pure.rug.nl/ws/files/495697175/8.0000435.pdf)

| File | Solvent | Matched grade and concentration combinations | Points |
|---|---|---:|---:|
| parisi_water.csv | Water | 15 | 283 |
| parisi_dmso.csv | DMSO | 6 | 101 |
| parisi_unresolved.csv | Water | Unresolved identity, exclude from grade-dependent models | 79 |

Use the solvent explicitly when concatenating files. Do not pool these rows as unrelated samples. All represent the same laboratory publication and all temperatures are 25 °C (298.15 K). All targets are steady-shear solution viscosities. The original y-axis is Pa s. Multiplying by 1000 yields cP. The concentrations are weight percent of polymer in the solution. Mn, Mw and dispersity are sample-level supplier values transcribed from Table I. Hydrolysis is approximate, and the authors discuss confirming it by NMR. Tables I and III are not flow-curve sources. No intrinsic viscosity, fit-generated prediction, low-shear plateau estimate or complex viscosity from an oscillatory test appears as a target in these CSVs.

## Point selection

The two source graphs are losslessly embedded RGB raster images. Each small panel has its own x-axis and a shared or panel-specific y-axis. Separate log10 affine transformations from labelled decades appear in the scripts. The coordinates x_px and y_px refer to pixel centres in the original embedded Figure 4 or 5. Points with independent outlines were picked by colour-mask connected components and checked against their published marker shapes. In dense bands, sharp separable tips or vertical extents were used when a symbol could be located. For the touching red-square bands in Figure 4A and 4C, just the two ends were kept: intermediate squares have ambiguous centres. No point was manufactured by assuming a regularly spaced shear-rate grid. Obscured, superposed and cropped symbols can therefore be missing in the tables.

| Water figure panel, concentration | 171k-80 | 144k-88 | 49k-75 | 24k-80 |
|---|---:|---:|---:|---:|
| 4A, 25 wt% | 2 | 8 | 21 | 0 |
| 4B, 17 wt% | 22 | 15 | 25 | 25 |
| 4C, 15 wt% | 2 | 35 | 31 | 18 |
| 4D, 4 wt% | 21 | 21 | 21 | 16 |

Figure 4’s blue triangle is labelled **48k-80**. Tables I and III list **48k-88**. Both strings are retained as published. The blue water records have empty grade metadata in parisi_unresolved.csv. Figure 5 independently labels its blue triangle **48k-88**. The grade match for that DMSO panel is unambiguous. The ambiguity in Figure 4 was preserved.

| DMSO figure panel, concentration | 171k-80 | 144k-88 | 49k-75 | 48k-88 |
|---|---:|---:|---:|---:|
| 5B, 15 wt% | 12 | 22 | 0 | 0 |
| 5C, 10 wt% | 16 | 21 | 15 | 15 |

Figure 5A is annotated 20 wt% in the image and 25 wt% in the caption. Panel D has 7 and 8 wt% annotations, with 7 wt% in the caption. Both panels were omitted. In 5B, the lower purple series is annotated 16.4 wt%, while the caption says 15 wt%; it was omitted. The two high-molar-mass grades in 5B are labelled 15 wt%, and Figure 7 independently identifies the two high-mass DMSO series at 15 wt%.

## Reproduction and digitization sensitivity

Install Python 3, numpy, scipy, Pillow and pymupdf. Obtain the original paper from the link above. The reference PDF SHA-256 is 20ca13fc37fe2ab077376545bc2165cbb9f4da14b519ec2e571c39f77575ff12. The extraction scripts verify pixel hashes, so a differently packaged publisher PDF with identical images can also work.

```sh
python digitize.py /path/to/paper.pdf
python digitize_dmso.py /path/to/paper.pdf
python check_extraction.py /path/to/paper.pdf
```

Add `--audit-image overlay.png` to either extraction command for local comparison. The overlays were visually checked against the symbol outlines for all selected panels. Independent dark-pixel tests verified printed decade ticks. A wider colour threshold and centre-position check gave a maximum shift of 1.0 pixel in Figure 4 and 0.5 pixel in Figure 5. A conservatively assigned ±2-pixel coordinate sensitivity corresponds to a factor of approximately 1.071 in viscosity for water and 1.057 and 1.043 for the DMSO top and bottom axes, respectively. This quantifies reading resolution only. The experimental repeatability, sample-to-sample variation and calibrated predictive uncertainty must be estimated separately. Decimal digits in machine-readable values serve reproducibility, rather than precision.

For one solvent contrast in the published 15 wt% series, grade 144k-88 has log–log interpolated viscosity at 10 s⁻¹ of about 9.7×10³ cP in water and 2.0×10⁴ cP in DMSO, at 298.15 K. These are approximate readings of the same-source curves. They provide a unit- and solvent-key check, not an independent validation of solvent transfer.
