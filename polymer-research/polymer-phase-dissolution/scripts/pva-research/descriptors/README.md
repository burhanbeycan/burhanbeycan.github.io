# PVA descriptors and polymer structure encoders

`pva_grade_descriptors.csv` pairs supplier-reported molar masses and nominal degrees of hydrolysis for five PVA grades with dry-chain site densities. Intrinsic viscosities in the final two columns refer to the same five grades in two solvents at 298.15 K. These are limiting solution measurements in mL/g, distinct from apparent viscosity in cP. Source: Parisi et al., *Journal of Rheology*, Tables I and III. [Primary paper](https://doi.org/10.1122/8.0000435)

## Computational chemistry

For the mole fraction of vinyl-alcohol units h = DH/100, the mean repeat mass is 44.053 h + 86.090 (1−h) g/mol. The topology has one backbone methylene per repeat, one OH on each alcohol unit and one acetate ester on each residual acetate unit. The nominal OH site density is 1000 h / M_repeat mmol/g. The calculation neglects ends, defects, impurities and bound water. Site densities count structural donor/acceptor opportunities. Neither contact occupancy nor hydrogen-bond lifetime follows from these counts. Simultaneously fitting DH, OH density, acetate density and repeat mass gives multiple representations of a single compositional degree of freedom. For a grade-level small-data model preselect one of these scalars. A physically distinct covariate is the molar mass distribution. Mₙ, M_w, and viscosity-average M_v have different definitions; M_w = Đ Mₙ can only be used when both operands are known. The two lower moments do not generally fix M_v. [IUPAC M_v](https://goldbook.iupac.org/terms/view/12220) · [Mark–Houwink equation](https://goldbook.iupac.org/terms/view/M03706)

`structure.py` checks atom valences with RDKit, the number and valence of polymer attachment atoms, and the declared repeat fractions. Unknown DH yields an unknown composition descriptor, even if Mₙ has been supplied. For polyelectrolytes the model distinguishes the mass of the covalent repeat from the declared dry counterion basis. `sample.json` is an example calculation using a known PVA grade, not a measured viscosity at the example concentration. Pure-solvent topology is separated from bulk solvent viscosity and solvent quality. Hydrolysis percentage and pH do not specify HPAM ionization under every salt condition. The CMC parent motif represents cellulose and requires a separate substitution descriptor. The optional number-average repeat count DPn = (Mn − Mends) / Mrepeat requires a stated combined end-group mass and mass basis, including any counterions. Open-chain attachment at O1/C4 followed by hydrolytic capping recovers a beta-D-glucose unit (PubChem CID 64689). Alginate M and G are beta-D-mannuronate and alpha-L-guluronate 1→4 anions; the two parent motifs differ at C5. Their relative abundance and block sequence remain grade metadata. [Cellulose motif stereochemistry](https://pubchem.ncbi.nlm.nih.gov/compound/64689) · [Mannuronic acid](https://pubchem.ncbi.nlm.nih.gov/compound/446145) · [Alginate M/G chemistry](https://doi.org/10.3390/md19050264)

The code calculates u = c[η] only when sample identity, solvent and temperature match and concentrations refer to solution volume. A mass fraction or wt% needs a measured or separately estimated *solution* density, with its uncertainty and conditions recorded by the caller. “g / 20 mL solvent initially charged” cannot automatically be converted to g/mL final solution. Concentration source, end-group mass, tacticity, blockiness, salt, aging and preshear remain study metadata. An independently measured relaxation time can be combined with shear rate as Wi = λγ̇. Fitting λ on a held-out flow curve before predicting that same curve leaks the target. External [η] may be used in a reference-assisted model; if it is predicted upstream, use out-of-fold upstream predictions for training the downstream model.

## Frozen transformer comparison

`polybert_embeddings.npz` contains one 600-component float32 vector per key in `repeat_catalog.json`. We computed these with the published polyBERT architecture. Inputs are the listed PSMILES strings. The official tokenizer inserts the special tokens; the encoder’s evaluation-mode final layer is attention-mask mean-pooled, including those tokens. Each saved vector represents a polymer motif. Concentration, solvent, Mₙ, relaxation time and processing history must be supplied as other covariates. Graph-equivalent inputs with alternate SMILES spelling map to a registered input; equivalent alternative choices of the repeat boundary need chemical review. Unknown polymer samples retain unknown composition parameters.

An embedding shared by an entire family supplies no within-family grade or Mn variation. The 11 motifs include hydroxyl/acetate/acrylamide/acrylate and alginate subunits and should not be confused with 11 independently observed polymer families. We obtained the published checkpoint format from the revision-pinned xushijie mirror. The archived PyTorch binary has the same SHA-256 in both available mirrors (469d1928fe140f3d8d93f69e69a11aebf540e7a2043b38a9a6758cf0bb640a9c). The inference script downloads the safetensors version, verifies its hash 21cfab4a3b0599852fb8ef529a67e41632246b9256ed834938ba9b46cb24148f, and checks the config, tokenizer, and pooling-config hashes. The checkpoint is a 12-layer DeBERTa-v2 polymer-language encoder, from the model originally released by Kuenneth & Ramprasad. Tested with transformers 4.49.0, torch 2.7.1+cpu, tokenizers 0.21.4 and RDKit 2025.03.5. The research package contains calculated fingerprints only. Consult the Georgia Tech/GTRC release license when obtaining and using the original model. [Publication](https://doi.org/10.1038/s41467-023-39868-6) · [Pinned model mirror](https://huggingface.co/xushijie/polyBERT/tree/e7dce434fb3eff37905dc114008660e5479ca9a8) · [License with matching model mirror](https://huggingface.co/cooperzvegintzov/polyBERT/blob/fe976ad0160fe446f5358d75794608719f19d7a1/LICENSE) · [Authors’ release](https://github.com/Ramprasad-Group/polyBERT)

POLYT5 by Sahu et al. uses an encoder–decoder transformer with polymer-specific SELFIES and a distinct tokenizer. It reports results for Tg, Tm, Td, band gap, dielectric constant, solubility and Tg-conditioned generation. Its published results do not establish viscosity error. Its study describes using private training data. Any future POLYT5 ablation must use a verifiable POLYT5 checkpoint with the authors' conversion/tokenization, add the same state inputs, and reuse the same family splits. [POLYT5 paper](https://doi.org/10.1038/s44387-026-00087-1)

## Execute

Basic descriptor and unit checks require Python, NumPy and RDKit.

```bash
python descriptors/structure.py descriptors/sample.json
python descriptors/test_structure.py
python descriptors/build_descriptor_reference.py
```

With torch, transformers, tokenizers and safetensors installed, reproduce the published-model features or encode a new linear repeat:

```bash
python descriptors/encode_polybert.py --model-dir ./polybert-checkpoint --download
python descriptors/encode_polybert.py --model-dir ./polybert-checkpoint --repeat '[*]CC(O)[*]' --out query.npz
```

A query matching PVA in an independent batch reproduced the saved PVA fingerprint within 9.54×10⁻⁷ maximum absolute element difference on the tested CPU.
