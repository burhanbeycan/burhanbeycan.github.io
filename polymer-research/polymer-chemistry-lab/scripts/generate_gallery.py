"""Reproduce the documented finite-chain gallery from the scientific kernel."""
from pathlib import Path

from rdkit.Chem import Draw
from polymer_chemistry_lab.model import TEMPLATES, build_oligomer


def main():
    directory = Path(__file__).resolve().parents[1] / "docs" / "assets"
    directory.mkdir(parents=True, exist_ok=True)
    chains = [build_oligomer(item["smiles"], "AAAA") for item in TEMPLATES.values()]
    legends = [
        f"{key.upper()} | 4 repeats | H/H caps\n"
        f"{chain.summary['formula']} | {chain.summary['average_molar_mass_g_mol']:.3f} g/mol"
        for key, chain in zip(TEMPLATES, chains)
    ]
    picture = Draw.MolsToGridImage(
        [chain.molecule for chain in chains], molsPerRow=2,
        subImgSize=(620, 320), legends=legends, useSVG=False,
    )
    picture.save(directory / "finite-chain-gallery.png")
    print(directory / "finite-chain-gallery.png")


if __name__ == "__main__":
    main()
