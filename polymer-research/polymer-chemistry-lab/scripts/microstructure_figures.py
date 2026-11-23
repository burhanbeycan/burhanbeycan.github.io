"""Render the two review figures from the Chemistry Lab models."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from polymer_chemistry_lab.plugin import calculate

DEST = Path(__file__).resolve().parents[1]/"docs"/"assets"

def draw_sequence_validation():
    r=calculate({"chemistry_mode":"pva_benchmark"})
    folds=r["pva_benchmark"]["sample_holdouts"]
    with plt.rc_context({"svg.fonttype":"none","font.size":10,"font.family":"DejaVu Sans"}):
        fig,ax=plt.subplots(figsize=(7.2,3.5),layout="constrained")
        x=range(len(folds))
        ax.axhline(1,color="#adb9bd",linestyle=":",label="Bernoulli reference")
        ax.scatter(x,[p[3] for p in folds],s=65,color="#0e7280",label="Measured NMR η",zorder=3)
        ax.scatter(x,[p[5] for p in folds],s=55,marker="x",color="#9f4c2c",label="Held-out prediction",zorder=4)
        ax.set_xticks(x,["r4","r16","r23","b9","b12","b15"])
        ax.set_xlabel("Acetylated PVA sample (r = solution, b = film)")
        ax.set_ylabel("Block character η"); ax.set_ylim(0,1.14)
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",color="#e0e6e9");ax.set_axisbelow(True)
        ax.legend(loc="lower left",bbox_to_anchor=(0,1.01),ncol=3,fontsize=8,frameon=False)
        fig.savefig(DEST/"pva-grade-holdout.svg",metadata={"Date":None});plt.close(fig)

def draw_charged_chain():
    r=calculate({"chemistry_mode":"petox_hydrolysis","ei_pattern":"EEIPEEIP"})
    mol=Chem.MolFromSmiles(r["tables"][0]["rows"][-1][5])
    colors={}
    for a in mol.GetAtoms():
        if a.GetFormalCharge()==1:colors[a.GetIdx()]=(.95,.63,.35)
        elif a.GetSymbol()=="N" and a.GetTotalNumHs()==1:colors[a.GetIdx()]=(.35,.77,.71)
    draw=rdMolDraw2D.MolDraw2DSVG(820,230)
    draw.drawOptions().legendFontSize=17
    draw.DrawMolecule(mol,highlightAtoms=list(colors),highlightAtomColors=colors,
                      legend="DP 8  |  4 EtOx, 2 neutral EI (teal), 2 protonated EI (amber)  |  2 chloride ions")
    draw.FinishDrawing()
    (DEST/"petox-ei-microstate.svg").write_text(draw.GetDrawingText(),encoding="utf-8")

if __name__=="__main__":
    draw_sequence_validation()
    draw_charged_chain()
