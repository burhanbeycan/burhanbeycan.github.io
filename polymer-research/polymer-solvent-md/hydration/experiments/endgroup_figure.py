"""Draw the exact H/H and Me/OH finite-chain graphs from the two native SDFs."""
import argparse
import io
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D


def depiction(sdf, capped):
    mol = Chem.MolFromMolFile(str(sdf),removeHs=True)
    mol.RemoveAllConformers()
    ns = mol.GetSubstructMatches(Chem.MolFromSmarts("[O;X1]=[C;X3]-[N;X3]"))
    if capped:
        methyl,first_n = mol.GetSubstructMatch(Chem.MolFromSmarts("[CH3]-[N;X3]"))
        _,_,last_c,oh = mol.GetSubstructMatch(Chem.MolFromSmarts("[N;X3]-[CH2]-[CH2]-[O;H1]"))
        start,end = methyl,oh
    else:
        (first_n,) = mol.GetSubstructMatch(Chem.MolFromSmarts("[N;X3;H1]"))
        _,_,last_c = mol.GetSubstructMatch(Chem.MolFromSmarts("[N;X3]-[CH2]-[CH3]"))
        start,end=first_n,last_c
    backbone=list(Chem.GetShortestPath(mol,first_n,last_c))
    if len(backbone)!=9 or len(ns)!=3:
        raise ValueError("Backbone does not match a finite trimer.")
    # Align the identical nine backbone heavy atoms in both chemical drawings.
    coords={idx:Point2D(k*1.30,0.0 if k%3==0 else -.75) for k,idx in enumerate(backbone)}
    if capped:
        coords[methyl]=Point2D(-1.30,-.75)
        coords[oh]=Point2D(9*1.30,0.0)
    rdDepictor.Compute2DCoords(mol,coordMap=coords,canonOrient=False)
    draw=rdMolDraw2D.MolDraw2DCairo(1670,400)
    op=draw.drawOptions()
    op.padding=.055
    op.minFontSize=26
    op.maxFontSize=36
    op.bondLineWidth=3.8
    op.continuousHighlight=False
    op.atomHighlightsAreCircles=True
    # Specify terminal CH3 glyphs explicitly even when the atom is a skeletal vertex.
    op.atomLabels[last_c]="CH2" if capped else "CH3"
    if capped:
        op.atomLabels[methyl]="CH3"
    colors={start:(0.96,0.78,0.49),end:(0.47,0.81,0.76)}
    draw.DrawMolecule(mol,highlightAtoms=[start,end],highlightBonds=[],highlightAtomColors=colors,
                      highlightAtomRadii={start:.51,end:.51})
    draw.FinishDrawing()
    img=np.asarray(Image.open(io.BytesIO(draw.GetDrawingText())).convert("RGB"))
    ys,xs=np.where(np.min(img,axis=2)<225)
    img=img[max(0,ys.min()-8):ys.max()+9,max(0,xs.min()-8):xs.max()+9]
    return img,rdMolDescriptors.CalcMolFormula(mol)


def render(base,case,out):
    hh,fh=depiction(base/"pilot/parameters/solute.sdf",False)
    cap,fc=depiction(case/"parameters/solute.sdf",True)
    fig=plt.figure(figsize=(12,6.6),facecolor="white",layout=None)
    fig.text(.055,.93,"Hangi uç grupları değişiyor?",size=26,weight="semibold",color="#142c38")
    fig.text(.055,.865,"Poli(2-etil-2-oksazolin)   ·   DP 3   ·   Aynı üç propiyonil yan grubu ve N–CH₂–CH₂ omurgası",
             size=13,color="#284959")
    for yy,img,title,formula,caption in [
        (.475,hh,"H/H",fh,"α ucunda sekonder amid N–H    ·    ω ucunda –CH₂CH₃"),
        (.111,cap,"α-Me / ω-OH",fc,"α ucunda tersiyer amid N–CH₃    ·    ω ucunda –CH₂CH₂OH")]:
        ax=fig.add_axes([.175,yy+.030,.68,.265]);ax.imshow(img);ax.axis('off')
        fig.text(.055,yy+.21,title,size=16,weight="bold",color="#173c43")
        fig.text(.055,yy+.16,formula,size=11,color="#465e69")
        fig.text(.22,yy,caption,size=12,color="#324e55")
    fig.add_artist(plt.Line2D([.055,.945],[.424,.424],transform=fig.transFigure,linewidth=.8,color="#cbdadb"))
    fig.text(.055,.039,"Sol uç α (kehribar), sağ uç ω (yeşil). İki uç aynı modelde birlikte değişir.",size=10.2,color="#496368")
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/"end_groups.png",dpi=180,bbox_inches=None)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base",type=Path,required=True)
    p.add_argument("--case",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    render(a.base,a.case,a.out)
