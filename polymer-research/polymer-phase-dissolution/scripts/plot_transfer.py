"""Plot the experimental checks from the exposed model results; no curve refitting."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def render(data, output):
    models = data["models"]
    points = data["coexistence"]["points"]
    fig, axes = plt.subplots(2, 1, figsize=(8, 8.9))
    fig.subplots_adjust(top=.865, bottom=.112, left=.127, right=.97, hspace=.60)
    fig.suptitle("Ozmotik ölçüm eklendiğinde", x=.127, y=.97, ha="left", fontsize=23, fontfamily="DejaVu Serif")
    fig.text(.127,.918,"Polistiren / siklohekzan   •   Kalibrasyona girmeyen ölçümler",color="#405c54",fontsize=11)
    styles = [("#ad582c", "--"), ("#126d64", "-")]
    labels = ["Yalnız bulutlanma", "Bulutlanma + ozmometri"]
    for ax in axes:
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y", alpha=.24)
        ax.set_axisbelow(True)
    ax=axes[0]
    ax.set_title("A   İki fazın bileşimi arasındaki fark", loc="left",fontsize=14,pad=9)
    points=sorted(points,key=lambda r:r["deltaT_K"])
    x=[r["deltaT_K"] for r in points]
    ax.set_xscale("log")
    for m,(color,ls),name in zip(models,styles,labels):
        preds=[r["predictions"][m["id"]] for r in points]
        y=[r["phi_rich"]-r["phi_lean"] if r["phases"]==2 else 0.0 for r in preds]
        ax.plot(x,y,color=color,linestyle=ls,linewidth=2.1,label=name)
        ax.scatter([t for t,r in zip(x,preds) if r["phases"]==1],[0 for r in preds if r["phases"]==1],facecolors="none",edgecolor=color,s=20)
    ax.scatter(x,[r["phi_rich"]-r["phi_lean"] for r in points],s=23,facecolor="#23312c",zorder=4,label="Deneysel eşfaz çifti")
    ax.set_ylabel("Hacim kesri farkı  Δφ")
    ax.set_xlabel("Ölçülen kritik sıcaklığa uzaklık  $T_{c,deney} - T$  (K)")
    ax.legend(frameon=False,loc="upper left",fontsize=9)
    ax.text(0,-.32,"M$_w$ = 1,56 MDa  ·  20 sıcaklık, bir PS numunesi  ·  Boş daire: model tek faz veriyor.",transform=ax.transAxes,fontsize=9,color="#405c54")
    ax=axes[1]
    ax.set_title("B   Seyreltik çözeltinin ikinci ozmotik virial katsayısı",loc="left",fontsize=14,pad=9)
    k=data["constants"]; t=np.linspace(295.5,313.8,180)
    for m,(color,ls),name in zip(models,styles,labels):
        y=1e4*k["vbar_ml_g"]**2/k["V_s_ml_mol"]*(.5-m["A"]-m["B"]/t+m["C"])
        ax.plot(t-273.15,y,color=color,linestyle=ls,linewidth=2,label=name)
    op=data["osmotic"]["holdoutPoints"]
    ax.errorbar([p["T_K"]-273.15 for p in op],[1e4*p["A2_mol_ml_g2"] for p in op],
                yerr=[1e4*p["reported_halfwidth_mol_ml_g2"] for p in op],fmt="o",color="#23312c",markersize=4.5,capsize=4,label="Bağımsız ozmometri")
    ax.axhline(0,color="#7a827e",linewidth=.7)
    ax.set_xlabel("Sıcaklık  (°C)")
    ax.set_ylabel("$A_2$  ($10^{-4}$ mL mol g$^{-2}$)")
    ax.legend(loc="lower right",frameon=False,fontsize=9)
    fig.text(.127,.045,"A: Nakata ve ark. (1978), Wohlfarth tablosu 1978NAK.   B: Striolo ve ark. (2000), Tablo 3.",fontsize=8.8,color="#405c54")
    fig.text(.127,.024,"Çizgiler aynı iki modelin hesaplarıdır. Alt panelde hata çubukları kaynakta verilen yarı genişliklerdir.",fontsize=8.6,color="#405c54")
    fig.savefig(output,dpi=185,bbox_inches="tight",facecolor="white")
    plt.close(fig)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    text=args.results.read_text()
    render(json.loads(text.split("=",1)[1].strip().removesuffix(";")),args.output)
