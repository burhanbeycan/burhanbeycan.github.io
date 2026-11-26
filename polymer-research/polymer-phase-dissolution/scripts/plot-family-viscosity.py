"""Plot release benchmark predictions, using no additional fitting."""
import argparse
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np


def load_hydro(path):
    return json.loads(subprocess.check_output(
        ["node", "-e", "console.log(JSON.stringify(require(process.argv[1])))", str(path.resolve())], text=True))


def render(web, output):
    theta = json.loads((web / "family-inputs/theta-family-benchmark.json").read_text())
    hydro = load_hydro(web / "hydro-family-benchmark.js")
    pullulan = next(r for r in hydro["folds"] if r["group"] == "pullulan-water")
    primary, learned, ink, muted = "#9B542B", "#16705D", "#172E32", "#4E6570"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10.3,"axes.labelcolor":ink,
                         "axes.edgecolor":"#ACB6BB","xtick.color":muted,"ytick.color":muted})
    fig = plt.figure(figsize=(7.4,8.15), layout="none")
    fig.subplots_adjust(left=.24,right=.93,top=.843,bottom=.119,hspace=.72)
    ax, bx = fig.subplots(2,1, gridspec_kw={"height_ratios":[0.85,1.15]})
    fig.text(.055,.962,"Viskozite tahmini yeni bir polimere aktarılabilir mi?",size=16.4,weight="bold",color=ink)
    fig.text(.055,.926,"İntrinsik viskozite [η]  ·  Aile veya numune grubunun tamamı eğitim dışında",size=10.0,color=muted)
    fig.text(.055,.878,"A   θ koşulunda aile aktarımı",size=12.5,weight="bold",color=ink)
    for y, row in enumerate(reversed(theta["families"])):
        v1,v2=(row["scores"][s]["mape"] for s in ["theta","ridge"])
        ax.plot([v1,v2],[y,y],color="#bbc9c5",lw=2,zorder=1)
        ax.plot(v1,y,marker="o",color=primary,ms=6,zorder=2)
        ax.plot(v2,y,marker="D",color=learned,ms=5,zorder=2)
    names={"polydiene":"Polidienler","polyolefin":"Poliizobütilen","vinylaromatic":"Vinil aromatikler"}
    ax.set_yticks(range(3), [names[r["family"]]+'\n(n = '+str(len(r["folds"]))+')' for r in reversed(theta["families"])])
    ax.set_ylim(-.4,2.42); ax.set_xlim(0,60); ax.set_xticks([0,10,20,30,40,50,60])
    ax.set_xlabel("Ortalama mutlak yüzde hata (MAPE, %)")
    ax.grid(axis="x",alpha=.21); ax.spines[["top","right","left"]].set_visible(False); ax.tick_params(axis="y",length=0)
    ax.legend([Line2D([],[],ls="none",marker="o",color=primary),Line2D([],[],ls="none",marker="D",color=learned)],
              ["Yalnız mol kütlesi","Mol kütlesi + omurga tanımlayıcısı"],loc="lower left",bbox_to_anchor=(-.245,-.49),
              frameon=False,ncol=2,handlelength=1,columnspacing=1.6,fontsize=9.5)
    fig.text(.055,.475,"B   Pullulanın sulu ortamda bağımsız sınaması",size=12.5,weight="bold",color=ink)
    rows=sorted(pullulan["rows"],key=lambda r:r["M_g_mol"])
    x=np.asarray([r["M_g_mol"]/1000 for r in rows])
    for k, color, marker, label in [("eta_mL_g",ink,"o","Deney"),("mass_prediction",primary,"^","Yalnız M"),
                                     ("size_prediction",learned,"D","M ve ölçülmüş Rɡ")]:
        bx.plot(x,[r[k] for r in rows],color=color,marker=marker,lw=1.35,ms=5.2,label=label,
                markerfacecolor="white" if k=="eta_mL_g" else color)
    bx.set_xscale("log");bx.set_yscale("log")
    bx.set_xticks([200,500,1000,2000]);bx.set_yticks([50,100,200,500,1000])
    bx.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f"{v:.0f}"))
    bx.yaxis.set_major_formatter(FuncFormatter(lambda v,pos:f"{v:.0f}"))
    bx.set_xlabel("Pullulanın $M_w$ değeri (kg/mol)");bx.set_ylabel("[η] (mL/g)")
    bx.grid(alpha=.22,which="major");bx.spines[["top","right"]].set_visible(False)
    bx.legend(loc="upper left",frameon=False,fontsize=9,handlelength=1.5)
    bx.set_xlim(x.min()*.88,x.max()*1.16)
    fig.text(.055,.045,"A · 56 gözlem, 5 kimya, 3 omurga grubu.  B · Eğitimde 22 hidrokarbon ölçümü, sınamada 7 pullulan fraksiyonu.",size=8.5,color=muted)
    fig.text(.055,.024,"B’de polimer, çözücü ve kaynak çalışma birlikte değişir. Noktaları birleştiren çizgiler yalnızca okumayı kolaylaştırır.",size=8.5,color=muted)
    fig.savefig(output,dpi=190,facecolor="white");plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--dist",type=Path,required=True);parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();render(args.dist,args.output)
