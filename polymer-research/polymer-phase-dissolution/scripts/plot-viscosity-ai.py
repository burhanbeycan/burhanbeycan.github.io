"""Render the released held-out predictions, without refitting or choosing test rows by their scores."""
import argparse
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


INK = "#202D38"
MUTED = "#4D606D"
GRID = "#DFE5E7"
BASE = "#915418"
GP = "#007E72"
QUAD = "#4E5BA2"


def locate(root, name):
    matches = list(root.rglob(name))
    if len(matches) != 1:
        raise ValueError(f"Expected one {name} under {root}, found {len(matches)}")
    return matches[0]


def tr(number, digits=2):
    return f"{number:.{digits}f}".replace(".", ",")


def axes_style(ax):
    ax.spines[["right","top"]].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=10, length=3)
    for side in ("bottom","left"):
        ax.spines[side].set_color("#8C9AA1")
        ax.spines[side].set_linewidth(.75)
    ax.set_axisbelow(True)


def load(root):
    ps = json.loads(locate(root, "viscosity-gp-validation.json").read_text(encoding="utf-8"))
    hpam_path = locate(root, "formulation-model.js")
    result = subprocess.run(["node", "-e",
        "process.stdout.write(JSON.stringify(require(process.argv[1])))",
        str(hpam_path.resolve())], check=True, capture_output=True, text=True)
    return ps, json.loads(result.stdout)


def plot(ps, hpam):
    plt.rcParams.update({"font.family":"DejaVu Sans","text.color":INK,"axes.labelcolor":INK})
    fig = plt.figure(figsize=(7.0, 8.55), facecolor="white")
    fig.text(.065,.951,"Viskozite modellerinin sınanması",size=19, weight="bold")
    fig.text(.068,.915,"Eğitim dışında tutulan numune ve formülasyonlarda tahmin hatası",size=11,color=MUTED)
    fig.text(.068,.853,"A   Polistirenin bağıl viskozitesi",size=13.2,weight="bold")
    fig.text(.068,.821,"Siklohekzanda (34,5 °C) bir fraksiyonun tümü, toluende en yüksek derişim ayrıldı.",
             size=9.3,color=MUTED)
    ax = fig.add_axes([.34,.564,.60,.228])
    axes_style(ax)
    ax.spines["left"].set_visible(False)
    ax.set_xscale("log")
    ax.tick_params(axis="y",length=0)
    ax.set_xlim(.25,20)
    ax.set_xticks([.5,1,2,5,10,20],["0,5","1","2","5","10","20"])
    ax.grid(axis="x", color=GRID, lw=.72)
    ax.set_xlabel("Bağıl hataların karelerinin ortalamasının karekökü (%)",size=9.7)
    specs = [
       ("luh1967-ps-cyclohexane-633000", "wholeGrade","Siklohekzan · 633 kg/mol"),
       ("luh1967-ps-cyclohexane-1770000","wholeGrade","Siklohekzan · 1770 kg/mol"),
       ("phywe-ps-toluene-25","highest","Toluen · 25 °C"),
       ("copenhaver-ps-toluene-30","highest","Toluen · 30 °C")
    ]
    labels = []
    for i,(identifier,strategy,label) in enumerate(specs):
        rows = [r for r in ps["within"]+ps["transfer"]
                if r["testSeriesId"]==identifier and r["strategy"]==strategy]
        if len(rows)!=1:
            raise ValueError(f"Missing or duplicate control for {identifier}")
        metrics=rows[0]["metrics"]
        labels.append(label+"\n"+f'n = {metrics["baseline"]["n"]}')
        for key, offset, color, marker in (("baseline",-.15,BASE,"o"),("gp",.15,GP,"D")):
            value=metrics[key]["rmse_pct"]
            ax.scatter(value,i+offset,marker=marker,color=color,s=30,zorder=3)
            ax.text(value*1.10,i+offset,tr(value),color=color,size=8.9,va="center")
    ax.set_yticks(range(4),labels,fontsize=9)
    ax.set_ylim(3.51,-.52)
    ax.legend(handles=[
        Line2D([],[],color=BASE,ls="",marker="o",label="Huggins"),
        Line2D([],[],color=GP,ls="",marker="D",label="Huggins + Gauss süreci")
    ],loc="upper left",bbox_to_anchor=(-.43,-.25),ncol=2,frameon=False,
       fontsize=9,handletextpad=.35,columnspacing=1.6)
    fig.text(.068,.430,"B   HPAM çözeltisinin görünür viskozitesi",size=13.2,weight="bold")
    fig.text(.068,.399,"24 °C · Her sınamada bir polimer–tuz formülasyonunun bütün kayma hızları ayrıldı.",
             size=9.3,color=MUTED)
    panel = fig.add_axes([.112,.105,.448,.258])
    axes_style(panel)
    panel.set_xlim(0,21); panel.set_ylim(0,21)
    panel.set_xticks([0,5,10,15,20]); panel.set_yticks([0,5,10,15,20])
    panel.set_xlabel("Ölçülen η (mPa·s)",size=9.6)
    panel.set_ylabel("Tahmin edilen η (mPa·s)",size=9.6)
    panel.grid(color=GRID,lw=.65)
    panel.plot([0,21],[0,21],ls=(0,(4,3)),c="#858F97",lw=1,zorder=2)
    m = hpam["test"]
    style=[("baseline","Kuvvet yasası",BASE,"o",28),
           ("ridge","Ridge · ikinci derece",GP,"D",28),
           ("quadratic","Cezasız ikinci derece",QUAD,"x",30)]
    fig.text(.594,.347,"Sınama RMSE'si",size=11.6,weight="bold")
    fig.text(.594,.321,"15 koşul birlikte · mPa·s",size=9.8,color=MUTED)
    for i,(key,label,color,marker,size) in enumerate(style):
        obs=np.array([p["y"] for p in m]); pred=np.array([p[key] for p in m])
        metric=np.sqrt(np.mean((pred-obs)**2))
        if not np.isclose(metric,hpam["scores"][key]["rmse"]):
            raise ValueError(f"Mismatch in released {key} validation metric")
        kwargs={"c":color} if marker=="x" else {"facecolors":"none","edgecolors":color}
        panel.scatter(obs,pred,marker=marker,s=size,lw=1.25,zorder=4,**kwargs)
        y=.282-i*.057
        fig.add_artist(Line2D([.601],[y+.004],transform=fig.transFigure,ls="",marker=marker,
                              color=color,markerfacecolor="none",markersize=6,markeredgewidth=1.4))
        fig.text(.624,y,label,size=9.7,color=color)
        fig.text(.624,y-.022,tr(metric,3),size=10.2,weight="bold",color=INK)
    fig.text(.594,.103,f'{hpam["nGroups"]} formülasyon · {hpam["nRuns"]} ölçüm',size=9.3,color=MUTED)
    fig.text(.065,.035,"HPAM merkezindeki altı tekrar tek koşul ortalamasıdır. Kesikli çizgi ölçümle eşitliği gösterir.",
             size=8.7,color=MUTED)
    return fig


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dist",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    fig=plot(*load(args.dist))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=200,facecolor="white")
    plt.close(fig)
