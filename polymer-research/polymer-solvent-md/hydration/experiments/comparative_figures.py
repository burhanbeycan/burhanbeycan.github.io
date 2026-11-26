"""Publication-style Turkish plots drawn from the measured replica arrays."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

INK = "#203947"
HH = "#157566"
CAP = "#A34731"
TINT = "#E1E9EA"
PALETTE = ["#13786F", "#DA8144", "#476C99"]
LABELS = {"HH":"H/H","MeOH":"Me/OH"}

def num(v, dec=2):
    return f"{v:.{dec}f}".replace(".",",")

def style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="y",linewidth=.65,color=TINT)
    ax.spines["left"].set_color("#8DA5A8")
    ax.spines["bottom"].set_color("#8DA5A8")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(colors=INK,width=.7,length=3)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x,p:f"{x:g}".replace(".",",")))

def points(ax, case, data, positions, color, j=0.0):
    for pos,d in zip(positions,data):
        x = pos+j
        vals = d["per_run_mean"]
        ax.scatter(x+np.linspace(-.045,.045,len(vals)), vals, s=23,color=color,alpha=.7,edgecolors="white",linewidths=.55,zorder=3)
        ax.errorbar(x,d["mean"],yerr=d["between_run_sd"],fmt="D",ms=6,capsize=4,
                    markeredgecolor="white",markeredgewidth=.65,color=color,linewidth=1.45,zorder=4)

def figures(sumry, trace, out):
    out=Path(out)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.labelcolor":INK,"text.color":INK,
                         "axes.titlecolor":INK,"axes.labelsize":10,"axes.titlesize":11,"figure.facecolor":"white"})
    fig,axs=plt.subplots(2,2,figsize=(10.5,7),gridspec_kw={"height_ratios":[1.05,1]})
    fig.subplots_adjust(top=.83,bottom=.16,hspace=.65,wspace=.32,left=.087,right=.96)
    fig.text(.084,.959,"PEtOx trimeri suda hangi zincir boyutlarını örnekledi?",fontsize=18,weight="medium")
    dur=sumry["matched_duration_ps"]
    fig.text(.085,.917,f"Her uç modeli için 3 koşu · {dur:g} ps üretim/koşu · 298,15 K · 500 TIP3P su molekülü",fontsize=11)
    legend=[Line2D([0],[0],color=c,lw=2,label=l) for c,l in zip(PALETTE,("Koşu 1","Koşu 2","Koşu 3"))]
    all_rg=[r["rg_nm"] for r in trace]
    pad=.06*(max(all_rg)-min(all_rg))+.002
    for ax,case in zip(axs[0],("HH","MeOH")):
        style(ax)
        for i,seed in enumerate(sumry["cases"][case]["seeds"]):
            rows=[r for r in trace if r["seed"]==seed and r["model"]==case]
            ax.plot([r["time_ps"] for r in rows],[r["rg_nm"] for r in rows],lw=1,alpha=.83,color=PALETTE[i])
        ax.set_title(f'{LABELS[case]} · Tüm zincirin jirasyon yarıçapı',loc="left")
        ax.set_ylabel("Rg (nm)")
        ax.set_xlabel("Yeni üretim süresi (ps)")
        ax.set_xlim(0,dur)
        ax.set_ylim(min(all_rg)-pad,max(all_rg)+pad)
    axs[0,1].legend(handles=legend,loc="upper right",ncol=3,fontsize=8,frameon=False,handlelength=1.2)
    ax=axs[1,0]; style(ax)
    sequences=[]
    if "archived_pilot" in sumry:
        sequences.append(("H/H\neski, 24 ps",sumry["archived_pilot"]["metrics"]["rg_nm"],"#809396"))
    sequences.extend([(LABELS[c]+"\n"+f"yeni, {dur:g} ps",sumry["cases"][c]["metrics"]["rg_nm"],color)
                      for c,color in [("HH",HH),("MeOH",CAP)]])
    for i,(name,stat,col) in enumerate(sequences):
        points(ax,name,[stat],[i],col)
    ax.set_xticks(range(len(sequences)),[e[0] for e in sequences],fontsize=9)
    ax.set_xlim(-.55,len(sequences)-.45)
    ax.set_ylabel("Rg (nm)");ax.set_title("Koşu ortalamaları ve aralarındaki dağılım",loc="left")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax=axs[1,1]; style(ax)
    for i,(case,col) in enumerate([("HH",HH),("MeOH",CAP)]):
        s=sumry["cases"][case]["metrics"]["rg_nm"]
        y=np.array(s["per_run_blocks"])
        for rep in y:
            ax.plot(np.arange(1,5),rep,color=col,alpha=.27,lw=.8)
        ax.plot(np.arange(1,5),y.mean(axis=0),color=col,marker="o",ms=4,lw=2,label=LABELS[case])
    ax.legend(frameon=False,ncol=2,fontsize=8)
    ax.set_title("Dört ardışık zaman bloğunda zincir boyutu",loc="left")
    ax.set_xlabel("Üretim bloğu (zaman sırasıyla)");ax.set_ylabel("Blok ortalaması Rg (nm)");ax.set_xticks([1,2,3,4])
    fig.text(.085,.055,"Noktalar tekil koşular, elmaslar üç koşunun ortalamasıdır. Hata çubukları koşu ortalamalarının örneklem SS’sini gösterir.",fontsize=9)
    fig.savefig(out/"convergence.png",dpi=180);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(10.5,4.9))
    fig.subplots_adjust(top=.73,bottom=.24,wspace=.29,left=.085,right=.97)
    fig.text(.08,.924,"Uç grupları değişince hangi karbonil çevresi değişiyor?",fontsize=17.5,weight="medium")
    fig.text(.081,.85,"α birimi başlatıcı ucuna, ω birimi hidroksillenebilen zincir ucuna komşudur.",fontsize=11)
    labels=["α-uç","Orta","ω-uç"]
    for ax,suffix,title,ylab in zip(axs,("_neighbors","_water_hb"),("Su oksijeni komşuları","Sudan karbonile yönlü hidrojen bağları"),
                                     ("Su O sayısı / karbonil O","H bağı sayısı / karbonil O")):
        style(ax)
        vb=[]
        for i,(case,col) in enumerate([("HH",HH),("MeOH",CAP)]):
            ms=[sumry["cases"][case]["metrics"][site+suffix] for site in ("alpha","middle","omega")]
            points(ax,case,ms,np.arange(3),col,j=(-.14,.14)[i])
            vb.extend([max(max(s["per_run_mean"]),s["mean"]+s["between_run_sd"]) for s in ms])
        ax.set_title(title,loc="left");ax.set_ylabel(ylab)
        ax.set_ylim(0,max(vb)*1.13)
        ax.set_xticks(np.arange(3),labels);ax.set_xlim(-.5,2.5)
    handles=[Line2D([0],[0],marker="D",markersize=7,color=col,lw=0,label=label) for col,label in [(HH,"H/H"),(CAP,"Me/OH")]]
    fig.legend(handles=handles,fontsize=9,ncol=2,loc="upper right",bbox_to_anchor=(.972,.88),frameon=False)
    fig.text(.081,.125,"Su komşuları için O···O ≤ 0,35 nm; H bağı için ayrıca O–H···O ≥ 150°. Elmaslar ortalama, çubuklar koşular arası SS’dir.",fontsize=9)
    fig.text(.081,.073,"İki uç birlikte değiştirilmiştir. Bölgesel farklar sonlu oligomerin her iki uçta yapılan değişime verdiği birleşik yanıttır.",fontsize=9)
    fig.savefig(out/"hydration.png",dpi=180);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(10.5,4.4),sharey=True)
    fig.subplots_adjust(top=.71,bottom=.26,wspace=.15,left=.095,right=.975)
    fig.text(.09,.925,"Karbonilin yanındaki aynı su ne kadar süre seçili kalıyor?",fontsize=17,weight="medium")
    fig.text(.091,.842,"0,35 nm karbonil–su komşuluğu · üç koşunun ortalaması · çözünmüş DP 3 PEtOx",fontsize=10.5)
    for ax,case in zip(axs,("HH","MeOH")):
        style(ax)
        pts=sumry["cases"][case]["water_exchange"]
        for site,col,label in zip(("alpha","middle","omega"),PALETTE,("α-uç","Orta","ω-uç")):
            rs=sorted((v for v in pts if v["site"]==site),key=lambda r:r["lag_ps"])
            x=[p["lag_ps"] for p in rs]
            ax.plot(x,[r["endpoint_mean"] for r in rs],marker="o",ms=3,lw=1.6,color=col,label=label)
            ax.plot(x,[r["all_sampled_mean"] for r in rs],lw=1.3,color=col,linestyle="--")
        ax.set_title(LABELS[case],loc="left")
        ax.set_xlabel("Karşılaştırılan kareler arasındaki süre (ps)")
        ax.set_xticks([1,2,5,10]);ax.set_ylim(0,1.02)
    axs[0].set_ylabel("Aynı suyun hâlâ komşu olma oranı")
    axs[1].legend(frameon=False,loc="upper right",ncol=3,fontsize=8)
    fig.text(.091,.13,"Düz çizgi yalnız ilk ve son karede aynı suyu arar. Kesikli çizgi aradaki kaydedilmiş tüm karelerde komşuluk arar.",fontsize=9)
    fig.text(.091,.077,"Her değer başlangıçtaki komşuluklara koşulludur. Kareler arasında gerçekleşen kısa ayrılma ve dönüşler çözülmez.",fontsize=9)
    fig.savefig(out/"water_exchange.png",dpi=180);plt.close(fig)
