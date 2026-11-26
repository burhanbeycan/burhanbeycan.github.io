"""Create the Turkish comparison figure from the exact fold predictions shown in the app."""
import argparse
import json
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def score(rows, key, metric, aggregation):
    def one(items):
        true = np.array([r['observed'] for r in items])
        estimate = np.array([r[key] for r in items])
        if metric == 'MAPE':
            return 100*np.mean(np.abs(estimate/true-1))
        return np.sqrt(np.mean((estimate-true)**2))
    if aggregation == 'pooled':
        return one(rows)
    return np.mean([one([r for r in rows if r['group']==g]) for g in dict.fromkeys(r['group'] for r in rows)])


def panel(ax, entries, xlabel):
    for i, (label, old, new) in enumerate(entries):
        ax.plot([old,new],[i,i],color='#9cabae',lw=2.4,zorder=1)
        ax.scatter(old,i,s=77,marker='o',color='#334b5a',zorder=2)
        ax.scatter(new,i,s=83,marker='D',color='#b05226',zorder=3)
        ax.annotate(f'{old:.2f}'.replace('.',','),(old,i),xytext=(0,10),textcoords='offset points',ha='center',color='#334b5a',fontsize=9)
        ax.annotate(f'{new:.2f}'.replace('.',','),(new,i),xytext=(0,-18),textcoords='offset points',ha='center',color='#a64622',fontsize=9)
    ax.set_yticks(range(len(entries)),[r[0] for r in entries],fontsize=10)
    ax.set_ylim(len(entries)-.1,-.7)
    ax.set_xlim(0,max(max(x[1:]) for x in entries)*1.14)
    ax.set_xlabel(xlabel,ha='right',x=1,labelpad=9)
    ax.grid(axis='x',color='#dce3e1',lw=.8)
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=10)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dist',type=Path,default=Path('dist'))
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    cmd="process.stdout.write(JSON.stringify(require(process.argv[1]).protocols))"
    reports=json.loads(subprocess.check_output(['node','-e',cmd,str((args.dist/'viscosity-improvement-data.js').resolve())]))
    records={p['id']:p for p in reports}
    def entry(id,label,group=None):
        p=records[id]
        rows=[r for r in p['rows'] if group is None or r['group']==group]
        return (label,*[score(rows,k,p['metric'],p['aggregation']) for k in ('baseline','candidate')])
    a=[entry('theta-one-reference','θ · kayıtlı referans (51)'),
       entry('theta-reference-weights','θ · referansa AI ağırlığı (51)'),
       entry('hydro-anchor','PI / PIB · başka yayından referans (18)'),
       entry('hydro-colby','Polibütadien · ayrı çalışma (8)'),
       entry('theta-zero-reference','θ · referanssız Ridge ayarı (56)')]
    b=[entry('pvoh-solvent-transfer','Su ve DMSO birlikte (10)'),
       entry('pvoh-solvent-transfer','Su (5)','Su'),entry('pvoh-solvent-transfer','DMSO (5)','DMSO')]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelcolor':'#31474b','text.color':'#173735'})
    fig=plt.figure(figsize=(10.8,8.6),facecolor='white')
    gs=fig.add_gridspec(2,1,height_ratios=[1.6,1],left=.35,right=.94,top=.78,bottom=.12,hspace=.67)
    ax=fig.add_subplot(gs[0])
    bx=fig.add_subplot(gs[1])
    panel(ax,a,'Ortalama mutlak yüzde hata (%)')
    panel(bx,b,'Ortalama kare hatanın karekökü (mL/g)')
    ax.set_title('Hidrokarbon polimerlerinde intrinsik viskozite',loc='left',pad=20,fontweight='bold',fontsize=13)
    bx.set_title('PVOH · ortak hidroliz etkisi ile çözücüye özgü iki eğim',loc='left',pad=20,fontweight='bold',fontsize=13)
    fig.text(.065,.955,'Viskozite tahminindeki hata nerede azaldı?',fontsize=23,weight='bold',va='top')
    fig.text(.065,.899,'Her satırda aynı saklı numuneler karşılaştırılır. Parantez içindeki sayılar numune–koşul sayısıdır.',fontsize=11)
    fig.legend(handles=[Line2D([],[],lw=0,marker='o',color='#334b5a',label='Karşılaştırma modeli'),Line2D([],[],lw=0,marker='D',color='#b05226',label='Sınanan yöntem')],loc='upper left',bbox_to_anchor=(.055,.877),ncol=2,frameon=False)
    fig.text(.065,.045,'θ karşılaştırmaları aile dengeli, PI / PIB karşılaştırması kimya dengelidir. Yalnız polibütadien içeren satır bir yayına dayanır.',fontsize=9,color='#405c5c')
    fig.savefig(args.output,dpi=160,bbox_inches='tight',pad_inches=.3)
    plt.close(fig)
if __name__=='__main__':
    main()
