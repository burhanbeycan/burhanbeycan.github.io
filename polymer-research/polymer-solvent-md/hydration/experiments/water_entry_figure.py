"""Draw a tracked water in two measured trajectory frames and export its geometric measurements."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from contact_pair import make_whole, mc


def angle_deg(x,y):
    return float(np.degrees(np.arccos(np.clip(np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y)),-1,1))))


def number(value,precision=3):
    return ('%.*f'%(precision,value)).replace('.',',')


def extract_pose(xyz,box,top,A,B,W):
    L=np.diag(box)
    a=make_whole(xyz,top['chain_a_indices'],top['bonds'],L)
    b=make_whole(xyz,top['chain_b_indices'],top['bonds'],L)
    # Present the same pair of closest imaged sites used by the contact selection.
    j=top['chain_b_indices'].index(B)
    b+=a[top['chain_a_indices'].index(A)]+mc(xyz[B]-xyz[A],L)-b[j]
    origin=a[top['chain_a_indices'].index(A)]
    aw=mc(xyz[W]-xyz[A],L)
    water=origin+aw+np.vstack([np.zeros(3),mc(xyz[W+1]-xyz[W],L),mc(xyz[W+2]-xyz[W],L)])
    return np.vstack([a,b,water])-origin


def measures(xyz,box,A,B,W,seed,t):
    L=np.diag(box)
    wa=mc(xyz[A]-xyz[W],L);wb=mc(xyz[B]-xyz[W],L)
    h=[mc(xyz[W+i]-xyz[W],L) for i in (1,2)]
    angles=[angle_deg(-v,wa-v) for v in h]
    da=float(np.linalg.norm(wa));db=float(np.linalg.norm(wb))
    ab=float(np.linalg.norm(mc(xyz[A]-xyz[B],L)))
    return dict(seed=seed,time_ps=t,atom_a=A,element_a='O',atom_b=B,element_b='C',water_o=W,element_w='O',
                a_b_nm=ab,a_ow_nm=da,b_ow_nm=db,a_ow_b_deg=angle_deg(wa,wb),
                water_donor_h= W+1+int(np.argmax(angles)),ow_h_a_deg=max(angles),
                satisfies_water_to_a_hbond=bool(da<=.35 and max(angles)>=150),
                pair_open_and_water_interposed=bool(ab>.45 and da<=.35 and db<=.35 and angle_deg(wa,wb)>=120))


def aligned_poses(x0,x1,ia,ib,iw):
    aa=x0[ia]-x0[ia].mean(axis=0);bb=x1[ia]-x1[ia].mean(axis=0)
    u,_,v=np.linalg.svd(bb.T@aa)
    sign=np.linalg.det(u@v)
    rot=u@np.diag([1,1,np.sign(sign)])@v
    x1=x1@rot
    # Site A is the origin; orient the camera perpendicular to the measured interaction plane at event time.
    ex=x1[ib]/np.linalg.norm(x1[ib])
    ey=x1[iw]-ex*np.dot(ex,x1[iw]);ey/=np.linalg.norm(ey)
    # A fixed tilt makes the covalent carbonyl fragment legible instead of projecting it end-on.
    ez=np.cross(ex,ey)
    tilt=np.deg2rad(-38)
    ey=ey*np.cos(tilt)+ez*np.sin(tilt)
    camera=np.column_stack([ex,ey,np.cross(ex,ey)])
    return x0@camera,x1@camera


def draw_pose(ax,xyz,top,A,B,W,limits,hbond,t):
    n=len(top['chain_a_indices'])+len(top['chain_b_indices'])
    heavy=np.array([i for i in range(n) if top['elements'][i]!='H'])
    # Only chemically connected two-bond fragments are shown around the tracked atoms.
    near={A,B}
    edges=[(i,j) for i,j in top['bonds'] if i in heavy and j in heavy]
    for _ in range(2):
        near |= {i for i,j in edges if j in near} | {j for i,j in edges if i in near}
    colorA='#425e82';colorB='#a26738';water='#008c86'
    ax.set_facecolor('#f2f6f4')
    for i,j in top['bonds']:
        if i in near and j in near:
            x,y=xyz[[i,j],:2].T
            col=colorA if i<n//2 else colorB
            if A in (i,j):
                norm=np.array([-(y[1]-y[0]),x[1]-x[0]]);norm*=.008/np.linalg.norm(norm)
                for sign in [-1,1]:
                    ax.plot(x+sign*norm[0],y+sign*norm[1],lw=1.7,color=col,solid_capstyle='round',zorder=2)
            else:
                ax.plot(x,y,lw=2.8,color=col,alpha=.69,solid_capstyle='round',zorder=2)
    for i in near:
        if i in (A,B): continue
        color=colorA if i<n//2 else colorB
        ax.scatter(*xyz[i,:2],s=25,c=color,alpha=.68,zorder=3)
    water_i=n
    for i in (n+1,n+2):
        xy=xyz[[n,i],:2]
        ax.plot(*xy.T,color=water,lw=3,zorder=4,solid_capstyle='round')
        ax.scatter(*xyz[i,:2],s=48,c='white',edgecolors=water,lw=1,zorder=5)
    ax.plot(*xyz[[A,B],:2].T,lw=1.7,ls=(0,(3,3)),color='#545d65',zorder=3.5)
    if hbond:
        for i in (A,B):
            ax.plot(*xyz[[water_i,i],:2].T,lw=1.7,ls=(0,(2,2)),color=water,zorder=3.6)
    focal=[(A,'A',colorA), (B,'B',colorB), (water_i,'W',water)]
    for i,letter,color in focal:
        ax.scatter(*xyz[i,:2],s=240,facecolor=color,edgecolor='white',lw=1.6,zorder=6)
        ax.annotate(letter,xy=xyz[i,:2],xytext=(0,-.3),textcoords='offset points',ha='center',va='center',
                    color='white',fontsize=8.4,weight='bold',zorder=7)
    ax.set_xlim(*limits[:2]);ax.set_ylim(*limits[2:]);ax.set_aspect('equal')
    x=limits[0]+.05;y=limits[2]+.055
    ax.plot([x,x+.20],[y,y],color='#314849',lw=2.4)
    ax.text(x,y+.014,'0,20 nm',fontsize=8,color='#314849',va='bottom')
    ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_color('#d9e1dc')
    ax.set_title(t,fontsize=12.3,weight='bold',loc='left',pad=8,color='#183c42')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',required=True,type=Path);ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--event-ps',type=float,default=10);ap.add_argument('--a',type=int,default=3)
    ap.add_argument('--b',type=int,default=61);ap.add_argument('--water',type=int,default=1243)
    args=ap.parse_args();folder=args.run;out=args.out;out.mkdir(parents=True,exist_ok=True)
    A,B,W=args.a,args.b,args.water
    top=json.loads((folder/'topology.json').read_text());meta=json.loads((folder/'protocol.json').read_text())
    with np.load(folder/'trajectory.npz',allow_pickle=False) as src:
        time=src['time_ps'];pos=src['positions_nm'];box=src['box_nm']
    if (A not in top['carbonyl_oxygens_a'] or B not in top['chain_b_indices'] or top['elements'][B]!='C'
            or W not in top['water_oxygens']): raise ValueError('Site selection does not match figure chemistry')
    sel=np.flatnonzero(time==args.event_ps)
    if len(sel)!=1: raise ValueError('Event frame unavailable')
    k=int(sel[0]);seed=meta['seed']
    rows=[measures(p,b,A,B,W,seed,float(t)) for p,b,t in zip(pos[:k+1],box[:k+1],time[:k+1],strict=True)]
    zero,last=rows[0],rows[-1]
    if zero['a_b_nm']>.45 or not last['pair_open_and_water_interposed']: raise ValueError('Selected frames do not satisfy event definition')
    csvpath=out/'water_entry.csv'
    with csvpath.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(zero));w.writeheader();w.writerows([zero,last])
    ia=top['chain_a_indices'];n=len(ia)+len(top['chain_b_indices'])
    heavy=[i for i in ia if top['elements'][i]!='H']
    poses=aligned_poses(extract_pose(pos[0],box[0],top,A,B,W),extract_pose(pos[k],box[k],top,A,B,W),heavy,B,n)
    focus=np.vstack([p[[A,B,n,n+1,n+2],:2] for p in poses]);x,y=focus.T
    mid=np.array([(min(x)+max(x))/2,(min(y)+max(y))/2])
    width=max(np.ptp(x)+.45,(np.ptp(y)+.38)*1.62,1.16);height=width/1.62
    limits=[mid[0]-width/2,mid[0]+width/2,mid[1]-height/2,mid[1]+height/2]
    fig=plt.figure(figsize=(11,8.1),facecolor='#fcfcf9')
    fig.text(.06,.94,'İki PEtOx zinciri arasına giren tek bir su',fontsize=23,color='#193b43',weight='bold',ha='left')
    fig.text(.06,.90,'H/H uçlu iki trimer · 700 TIP3P su · 298,15 K · 1 bar · bağımsız koşu '+str(seed),
             fontsize=11.5,color='#526969')
    def leg(x,c,label):
        fig.text(x,.857,'●',fontsize=13,color=c,va='center');fig.text(x+.022,.857,label,fontsize=10,color='#263d42',va='center')
    leg(.06,'#425e82',f'A · Karbonil O[{A}]');leg(.31,'#a26738',f'B · Alifatik C[{B}]');leg(.55,'#008c86',f'W · Aynı suyun O[{W}] atomu')
    ax0=fig.add_axes([.06,.495,.419,.31]);ax1=fig.add_axes([.515,.495,.419,.31])
    draw_pose(ax0,poses[0],top,A,B,W,limits,False,'0 ps  ·  İki zincir temas halinde')
    draw_pose(ax1,poses[1],top,A,B,W,limits,True,f'{args.event_ps:g} ps  ·  Aynı su arada')
    fig.text(.06,.445,'A–B   '+number(zero['a_b_nm'])+' nm',fontsize=12,weight='bold',color='#243b43')
    fig.text(.515,.445,'A–B   '+number(last['a_b_nm'])+' nm       A···W···B   '+number(last['a_ow_b_deg'],1)+'°',fontsize=12,weight='bold',color='#243b43')
    for start,r in ((.06,zero),(.515,last)):
        fig.text(start,.415,'A···W   '+number(r['a_ow_nm'])+' nm       B···W   '+number(r['b_ow_nm'])+' nm',fontsize=11,color='#405556')
    ax=fig.add_axes([.1,.174,.81,.151]);ax.set_facecolor('#fcfcf9')
    ax.plot(time[:k+1],[r['a_b_nm'] for r in rows],c='#60545c',marker='o',ms=3,label='A–B')
    ax.plot(time[:k+1],[r['a_ow_nm'] for r in rows],c='#206a95',marker='o',ms=3,label='A···W')
    ax.plot(time[:k+1],[r['b_ow_nm'] for r in rows],c='#ae7751',marker='o',ms=3,label='B···W')
    ax.axhline(.35,c='#58938e',lw=.9,ls=':',label='0,35 nm su komşuluğu')
    ax.set_xlabel('Serbest MD süresi (ps)',fontsize=10);ax.set_ylabel('Uzaklık (nm)',fontsize=10)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.18);ax.tick_params(labelsize=9)
    ax.legend(ncol=4,loc='lower center',bbox_to_anchor=(.5,1.015),frameon=False,fontsize=9)
    fig.text(.06,.079,
        'A ve B çevresinde iki bağlık zincir parçaları aynı kamera yönünde gösterildi. Uzaklıklar üç boyutta periyodik sınırlarla ölçüldü.',
        fontsize=10,color='#455d60')
    fig.text(.06,.049,
        f'{args.event_ps:g} ps’de W–H[{last["water_donor_h"]}]···A açısı {number(last["ow_h_a_deg"],1)}°. B···W alifatik karbona yakınlığı gösterir.',
        fontsize=10,color='#455d60')
    fig.savefig(out/'water_entry.png',dpi=180,facecolor=fig.get_facecolor());plt.close(fig)

if __name__=='__main__':
    main()
