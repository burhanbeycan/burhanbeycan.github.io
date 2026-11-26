"""Render saved discrete MD samples, using a fixed camera per run and minimum images."""
from pathlib import Path
import argparse
import csv
import json
import math
import subprocess
import numpy as np
from PIL import Image, ImageDraw, ImageFont

APP = Path(__file__).resolve().parents[1]
DEFAULT_TRAJECTORIES = APP.parent / 'polymer-solvent-md' / 'hydration' / 'data' / 'pair'
W, H = 920, 1100
BG = (244,246,245)
INK = (21,43,52)
MUTED = (66,91,99)
TEAL = (54,204,184)
LILAC = (177,150,247)
OXY = (254,124,103)
NIT = (106,164,252)
HYD = (229,239,247)
GOLD = (255,201,97)
HB = (145,237,135)
COLOUR_A = (21,124,111)
COLOUR_B = (109,70,169)
fontdir=Path('/usr/share/fonts/truetype/dejavu')
def F(size, bold=False):
    return ImageFont.truetype(str(fontdir / ('DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf')), size)
fonts={(s,b):F(s,b) for s in [15,16,17,18,19,20,21,22,23,24,26,27,28,30,32,34,35,40,43] for b in [False,True]}
def txt(d,xy,s,size,colour=INK,b=False,anchor='la'):
    d.text(xy,s,fill=colour,font=fonts[(size,b)],anchor=anchor)
def tr(v,n):
    return f'{v:.{n}f}'.replace('.',',')
def angular(u,v):
    return np.degrees(np.arccos(np.clip(u@v/(np.linalg.norm(u)*np.linalg.norm(v)),-1,1)))
def dashed(d,a,b,colour,width=3,spacing=11):
    a,b=np.array(a),np.array(b)
    n=max(1,math.ceil(np.linalg.norm(b-a)/spacing))
    for i in range(n):
        d.line([tuple(a+(b-a)*i/n), tuple(a+(b-a)*min((i+.47)/n,1))],fill=colour,width=width)
def ball(d,xy,r,c):
    x,y=xy
    for i,m in [(0,.60),(1,.77),(2,.91),(3,1.0)]:
        rr=r*(1-i*.115)
        dx,dy=-r*.07*i,-r*.07*i
        d.ellipse((x+dx-rr,y+dy-rr,x+dx+rr,y+dy+rr),fill=tuple(int(a*m) for a in c))
def unit(x):
    return x/np.linalg.norm(x)
class Run:
    def __init__(self,seed,chapter,trajectory_directory):
        p=trajectory_directory/f'seed_{seed}'
        self.seed,self.chapter=seed,chapter
        top=json.loads((p/'topology.json').read_text())
        traj=np.load(p/'trajectory.npz')
        self.xyz=traj['positions_nm'].astype('float64')
        self.boxes=traj['box_nm'].astype('float64')
        self.times=traj['time_ps']
        self.elements=top['elements']
        self.bonds=top['bonds']
        self.carbonyl=top['carbonyl_oxygens_a']+top['carbonyl_oxygens_b']
        with (p/'interposed_waters.csv').open() as f:
            row=next(csv.DictReader(f))
        self.time=int(float(row['time_ps']))
        self.a,self.b,self.w=[int(row[k]) for k in ['atom_a','atom_b','water_oxygen']]
        self.initial=float(row['initial_ab_nm'])
        self.water=top['water_oxygens'].index(self.w)+1
        self.stop=self.time+5
        el=self.elements
        show=[i for i in top['chain_a_indices']+top['chain_b_indices'] if el[i]!='H']
        self.show=show+[self.w,self.w+1,self.w+2]
        self.lines=[(i,j) for i,j in self.bonds if i in self.show and j in self.show]
        self.adj={i:[] for i in self.show}
        for i,j in self.lines:
            self.adj[i].append(j)
            self.adj[j].append(i)
        self.qs=[self.geometry(t) for t in range(self.stop+1)]
        assert self.qs[self.time]['hit'] and not any(q['hit'] for q in self.qs[:self.time])
        q=self.qs[self.time]
        assert max(abs(q[key]-float(row[rk])) for key,rk in [('ao','a_ow_nm'),('bo','b_ow_nm'),('ab','current_ab_nm'),('angle','angle_a_ow_b_deg')])<0.001
        # Orient the x–y plane along the event triangle. Keep this rotation for the entire chapter.
        v=q['unwrapped']; x=unit(v[self.b]-v[self.a]); plane_v=v[self.w]-(v[self.a]+v[self.b])/2; y=unit(plane_v-x*(plane_v@x)); z=np.cross(x,y)
        self.rot=np.stack((x,-y,z),axis=1)
        self.views=[]
        for q in self.qs:
            xyz=q['unwrapped']; mid=(xyz[self.a]+xyz[self.b])/2; self.views.append({i:(v-mid)@self.rot for i,v in xyz.items()})
        self.labels=['A · '+self.elements[self.a]+str(self.a),'B · '+self.elements[self.b]+str(self.b),'W'+str(self.water)]
        self.colour=[TEAL,LILAC,GOLD]
    def geometry(self,t):
        pos=self.xyz[t]; box=self.boxes[t]; inv=np.linalg.inv(box)
        def mic(v):
            frac=v@inv
            return (frac-np.rint(frac))@box
        a,b,w=self.a,self.b,self.w
        un={}
        for root,base in [(a,np.zeros(3)),(b,mic(pos[b]-pos[a])),(w,mic(pos[w]-pos[a]))]:
            un[root]=base; stack=[root]
            while stack:
                i=stack.pop()
                for j in self.adj[i]:
                    if j not in un:
                        un[j]=un[i]+mic(pos[j]-pos[i]);stack.append(j)
        ao,bo=mic(pos[a]-pos[w]),mic(pos[b]-pos[w])
        ab=np.linalg.norm(mic(pos[a]-pos[b]))
        dista,distb=np.linalg.norm(ao),np.linalg.norm(bo)
        angle=angular(ao,bo)
        hbs=[]
        for i,v in ((a,dista),(b,distb)):
            if v<=.35 and i in self.carbonyl:
                for h in (w+1,w+2):
                    ang=angular(mic(pos[w]-pos[h]),mic(pos[i]-pos[h]))
                    if ang>=150: hbs.append((h,i,ang))
        return dict(unwrapped=un,ao=dista,bo=distb,ab=ab,angle=angle,hbs=hbs,hit=self.initial<=.45 and ab>.45 and max(dista,distb)<=.35 and angle>=120)
def pick_label_pos(d,xy,text,color,occupied,points):
    x,y=xy; bx=d.textbbox((0,0),text,font=fonts[18,True]); lw=bx[2]-bx[0]+14; lh=29
    candidates=[]
    for r in [49,67,91]:
        for ang in [-45,45,-135,135,-90,90,0,180]:
            a=math.radians(ang);cx=x+r*math.cos(a);cy=y+r*math.sin(a)
            box=(cx-lw/2,cy-lh/2,cx+lw/2,cy+lh/2)
            over=sum(max(0,min(box[2],p[2])-max(box[0],p[0]))*max(0,min(box[3],p[3])-max(box[1],p[1])) for p in occupied)
            dim=sum(box[0]-9 < p[0] < box[2]+9 and box[1]-9 < p[1] < box[3]+9 for p in points)
            penalty=(max(0,457-box[0])+max(0,box[2]-868)+max(0,269-box[1])+max(0,box[3]-552))*200
            candidates.append((over+dim*95+penalty+r*.3,box,cx,cy))
    _,box,cx,cy=min(candidates,key=lambda x:x[0])
    d.line((x,y,cx,cy),fill=color,width=2)
    d.rounded_rectangle(box,radius=6,fill=(19,39,52),outline=color,width=1)
    txt(d,(cx,cy),text,18,color,True,'mm')
    occupied.append(box)
def draw_scene(d,run,t):
    q=run.qs[t]; a,b,w=run.a,run.b,run.w
    focus=[a,b,w,w+1,w+2]
    def projection(inds,centre,scl):
        v=run.views[t]
        p=np.array([v[i] for i in inds]); mid=(p.max(0)+p.min(0))/2
        return {i:np.array([centre[0]+(v[i][0]-mid[0])*scl,centre[1]+(v[i][1]-mid[1])*scl,v[i][2]]) for i in inds}
    def extent_scale(inds,dim):
        all_spans=[np.ptp(np.array([v[i] for i in inds]),axis=0) for v in run.views]
        ext=np.max(all_spans,axis=0)
        return min(dim[0]/ext[0],dim[1]/ext[1])
    pp=projection(run.show,(240,403),extent_scale(run.show,(340,250)))
    for i,j in sorted(run.lines,key=lambda b:(pp[b[0]][2]+pp[b[1]][2])/2):
        x,y=pp[i][:2],pp[j][:2]; c=(47,124,124) if i<50 else ((115,97,155) if i<100 else (202,194,158))
        d.line((tuple(x),tuple(y)),fill=c,width=5)
    for i in sorted(pp,key=lambda i:pp[i][2]):
        el=run.elements[i]; cc={'C':TEAL if i<50 else LILAC,'O':OXY,'N':NIT,'H':HYD}[el]
        ball(d,pp[i][:2],8 if el in ('O','N') else (4 if el=='H' else 7),cc)
    for j,cc in [(a,TEAL),(b,LILAC),(w,GOLD)]:
        x,y=pp[j][:2]
        d.ellipse((x-12,y-12,x+12,y+12),outline=cc,width=2)
    # Independent, fixed magnification of the same projected event. Labels maintain element identities.
    d.line((444,239,444,540), fill=(54,76,88),width=1)
    txt(d,(49,241),'Zincirler ve izlenen su',19,(176,198,211))
    txt(d,(463,241),'Yakın plan',19,(176,198,211))
    pp=projection(focus,(663,406),extent_scale(focus,(318,225)))
    for i in (a,b):
        dashed(d,pp[i][:2],pp[w][:2],(136,162,172),3,11)
    for h,j,ang in q['hbs']:
        dashed(d,pp[h][:2],pp[j][:2],HB,5,12)
    for h in (w+1,w+2):
        d.line([tuple(pp[w][:2]),tuple(pp[h][:2])],fill=(218,211,189),width=8)
    for j,cc in [(a,TEAL),(b,LILAC),(w,GOLD)]:
        x,y=pp[j][:2]
        d.ellipse((x-17,y-17,x+17,y+17),outline=cc,width=3)
    for j in (a,b,w,w+1,w+2):
        el=run.elements[j]
        cc=OXY if el=='O' else (HYD if el=='H' else LILAC)
        ball(d,pp[j][:2],13 if el=='O' else (6 if el=='H' else 12),cc)
    occupied=[tuple([p[0]-15,p[1]-15,p[0]+15,p[1]+15]) for p in [pp[i] for i in focus]]
    for k,j in enumerate((a,b,w)):
        pick_label_pos(d,pp[j][:2],run.labels[k],run.colour[k],occupied,[p[:2] for p in pp.values()])
def chip(d,x,y,s,fg,fill):
    width=d.textlength(s,font=fonts[18,True])+22
    d.rounded_rectangle((x,y,x+width,y+31),15,fill=fill);txt(d,(x+11,y+15),s,18,fg,True,'lm')
    return width
def spark(d,r,t):
    xl,xr,yt,yb=112,854,776,869
    mx=max(.75,max(max(q['ao'],q['bo']) for q in r.qs)*1.03)
    def X(i):return xl+(xr-xl)*i/r.stop
    def Y(v):return yb-(yb-yt)*v/mx
    for val in [0,.35,.70]:
        y=Y(val);d.line((xl,y,xr,y),fill=(206,215,216),width=1)
        txt(d,(xl-13,y),tr(val,2),16,MUTED,False,'rm')
    dashed(d,(xl,Y(.35)),(xr,Y(.35)),(136,151,157),2,14)
    for key,colour in [('ao',COLOUR_A),('bo',COLOUR_B)]:
        pts=[(X(k),Y(q[key])) for k,q in enumerate(r.qs[:t+1])]
        if len(pts)>1:d.line(pts,fill=colour,width=4)
        x,y=pts[-1];d.ellipse((x-5,y-5,x+5,y+5),fill=colour)
    for k,tag,an in [(0,'0 ps','la'),(r.time,str(r.time)+' ps','ma'),(r.stop,str(r.stop)+' ps','ra')]:
        txt(d,(X(k),yb+6),tag,16,MUTED,False,an)
    x=X(t);d.line((x,yt,x,yb),fill=(139,147,150),width=1)
def paint(r,t):
    im=Image.new('RGB',(W,H),BG); d=ImageDraw.Draw(im); q=r.qs[t]
    d.rectangle((0,0,W,7),fill=(27,108,107))
    txt(d,(36,26),'POLYMER RESEARCH STUDIO',19,COLOUR_A,True)
    txt(d,(882,27),'MOLEKÜLER DİNAMİK',18,MUTED,False,'ra')
    txt(d,(34,62),'İki zincir arasında bir su',43,INK,True)
    txt(d,(36,127),'PEtOx  •  İki trimer (DP = 3)  •  700 su molekülü  •  298,15 K',23,MUTED)
    d.rounded_rectangle((28,174,892,605),19,fill=(19,38,51))
    txt(d,(48,192),f'KOŞU {r.chapter} / 2     •     Tohum {r.seed}',19,(175,201,207),True)
    txt(d,(870,185),f'{t:g} ps',32,HYD,True,'ra')
    draw_scene(d,r,t)
    txt(d,(51,565),'A zinciri',19,TEAL,True);txt(d,(163,565),'B zinciri',19,LILAC,True)
    txt(d,(278,565),'O',19,OXY,True);txt(d,(311,565),'N',19,NIT,True);txt(d,(343,565),'H',19,HYD,True)
    dashed(d,(410,578),(446,578),HB,4,9);txt(d,(461,565),'H bağı',19,HYD)
    txt(d,(866,565),'Hedef su ve zincir ağır atomları',17,(166,193,204),False,'ra')
    d.rounded_rectangle((28,620,892,719),13,fill=(255,255,255))
    typ1='O';typ2=r.elements[r.b]; labels=[f'A({typ1})–O(su)',f'B({typ2})–O(su)','A–B','∠ A–O(su)–B']
    vals=[tr(q['ao'],3)+' nm',tr(q['bo'],3)+' nm',tr(q['ab'],3)+' nm',tr(q['angle'],1)+'°']
    for k,(label,val) in enumerate(zip(labels,vals)):
        x=49+k*213
        if k:d.line((x-13,636,x-13,704),fill=(219,224,225),width=1)
        txt(d,(x,632),label,18,[COLOUR_A,COLOUR_B,MUTED,MUTED][k],True)
        txt(d,(x,664),val,27,INK,True)
    txt(d,(40,734),'Su oksijenine uzaklık (nm)',21,INK,True)
    txt(d,(873,737),'Kesikli eşik  0,35 nm',18,MUTED,False,'ra')
    spark(d,r,t)
    hit=q['hit']; first=t==r.time
    current=('İLK ARA KONUM' if first else ('ARA KONUM' if hit else 'İZLENEN SU'))
    chip(d,39,908,current,(24,91,75) if hit else MUTED,(206,239,216) if hit else (225,232,232))
    hbs=q['hbs']
    if hbs:
        which=' ve '.join(['A' if ix==r.a else 'B' for _,ix,_ in hbs])
        note=f'W{r.water}, {which} karboniline H bağı veriyor.'
    else:note=f'W{r.water} ile seçilen atomlar arasında H bağı gözlenmedi.'
    txt(d,(39,952),note,23,INK,True)
    if first:
        detail='A tarafında H-bağı açısı 170,4°. B tarafında alifatik karbon var.' if r.chapter==1 else 'B karbonili 150,3° ile eşiği geçiyor. A tarafındaki en iyi açı 133,9°.'
    else:detail=f'W{r.water} her karede aynı su. Gri çizgiler iki uzaklığı gösterir.'
    txt(d,(39,990),detail,19,MUTED)
    d.line((37,1032,883,1032),fill=(204,214,215))
    txt(d,(39,1051),'OpenMM  •  GAFF2 / AM1-BCC  •  TIP3P',17,MUTED)
    txt(d,(881,1051),'Kayıt aralığı 1 ps',17,MUTED,False,'ra')
    return im
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trajectories', type=Path, default=DEFAULT_TRAJECTORIES)
    parser.add_argument('--output', type=Path, default=APP)
    options = parser.parse_args()
    out = options.output
    out.mkdir(parents=True, exist_ok=True)
    rs = [Run(73111, 1, options.trajectories), Run(86201, 2, options.trajectories)]
    for r in rs:
        print('Seed', r.seed, 'event frame', r.time, 'H bonds', r.qs[r.time]['hbs'], flush=True)
    fps=10
    args=['ffmpeg','-y','-v','error','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(fps),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','21','-pix_fmt','yuv420p','-movflags','+faststart',str(out/'PEtOx_Iki_Temas.mp4')]
    with subprocess.Popen(args,stdin=subprocess.PIPE) as p:
        for r in rs:
            for t in range(r.stop+1):
                im=paint(r,t)
                rep=33 if t==r.time else (13 if t in [0,r.stop] else 4)
                p.stdin.write(im.tobytes()*rep)
        p.stdin.close()
        if p.wait()!=0:raise RuntimeError('ffmpeg encode failed')
    paint(rs[0],rs[0].time).save(out/'PEtOx_Iki_Temas.png')
if __name__=='__main__':main()
