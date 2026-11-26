"""Digitize visible Figure 4 steady-shear symbols from Parisi et al., J Rheol 2022.
Run python digitize.py path/to/original-paper.pdf. Coordinates below refer to the embedded
1000x833 pixel Figure 4 (not a rendered page). Only complete, assignable symbols are used.
"""
from pathlib import Path
from io import BytesIO
import argparse,csv,hashlib
import pymupdf
from PIL import Image, ImageDraw
import numpy as np
from scipy import ndimage
from scipy.signal import find_peaks
D=Path(__file__).resolve().parent
COL={'171k-80':(255,0,0),'144k-88':(255,132,0),'49k-75':(132,0,255),'24k-80':(0,132,0),'48k-80':(0,0,255)}
META={'171k-80':(112300,171000,1.52,80),'144k-88':(90900,144500,1.59,88),'49k-75':(28000,49500,1.77,75),'24k-80':(15100,24500,1.63,80)}
PANEL={'A':(25,(85,515,39,355),[(-2,123.5),(-1,199.6),(0,275.9),(1,352.0),(2,428.1),(3,504.0)]),'B':(17,(553,975,39,351),[(0,578.3),(1,705.0),(2,831.8),(3,958.5)]),'C':(15,(88,524,429,680),[(-1,113.6),(0,212.7),(1,312.0),(2,411.3),(3,510.4)]),'D':(4,(550,976,638,746),[(1,569.2),(2,755.0),(3,940.8)])}
Y_TOP=[(3,33),(2,100.1),(1,167.2),(0,234.3),(-1,301.4),(-2,368.5)]
Y_BOTTOM=[(2,416.4),(1,483.6),(0,550.7),(-1,617.9),(-2,685.1),(-3,752.2)]
def line(ticks):
 logs,pixels=zip(*ticks)
 m,b=np.polyfit(pixels,logs,1)
 return m,b
def extract_image(pdf):
 with pymupdf.open(pdf) as doc:
  matches=[q for page in doc for q in page.get_images() if q[2:4]==(1000,833)]
  for q in matches:
   im=Image.open(BytesIO(doc.extract_image(q[0])['image'])).convert('RGB')
   if hashlib.sha256(im.tobytes()).hexdigest()=='a11355c4ed1723fcf014110b9c28f426dd5cfcf88b9a4fe0974852b1da6db0b6':return im
 raise ValueError('Matching original Figure 4 RGB image not found.')
def detect_separate(mask,lo,hi):
 labels,_=ndimage.label(mask)
 found=[]
 for idx,sl in enumerate(ndimage.find_objects(labels),1):
  if sl is None:continue
  ys,xs=sl;w=xs.stop-xs.start;h=ys.stop-ys.start;area=np.sum(labels[sl]==idx)
  if not (lo[0]<=w<=hi[0] and lo[1]<=h<=hi[1] and lo[2]<=area<=hi[2]): continue
  found.append(((xs.start+xs.stop-1)/2,(ys.start+ys.stop-1)/2,'separate_complete'))
 return sorted(found)
def identify(mask,panel,grade):
 bounds={'171k-80':((8,8,70),(14,14,135)), '144k-88':((8,10,70),(14,16,135)),
         '49k-75':((10,10,70),(16,16,135)), '24k-80':((8,8,63),(14,14,115)), '48k-80':((10,10,62),(16,15,110))}
 if panel in ('B','D'):
  return detect_separate(mask,*bounds[grade])
 # In the two crowded panels the touching filled squares form a solid band.
 # Only its two terminal squares have an observable left or right boundary.
 active=np.flatnonzero(mask.any(axis=0))
 if not len(active):return []
 if grade=='171k-80':
  xs=(active[0]+5,active[-1]-5)
  method='terminal_square'
 else:
  counts=mask.sum(axis=0)
  minimum_height={'144k-88':11,'49k-75':13,'24k-80':10,'48k-80':11}[grade]
  peaks,props=find_peaks(counts,prominence=2,distance=6,plateau_size=(1,5))
  xs=[(l+r)/2 for i,l,r in zip(peaks,props['left_edges'],props['right_edges']) if counts[i]>=minimum_height]
  method='resolved_outline_peak'
 coords=[]
 for x in xs:
  ys=np.flatnonzero(mask[:,round(x)])
  if len(ys)==0:continue
  coords.append((float(x),float((ys[0]+ys[-1])/2),method))
 return coords
def run(pdf,audit_image=None):
 im=extract_image(pdf);rgb=np.array(im,dtype='int16');marked=im.copy();draw=ImageDraw.Draw(marked)
 observations=[];unresolved=[]
 for panel,(wt,roi,xticks) in PANEL.items():
  x0,x1,y0,y1=roi;ym,yb=line(Y_TOP if panel in 'AB' else Y_BOTTOM);xm,xb=line(xticks)
  for grade,color in COL.items():
   mask=np.max(np.abs(rgb[y0:y1+1,x0:x1+1]-color),axis=2)<9
   for lx,ly,method in identify(mask,panel,grade):
    x=lx+x0;y=ly+y0;rate=10**(xm*x+xb);visc=1000*10**(ym*y+yb);mn,mw,dd,dh=META.get(grade,('','','',''))
    row=dict(study_id='Parisi2022',doi='10.1122/8.0000435',source_figure=4,panel=panel,grade_id=grade,mn_g_mol=mn,mw_g_mol=mw,dispersity=dd,dh_mol_pct=dh,temperature_K=298.15,concentration_wt_pct=wt,shear_rate_s_inv=round(rate,5),viscosity_cP=round(visc,5),x_px=x,y_px=y,pixel_bound_px=2,extraction_method=method)
    (observations if grade in META else unresolved).append(row)
    draw.ellipse((x-2,y-2,x+2,y+2),fill='black')
 for rows,name in [(observations,'parisi_water.csv'),(unresolved,'parisi_unresolved.csv')]:
  rows.sort(key=lambda r:(r['panel'],r['grade_id'],r['x_px']))
  with (D/name).open('w',newline='') as h:
   w=csv.DictWriter(h,rows[0].keys());w.writeheader();w.writerows(rows)
 if audit_image:marked.save(audit_image)
 return observations,unresolved
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('pdf',type=Path);parser.add_argument('--audit-image',type=Path);args=parser.parse_args()
 print('Joined / unresolved observations:',*[len(r) for r in run(args.pdf,args.audit_image)])
