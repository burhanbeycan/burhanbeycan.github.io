"""Reproduce validated B and C panel subsets of the Parisi DMSO Figure 5."""
from pathlib import Path
import argparse,csv,hashlib
from io import BytesIO
import numpy as np,pymupdf
from PIL import Image,ImageDraw
from scipy.signal import find_peaks
from digitize import line, detect_separate, META, COL
D=Path(__file__).resolve().parent
REF_HASH='18ab6345c8c92abee374396c0fefb0336bd98284b8ceeac8f31c2bb6e2b1ecc4'
def extract_image(pdf):
 with pymupdf.open(pdf) as doc:
  for page in doc:
   for entry in page.get_images():
    if entry[2:4]!=(1000,818):continue
    im=Image.open(BytesIO(doc.extract_image(entry[0])['image'])).convert('RGB')
    if hashlib.sha256(im.tobytes()).hexdigest()==REF_HASH:return im
 raise ValueError('Cannot locate the published Figure 5.')
def run(pdf,audit=None):
 im=extract_image(pdf);rgb=np.array(im,dtype='int16');marked=im.copy();draw=ImageDraw.Draw(marked);rows=[]
 for panel,wt,xticks,yticks in [
    ('B',15,[(0,577),(1,704.4),(2,831.7),(3,959)],[(3,34.3),(2,116.8),(1,199.3),(0,281.7),(-1,364.2)]),
    ('C',10,[(-1,136.1),(0,228.2),(1,320.2),(2,412.3),(3,504.3)],[(1,411),(0,521.1),(-1,631.1),(-2,741.1)])]:
  xm,xb=line(xticks);ym,yb=line(yticks)
  # B purple is annotated 16.4% in the figure but 15% in the caption. Exclude.
  for grade in (['171k-80','144k-88'] if panel=='B' else ['171k-80','144k-88','49k-75','48k-88']):
   color=COL['48k-80' if grade=='48k-88' else grade]
   if panel=='B':x0,x1,y0,y1=556,969,150,208
   elif grade in ('171k-80','144k-88'):x0,x1,y0,y1=112,525,448,530
   elif grade=='49k-75':x0,x1,y0,y1=300,520,630,652
   else:x0,x1,y0,y1=300,520,610,631
   mask=np.max(np.abs(rgb[y0:y1+1,x0:x1+1]-color),axis=2)<31
   if grade in ('171k-80','144k-88'):
    # Require a nearly complete bounding rectangle. Partially hidden red outlines are omitted.
    coords=detect_separate(mask,(9,9,52),(15,16,120))
   else:
    # Each open diamond/triangle has an isolated lower tip even where outlines touch laterally.
    edge=np.max(np.where(mask,np.arange(mask.shape[0])[:,None],-1),axis=0)
    peaks,pr=find_peaks(edge,plateau_size=(1,5),prominence=3,distance=7);coords=[]
    for k,l,r in zip(peaks,pr['left_edges'],pr['right_edges']):
     xx=(l+r)/2;ys=np.flatnonzero(mask[:,round(xx)]);h=ys[-1]-ys[0]+1
     if 10<=h<=18:coords.append((xx,(ys[-1]+ys[0])/2,'resolved_lower_tip'))
   for lx,ly,how in coords:
    x=float(x0+lx);y=float(y0+ly);mn,mw,dd,dh=META.get(grade,(30000,48100,1.61,88))
    rows.append(dict(study_id='Parisi2022',doi='10.1122/8.0000435',source_figure=5,panel=panel,grade_id=grade,mn_g_mol=mn,mw_g_mol=mw,dispersity=dd,dh_mol_pct=dh,temperature_K=298.15,concentration_wt_pct=wt,shear_rate_s_inv=round(10**(xm*x+xb),5),viscosity_cP=round(1000*10**(ym*y+yb),5),x_px=x,y_px=y,pixel_bound_px=2,extraction_method=how))
    draw.ellipse((x-2,y-2,x+2,y+2),fill='black')
 rows.sort(key=lambda r:(r['panel'],r['grade_id'],r['x_px']))
 with (D/'parisi_dmso.csv').open('w',newline='') as h:
  w=csv.DictWriter(h,rows[0].keys());w.writeheader();w.writerows(rows)
 if audit:marked.save(audit)
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('pdf',type=Path);p.add_argument('--audit-image',type=Path);args=p.parse_args()
 print('Joined DMSO observations:',len(run(args.pdf,args.audit_image)))
