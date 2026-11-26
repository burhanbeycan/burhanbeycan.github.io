"""Independent marker- and axis-level checks, without using pixel masks from digitize.py."""
from collections import Counter
import csv, math,sys
from pathlib import Path
import numpy as np
from digitize import COL, extract_image
D=Path(__file__).resolve().parent
rows=list(csv.DictReader((D/'parisi_water.csv').open()))
extra=list(csv.DictReader((D/'parisi_unresolved.csv').open()))
assert len(rows)==283 and len(extra)==79
assert len({(q['panel'],q['grade_id'],q['x_px']) for q in rows+extra})==len(rows+extra)
assert all(q['grade_id']!='48k-80' and q['mn_g_mol'] for q in rows)
assert all(q['grade_id']=='48k-80' and not q['mn_g_mol'] for q in extra)
assert len(set((q['panel'],q['grade_id']) for q in rows))==15
im=np.asarray(extract_image(Path(sys.argv[1])),dtype=int)
# Separate reading of labelled major ticks from the original pixel grid.
# All selected anchors must meet a visible dark tick extending from the axis line.
for x,y in [(199,37),(276,37),(352,37),(428,37),(504,37),(705,37),(958,37),
            (113,756),(213,756),(312,756),(411,756),(510,756),(569,756),(755,756),(941,756)]:
 assert im[y-1:y+2,x-1:x+2].mean(axis=2).min()<160,(x,y)
for y in [33,100,167,234,301,368,416,484,551,618,685,752]:
 assert im[y-1:y+2,79:84].mean(axis=2).min()<180,y
# Verification method independent from a marker's original reconstruction:
# rectangular windows, independently broadened RGB-class threshold and the observed
# top and bottom edge of the centre-column pixels (including antialiasing).
max_vertical_offset=0
for q in rows:
 x,y=(round(float(q[k])) for k in ['x_px','y_px'])
 patch=im[y-9:y+10,x,:];color=np.asarray(COL[q['grade_id']])
 ids=np.flatnonzero(np.max(abs(patch-color),axis=1)<32)
 assert len(ids)>=6,q
 centre=y-9+(ids[0]+ids[-1])/2
 max_vertical_offset=max(max_vertical_offset,abs(float(q['y_px'])-centre))
assert max_vertical_offset<=1.5,max_vertical_offset
print('15 named grade/condition groups, 283 joined, 79 unjoined, all independent tick checks pass.')
print('Maximum centre-column displacement under wider colour threshold (pixels):',max_vertical_offset)
print(Counter((q['panel'],q['grade_id']) for q in rows))
from digitize_dmso import extract_image as extract_dmso
dmso=list(csv.DictReader((D/'parisi_dmso.csv').open())); im=np.array(extract_dmso(Path(sys.argv[1])),dtype=int)
assert len(dmso)==101 and len({(r['panel'],r['grade_id']) for r in dmso})==6
assert {r['panel'] for r in dmso}=={'B','C'}
assert all(r['grade_id'] in ('144k-88','171k-80') for r in dmso if r['panel']=='B')
assert len({(r['panel'],r['grade_id'],r['x_px']) for r in dmso})==len(dmso)
for x,y in [(577,369),(704,369),(832,369),(959,369),(136,745),(228,745),(320,745),(412,745),(504,745)]:
 assert im[y-1:y+2,x-1:x+2].mean(axis=2).min()<170,(x,y)
for y in (34,117,199,282,364,411,521,631,741):
 assert im[y-1:y+2,94:100].mean(axis=2).min()<180,y
offset=0
for q in dmso:
 x,y=(round(float(q[k])) for k in ['x_px','y_px'])
 col=COL['48k-80' if q['grade_id']=='48k-88' else q['grade_id']]
 # Check the bounding-box centre of a wider colour mask in the local marker window.
 yy,xx=np.where(np.max(np.abs(im[y-9:y+10,x-6:x+7]-col),2)<42)
 assert len(xx)>15,q
 mid=(max(yy)+min(yy))/2+y-9
 offset=max(offset,abs(mid-float(q['y_px'])))
assert offset<=2.0,offset
print('6 DMSO grade/condition groups, 101 joined; max alternate centre offset (pixels):',offset)
