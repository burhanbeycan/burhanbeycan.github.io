"""Read experimental marker centers in Hager et al., arXiv:cond-mat/0206464v1, Fig. 2.

Usage: python extract_hager_fig2.py /path/to/hager.pdf /output/directory
Dependencies: PyMuPDF, Python 3. Source: https://arxiv.org/pdf/cond-mat/0206464v1 .
Indices refer to PDF vector paths. Hash pin prevents accidental extraction of another version.
"""
from collections import defaultdict
import csv
from hashlib import sha256
from pathlib import Path
from statistics import mean
import sys
import fitz

def calibrate(positions, values):
    mx, my = mean(positions), mean(values)
    slope = sum((x-mx)*(y-my) for x,y in zip(positions,values))/sum((x-mx)**2 for x in positions)
    return lambda p: my + slope*(p-mx)

def extract(pdf, output_dir):
    if sha256(pdf.read_bytes()).hexdigest()!='2572457e1f46be4a9bfc64856d093b625cebca1a01330002cec82f802abb6d31':
        raise ValueError('Use the v1 PDF for arXiv:cond-mat/0206464')
    doc=fitz.open(pdf)
    paths=doc[10].get_drawings()
    # Page 11 labels and major tick paths agree (x 0–.28 in steps of .04, y 290–310 K in steps of 5).
    xaxis=calibrate([paths[i]['rect'].x0 for i in range(238,253,2)], [i*.04 for i in range(8)])
    yaxis=calibrate([paths[i]['rect'].y0 for i in range(255,264,2)], [290+i*5 for i in range(5)])
    for mw, first, last, doi, name in [(1560000,1,40,'10.1103/PhysRevA.18.2683','fig2_crosscheck.csv'),
                                        (200000,42,79,'10.1063/1.430349','tie_lines_200k_aux.csv')]:
        markers=defaultdict(list)
        for i in range(first,last+1):
            path=paths[i]
            rect=path['rect']
            assert path['type']=='fs' and len(path['items'])==4 and all(s[0]=='c' for s in path['items'])
            assert 3.98 < rect.width < 4.02 and 3.98 < rect.height < 4.02
            y=(rect.y0+rect.y1)/2
            for grouped_y in markers:
                if abs(grouped_y-y)<.002:
                    break
            else:
                grouped_y=y
            markers[grouped_y].append((i,(rect.x0+rect.x1)/2))
        with (output_dir/name).open('w',newline='') as fh:
            writer=csv.writer(fh)
            writer.writerow(['T_K','phi_lean','phi_rich','Mw_g_mol','Mn_g_mol','dispersity',
                             'lower_path_id','upper_path_id','original_doi','figure_source','method'])
            lowcut=(first+last+1)//2
            for y,points in sorted(markers.items(),reverse=True):
                if len(points)!=2: continue  # zoom cannot recover temperatures merged on a vector-plot row
                (i,a),(j,b)=sorted(points,key=lambda p:p[1])
                if not i<lowcut<=j: continue   # require one point from each arm
                writer.writerow([f'{yaxis(y):.3f}',f'{xaxis(a):.5f}',f'{xaxis(b):.5f}', mw,
                                 1515000 if mw==1560000 else '',f'{1560000/1515000:.8f}' if mw==1560000 else '',
                                 i,j,doi,'arXiv:cond-mat/0206464v1 Figure 2','digitized_matched_experimental_symbol_centres'])
if __name__=='__main__':
    extract(Path(sys.argv[1]),Path(sys.argv[2]))
