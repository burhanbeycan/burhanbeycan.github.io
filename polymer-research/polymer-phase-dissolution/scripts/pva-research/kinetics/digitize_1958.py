"""Recover Figure 1 measurements from manually located marker centers in the primary scan.

For optional reproduction/inspection of the source crop:
python digitize_1958.py --pdf article.pdf --overlay /tmp/markers.png
PyMuPDF and Pillow are needed only for the overlay. This script otherwise uses the standard library.
"""
from pathlib import Path
import argparse
import csv

ROOT = Path(__file__).resolve().parent

# Page 1, PyMuPDF crop Rect(278,449,454,643), render scale 5; image 880 x 970 px.
# Axis 0,10,20,30 h ticks near x=77,269,458,648 respectively.
# At the left axis, 50,70,80,90,100,110 P tick centers are y=696,503,407,312,215,120.
# Small scan rotation lowers page-y by ~4 px along 756 px of the horizontal axis.
TICKS = {"x0": 77.0, "x30": 648.0, "y50_at_x0": 696.0, "y110_at_x0": 120.0,
         "y_slope_with_x": -4 / 756}
POINTS = [
    ("L1490", 1490, 16.04, "open", [(77,601),(125,585),(169,577),(233,549),
        (278,535),(366,500),(495,349),(527,249),(548,123)]),
    ("L1705", 1705, 14.63, "filled", [(77,653),(121,641),(145,643),(191,624),
        (277,611),(339,597),(479,560),(637,504),(697,444),(741,367)])
]

def measured_rows():
    for sample, dp, c, symbol, points in POINTS:
        for x,y in points:
            dx = x - TICKS["x0"]
            time_h = 30 * dx / (TICKS["x30"] - TICKS["x0"])
            y_left = y - dx * TICKS["y_slope_with_x"]
            eta_poise = 50 + (TICKS["y50_at_x0"] - y_left) * 60 / (TICKS["y50_at_x0"] - TICKS["y110_at_x0"])
            yield dict(material_id=sample, reported_dp_index=dp, wt_percent=c,
                       storage_T_K=303.15, shear_rate_s_inv="", age_h=round(time_h,2),
                       apparent_eta_cP=round(eta_poise*100,1), marker=symbol,
                       scan_center_x=x,scan_center_y=y)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf",type=Path)
    parser.add_argument("--overlay",type=Path)
    args=parser.parse_args()
    rows=list(measured_rows())
    with (ROOT/"aging_fig1.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    if args.overlay:
        if not args.pdf: parser.error("--overlay requires --pdf")
        import fitz
        from PIL import Image, ImageDraw
        with fitz.open(args.pdf) as doc:
            pm=doc[0].get_pixmap(matrix=fitz.Matrix(5,5),clip=fitz.Rect(278,449,454,643),alpha=False)
            im=Image.frombytes("RGB",[pm.width,pm.height],pm.samples)
        draw=ImageDraw.Draw(im)
        for r in rows:
            x,y=r["scan_center_x"],r["scan_center_y"]
            draw.ellipse((x-10,y-10,x+10,y+10),outline=(185,28,0),width=2)
            draw.line((x-2,y,x+2,y),fill=(25,115,210))
        im.save(args.overlay)
if __name__=="__main__": main()
