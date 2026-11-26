"""Render a rigid rotation of one sampled theta-reference chain. Requires Pillow, ffmpeg, Node."""
import json
import math
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

root = Path(__file__).resolve().parents[1]
web = root / 'dist' if (root / 'dist').is_dir() else root
js = """global.window=globalThis;require('./chain-data.js');const row=CHAIN_DATA.rows.find(r=>r.id==='IK1500-1');const g=require('./chain-engine.js').gaussian(row.Rg_nm,17,64,384);process.stdout.write(JSON.stringify({row,g}));"""
info = json.loads(subprocess.check_output(['node', '-e', js], cwd=web))
g, sample = info['g'], info['row']
w,h=1024,560
points=g['points']
scale=172/max(math.sqrt(x*x+y*y+z*z) for x,y,z in points)
fonts='/usr/share/fonts/truetype/dejavu/'
def font(size, name='DejaVuSans.ttf'):
    return ImageFont.truetype(fonts+name, size)
def tr(number, digits=1):
    return f'{number:.{digits}f}'.replace('.',',')
title=font(35,'DejaVuSerif.ttf'); body=font(19); small=font(15); number=font(27,'DejaVuSerif.ttf')
fps=18
proc=subprocess.Popen(['ffmpeg','-hide_banner','-loglevel','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{w}x{h}','-r',str(fps),'-i','-','-an','-c:v','libx264','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(web/'chain-rotation.mp4')],stdin=subprocess.PIPE)
for i in range(8*fps):
    im=Image.new('RGB',(w,h),'#f4f7f3'); d=ImageDraw.Draw(im)
    d.rectangle([0,0,w,59],fill='#173e36'); d.text((28,17),'Polymer Research Studio',font=body,fill='white')
    d.text((29,79),'θ koşulunda bir polistiren bobini',font=title,fill='#18392e')
    d.text((31,135),'IK1500-1  ·  siklohekzan  ·  34,5 °C',font=body,fill='#495f53')
    d.rounded_rectangle([20,181,575,491],12,fill='#ffffff',outline='#d0dcd4')
    d.rounded_rectangle([593,181,1001,491],12,fill='#e5eee8')
    ax=2*math.pi*i/(8*fps)
    cs,sn=math.cos(ax),math.sin(ax)
    coords=[]
    for x,y,z in points:
        t=sn*x+cs*z
        coords.append((300+scale*(cs*x-sn*z),321+scale*(-.92*y+.39*t),.39*y+.92*t))
    segments=sorted(range(len(coords)-1), key=lambda j:(coords[j][2]+coords[j+1][2])/2)
    for j in segments:
        a,b=coords[j],coords[j+1]
        mix=j/(len(segments)-1)
        rgb=tuple(round(u+(v-u)*mix) for u,v in zip((75,105,143),(25,127,93)))
        d.line([a[:2],b[:2]],fill=rgb,width=4)
    for j,clr in [(0,'#405e88'),(-1,'#af5b3b')]:
        x,y,_=coords[j]; d.ellipse((x-4,y-4,x+4,y+4),fill=clr)
    length=200*scale
    d.line((41,464,41+length,464),fill='#253f35',width=3);d.text((41,436),'200 nm',font=small,fill='#253f35')
    for y,label,val in [(203,'Saçılma ölçümü · √⟨Rg²⟩z',str(int(sample['Rg_nm']))+' nm'),(293,'Gösterilen tek örneğin Rg değeri',tr(g['sampleRg_nm'])+' nm'),(382,'384 bağımsız örnek · √⟨Rg²⟩',tr(g['ensembleRmsRg_nm'])+' nm')]:
        d.text((616,y),label,font=small,fill='#496759');d.text((616,y+24),val,font=number,fill='#143c30')
    d.text((29,507),'Sabit bir Gauss zinciri kamera karşısında döner. 64 istatistiksel adım; renk zincir sırasını gösterir.',font=small,fill='#334f41')
    d.text((29,530),'Uzunluk ölçeği: Miyaki, Einaga ve Fujita (1978), Tablo I. Video süresi oynatma süresidir.',font=small,fill='#435a4c')
    proc.stdin.write(im.tobytes())
proc.stdin.close()
if proc.wait():
    raise SystemExit('ffmpeg failed')
