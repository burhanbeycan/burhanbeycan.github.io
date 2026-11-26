"""Export the actual pair trajectory through the offline canvas viewer (Playwright, Chromium, ffmpeg)."""
import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from playwright.sync_api import sync_playwright


def export(html, movie, poster, run_index=0, fps=5):
    with sync_playwright() as pw, tempfile.TemporaryDirectory() as tmp:
        browser = pw.chromium.launch(executable_path=os.environ.get("CHROMIUM_EXECUTABLE") or shutil.which("chromium") or shutil.which("google-chrome"),
                                     headless=True,args=["--no-sandbox","--disable-gpu","--disable-dev-shm-usage"])
        page = browser.new_page(viewport={"width":1400,"height":830},device_scale_factor=1)
        page.goto(html.resolve().as_uri()+"?mode=movie")
        page.evaluate("(i)=>load(i)",run_index)
        page.wait_for_function("run && run.frames.length>0")
        n=page.evaluate("run.frames.length")
        # Crop on even boundaries for yuv420p. All movie frames show the same field of view.
        height=page.evaluate("Math.ceil((document.querySelector('.grid').getBoundingClientRect().bottom+18)/2)*2")
        page.set_viewport_size({"width":1400,"height":height})
        for i in range(n):
            page.evaluate("(i)=>show(i)",i)
            page.screenshot(path=str(Path(tmp)/f"frame_{i:04d}.png"), animations="disabled")
        pick=page.evaluate("run.events.length ? run.frames.findIndex(f=>f.t===run.events[0].t) : run.frames.reduce((best,f,i,a)=>f.s>a[best].s?i:best,0)")
        if pick>=0:
            page.evaluate("(i)=>show(i)",pick)
        page.screenshot(path=str(poster), animations="disabled")
        browser.close()
        subprocess.run(["ffmpeg","-hide_banner","-loglevel","error","-y","-framerate",str(fps),"-i",str(Path(tmp)/"frame_%04d.png"),
                        "-c:v","libx264","-threads","1","-vf","fps=25","-pix_fmt","yuv420p","-crf","20","-movflags","+faststart",str(movie)],check=True)
        return n,height,pick


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("html","mp4","poster"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--run",type=int,default=0)
    args=p.parse_args()
    print(export(args.html,args.mp4,args.poster,args.run))
