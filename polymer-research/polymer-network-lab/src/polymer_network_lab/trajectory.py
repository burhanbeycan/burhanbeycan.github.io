"""Read an orthorhombic GROMACS multi-frame GRO in atom order and analyze NH contacts."""

import argparse
import json
from pathlib import Path
import re

from .ldi import contact_analysis


def gro_frames(stream, spacing_ps=None):
    if spacing_ps is not None and spacing_ps <= 0:
        raise ValueError("Frame spacing must be positive.")
    frame_index=0
    atom_records=None
    while True:
        title = stream.readline()
        if not title:
            return
        if not title.strip():
            continue
        match = re.search(r"(?:^|\s)t\s*=\s*([-+.\dEe]+)", title)
        if not match and spacing_ps is None:
            raise ValueError("GRO frame lacks t=... . Supply --spacing-ps for a regularly sampled trajectory.")
        t = float(match.group(1)) if match else frame_index*spacing_ps
        try:
            n = int(stream.readline())
            if not 3 <= n <= 250000:
                raise ValueError()
            xyz = []
            identities=[]
            for _ in range(n):
                line = stream.readline()
                identities.append(line[:20])
                # Coordinate width follows GROMACS' n+5 format, n being precision.
                periods=[i for i,c in enumerate(line[20:]) if c=="."]
                if len(periods)<3 or periods[1]-periods[0] != periods[2]-periods[1]:
                    raise ValueError()
                width=periods[1]-periods[0]
                xyz.append([float(line[20+j*width:20+(j+1)*width]) for j in range(3)])
            box = [float(v) for v in stream.readline().split()]
        except (TypeError, ValueError) as exc:
            raise ValueError("Malformed GRO frame or unsupported coordinate precision.") from exc
        if atom_records is None:
            atom_records=identities
        elif identities!=atom_records:
            raise ValueError("GRO atom identifiers or ordering changed between frames.")
        if len(box)==9 and all(v==0 for v in box[3:]):
            box=box[:3]
        if len(box)!=3:
            raise ValueError("Convert the trajectory to an orthorhombic representation or use a triclinic-aware analyzer.")
        yield {"time_ps":t, "xyz_nm":xyz, "box_nm":box}
        frame_index+=1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("topology_json",type=Path,help="Zero-based donor/acceptor and chain index JSON.")
    parser.add_argument("trajectory_gro",type=Path,help="Concatenated frames in system atom order and nm.")
    parser.add_argument("--spacing-ps",type=float,help="Used only if timestamps are absent.")
    parser.add_argument("--cutoff-nm",type=float,default=.35)
    parser.add_argument("--angle-deg",type=float,default=150.)
    args=parser.parse_args()
    with args.topology_json.open(encoding="utf-8") as src:
        meta=json.load(src)
    with args.trajectory_gro.open(encoding="utf-8") as coords:
        out=contact_analysis(meta,gro_frames(coords,args.spacing_ps),args.cutoff_nm,args.angle_deg)
    print(json.dumps(out,allow_nan=False,indent=2))


if __name__=="__main__":
    main()
