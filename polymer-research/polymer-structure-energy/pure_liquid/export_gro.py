"""Export one verified final TIP3P state as a GROMACS coordinate file."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parent


def export(directory, seed, destination):
    doc = json.loads((directory / "water_npt.json").read_text("utf-8"))
    runs = {row["seed"]: row for row in doc["replicates"]}
    if seed not in runs:
        raise ValueError("Kayıtlı bir su koşusunun tohumunu seçin.")
    row = runs[seed]
    raw = (directory / row["final_state"]).read_bytes()
    if hashlib.sha256(raw).hexdigest() != row["final_state_sha256"]:
        raise ValueError("Kaynak koordinatlarının dosya özeti uyuşmuyor.")
    state = ET.fromstring(gzip.decompress(raw))
    positions = state.findall("./Positions/Position")
    velocities = state.findall("./Velocities/Velocity")
    if len(positions) != len(velocities) or len(positions) != 3 * doc["protocol"]["n_water"]:
        raise ValueError("TIP3P atom sayısı veya hız kaydı eksik.")
    vectors = [tuple(float(node.attrib[a]) for a in "xyz") for node in state.find("PeriodicBoxVectors")]
    if len(vectors) != 3 or any(abs(v[a]) > 1e-10 for i,v in enumerate(vectors) for a in range(3) if i != a):
        raise ValueError("Bu dışa aktarım dik periyodik su kutularına özeldir.")
    names = ("OW", "HW1", "HW2")  # OpenMM tip3p.pdb atom order is O, H1, H2.
    lines = [f"OpenMM TIP3P water seed {seed}  t= {float(state.attrib['time']):.4f}", str(len(positions))]
    for i,(pt,velocity) in enumerate(zip(positions,velocities)):
        prefix = f"{i//3+1:5d}{'SOL':<5}{names[i%3]:>5}{i+1:5d}"
        xyz = "".join(f"{float(pt.attrib[a]):8.3f}" for a in "xyz")
        speed = "".join(f"{float(velocity.attrib[a]):8.4f}" for a in "xyz")
        lines.append(prefix + xyz + speed)
    lines.append(" ".join(f"{v[i]:.9f}" for i,v in enumerate(vectors)))
    destination.write_text("\n".join(lines)+"\n", encoding="ascii")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--directory", type=Path, default=BASE)
    args = parser.parse_args()
    export(args.directory,args.seed,args.out)
