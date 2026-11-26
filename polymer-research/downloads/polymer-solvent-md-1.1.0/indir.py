"""Download the complete Polymer–Solvent MD 1.1.0 ZIP using Python 3."""
import argparse
import hashlib
import json
from pathlib import Path
import os
import tempfile
from urllib.request import urlopen

BASE_URL = "https://raw.githubusercontent.com/burhanbeycan/burhanbeycan.github.io/polymer-research-review-v0.2.0/polymer-research/downloads/polymer-solvent-md-1.1.0/"
FILENAME = 'polymer-solvent-md-1.1.0.zip'
SHA256 = 'f82ec6f044de90173d88c15af95184072af2fe956e2a2018efd4b68badb8c6e6'
PARTS = [
  {
    "name": "polymer-solvent-md-1.1.0.zip.001",
    "size": 6000000,
    "sha256": "6315e898d8746e91b8cacb511cd39e1f83c3a201c28594e7975065aaaff2efad"
  },
  {
    "name": "polymer-solvent-md-1.1.0.zip.002",
    "size": 6000000,
    "sha256": "de48b12d894168d806601192f6b70dd7c766ea244e66c4463ab0e13f027e6cf8"
  },
  {
    "name": "polymer-solvent-md-1.1.0.zip.003",
    "size": 6000000,
    "sha256": "c25fd8547ee9a365cbffd8f77d13d13c35e951f6c261f408c735302a4c4550c6"
  },
  {
    "name": "polymer-solvent-md-1.1.0.zip.004",
    "size": 6000000,
    "sha256": "fe05332bb15605db4f3bfc720baa1137b0fa1c4401670fbce2bc461404f78dda"
  },
  {
    "name": "polymer-solvent-md-1.1.0.zip.005",
    "size": 6000000,
    "sha256": "7863c152e0af7b16a6d536c9dafbc6c78365e3ba55604b1384195a4d5ccec588"
  },
  {
    "name": "polymer-solvent-md-1.1.0.zip.006",
    "size": 3394420,
    "sha256": "03ef20d4b2c7165a1a936f87adbebd0dc03184f774f23ec263ed4f0e9ac01998"
  }
]


def verify(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(65536), b""):
            result.update(chunk)
    return result.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=Path(FILENAME))
    ap.add_argument("--parts", type=Path, default=Path(__file__).resolve().parent,
                    help="folder with previously downloaded numbered parts")
    args = ap.parse_args()
    dest = args.output
    if dest.is_file():
        if verify(dest) != SHA256:
            raise SystemExit("Existing destination file differs. Choose another --output.")
        print("Verified", dest)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp = None
    try:
        with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as target:
            temp = Path(target.name)
            for part in PARTS:
                source_path = args.parts / part["name"]
                print("Reading", part["name"], flush=True)
                digest = hashlib.sha256()
                count = 0
                with (source_path.open("rb") if source_path.is_file() else
                      urlopen(BASE_URL + part["name"], timeout=60)) as source:
                    for chunk in iter(lambda: source.read(65536), b""):
                        digest.update(chunk)
                        target.write(chunk)
                        count += len(chunk)
                if count != part["size"] or digest.hexdigest() != part["sha256"]:
                    raise ValueError("Incorrect bytes in " + part["name"])
        if verify(temp) != SHA256:
            raise ValueError("Complete ZIP hash mismatch")
        os.replace(temp, dest)
        print("Verified", dest, SHA256)
    finally:
        if temp and temp.exists():
            temp.unlink()


if __name__ == "__main__":
    main()
