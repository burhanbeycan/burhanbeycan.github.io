"""Project an existing single-chain water trajectory into chemically annotated frames."""
import json
from pathlib import Path

import numpy as np

from trajectory_analysis import (
    HBOND_MIN_ANGLE_DEG, NEIGHBOR_RADIUS_NM, Selections, _snapshot, box_lengths,
    file_hash, minimum_image, whole_chain,
)


def hydrogen_links(xyz, box, selection, snapshot):
    """Return donor O, donor H, acceptor O as trajectory (zero-based) atom indices."""
    links = []
    threshold = np.cos(np.deg2rad(HBOND_MIN_ANGLE_DEG))
    for atom in snapshot["atoms"]:
        if atom["group"] != "water" or atom["element"] != "O":
            continue
        donor = atom["topology_atom_index"]
        hydrogens = selection.water_h[donor]
        oh = minimum_image(xyz[hydrogens] - xyz[donor], box)
        for acceptor in selection.carbonyl:
            oa = minimum_image(xyz[acceptor] - xyz[donor], box)
            if oa @ oa > NEIGHBOR_RADIUS_NM**2:
                continue
            ha = oa - oh
            denom = np.linalg.norm(oh,axis=1) * np.linalg.norm(ha,axis=1)
            if np.any(denom <= 0):
                raise ValueError("Su–karbonil açısı için çakışan atomlar.")
            angles = np.einsum("ij,ij->i", -oh, ha) / denom
            for h, accepted in zip(hydrogens, angles <= threshold):
                if accepted:
                    links.append([donor, int(h), int(acceptor)])
    return links


def trajectory_frames(folder, record, seed, limit=160):
    folder = Path(folder)
    if not 2 <= limit <= 250 or seed not in record["method"]["seeds"]:
        raise ValueError("Kayıtlı bir MD koşusu seçin.")
    topology, trajectory = "run/topology.json", f"run/seed_{seed}.npz"
    for name in (topology, trajectory):
        if file_hash(folder/name) != record["files_sha256"][name]:
            raise ValueError("Yörünge veya topoloji içerik özeti uyuşmuyor.")
    selection = Selections(json.loads((folder/topology).read_bytes()))
    frames = []
    with np.load(folder/trajectory, allow_pickle=False) as data:
        times = data["time_ps"]
        boxes = box_lengths(data["box_nm"])
        positions = data["positions_nm"]
        if len(times) != len(positions) or len(times) != len(boxes) or not len(times):
            raise ValueError("Yörünge dizileri uyumsuz.")
        chosen = np.linspace(0,len(times)-1,num=min(limit,len(times)),dtype=int)
        for index in chosen:
            xyz = positions[index].astype(float)
            box = boxes[index]
            whole = whole_chain(xyz[selection.sol], box, selection)
            snapshot = _snapshot(xyz, whole, box, selection, seed, times[index])
            snapshot["frame_selection"] = "Saved production frame. Coordinates translated to the solute centre of mass."
            snapshot["hydrogen_links"] = hydrogen_links(xyz, box, selection, snapshot)
            frames.append(snapshot)
    return {"study_id": record["study_id"], "seed":seed, "frames":frames, "source_frames":len(times),
            "trajectory_sha256":record["files_sha256"][trajectory],
            "selection_cutoff_nm":NEIGHBOR_RADIUS_NM,"min_angle_deg":HBOND_MIN_ANGLE_DEG}
