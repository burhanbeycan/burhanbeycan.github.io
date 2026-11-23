"""Oxygen-resolved coordination and continuous contact survival in PEO/LiTFSI."""
from __future__ import annotations

import json
from collections import defaultdict
from numbers import Integral

import numpy as np

from .model import finite


def object_record(value, label):
    if isinstance(value, str):
        if len(value) > 1_000_000:
            raise ValueError(f"{label} JSON exceeds 1 MB.")
        try:
            value = json.loads(value)
        except ValueError as exc:
            raise ValueError(f"{label} must be valid JSON.") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object.")
    return value


def integer(value, name, low, high):
    if isinstance(value, bool) or not isinstance(value, Integral) or not low <= value <= high:
        raise ValueError(f"{name} must be an integer from {low} to {high}.")
    return int(value)


def text_field(record, key):
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Provide a nonempty {key}.")
    return value.strip()


def positive(value, name):
    return finite(value, name, lower=0, strict=True)


def pairset(values, n_li, n_oxygen):
    if not isinstance(values, list):
        raise ValueError("List each Li–oxygen atom pair for every frame.")
    result = set()
    for pair in values:
        if not isinstance(pair, list) or len(pair) != 2:
            raise ValueError("Contact pairs use [Li index, oxygen atom index].")
        key = (integer(pair[0], "Li index", 0, n_li-1),
               integer(pair[1], "oxygen index", 0, n_oxygen-1))
        if key in result:
            raise ValueError("Duplicate pair in the same frame.")
        result.add(key)
    return result


def survival(frames, max_lag):
    """Origins having contact pairs contribute to the denominator if the entire lag is observed."""
    kept = np.zeros(max_lag + 1)
    starts = np.zeros(max_lag + 1)
    for origin, initial in enumerate(frames):
        common = initial.copy()
        last_lag = min(max_lag, len(frames)-1-origin)
        starts[:last_lag+1] += len(initial)
        for lag in range(last_lag+1):
            common.intersection_update(frames[origin+lag])
            kept[lag] += len(common)
            if not common:
                break
    return [float(k/n) if n else None for k, n in zip(kept, starts)]


def shell_runs(record):
    record = object_record(record, "Contact record")
    state = text_field(record, "state_id")
    positive(record.get("temperature_K"), "temperature_K")
    n_li = integer(record.get("n_li"), "n_li", 1, 10**5)
    n_ether = integer(record.get("n_ether_o"), "n_ether_o", 1, 10**6)
    dt = positive(record.get("dt_ps"), "dt_ps")
    cutoff = positive(record.get("cutoff_nm"), "cutoff_nm")
    lag = integer(record.get("max_lag"), "max_lag", 1, 64)
    runs = record.get("runs")
    if not isinstance(runs, list) or not 1 <= len(runs) <= 64:
        raise ValueError("Provide 1–64 independently prepared runs.")
    seen = set()
    model_info = {}
    out = []
    for run in runs:
        if not isinstance(run, dict) or run.get("state_id") != state:
            raise ValueError("All runs must identify the same material and measurement state.")
        model, replica = text_field(run, "model"), text_field(run, "replica")
        source = text_field(run, "force_field_source")
        text_field(run, "trajectory_source")
        if model in model_info and model_info[model] != source:
            raise ValueError("A model label has conflicting force-field sources.")
        model_info[model] = source
        if (model, replica) in seen:
            raise ValueError("Duplicate model and replica ID.")
        seen.add((model, replica))
        raw_frames = run.get("frames")
        if not isinstance(raw_frames, list) or not 2*lag+1 <= len(raw_frames) <= 10000:
            raise ValueError("Supply at least 2×max_lag+1 complete frames and at most 10000.")
        if any(not isinstance(f, dict) or set(f) != {"ether", "tfsi"} for f in raw_frames):
            raise ValueError("Explicit ether and tfsi pair lists are required in each frame, including empty frames.")
        for kind, count in (("ether", n_ether), ("tfsi", 4*n_li)):
            frames = [pairset(f[kind], n_li, count) for f in raw_frames]
            cn = sum(map(len, frames))/(len(frames)*n_li)
            curve = survival(frames, lag)
            tau = (dt * (sum(curve) - (curve[0]+curve[-1])/2)
                   if all(s is not None for s in curve) else None)
            out.append(dict(model=model, replica=replica, area_ps=tau, curve=curve,
                            cn=cn, species=kind, frames=len(frames)))
    return out, dt, cutoff, n_ether/n_li


def table(title, columns, rows):
    return dict(title=title, columns=columns, rows=rows)


def plot(title, x_label, y_label, series):
    return dict(title=title, x_label=x_label, y_label=y_label,
                x_scale="linear", y_scale="linear", series=series)


def result(tables, plots=(), notes=(), warnings=(), metrics=()):
    return dict(tables=tables, plots=list(plots), notes=list(notes), warnings=list(warnings),
                metrics=list(metrics), handoffs=[], material_updates={})


def metric(key, label, value, unit):
    return dict(key=key, label=label, value=value, unit=unit)


def aggregate(values):
    return float(np.mean(values)), float(np.std(values, ddof=1)) if len(values)>1 else None


def summarize_shell(record):
    out, dt, cutoff, eo_li = shell_runs(record)
    groups = defaultdict(list)
    rows, series = [], []
    for r in out:
        groups[r["model"], r["species"]].append(r)
        rows.append([r["model"], r["replica"], r["species"], r["frames"], r["cn"], r["area_ps"]])
        present = [(i*dt, s) for i, s in enumerate(r["curve"]) if s is not None]
        if present:
            series.append(dict(name=f'{r["model"]} / {r["replica"]} / {r["species"]}',
                               x=[x for x,_ in present], y=[y for _,y in present]))
    summary = []
    for (model, species), runs in groups.items():
        mean, sd = aggregate([r["cn"] for r in runs])
        areas = [r["area_ps"] for r in runs if r["area_ps"] is not None]
        area, area_sd = aggregate(areas) if areas else (None, None)
        summary.append([model, species, len(runs), mean, sd, len(areas), area, area_sd])
    label = "Continuous survival (truncated integral)"
    return result([
        table("Independent-run summaries", ["Force field","O donor","Runs","Mean O/Li","Run SD","Area runs","Mean area (ps)","Area SD (ps)"], summary),
        table(label, ["Force field","Replica","O donor","Frames","O/Li","Area (ps)"], rows)],
        [plot("Continuous survival by oxygen atom", "Lag (ps)", "Surviving fraction", series)] if series else [],
        [f"EO:Li = {eo_li:g}; distance cutoff = {cutoff:g} nm; sampling interval = {dt:g} ps.",
         "A broken contact ends continuous survival even if the same oxygen returns later. "
         "The finite-lag area measures residual persistence over the chosen window; it is not a fitted full residence time.",
         "Means use equal weights per independently started trajectory. SD describes variation between runs. "
         "A single run has no between-run uncertainty estimate.",
         "Each TFSI has four sulfonyl oxygens. Multiple oxygens on one anion count separately."])


def orthorhombic_contacts(li_nm, ether_o_nm, tfsi_o_nm, box_nm, cutoff_nm):
    """Time series of contacts from nm coordinates and orthorhombic box lengths."""
    li, ether, tfsi, box = [np.asarray(x,dtype=float) for x in (li_nm,ether_o_nm,tfsi_o_nm,box_nm)]
    cutoff = positive(cutoff_nm,"cutoff_nm")
    if any(a.ndim!=3 or a.shape[2]!=3 or not a.shape[1] for a in (li,ether,tfsi)):
        raise ValueError("Use atom coordinates shaped (frames, atom count, 3).")
    if not (len(li)==len(ether)==len(tfsi) and len(li)>0 and tfsi.shape[1]==4*li.shape[1]):
        raise ValueError("Match frame counts and supply four sulfonyl O per Li for neutral LiTFSI.")
    if box.shape==(3,):
        box = np.broadcast_to(box,(len(li),3))
    if box.shape!=(len(li),3) or any(not np.isfinite(a).all() for a in (li,ether,tfsi,box)) or np.min(box)<=2*cutoff:
        raise ValueError("Supply positive orthorhombic box lengths exceeding twice the cutoff in every frame.")
    frames=[]
    for i, lengths in enumerate(box):
        frame={}
        for name, oxygens in (("ether",ether[i]),("tfsi",tfsi[i])):
            pairs=[]
            for ion, xyz in enumerate(li[i]):
                d = oxygens-xyz
                d -= lengths*np.rint(d/lengths)
                selected=np.flatnonzero(np.sum(d*d,axis=1)<=cutoff**2)
                pairs.extend([ion,int(o)] for o in selected)
            frame[name]=pairs
        frames.append(frame)
    return frames
