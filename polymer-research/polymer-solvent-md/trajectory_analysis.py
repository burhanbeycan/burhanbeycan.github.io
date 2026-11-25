"""Reuse explicit-water MD production frames for chain dimensions and hydration.

Distances in nm, time in ps. RDF is a frame-wise bulk-water-normalized
cross radial distribution with exact spherical shell volumes in an orthorhombic
periodic cell. Solute atoms are first made whole using the covalent bond graph.
"""
import argparse
from collections import deque
import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile

import numpy as np


SCHEMA = "polymer-water-observables/1"
BIN_WIDTH_NM = 0.02
NEIGHBOR_RADIUS_NM = 0.35
HBOND_MIN_ANGLE_DEG = 150.0


def file_hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def minimum_image(displacement, lengths):
    return displacement - lengths * np.rint(displacement / lengths)


def box_lengths(box):
    box = np.asarray(box, dtype=float)
    if box.ndim == 3 and box.shape[1:] == (3, 3):
        lengths = np.diagonal(box, axis1=1, axis2=2)
        axes = np.zeros_like(box)
        axes[:, np.arange(3), np.arange(3)] = lengths
        if not np.allclose(box, axes, rtol=0, atol=1e-7):
            raise ValueError("Bu sürüm dik açılı periyodik hücreler için doğrulandı.")
    elif box.ndim == 2 and box.shape[1] == 3:
        lengths = box
    else:
        raise ValueError("Kutu vektörleri (kare, 3, 3) veya uzunlukları (kare, 3) olmalıdır.")
    if not lengths.size or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("Kutu uzunlukları sonlu ve pozitif olmalıdır.")
    return lengths


def _indices(values, n, label):
    a = np.asarray(values)
    if a.ndim != 1 or a.size and (not np.issubdtype(a.dtype, np.integer)
                                  or np.min(a) < 0 or np.max(a) >= n):
        raise ValueError(f"Geçersiz atom seçimi: {label}")
    if len(np.unique(a)) != a.size:
        raise ValueError(f"Yinelenmiş atom indeksi: {label}")
    return a.astype(int)


class Selections:
    def __init__(self, top):
        self.elements = list(top["elements"])
        n = len(self.elements)
        mass_all = np.asarray(top["mass_da"], dtype=float)
        if mass_all.shape != (n,) or not np.isfinite(mass_all).all() or np.any(mass_all < 0):
            raise ValueError("Atom kütleleri geçersiz.")
        self.sol = _indices(top["solute_indices"], n, "zincir")
        self.carbonyl = _indices(top["carbonyl_oxygens"], n, "amid karbonil oksijenleri")
        self.water = _indices(top["water_oxygens"], n, "su oksijenleri")
        if not len(self.sol) or not len(self.water):
            raise ValueError("Zincir ve su seçimleri gereklidir.")
        solset = set(self.sol)
        if not set(self.carbonyl).issubset(solset) or solset.intersection(self.water):
            raise ValueError("Karbonil ve su seçimi polimer topolojisiyle uyuşmuyor.")
        if any(self.elements[i] != "O" for i in (*self.carbonyl, *self.water)):
            raise ValueError("Karbonil ve su seçimlerinde yalnız oksijen bulunabilir.")
        self.mass = mass_all[self.sol]
        self.heavy = np.array([self.elements[i] != "H" for i in self.sol])
        self.local = {i: j for j, i in enumerate(self.sol)}
        self.refs = [self.local[i] for i in self.carbonyl]
        if not np.all(self.mass > 0) or not np.any(self.heavy):
            raise ValueError("Zincirin fiziksel atomları pozitif kütleye sahip olmalıdır.")
        neighbours = [[] for _ in range(n)]
        self.bonds = []
        # Optional explicit bond orders supplied by the simulation producer.
        supplied_orders = top.get("bond_orders")
        if supplied_orders is not None and len(supplied_orders) != len(top["bonds"]):
            raise ValueError("Bağ dereceleri ve bağ listesi farklı uzunluklarda.")
        for k, entry in enumerate(top["bonds"]):
            a, b = map(int, entry[:2])
            if min(a, b) < 0 or max(a, b) >= n or a == b:
                raise ValueError("Topolojide geçersiz bağ.")
            neighbours[a].append(b)
            neighbours[b].append(a)
            if a in solset and b in solset:
                order = supplied_orders[k] if supplied_orders is not None else (entry[2] if len(entry) > 2 else None)
                self.bonds.append((self.local[a], self.local[b], order))
            elif (a in solset) != (b in solset):
                raise ValueError("Çözünen atom seçimi bir kovalent bağı kesiyor.")
        self.water_h = {}
        for i in self.water:
            hs = [j for j in neighbours[i] if self.elements[j] == "H"]
            if len(hs) != 2:
                raise ValueError("Görselleştirme için her su oksijeninin iki H komşusu gereklidir.")
            self.water_h[i] = hs
        # Ordered BFS edges define which already-unwrapped atom anchors its child.
        queue = deque([int(self.sol[0])])
        seen = {queue[0]}
        self.tree = []
        while queue:
            a = queue.popleft()
            for b in neighbours[a]:
                if b in solset and b not in seen:
                    seen.add(b)
                    self.tree.append((self.local[a], self.local[b]))
                    queue.append(b)
        if seen != solset:
            raise ValueError("Jirasyon analizi tek bir bağlı zincir gerektirir.")
        self.edges = np.array([[a, b] for a, b, _ in self.bonds], dtype=int).reshape(-1, 2)
        self.sites = amide_site_identities(self.carbonyl, self.elements, neighbours)


def amide_site_identities(oxygens, elements, neighbours):
    """Annotate the preselected amide carbonyls from the explicit bonded graph."""
    sites = []
    anchors = set()
    for oxygen in oxygens:
        oxygen = int(oxygen)
        attached = neighbours[oxygen]
        carbon = attached[0] if len(attached) == 1 and elements[attached[0]] == "C" else None
        nitrogens = [i for i in neighbours[carbon] if elements[i] == "N"] if carbon is not None else []
        nitrogen = nitrogens[0] if len(nitrogens) == 1 else None
        bonded_h, kind = None, "unclassified"
        if nitrogen is not None:
            ligands = [elements[i] for i in neighbours[nitrogen]]
            bonded_h = ligands.count("H")
            if len(ligands) == 3 and ligands.count("C") + bonded_h == 3:
                kind = {2:"primary_amide", 1:"secondary_amide", 0:"tertiary_amide"}.get(bonded_h, "unclassified")
            if kind == "secondary_amide":
                anchors.add(nitrogen)
        sites.append(dict(site_id=f"O[{oxygen}]",oxygen_atom_index=oxygen,carbonyl_c_atom_index=carbon,
                          amide_n_atom_index=nitrogen,amide_class=kind,n_bonded_hydrogens=bonded_h))
    distances = {i:0 for i in anchors}
    queue = deque(anchors)
    while queue:
        a = queue.popleft()
        for b in neighbours[a]:
            if elements[b] != "H" and b not in distances:
                distances[b] = distances[a] + 1
                queue.append(b)
    for site in sites:
        site["n_bonds_to_nearest_secondary_amide_n"] = distances.get(site["amide_n_atom_index"])
    return sites


def whole_chain(wrapped, lengths, selection):
    xyz = np.empty_like(wrapped, dtype=float)
    xyz[0] = wrapped[0]
    for a, b in selection.tree:
        xyz[b] = xyz[a] + minimum_image(wrapped[b] - wrapped[a], lengths)
    edges = selection.edges
    if len(edges):
        a, b = edges.T
        if not np.allclose(xyz[a] - xyz[b], minimum_image(wrapped[a] - wrapped[b], lengths),
                           atol=1e-5, rtol=0):
            raise ValueError("Kovalent çevrim periyodik hücrede tutarlı biçimde açılamadı.")
    return xyz


def nearest_self_image(xyz, lengths):
    """Minimum heavy-atom separation to *any nonzero periodic image* of same chain."""
    delta = xyz[:, None, :] - xyz[None, :, :]
    translations = np.rint(delta / lengths)
    close = delta - lengths * translations
    sq = np.sum(close**2, axis=-1)
    zero = np.all(translations == 0, axis=-1)
    # When the closest image is the unshifted cell, translate along its cheapest axis.
    sq[zero] += np.min(lengths**2 - 2 * np.abs(close[zero]) * lengths, axis=-1)
    return float(np.sqrt(max(0, sq.min())))


def accepted_water_hbonds(frame, lengths, selection, pair_vectors, pair_distances):
    """Count directional water D-H···O(carbonyl) triplets among close O-O pairs."""
    ref_rows, water_cols = np.where(pair_distances <= NEIGHBOR_RADIUS_NM)
    if not len(ref_rows):
        return np.zeros(len(selection.carbonyl), dtype=int)
    donors = selection.water[water_cols]
    hydrogens = np.asarray([selection.water_h[d] for d in donors])
    donor_to_h = minimum_image(frame[hydrogens] - frame[donors, None, :], lengths)
    donor_to_acceptor = pair_vectors[ref_rows, water_cols]
    h_to_acceptor = donor_to_acceptor[:, None, :] - donor_to_h
    denominators = np.linalg.norm(donor_to_h, axis=-1) * np.linalg.norm(h_to_acceptor, axis=-1)
    if np.any(denominators < 1e-10):
        raise ValueError("Hidrojen bağı açısı için çakışan atomlar bulundu.")
    cos_dha = np.einsum("ijk,ijk->ij", -donor_to_h, h_to_acceptor) / denominators
    counts = np.count_nonzero(cos_dha <= np.cos(np.deg2rad(HBOND_MIN_ANGLE_DEG)),axis=1)
    return np.bincount(ref_rows,weights=counts,minlength=len(selection.carbonyl)).astype(int)


def _snapshot(frame, whole, lengths, selection, seed, time):
    com = np.average(whole, axis=0, weights=selection.mass)
    atoms = []
    bonds = []
    def add_atom(atom_index, pos, group):
        p = pos - com
        atoms.append(dict(element=selection.elements[atom_index], x=float(p[0]), y=float(p[1]),
                          z=float(p[2]), group=group, topology_atom_index=int(atom_index)))
    for i, p in zip(selection.sol, whole):
        add_atom(i, p, "polymer")
    for a, b, order in selection.bonds:
        bonds.append(dict(a=a, b=b, order=order))
    water_count = 0
    refs = whole[selection.refs]
    for o in selection.water:
        if not len(refs):
            break
        displacement = minimum_image(frame[o] - refs, lengths)
        rr = np.einsum("ij,ij->i", displacement, displacement)
        j = np.argmin(rr)
        if rr[j] > NEIGHBOR_RADIUS_NM**2:
            continue
        o_pos = refs[j] + displacement[j]
        o_id = len(atoms)
        add_atom(o, o_pos, "water")
        for h in selection.water_h[o]:
            h_pos = o_pos + minimum_image(frame[h] - frame[o], lengths)
            bonds.append(dict(a=o_id, b=len(atoms), order=1.0))
            add_atom(h, h_pos, "water")
        water_count += 1
    return dict(source_seed=seed, production_time_ps=float(time), unit="nm", selected_water_count=water_count,
                selection_cutoff_nm=NEIGHBOR_RADIUS_NM, atoms=atoms, bonds=bonds,
                water_selection="Distinct waters whose O is within 0.35 nm of any selected carbonyl O; distance criterion only.",
                frame_selection="Final saved production frame of lowest-numbered seed; not an ensemble average.",
                bond_note="null bond order denotes an untyped covalent connection.")


def one_run(source, lengths, sel, edges, seed, want_snapshot):
    with np.load(source, allow_pickle=False) as run:
        time = np.asarray(run["time_ps"], dtype=float)
        positions = np.asarray(run["positions_nm"], dtype=float)
    if positions.shape != (len(lengths), len(sel.elements), 3) or time.shape != (len(lengths),):
        raise ValueError("Zaman, kutu, koordinat ve topoloji boyutları uyumsuz.")
    if len(time) < 4 or not np.isfinite(positions).all() or not np.isfinite(time).all() or np.any(time < 0):
        raise ValueError("En az dört sonlu üretim karesi gereklidir.")
    steps = np.diff(time)
    if np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-5, atol=1e-6):
        raise ValueError("Üretim örnekleri eşit aralıklı ve zaman damgaları artan sırada olmalıdır.")
    rg, box_approaches = [], []
    n_ref, n_water = len(sel.carbonyl), len(sel.water)
    shell = 4 * np.pi * np.diff(edges**3) / 3
    rdf_sum = np.zeros_like(shell)
    pair_counts = np.zeros(len(shell), dtype=np.int64)
    local_035 = []
    hydrogen_bonds = []
    site_neighbors = np.zeros(n_ref, dtype=int)
    site_hbonds = np.zeros(n_ref, dtype=int)
    for xyz, L in zip(positions, lengths):
        whole = whole_chain(xyz[sel.sol], L, sel)
        com = np.average(whole, axis=0, weights=sel.mass)
        rg.append(float(np.sqrt(np.average(((whole - com)**2).sum(axis=1), weights=sel.mass))))
        box_approaches.append(nearest_self_image(whole[sel.heavy], L))
        if n_ref:
            diffs = minimum_image(whole[sel.refs, None, :] - xyz[sel.water][None, :, :], L)
            distances = np.sqrt(np.einsum("ijk,ijk->ij", diffs, diffs))
            per_site_hbonds = accepted_water_hbonds(xyz, L, sel, diffs, distances)
            per_site_neighbors = np.count_nonzero(distances <= NEIGHBOR_RADIUS_NM, axis=1)
            site_hbonds += per_site_hbonds
            site_neighbors += per_site_neighbors
            hydrogen_bonds.append(float(per_site_hbonds.mean()))
            counts = np.histogram(distances, bins=edges)[0]
            pair_counts += counts
            rdf_sum += counts * np.prod(L) / (shell * n_water * n_ref)
            local_035.append(float(per_site_neighbors.mean()))
    blocks = [float(np.mean(values)) for values in np.array_split(rg, 4)]
    out = dict(seed=seed, frames=len(time), first_saved_ps=float(time[0]), last_saved_ps=float(time[-1]),
               duration_ps=float(time[-1]), sample_interval_ps=float(steps[0]), time_ps=time.tolist(), rg_nm=rg,
               rg_mean_nm=float(np.mean(rg)), rg_block_means_nm=blocks,
               min_solute_image_heavy_distance_nm=min(box_approaches),
               min_cube_edge_nm=float(lengths.min()),
               rdf_g=(rdf_sum/len(time)).tolist() if n_ref else [],
               coordination_to_outer_edge=(np.cumsum(pair_counts)/(n_ref*len(time))).tolist() if n_ref else [],
               coordination_mean_at_035_nm=float(np.mean(local_035)) if n_ref else None,
               site_mean_neighbors=(site_neighbors/len(time)).tolist(),
               site_mean_hbonds=(site_hbonds/len(time)).tolist(),
               hydrogen_bonds_per_carbonyl=hydrogen_bonds,
               hydrogen_bond_mean_per_carbonyl=float(np.mean(hydrogen_bonds)) if n_ref else None)
    snap = _snapshot(positions[-1], whole, lengths[-1], sel, seed, time[-1]) if want_snapshot else None
    return out, snap


def _write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(header)
        writer.writerows(rows)


def _json(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n",
                    encoding="utf-8")


def _simulation_records(sources, topology_sha256, trajectory_hashes):
    """If simulator manifests accompany a trajectory set, verify the entire set."""
    manifests = [path.with_suffix(".json") for _, path, _ in sources]
    if not any(p.exists() for p in manifests):
        return {}
    if not all(p.is_file() for p in manifests):
        raise ValueError("Her yörünge için tamamlanmış simülasyon manifesti gereklidir.")
    varying = {"seed", "thermostat_seed", "barostat_seed", "velocity_seed",
               "initial_geometry", "initial_coordinates_sha256"}
    shared_protocol = None
    ids = {}
    for (seed, path, _), manifest in zip(sources, manifests):
        record = json.loads(manifest.read_text("utf-8"))
        params = record["parameters"]
        if (record["trajectory_sha256"] != trajectory_hashes[path.name]
                or record["shared_topology_sha256"] != topology_sha256
                or params["seed"] != seed or record["trajectory"] != path.name):
            raise ValueError("Yörünge, tohum veya topoloji simülasyon manifestiyle uyuşmuyor.")
        common = {key: value for key, value in params.items() if key not in varying}
        if shared_protocol is not None and common != shared_protocol:
            raise ValueError("Tek bir ortalama için her replikanın fiziksel protokolü aynı olmalıdır.")
        shared_protocol = common
        ids[path.name] = record["run_id"]
    return ids


def analyze_replicas(topology_path, trajectories, output_dir, *, bin_width_nm=BIN_WIDTH_NM):
    topology_path, output_dir = Path(topology_path), Path(output_dir)
    if not np.isfinite(bin_width_nm) or bin_width_nm <= 0:
        raise ValueError("Pozitif histogram aralığı seçin.")
    top = json.loads(topology_path.read_text("utf-8"))
    selection = Selections(top)
    if not trajectories:
        raise ValueError("Üretim yörüngesi bulunamadı.")
    sources = []
    for value in trajectories:
        file = Path(value)
        m = re.fullmatch(r"seed_(\d+)\.npz", file.name)
        if not m:
            raise ValueError("Koşu dosyası seed_<tamsayı>.npz biçiminde olmalıdır.")
        with np.load(file, allow_pickle=False) as inp:
            lengths = box_lengths(inp["box_nm"])
        sources.append((int(m.group(1)), file, lengths))
    sources.sort(key=lambda run: run[0])
    if len(sources) != len({s for s, _, _ in sources}):
        raise ValueError("Aynı tohum iki kez bulundu.")
    common_half_box = min(L.min() for _, _, L in sources) / 2
    if common_half_box < NEIGHBOR_RADIUS_NM:
        raise ValueError("Hücre, hidratasyon sayım yarıçapı için çok küçüktür.")
    nbins = int(np.floor(common_half_box / bin_width_nm))
    if not 1 <= nbins <= 2500:
        raise ValueError("Histogram aralığı kutu boyutuna uygun olmalıdır.")
    edges = np.arange(nbins + 1) * bin_width_nm
    topology_sha256 = file_hash(topology_path)
    hashes = {path.name: file_hash(path) for _, path, _ in sources}
    run_ids = _simulation_records(sources, topology_sha256, hashes)
    runs, snapshot = [], None
    for k,(seed,path,lengths) in enumerate(sources):
        r, snap = one_run(path, lengths, selection, edges, seed, k == 0)
        runs.append(r)
        if snap is not None:
            snapshot = snap
    rg_all = [v for run in runs for v in run["rg_nm"]]
    rg_range = None
    if np.ptp(rg_all) < 1e-12:
        lower = max(0, rg_all[0] - 0.0005)
        rg_range = (lower, lower + 0.001)
    rg_edges = np.histogram_bin_edges(rg_all, bins=min(20, max(4, int(np.sqrt(len(rg_all))))),
                                      range=rg_range)
    rg_pdf = [np.histogram(r["rg_nm"], bins=rg_edges, density=True)[0] for r in runs]
    rg_means = [r["rg_mean_nm"] for r in runs]
    nrep = len(runs)
    hydrated = bool(len(selection.carbonyl))
    rs = np.array([r["rdf_g"] for r in runs])
    nc = np.array([r["coordination_to_outer_edge"] for r in runs])
    cn_fixed = [r["coordination_mean_at_035_nm"] for r in runs]
    rdf_sd = np.std(rs,axis=0,ddof=1).tolist() if nrep > 1 and hydrated else []
    mean_r = np.mean(rs,axis=0).tolist() if hydrated else []
    mean_nc = np.mean(nc,axis=0).tolist() if hydrated else []
    meta = dict(n_replicas=nrep, n_frames=len(rg_all), solute_atoms=len(selection.sol),
                carbonyl_oxygens=len(selection.carbonyl), water_oxygens=len(selection.water),
                phase="explicit_solvent", secondary_scope_note="Short model trajectories; between-replica SD is descriptive, not a confidence interval.")
    rg = dict(mean_nm=float(np.mean(rg_means)), between_replica_sd_nm=float(np.std(rg_means,ddof=1)) if nrep>1 else None,
              bin_edges_nm=rg_edges.tolist(), bin_mid_nm=((rg_edges[1:]+rg_edges[:-1])/2).tolist(),
              density_per_nm=np.mean(rg_pdf,axis=0).tolist(),
              definition="Mass-weighted radius of gyration of the complete bonded solute, including explicit H. Weights are equal across replica means.")
    rdf = dict(status="calculated" if hydrated else "no_matching_solute_sites",
               r_nm=((edges[1:]+edges[:-1])/2).tolist() if hydrated else [],
               r_outer_nm=edges[1:].tolist() if hydrated else [], mean_g=mean_r,
               between_replica_sd_g=rdf_sd, mean_coordination=mean_nc, bin_width_nm=bin_width_nm,
               r_max_nm=float(edges[-1]), selection="Amide carbonyl oxygens of the solute versus oxygen sites of all water molecules.",
               normalization="Within each sampled frame, g_k=H_k V/(N_ref N_water (4π/3)(r_out³−r_in³)). Arithmetic time average per replica, then equal weight per replica; PBC minimum-image distances. Density uses the full cell including the solute excluded volume.")
    coord = dict(cutoff_nm=NEIGHBOR_RADIUS_NM, mean_neighbors=float(np.mean(cn_fixed)) if hydrated else None,
                 between_replica_sd=float(np.std(cn_fixed,ddof=1)) if hydrated and nrep>1 else None,
                 definition="Direct count of water O ≤0.35 nm, averaged per amide carbonyl O and over saved frames. A water neighboring two references contributes to each site. Fixed distance is not a fitted first minimum or an H-bond definition.")
    hb_means = [r["hydrogen_bond_mean_per_carbonyl"] for r in runs]
    hbond = dict(donor="water oxygen bonded to its explicit hydrogen", acceptor="solute amide carbonyl oxygen",
                 max_donor_acceptor_nm=NEIGHBOR_RADIUS_NM, min_DHA_angle_deg=HBOND_MIN_ANGLE_DEG,
                 mean_per_carbonyl=float(np.mean(hb_means)) if hydrated else None,
                 between_replica_sd=float(np.std(hb_means,ddof=1)) if hydrated and nrep>1 else None,
                 definition="Geometric count of water-donor D–H···O(carbonyl) triplets per selected carbonyl O over production frames; minimum-image vectors. Excludes terminal N–H donors; no hydrogen-bond lifetime or energy inferred.")
    site_hydration = []
    for i, site in enumerate(selection.sites):
        neighbors = [r["site_mean_neighbors"][i] for r in runs]
        hbonds = [r["site_mean_hbonds"][i] for r in runs]
        site_hydration.append(dict(**site, mean_neighbors=float(np.mean(neighbors)),mean_hbonds=float(np.mean(hbonds)),
                              between_replica_sd_neighbors=float(np.std(neighbors,ddof=1)) if nrep>1 else None,
                              between_replica_sd_hbonds=float(np.std(hbonds,ddof=1)) if nrep>1 else None,
                              replicas=[dict(seed=r["seed"],mean_neighbors=n,mean_hbonds=h)
                                        for r,n,h in zip(runs,neighbors,hbonds)]))
    summary = dict(schema_version=SCHEMA,status="pilot",meta=meta,rg=rg,rdf=rdf,coordination=coord,
                   hydrogen_bonds=hbond,site_hydration=site_hydration,replicas=runs,
                   provenance=dict(topology_sha256=topology_sha256,trajectories_sha256=hashes,simulation_run_ids=run_ids,
                                   code_sha256=file_hash(__file__),numpy_version=np.__version__,
                                   canonical_smiles=top.get("canonical_smiles")))
    output_dir.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output_dir) as tmp:
        stage = Path(tmp)
        _write_csv(stage/"rg_trace.csv", ["seed","time_ps","rg_nm","water_donor_hbonds_per_carbonyl"],
                   ([r["seed"],t,r["rg_nm"][i],r["hydrogen_bonds_per_carbonyl"][i] if hydrated else ""]
                    for r in runs for i,t in enumerate(r["time_ps"])))
        _write_csv(stage/"rg_distribution.csv",["r_low_nm","r_high_nm","density_per_nm"],zip(rg_edges[:-1],rg_edges[1:],rg["density_per_nm"]))
        _write_csv(stage/"rdf.csv",["r_mid_nm","r_outer_nm","mean_g","sd_between_replica_g","mean_coordination"]+[f'g_seed_{r["seed"]}' for r in runs],
                   ([r,outer,mean_r[i],rdf_sd[i] if rdf_sd else "",mean_nc[i],*[run["rdf_g"][i] for run in runs]]
                    for i,(r,outer) in enumerate(zip(rdf["r_nm"],rdf["r_outer_nm"]))))
        id_columns = ["site_id","oxygen_atom_index","carbonyl_c_atom_index","amide_n_atom_index","amide_class",
                      "n_bonded_hydrogens","n_bonds_to_nearest_secondary_amide_n"]
        site_columns = id_columns+["mean_neighbors","between_replica_sd_neighbors","mean_hbonds","between_replica_sd_hbonds"]
        replica_columns = [f'{metric}_seed_{r["seed"]}' for r in runs for metric in ("mean_neighbors","mean_hbonds")]
        _write_csv(stage/"site_hydration.csv",site_columns+replica_columns,
                   ([site[k] for k in site_columns]+[r[k] for r in site["replicas"] for k in ("mean_neighbors","mean_hbonds")]
                    for site in site_hydration))
        _json(stage/"snapshot.json",snapshot)
        summary["provenance"]["outputs_sha256"] = {p.name:file_hash(p) for p in sorted(stage.iterdir())}
        _json(stage/"observables.json",summary)
        for name in ["rg_trace.csv","rg_distribution.csv","rdf.csv","site_hydration.csv","snapshot.json","observables.json"]:
            (stage/name).replace(output_dir/name)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("trajectories",nargs="+",type=Path)
    args = parser.parse_args()
    result=analyze_replicas(args.topology,args.trajectories,args.out)
    print(f"{result['meta']['n_replicas']} koşu; toplam {result['meta']['n_frames']} üretim karesi.")
