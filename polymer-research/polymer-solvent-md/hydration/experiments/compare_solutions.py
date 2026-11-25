"""Compare the explicitly hydrated DP3 PEtOx H/H and Me/OH models.

Usage (run from the study directory):
    python experiments/compare_solutions.py --hh hh-long --caps caps/run \
        --pilot ../polymer-solvent-md/pilot/run --out results --plots

The primary window starts at the first new production sample. All production
replicates must share a cadence. Intersect saved times and weight each
replicate equally. No transient frames are excluded after looking at results.
"""
from collections import deque
import argparse
import csv
import hashlib
import json
from pathlib import Path
from statistics import NormalDist

import numpy as np

RADIUS = 0.35
MIN_ANGLE = 150.0
SHORT_RANGE_CUTOFF = 0.80
SITES = ("alpha", "middle", "omega")
PRIMARY = ("rg_nm", "core_rg_nm", "carbonyl_neighbors", "water_to_carbonyl_hb")
SITE_NAMES = {"alpha": "α-uç", "middle": "Orta", "omega": "ω-uç"}


def file_sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(Path(path).read_bytes())


def shortest(graph, root):
    distance, edges = {root: 0}, []
    queue = deque([root])
    while queue:
        a = queue.popleft()
        for b in sorted(graph[a]):
            if b not in distance:
                edges.append((a, b))
                distance[b] = distance[a] + 1
                queue.append(b)
    return distance, edges


class Model:
    def __init__(self, topology, cap):
        self.cap = cap
        self.elem = topology["elements"]
        sol = set(topology["solute_indices"])
        self.sol = sorted(sol)
        self.mass = np.asarray(topology["mass_da"], float)
        self.heavy = sorted(i for i in sol if self.elem[i] != "H")
        self.water = np.asarray(topology["water_oxygens"], int)
        graph = {i: set() for i in range(len(self.elem))}
        for a, b in topology["bonds"]:
            graph[a].add(b)
            graph[b].add(a)
        self.graph = graph
        self.h_water = np.array([sorted(b for b in graph[w] if self.elem[b] == "H") for w in self.water])
        if self.h_water.shape != (len(self.water), 2) or len(self.water) != 500:
            raise ValueError("Expected 500 explicit three-site waters per single solute.")
        self.o_to_n = {}
        for o in topology["carbonyl_oxygens"]:
            c, = graph[o]
            n, = [a for a in graph[c] if self.elem[a] == "N"]
            self.o_to_n[o] = n
        if len(self.o_to_n) != 3:
            raise ValueError("Expected exactly three amide carbonyls.")
        def neighbors(a, el):
            return sorted(b for b in graph[a] if self.elem[b] == el)
        caps = set()
        if cap == "HH":
            n_alpha, = [n for n in self.o_to_n.values() if len(neighbors(n, "H")) == 1]
            self.terminal_donor = n_alpha
            self.terminal_acceptor = None
        elif cap == "MeOH":
            n_alpha, methyl = next((n,b) for n in self.o_to_n.values() for b in neighbors(n,"C")
                                   if len(neighbors(b,"H")) == 3)
            if any(neighbors(n, "H") for n in self.o_to_n.values()):
                raise ValueError("Me/OH must contain three tertiary amides.")
            oh, = [i for i in sol if self.elem[i] == "O" and len(neighbors(i, "H")) == 1]
            caps = {methyl, oh}
            self.terminal_acceptor = oh
            self.terminal_donor = oh
        else:
            raise ValueError("Unknown end chemistry.")
        self.donor_h, = neighbors(self.terminal_donor,"H")
        d, self.tree = shortest(graph, n_alpha)
        if set(d) != sol:
            raise ValueError("Solute graph is disconnected or contains a solvent bond.")
        self.anchor = n_alpha
        os = sorted(self.o_to_n, key=lambda o:d[self.o_to_n[o]])
        if [d[self.o_to_n[o]] for o in os] != [0,3,6]:
            raise ValueError("Unexpected DP3 backbone graph.")
        self.refs = np.array(os)
        if self.terminal_acceptor is not None and d[self.terminal_acceptor] != 9:
            raise ValueError("Hydroxyl assignment is inconsistent with the opposite backbone end.")
        self.core = sorted(set(self.heavy) - caps)
        from collections import Counter
        if Counter(self.elem[i] for i in self.core) != {"C":15,"N":3,"O":3}:
            raise ValueError("Common heavy-atom scaffold must be C15N3O3.")
        if not np.all(np.isfinite(self.mass)) or np.any(self.mass <= 0):
            raise ValueError("Invalid masses.")

    def selection(self):
        return {"shared_core_indices":self.core,
                "terminal_donor":self.terminal_donor,"terminal_hydrogen":self.donor_h,
                "terminal_acceptor":self.terminal_acceptor,
                "carbonyl_oxygens":dict(zip(SITES,map(int,self.refs))),
                "carbonyl_amide_nitrogens":dict(zip(SITES,[self.o_to_n[o] for o in self.refs]))}


def mic(vector, L):
    return vector - L*np.rint(vector / L)


def water_to_acceptors(xyz, L, O, H, acceptors, *, keep_occupancy=False):
    """Count D–H···A water-donor triplets with per-reference neighbour counts."""
    da = mic(xyz[O][None,:,:] - xyz[acceptors][:,None,:], L)
    near = np.sum(da*da, axis=-1) <= RADIUS**2
    hb = np.zeros(len(acceptors), int)
    a, w = np.where(near)
    if a.size:
        wh = xyz[H[w]]
        hd = mic(xyz[O[w],None,:] - wh, L)
        ha = mic(xyz[acceptors[a],None,:] - wh, L)
        norm = np.linalg.norm(hd,axis=-1)*np.linalg.norm(ha,axis=-1)
        if np.any(norm < 1e-9):
            raise ValueError("Degenerate hydrogen-bond geometry.")
        cosine = np.sum(hd*ha,axis=-1)/norm
        hb += np.bincount(a, weights=np.sum(cosine <= np.cos(np.deg2rad(MIN_ANGLE)), axis=1),
                          minlength=len(acceptors)).astype(int)
    counts = np.sum(near,axis=1)
    return (counts, hb, near) if keep_occupancy else (counts, hb)


def solute_donor_to_water(xyz, L, O, donor, h):
    mask = np.linalg.norm(mic(xyz[O]-xyz[donor],L),axis=1) <= RADIUS
    if not np.any(mask):
        return 0
    hd = mic(xyz[donor]-xyz[h],L)
    ha = mic(xyz[O[mask]]-xyz[h],L)
    cosine = (ha@hd)/(np.linalg.norm(ha,axis=1)*np.linalg.norm(hd))
    return int(np.sum(cosine <= np.cos(np.deg2rad(MIN_ANGLE))))


def mass_rg(coords, masses):
    center = np.average(coords,axis=0,weights=masses)
    return float(np.sqrt(np.average(np.sum((coords-center)**2,axis=1),weights=masses)))


def frame(xyz, L, model):
    # Connected, whole coordinates with the exact published solute atom indexing.
    whole = np.empty_like(xyz)
    whole[model.anchor] = xyz[model.anchor]
    for a,b in model.tree:
        whole[b] = whole[a] + mic(xyz[b]-xyz[a],L)
    rg = mass_rg(whole[model.sol],model.mass[model.sol])
    core = mass_rg(whole[model.core],model.mass[model.core])
    refs = model.refs
    if model.terminal_acceptor is not None:
        refs = np.append(refs,model.terminal_acceptor)
    cn,hb,occupancy = water_to_acceptors(xyz,L,model.water,model.h_water,refs,keep_occupancy=True)
    record = {"rg_nm": rg, "core_rg_nm":core, "carbonyl_neighbors": float(np.mean(cn[:3])),
              "water_to_carbonyl_hb": float(np.mean(hb[:3])),
              "terminal_to_water_hb":solute_donor_to_water(xyz,L,model.water,model.terminal_donor,model.donor_h)}
    if model.terminal_acceptor is not None:
        record.update(water_to_terminal_OH_hb=int(hb[3]), terminal_OH_neighbors=int(cn[3]))
    for i,site in enumerate(SITES):
        record[f"{site}_neighbors"] = int(cn[i])
        record[f"{site}_water_hb"] = int(hb[i])
    return record, occupancy[:3]


def load_case(folder, cap):
    folder = Path(folder)
    topfile = folder/"topology.json"
    if not topfile.is_file() and (folder/"run"/"topology.json").is_file():
        folder = folder/"run"
        topfile = folder/"topology.json"
    model = Model(read(topfile),cap)
    entries=[]
    topology_digest=file_sha(topfile)
    for file in sorted(folder.glob("seed_*.npz")):
        meta_path = file.with_suffix(".json")
        if not meta_path.is_file():
            raise ValueError(f"Missing run manifest for {file.name}")
        meta=read(meta_path)
        actual=file_sha(file)
        registered=meta.get("trajectory_sha256")
        if not registered or registered != actual:
            raise ValueError(f"Trajectory checksum missing or mismatched for {file.name}.")
        param=meta.get("parameters",{})
        if (meta.get("shared_topology_sha256")!=topology_digest or meta.get("trajectory")!=file.name
                or param.get("seed")!=int(file.stem.split("_")[1])
                or param.get("n_water")!=len(model.water) or param.get("ensemble")!="NPT"
                or param.get("solute_model")!="GAFF2 + AM1-BCC" or param.get("water_model")!="TIP3P"
                or not np.isclose(param.get("temperature_K",-1),298.15)
                or not np.isclose(param.get("pressure_bar",-1),1.)
                or not np.isclose(param.get("cutoff_nm",-1),SHORT_RANGE_CUTOFF)):
            raise ValueError("Run manifest does not match solute/solvent method, topology, seed or registered conditions.")
        raw=np.load(file,allow_pickle=False)
        with raw:
            xyz = raw["positions_nm"]
            t = np.array(raw["time_ps"],float)
            boxes = np.array(raw["box_nm"],float)
            if boxes.ndim == 2 and boxes.shape == (len(t),3):
                L=boxes
            elif boxes.shape==(len(t),3,3):
                L=np.diagonal(boxes,axis1=1,axis2=2)
                if np.max(np.abs(boxes-(np.eye(3)[None]*L[:,None,:])))>1e-5:
                    raise ValueError("Analysis expects orthorhombic periodic vectors.")
            else:
                raise ValueError("Invalid periodic cell arrays.")
            if xyz.shape[:2] != (len(t),len(model.elem)) or len(t)<8 or xyz.shape[-1]!=3 or not np.isfinite(xyz).all():
                raise ValueError("Invalid trajectory arrays.")
            if np.any(~np.isfinite(L)) or np.min(L)<=2*SHORT_RANGE_CUTOFF:
                raise ValueError("Cells are too small for the registered direct interaction cutoff.")
            dt=float(np.median(np.diff(t)))
            if (meta.get("frames")!=len(t) or not np.isclose(param.get("sample_ps",-1),dt)
                    or not np.isclose(param.get("production_ps",-1),t[-1])):
                raise ValueError("Duration or stride reported by the engine disagrees with the stored coordinates.")
            if dt<=0 or not np.allclose(np.diff(t),dt,atol=0.0001):
                raise ValueError("Production time axis must be evenly sampled.")
            vals=[frame(np.asarray(p,dtype=float),l,model) for p,l in zip(xyz,L)]
            arr, occupancy = zip(*vals)
            occupancy = np.asarray(occupancy, dtype=bool)
            for key in ("density_g_ml","temperature_K"):
                if key in raw:
                    for row,val in zip(arr,raw[key],strict=True):
                        row[key]=float(val)
        entries.append({"name":file.stem,"seed":int(file.stem.split("_")[1]),"run_id":meta.get("run_id"),
                        "time_ps":t,"stride_ps":dt,"frames":arr,"shell_occupancy":occupancy,"trajectory_sha256":actual,
                        "run_manifest_sha256":file_sha(meta_path),"topology_sha256":topology_digest,
                        "minimum_self_image_distance_nm":meta["minimum_self_image_distance_nm"],
                        "sample_frames_with_self_image_within_cutoff":meta["sample_frames_with_self_image_within_cutoff"],
                        "solute_mass_fraction":meta["solute_mass_fraction"],
                        "protocol":meta.get("parameters",meta.get("protocol",{}))})
    if len(entries)!=3:
        raise ValueError(f"{cap} requires exactly three completed independent trajectories (found {len(entries)}).")
    return model,entries


def ranks_as_normal(values):
    flat=np.asarray(values).ravel()
    _, inv, cnt=np.unique(flat,return_inverse=True,return_counts=True)
    offset=np.cumsum(cnt)-cnt
    ranks=offset+(cnt+1)/2
    inv_cdf=NormalDist().inv_cdf
    q=np.array([inv_cdf((x-3/8)/(len(flat)+1/4)) for x in ranks])
    return q[inv].reshape(np.shape(values))


def rhat_stat(x):
    m,n=x.shape
    W=np.mean(np.var(x,axis=1,ddof=1))
    if W<=1e-15:
        return None
    B=n*np.var(x.mean(axis=1),ddof=1)
    return float(np.sqrt(((n-1)*W/n+B/n)/W))


def rank_folded_split_rhat(values):
    v=np.array(values,float)
    n=v.shape[1]//2
    if n<4:
        return None
    split=np.concatenate([v[:,:n],v[:,-n:]],axis=0)
    bulk=rhat_stat(ranks_as_normal(split))
    tail=rhat_stat(ranks_as_normal(np.abs(split-np.median(split))))
    a=[r for r in (bulk,tail) if r is not None]
    return max(a) if a else None


def acf_diagnostic(x,dt):
    """Finite-series initial-positive/monotone autocorrelation pair sums."""
    x=np.asarray(x,float)-np.mean(x)
    n=len(x)
    if n<8 or np.dot(x,x)<1e-15:
        return {"g":None,"tau_int_ps":None,"ess":None}
    rho=np.correlate(x,x,"full")[n-1:]/np.dot(x,x)
    pairs=[]
    for a,b in zip(rho[::2],rho[1::2]):
        p=float(a+b)
        if p<=0:
            break
        pairs.append(min(p,pairs[-1]) if pairs else p)
    g=float(np.clip(2*sum(pairs)-1,1,n))
    return {"g":g,"tau_int_ps":float(dt*g/2),"ess":float(n/g)}


def summarize(x,dt):
    x=np.asarray(x,float)
    means=np.mean(x,axis=1)
    blocks=[np.mean(b,axis=1).tolist() for b in np.array_split(x,4,axis=1)]
    n=x.shape[1]//2
    delta=x[:,-n:].mean(axis=1)-x[:,:n].mean(axis=1)
    # Equal-size empirical distributions: 1D Wasserstein-1 = mean absolute
    # difference of their ordered samples, in the original observable units.
    ordered=np.sort(x,axis=1)
    w1=[float(np.mean(np.abs(ordered[i]-ordered[j]))) for i in range(len(ordered)) for j in range(i)]
    return {"replica_mean_range":float(np.ptp(means)),"maximum_pairwise_w1":max(w1),"mean":float(np.mean(means)),"between_run_sd":float(np.std(means,ddof=1)),
            "per_run_mean":means.tolist(),"per_run_blocks":np.array(blocks).T.tolist(),
            "per_run_last_minus_first_half":delta.tolist(),"split_rank_rhat":rank_folded_split_rhat(x),
            "per_run_acf":[acf_diagnostic(z,dt) for z in x]}


def retention(occupancy, lag):
    """Fraction of initially occupied (time, site, water) pairs by site.

    Water identity is fixed within a replica. Intermittent retention allows
    departures/returns between endpoints. Sampled-continuous retention
    requires all recorded frames; brief unsampled exits remain unresolved.
    """
    if not 1 <= lag < len(occupancy):
        raise ValueError("Positive lag must fit inside the production window.")
    origins=occupancy[:-lag]
    trials=origins.sum(axis=(0,2))
    endpoint=origins & occupancy[lag:]
    uninterrupted=endpoint.copy()
    for step in range(1,lag):
        uninterrupted &= occupancy[step:len(occupancy)-lag+step]
    return trials, endpoint.sum(axis=(0,2)), uninterrupted.sum(axis=(0,2))


def write_csv(path,rows,keys=None):
    if not rows:
        raise ValueError(f"Empty CSV {path}")
    keys=keys or list(rows[0])
    with Path(path).open("w",newline="",encoding="utf-8") as fp:
        w=csv.DictWriter(fp,fieldnames=keys,extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def compare(hh, caps, out, pilot=None):
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    cases={"HH":load_case(hh,"HH"),"MeOH":load_case(caps,"MeOH")}
    runs=[r for _,reps in cases.values() for r in reps]
    dt=runs[0]["stride_ps"]
    if any(not np.isclose(dt,r["stride_ps"]) for r in runs):
        raise ValueError("Matched end-group study requires a common stored cadence.")
    keys=("ensemble","temperature_K","pressure_bar","sample_ps","equil_ps","production_ps","timestep_ps",
          "cutoff_nm","pme_relative_tolerance","rigid_water","solute_constraints","water_parameters_sha256",
          "thermostat","barostat","barostat_interval_steps","friction_per_ps","solute_model","engine_version")
    if any(r["protocol"][k]!=runs[0]["protocol"][k] for r in runs for k in keys):
        raise ValueError("Matched production protocols differ between replicas or models.")
    # Compare registered initial segment only. Keep each run's absolute production time.
    if any(not np.isclose(r["time_ps"][0],dt) for r in runs):
        raise ValueError("Require first production frame one stride after t=0.")
    ns=min(len(r["frames"]) for r in runs)
    if not all(np.allclose(r["time_ps"][:ns],runs[0]["time_ps"][:ns]) for r in runs):
        raise ValueError("New production frames do not have common physical times.")
    result={"schema":"petox-end-group-comparison/1","matched_frames_per_run":ns,
            "sample_interval_ps":dt,"matched_interval_ps":[float(runs[0]["time_ps"][0]),float(runs[0]["time_ps"][ns-1])],
            "matched_duration_ps":float(ns*dt),"runs_per_model":3,"cases":{},"contrast":{},
            "analysis_sha256":file_sha(__file__), "analysis_numpy":np.__version__,
            "common_protocol":{k:runs[0]["protocol"][k] for k in keys},
            "definitions":{"CN":"water O within 0.35 nm of carbonyl O. Site counts may share a water.",
              "HB":"Water-donor D–H···carbonyl O, D···A <=0.35 nm, D–H···A >=150°. Triplets per carbonyl.",
              "core_rg":"Mass-weighted common repeat heavy atoms C15N3O3, with both models excluding H and Me/OH model also excluding cap methyl C and terminal O.",
              "rhat":"Rank-normalized split R-hat and folded rank-normalized split R-hat, upper value. Exploratory diagnostic for MD replicas.",
              "ess":"Per-run finite-lag, positive/monotone-pair autocorrelation estimate. tau_int=0.5*g*dt. Short records can understate slow modes.",
              "W1":"Maximum pairwise Wasserstein-1 distance among empirical replica distributions, reported in each observable’s native unit; descriptive, without an independent-frame p-value.",
              "SD":"Descriptive sample SD across three independently initialized run means; not a confidence interval.",
              "contrast":"Mean(Me/OH) - mean(H/H). Replica means weighted equally, and whole-cap change is joint. No pairing across chemistries.",
              "blocks":"Four time blocks in order; arrays ordered by increasing numeric replica seed."}}
    trace=[];sites_rows=[];stat_rows=[];exchange_rows=[]
    for cap,(model,reps) in cases.items():
        names=list(reps[0]["frames"][0])
        metrics={}
        for k in names:
            vals=np.array([[f[k] for f in r["frames"][:ns]] for r in reps])
            metrics[k]=summarize(vals,dt)
            s=metrics[k]
            for i,rep in enumerate(reps):
                stat_rows.append({"model":cap,"seed":rep["seed"],"metric":k,"mean":s["per_run_mean"][i],
                                  "early_to_late_change":s["per_run_last_minus_first_half"][i],
                                  **{"block_"+str(b+1):y for b,y in enumerate(s["per_run_blocks"][i])},
                                  **s["per_run_acf"][i], "model_mean":s["mean"], "between_run_sd":s["between_run_sd"],
                                  "rank_folded_split_rhat":s["split_rank_rhat"]})
        for r in reps:
            for t,f in zip(r["time_ps"][:ns],r["frames"][:ns]):
                trace.append({"model":cap,"seed":r["seed"],"time_ps":float(t),**f})
        for i,name in enumerate(SITES):
            c=metrics[name+"_neighbors"];h=metrics[name+"_water_hb"]
            sites_rows.append({"model":cap,"site":name,"carbonyl_O_index":int(model.refs[i]),
                               "amide_N_index":model.o_to_n[model.refs[i]],"neighbors":c["mean"],
                               "neighbors_between_run_sd":c["between_run_sd"],"water_donor_hb":h["mean"],
                               "hb_between_run_sd":h["between_run_sd"]})
        result["cases"][cap]={"seeds":[r["seed"] for r in reps],"selection":model.selection(),
                               "metrics":metrics,"runs":[{k:r[k] for k in ("seed","run_id","trajectory_sha256",
                                   "run_manifest_sha256","topology_sha256","minimum_self_image_distance_nm",
                                   "sample_frames_with_self_image_within_cutoff","solute_mass_fraction")} for r in reps]}
        curves=[]
        for lag_ps in (1.0,2.0,5.0,10.0):
            lag=round(lag_ps/dt)
            if not 1<=lag<ns or not np.isclose(lag*dt,lag_ps):
                continue
            counts=[retention(r["shell_occupancy"][:ns],lag) for r in reps]
            for i,site in enumerate(SITES):
                endpoints=[];continues=[]
                for rep,(den, ep, co) in zip(reps,counts):
                    if den[i]==0:
                        raise ValueError("No eligible hydration-shell time origins for conditional retention.")
                    epfrac=float(ep[i]/den[i]); cofrac=float(co[i]/den[i])
                    endpoints.append(epfrac);continues.append(cofrac)
                    exchange_rows.append({"model":cap,"seed":rep["seed"],"site":site,"lag_ps":lag_ps,
                        "occupancy_origins":int(den[i]),"endpoint_matches":int(ep[i]),
                        "all_sampled_matches":int(co[i]),"endpoint_fraction":epfrac,"all_sampled_fraction":cofrac})
                curves.append({"site":site,"lag_ps":lag_ps,"endpoint_mean":float(np.mean(endpoints)),
                               "endpoint_between_run_sd":float(np.std(endpoints,ddof=1)),
                               "all_sampled_mean":float(np.mean(continues)),"all_sampled_between_run_sd":float(np.std(continues,ddof=1))})
        result["cases"][cap]["water_exchange"]=curves
    for k in (set(result["cases"]["HH"]["metrics"]) & set(result["cases"]["MeOH"]["metrics"])):
        a=result["cases"]["HH"]["metrics"][k]
        b=result["cases"]["MeOH"]["metrics"][k]
        result["contrast"][k]={"mean_difference":b["mean"]-a["mean"],"HH_mean":a["mean"],"MeOH_mean":b["mean"]}
    if pilot is not None:
        _,old=load_case(pilot,"HH")
        olddt=old[0]["stride_ps"]
        metrics={}
        for k in PRIMARY:
            metrics[k]=summarize([[r0[k] for r0 in r["frames"]] for r in old],olddt)
        result["archived_pilot"]={"seeds":[r["seed"] for r in old],"duration_ps":float(old[0]["time_ps"][-1]),
            "stride_ps":olddt,"metrics":metrics,"trajectory_sha256":[r["trajectory_sha256"] for r in old]}
        # H/H lineage supports a longitudinal comparison. Chemistry contrasts
        # are separate, with no atomwise pairing across differently sized solutes.
        old_by_seed={r["seed"]:r for r in old}
        current=cases["HH"][1]
        has_parent=[r["protocol"]["initial_geometry"].get("source_trajectory_sha256")
                    ==old_by_seed[r["seed"]]["trajectory_sha256"] for r in current]
        if all(has_parent):
            old_duration=result["archived_pilot"]["duration_ps"]
            nwin=int(np.floor((old_duration+1e-7)/dt))
            paired=[]
            for r in current:
                prior=old_by_seed[r["seed"]]
                if nwin<1 or nwin*2>ns:
                    raise ValueError("Extension must supply separate early and late windows matching the archived pilot.")
                for k in PRIMARY:
                    x=np.array([f[k] for f in r["frames"][:ns]])
                    y=np.mean([f[k] for f in prior["frames"]])
                    paired.append({"seed":r["seed"],"metric":k,"original_pilot":float(y),
                                   "new_first_window":float(np.mean(x[:nwin])),
                                   "new_last_window":float(np.mean(x[-nwin:])),
                                   "new_overall":float(np.mean(x)),"new_overall_minus_pilot":float(np.mean(x)-y)})
            result["longitudinal_hh"]={"window_ps":nwin*dt,"rows":paired,
                "parent_trajectory_verified":True,
                "scope":"Initial and final new windows and archived pilot are separately summarized. Velocities were resampled before the new production."}
            write_csv(out/"hh_lineages.csv",paired)
    json_text=json.dumps(result,ensure_ascii=False,allow_nan=False,sort_keys=True,indent=2)
    (out/"comparison.json").write_text(json_text+"\n",encoding="utf-8")
    ordered={k:None for row in trace for k in row}
    write_csv(out/"replica_trace.csv",trace,ordered)
    write_csv(out/"replica_statistics.csv",stat_rows)
    write_csv(out/"site_hydration.csv",sites_rows)
    if exchange_rows:
        write_csv(out/"water_exchange.csv",exchange_rows)
    return result,trace


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hh",required=True,type=Path)
    parser.add_argument("--caps",required=True,type=Path)
    parser.add_argument("--pilot",type=Path)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--plots",action="store_true")
    args=parser.parse_args()
    result,trace=compare(args.hh,args.caps,args.out,args.pilot)
    if args.plots:
        from comparative_figures import figures
        figures(result,trace,args.out)
    print(result["matched_duration_ps"],{c:result["cases"][c]["metrics"]["rg_nm"]["mean"] for c in result["cases"]})
