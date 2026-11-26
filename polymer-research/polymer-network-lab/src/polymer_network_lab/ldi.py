"""LDI polyurethane formulation scenarios and group-specific trajectory contacts."""

from __future__ import annotations

import json
from math import isfinite

import numpy as np

from .wbpu import label, numeric, read_rows


BEYCAN = "10.1021/acs.iecr.4c00759"
HAO = "10.1039/C6BM00588H"
CLASSES = ("PEG_ether", "PCL_ester", "LDI_ester", "urethane_carbonyl", "urea_carbonyl")
FIXTURE_NOTE = "Constructed coordinate frames for checking contact geometry; no force field or molecular dynamics."
SERIES_NOTE = "Illustrative formulation series with constructed contact descriptors and FTIR targets; no experimental accuracy is estimated."
SERIES_HEADER = ("formulation,specimen,protocol,peg_soft_mass_fraction,nco_oh,ether_enrichment,"
                 "carbonyl_area_fraction,trajectory_id,spectrum_id")

# Geometry intentionally places one donor over an orthorhombic periodic boundary.
def fixture():
    frames = []
    for i, (x_ether, x_pcl, y_urethane) in enumerate(
            [(0.12, 1.1, 0.18), (0.12, 1.1, 1.2), (1.1, 0.13, 1.2), (0.12, 0.13, 1.2)]):
        frames.append({"time_ps": float(i), "box_nm": [2, 2, 2],
                       "xyz_nm": [[1.85, .5, .5], [1.95, .5, .5], [x_ether, .5, .5],
                                  [1.2, 1, 1], [x_pcl, .5, .5], [0.12, y_urethane, .5],
                                  [1.99, .5, .5]]})
    return {"atom_chain": ["A","A","B","B","C","C","A"], "donor_NH": [[0, 1]],
            "acceptors": [[2,"PEG_ether"], [3,"PEG_ether"], [4,"PCL_ester"],
                          [5,"urethane_carbonyl"], [6,"urethane_carbonyl"]],
            "frames": frames}


DEMO_TRAJECTORY = json.dumps(fixture(), separators=(",", ":"))
DEMO_SERIES = SERIES_HEADER + """
F1,F1-a,illustrative dry films,0.15,1.0,0.65,0.35,sim1,example1
F2,F2-a,illustrative dry films,0.25,1.1,1.36,0.50,sim2,example2
F3,F3-a,illustrative dry films,0.35,1.0,0.88,0.42,sim3,example3
F4,F4-a,illustrative dry films,0.45,1.1,1.17,0.49,sim4,example4
F5,F5-a,illustrative dry films,0.55,1.0,0.75,0.42,sim5,example5
F6,F6-a,illustrative dry films,0.65,1.1,1.42,0.57,sim6,example6
F7,F7-a,illustrative dry films,0.75,1.0,1.04,0.51,sim7,example7
F8,F8-a,illustrative dry films,0.85,1.1,0.80,0.47,sim8,example8"""


def feed_scenario(peg_mmol, pcl_mmol, ldi_mmol, extender_mmol, extender, residual_nco, lysine_mmol):
    """Nominal, dry reactant amounts in mmol; optional measured NCO at dispersion."""
    amounts = [numeric(v, k, 0, 1e8, include_low=True) for k, v in
               (("PEG",peg_mmol),("PCL",pcl_mmol),("LDI",ldi_mmol),
                ("extender",extender_mmol),("L-lysine",lysine_mmol))]
    peg, pcl, ldi, ext, lys = amounts
    if extender not in ("BDO", "PDO", "glycerol"):
        raise ValueError("Choose BDO, PDO or glycerol.")
    oh_ext = ext * (3 if extender == "glycerol" else 2)
    nco, oh = 2*ldi, 2*(peg+pcl) + oh_ext
    if not peg+pcl or not nco:
        raise ValueError("Supply a macrodiol and LDI.")
    measured = None if residual_nco is None or str(residual_nco).strip() == "" else numeric(
        residual_nco, "residual NCO mmol before dispersion", 0, nco, include_low=True)
    # One primary amine site per reaction event; two primary amine sites per L-lysine molecule.
    upper = min(measured, 2*lys) if measured is not None else None
    return {"index": 100*nco/oh, "oh": oh, "nco": nco, "measured": measured, "urea_upper": upper,
            "triol": ext if extender=="glycerol" else 0,
            "rows":[["PEG diol",peg,2*peg,"OH"],["PCL diol",pcl,2*pcl,"OH"],["LDI",ldi,nco,"NCO"],
                    [extender,ext,oh_ext,"OH"],["L-lysine (aqueous stage)",lys,2*lys,"primary amine sites"]]}


def index(value, count):
    if type(value) is not int or value<0 or value>=count:
        raise ValueError("Atom indices must be distinct, zero-based integers within the topology.")
    return value


def topology(meta):
    if not isinstance(meta,dict):
        raise ValueError("Supply a JSON object with atom_chain, donor_NH, acceptors and frames.")
    chain = meta.get("atom_chain")
    if not isinstance(chain,list) or not 3<=len(chain)<=250000 or any(not isinstance(c,str) or not c.strip() for c in chain):
        raise ValueError("Assign a molecular chain ID to every atom.")
    pairs = meta.get("donor_NH")
    acceptors = meta.get("acceptors")
    if not isinstance(pairs,list) or not pairs or not isinstance(acceptors,list) or not acceptors:
        raise ValueError("Supply urethane [N,H] atom pairs and [oxygen,class] acceptors.")
    if any(not isinstance(p,list) or len(p)!=2 for p in [*pairs,*acceptors]):
        raise ValueError("A donor or acceptor record must have exactly two entries.")
    donors, hydrogen = zip(*[(index(n,len(chain)),index(h,len(chain))) for n,h in pairs])
    oxygen = [index(o,len(chain)) for o,k in acceptors]
    if len(set([*donors,*hydrogen,*oxygen])) != 2*len(donors)+len(oxygen):
        raise ValueError("Use one unique N,H pair per urethane donor and one group per distinct acceptor oxygen.")
    if any(chain[d]!=chain[h] for d,h in zip(donors,hydrogen)):
        raise ValueError("A bonded N and H must belong to the same chain.")
    kinds=[k for _,k in acceptors]
    if any(k not in CLASSES for k in kinds):
        raise ValueError("Unknown acceptor class: choose "+", ".join(CLASSES))
    return np.array(donors),np.array(hydrogen),np.array(oxygen),np.array(kinds),np.array(chain)


def contact_analysis(meta, frames, cutoff_nm=0.35, min_angle_deg=150.):
    """Count interchain D–H···A triplets; frames are equally spaced and lengths are nm."""
    cutoff=numeric(cutoff_nm,"N···O cutoff (nm)",0,0.6)
    angle=numeric(min_angle_deg,"minimum N–H···O angle",0,180)
    d,h,o,kinds,chain=topology(meta)
    possible=np.array([sum((chain[o][kinds==k]!=chain[di]).sum() for di in d) for k in CLASSES])
    if not possible.sum():
        raise ValueError("Assign at least one acceptor on another chain.")
    totals=np.zeros(len(CLASSES),int)
    by_frame=[]
    previous=None
    stride=None
    previous_pairs=set()
    class_by_oxygen=dict(zip(o,kinds))
    carried=np.zeros(len(CLASSES),int)
    total_preceding=np.zeros(len(CLASSES),int)
    starts={}
    longest=np.zeros(len(CLASSES),int)
    complete_runs=[[] for _ in CLASSES]
    for frame in frames:
        if not isinstance(frame,dict):
            raise ValueError("Each frame must contain time_ps, box_nm and xyz_nm.")
        try:
            xyz=np.asarray(frame["xyz_nm"],dtype=float)
            box=np.asarray(frame["box_nm"],dtype=float)
            t=float(frame["time_ps"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Invalid frame coordinates, time or orthorhombic box.") from exc
        if xyz.shape!=(len(chain),3) or not np.isfinite(xyz).all() or box.shape!=(3,) or not np.isfinite(box).all() or np.any(box<=2*cutoff) or not isfinite(t):
            raise ValueError("Finite xyz in nm and an orthorhombic box longer than twice the cutoff are required.")
        if previous is not None:
            dt=t-previous
            if dt<=0 or (stride is not None and not np.isclose(dt,stride,rtol=1e-4,atol=1e-8)):
                raise ValueError("Provide monotonically increasing, equally spaced trajectory frames.")
            if stride is None:
                stride=dt
        previous=t
        def image(v):
            return v-box*np.rint(v/box)
        hd=image(xyz[d]-xyz[h])
        nh_len=np.linalg.norm(hd,axis=1)
        if np.any((nh_len<0.06)|(nh_len>0.15)):
            raise ValueError("N–H bond outside 0.06–0.15 nm; check atom selection and units.")
        events=set()
        counts=np.zeros(len(CLASSES),int)
        for j,(n,hy) in enumerate(zip(d,h)):
            da=image(xyz[o]-xyz[n])
            cand=np.flatnonzero((chain[o]!=chain[n]) & (np.einsum("ij,ij->i",da,da)<=cutoff**2))
            if not len(cand):
                continue
            ha=image(xyz[o[cand]]-xyz[hy])
            dist=np.linalg.norm(ha,axis=1)
            if np.any(dist<1e-8):
                raise ValueError("Hydrogen and acceptor coordinates coincide.")
            cosine=ha@hd[j]/(dist*nh_len[j])
            for ix in cand[cosine<=np.cos(np.deg2rad(angle))+1e-12]:
                counts[CLASSES.index(kinds[ix])]+=1
                events.add((int(n),int(o[ix])))
        frame_number=len(by_frame)
        for pair in previous_pairs:
            group=CLASSES.index(class_by_oxygen[pair[1]])
            total_preceding[group]+=1
            if pair in events:
                carried[group]+=1
            else:
                first=starts.pop(pair)
                if first>0:  # Both sides of this episode were observed unbonded.
                    complete_runs[group].append(frame_number-1-first)
        for pair in events:
            first=starts.setdefault(pair,frame_number)
            group=CLASSES.index(class_by_oxygen[pair[1]])
            longest[group]=max(longest[group],frame_number-first)
        previous_pairs=events
        totals+=counts
        by_frame.append([t,*counts.tolist()])
    if not by_frame:
        raise ValueError("Supply at least one frame.")
    rows=[]
    for k,possible_pairs,n_events in zip(CLASSES,possible,totals):
        if not np.any(kinds==k):
            continue
        expected=possible_pairs/possible.sum()
        share=n_events/totals.sum() if totals.sum() else None
        enrich=float(share/expected) if share is not None and expected else None
        rows.append([k,int((kinds==k).sum()),int(possible_pairs),int(n_events),
                     float(n_events/(len(d)*len(by_frame))),float(expected),enrich])
    persistence=[]
    for i,k in enumerate(CLASSES):
        if not np.any(kinds==k):
            continue
        runs=complete_runs[i]
        lag=float(carried[i]/total_preceding[i]) if total_preceding[i] else None
        max_span=float(longest[i]*stride) if totals[i] and stride is not None else None
        mean_span=float(np.mean(runs)*stride) if runs and stride is not None else None
        persistence.append([k,lag,max_span,len(runs),mean_span])
    return {"rows":rows,"persistence":persistence,"by_frame":by_frame,"frames":len(by_frame),"donors":len(d),
            "total":int(totals.sum()),"stride":stride,
            "lag_survival":float(carried.sum()/total_preceding.sum()) if total_preceding.sum() else None}


def contact_folds(csv_text):
    """Fixed ridge (alpha=1), standardized within each training fold, group means."""
    groups={}
    ids=set()
    ownership={"trajectory_id":{}, "spectrum_id":{}}
    protocols=set()
    for row in read_rows(csv_text, SERIES_HEADER):
        f=label(row,"formulation")
        specimen=label(row,"specimen")
        if specimen in ids:
            raise ValueError("Each specimen row must be identified once. Aggregate spectra and MD seeds per specimen first.")
        ids.add(specimen)
        protocols.add(label(row,"protocol"))
        for key, owner in ownership.items():
            origin=label(row,key)
            if origin in owner and owner[origin]!=f:
                raise ValueError(f"{key} is assigned to more than one formulation. Resolve shared observations before fitting.")
            owner[origin]=f
        peg=numeric(row["peg_soft_mass_fraction"],"PEG fraction",0,1,include_low=True)
        ratio=numeric(row["nco_oh"],"NCO/OH",0,20)
        enrichment=numeric(row["ether_enrichment"],"ether enrichment",0,1000,include_low=True)
        y=numeric(row["carbonyl_area_fraction"],"fitted carbonyl area fraction",0,1,include_low=True)
        values=groups.setdefault(f,[])
        if values and (peg,ratio)!=tuple(values[0][:2]):
            raise ValueError("Each formulation must have one consistent PEG mass fraction and NCO/OH.")
        values.append([peg,ratio,enrichment,y])
    if len(groups)<5 or len(protocols)!=1:
        raise ValueError("Supply at least five different formulations with one matched FTIR fitting and conditioning protocol.")
    means=np.array([np.mean(v,axis=0) for v in groups.values()])
    output=[]
    for i, group in enumerate(groups):
        train=np.arange(len(means))!=i
        y=means[train,3]
        preds=[]
        for nfeatures in (2,3):
            x=means[train,:nfeatures]
            avg=x.mean(axis=0)
            std=x.std(axis=0)
            std[std<1e-12]=1.0
            z=(x-avg)/std
            coeff=np.linalg.solve(z.T@z+np.eye(nfeatures),z.T@(y-y.mean()))
            preds.append(float(y.mean()+((means[i,:nfeatures]-avg)/std)@coeff))
        output.append([group,len(groups[group]),float(means[i,3]),*preds])
    err=np.abs(np.array([r[3:] for r in output])-means[:,3:4])
    return {"rows":output,"mae_pp":(err.mean(axis=0)*100).tolist(),"protocol":next(iter(protocols))}
