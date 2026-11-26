"""Audit contacting, periodically solvated trimer pairs and resolve interfacial water identities.

Input directories are produced by contact_pair.py. Summary statistics are descriptive;
short, correlated frames are never treated as independent samples in an error estimate.
"""
import argparse
import csv
import itertools
from pathlib import Path

import numpy as np

from contact_pair import checksum, make_whole, mc, metrics, pair_dist, save_json

SHIFTS=np.array([x for x in itertools.product((-1,0,1),repeat=3) if x!=(0,0,0)])


def min_between(x,y):
    d=x[:,None,:]-y[None,:,:]
    return float(np.sqrt(np.min(np.einsum('ijk,ijk->ij',d,d))))


def radii_and_image_clearance(xyz,box,top):
    """Evaluate whole-chain geometry independently of stored periodic images."""
    L=np.diag(box)
    sol=[]
    centers=[]
    for group in [top['chain_a_indices'],top['chain_b_indices']]:
        x=make_whole(xyz,group,top['bonds'],L)
        centers.append(np.average(x,axis=0,weights=np.array(top['mass_da'])[group]))
        sol.append(x[[i for i,j in enumerate(group) if top['elements'][j]!='H']])
    a,b=sol
    # Coordinates in this image are those of the nearest COMs; all other images
    # also count in an atomic minimum-image analysis when two slender chains span the cell.
    b=b-np.rint((centers[1]-centers[0])/L)*L
    self_gap=min(min_between(x,x+shift*L) for x in (a,b) for shift in SHIFTS)
    other_gap=min(min_between(a,b+shift*L) for shift in SHIFTS)
    gaps=[min_between(a,b)]+[min_between(a,b+shift*L) for shift in SHIFTS]
    return min(gaps),self_gap,other_gap,sum(d<=.45 for d in gaps)


EVENT_FIELDS=['seed','time_ps','atom_a','atom_b','water_oxygen','initial_ab_nm','current_ab_nm','a_ow_nm','b_ow_nm','angle_a_ow_b_deg']


def interposed_waters(xyz,box,initial_pairs,waters,t,seed):
    """Geometric interposition: both original sites opposite the same water oxygen."""
    L=np.diag(box);ow=xyz[waters];rows=[]
    for ai,bi,d0 in initial_pairs:
        ab=mc(xyz[ai]-xyz[bi],L)
        distance=float(np.linalg.norm(ab))
        if distance<=.45:
            continue
        wa=mc(xyz[ai]-ow,L);wb=mc(xyz[bi]-ow,L)
        da=np.linalg.norm(wa,axis=1);db=np.linalg.norm(wb,axis=1)
        idx=np.flatnonzero((da<=.35)&(db<=.35))
        if not len(idx): continue
        cosine=np.einsum('ij,ij->i',wa[idx],wb[idx])/(da[idx]*db[idx])
        for j,c in zip(idx,cosine):
            if c<=-.5:
                rows.append(dict(seed=seed,time_ps=float(t),atom_a=int(ai),atom_b=int(bi),water_oxygen=int(waters[j]),
                                 initial_ab_nm=float(d0),current_ab_nm=distance,a_ow_nm=float(da[j]),b_ow_nm=float(db[j]),
                                 angle_a_ow_b_deg=float(np.degrees(np.arccos(np.clip(c,-1,1))))))
    return rows


def interposition_refs(xyz,box,A,B):
    ab=pair_dist(xyz[A],xyz[B],np.diag(box))
    return [(A[i],B[j],ab[i,j]) for i,j in zip(*np.nonzero(ab<=.45))]


def event_lists(times,flags,waters):
    episodes=[]
    for j,atom in enumerate(waters):
        mask=flags[:,j]
        if not mask.any():
            continue
        starts=np.flatnonzero(mask&(~np.r_[False,mask[:-1]]))
        ends=np.flatnonzero(mask&(~np.r_[mask[1:],False]))
        for i,k in zip(starts,ends):
            episodes.append(dict(water_oxygen_index=atom,first_frame=int(i),last_frame=int(k),first_sample_ps=float(times[i]),
                                 last_sample_ps=float(times[k]),sampled_frames=int(k-i+1),
                                 left_censored=bool(i==0),right_censored=bool(k==len(times)-1)))
    return sorted(episodes,key=lambda r:(r['first_frame'],r['water_oxygen_index']))


def csv_write(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def analyze(folder):
    import json
    raw=json.loads((folder/'protocol.json').read_text())
    if not raw.get('complete'):
        raise ValueError('Incomplete simulation '+str(folder))
    for k,v in raw['files_sha256'].items():
        if checksum(folder/k)!=v: raise ValueError('Checksum mismatch: '+str(folder/k))
    top=json.loads((folder/'topology.json').read_text())
    saved=list(csv.DictReader((folder/'metrics.csv').open()))
    wo=top['water_oxygens']
    groups=[[i for i in top[key] if top['elements'][i]!='H'] for key in ('chain_a_indices','chain_b_indices')]
    with np.load(folder/'trajectory.npz',allow_pickle=False) as src:
        times=src['time_ps'];pos=src['positions_nm'];boxes=src['box_nm']
    if len(times)!=raw['frames'] or len(times)!=len(saved) or not np.array_equal(np.diff(times),np.full(len(times)-1,raw['frame_stride_ps'])):
        raise ValueError('Frame metadata mismatch')
    refs=interposition_refs(pos[0],boxes[0],*groups)
    flags=[];contacts=[];water_counts=[];dist=[];self_gap=[];alt=[];nimgs=[];seps=[];hbs=[];events=[];thresholds=[]
    for k,(xyz,box) in enumerate(zip(pos,boxes,strict=True)):
        m=metrics(xyz,box,top)
        if not all(np.isclose(v,float(saved[k][n]),rtol=1e-8,atol=1e-8) for n,v in m.items()):
            raise ValueError('Metric not reproducible at frame '+str(k))
        events.extend(interposed_waters(xyz,box,refs,wo,times[k],raw['seed']))
        ab=pair_dist(xyz[groups[0]],xyz[groups[1]],np.diag(box))
        thresholds.append([int(np.count_nonzero(ab<=r)) for r in [.4,.45,.5]])
        hbs.append((m['water_donor_hbonds_a']+m['water_donor_hbonds_b']) /
                   (len(top['carbonyl_oxygens_a'])+len(top['carbonyl_oxygens_b'])))
        A,B=[pair_dist(xyz[ix],xyz[wo],np.diag(box)).min(axis=0)<=.35 for ix in groups]
        shared=A&B
        if np.count_nonzero(shared)!=m['shared_water_035']: raise ValueError('Shared-water count mismatch')
        nearest,own,other,nim=radii_and_image_clearance(xyz,box,top)
        if not np.isclose(nearest,m['min_heavy_distance_nm'],atol=2e-5):
            raise ValueError('Explicit image enumeration disagrees with atomic minimum image')
        dist.append(nearest);self_gap.append(own);alt.append(other);nimgs.append(nim)
        contacts.append(m['heavy_pair_contacts_045']);water_counts.append(int(shared.sum()));flags.append(shared);seps.append(m['com_separation_nm'])
    flags=np.asarray(flags)
    ep=event_lists(times,flags,wo)
    ecsv=folder/'shared_water_events.csv'
    csv_write(ecsv,ep,['water_oxygen_index','first_frame','last_frame','first_sample_ps','last_sample_ps','sampled_frames','left_censored','right_censored'])
    eventcsv=folder/'interposed_waters.csv'
    csv_write(eventcsv,events,EVENT_FIELDS)
    window=(times[-1]-times[0])/4
    def window_means(mask):
        return dict(first_ps=float(times[mask][0]),last_ps=float(times[mask][-1]),n_frames=int(mask.sum()),
                    contact_pairs=float(np.mean(np.array(contacts)[mask])),shared_waters=float(np.mean(np.array(water_counts)[mask])),
                    min_heavy_nm=float(np.mean(np.array(dist)[mask])),water_donor_hbs_per_carbonyl=float(np.mean(np.array(hbs)[mask])))
    windows=dict(first_quarter=window_means(times<times[0]+window),last_quarter=window_means(times>times[-1]-window))
    cz=np.asarray(contacts)==0
    zero_run=max((sum(1 for _ in g) for val,g in itertools.groupby(cz) if val),default=0)
    recontact=bool(cz.any() and (~cz[np.flatnonzero(cz)[0]:]).any())
    k=int(np.argmax(water_counts));j=int(np.argmin(contacts))
    after0=set(wo[i] for i in np.flatnonzero(flags[1:].any(axis=0)&~flags[0]))
    zero=times[np.array(contacts)==0]
    dens=np.array([float(i['density_g_ml']) for i in saved])
    volumes=np.linalg.det(boxes.astype(float));concentrations=2/(.602214076*volumes)
    return dict(seed=raw['seed'],free_frames=len(times),duration_ps=float(times[-1]-times[0]),
                stride_ps=float(times[1]-times[0]),water_count=raw['water_molecules'],
                contact_pairs_t0=contacts[0],contact_pairs_last=contacts[-1],contact_pairs_min=min(contacts),
                contact_pairs_max=max(contacts),contact_pairs_mean=float(np.mean(contacts)),
                first_minimum_contact_sample_ps=float(times[j]),closest_heavy_distance_range_nm=[min(dist),max(dist)],
                com_t0_and_final_nm=[seps[0],seps[-1]],com_range_nm=[min(seps),max(seps)],
                merged_contact_image_count_max=int(max(nimgs)),closest_heavy_own_image_nm=min(self_gap),
                closest_heavy_alternate_cross_image_nm=min(alt),samples_with_zero_contact=len(zero),
                windows=windows,recontact_after_first_zero=recontact,longest_consecutive_zero_samples=int(zero_run),
                new_interposed_water_count=len({r['water_oxygen'] for r in events}),
                initial_pairs_with_interposition=len({(r['atom_a'],r['atom_b']) for r in events}),
                interposition_first_ps=events[0]['time_ps'] if events else None,
                frames_with_interposition=len({r['time_ps'] for r in events}),
                shared_interposition_event_csv_sha256=checksum(eventcsv),
                cutoff_sensitivity={str(r):dict(contact_free_sample_count=int(np.sum(np.asarray(thresholds)[:,j]==0)),
                    mean_contact_pairs=float(np.mean(np.asarray(thresholds)[:,j]))) for j,r in enumerate([.4,.45,.5])},
                first_zero_contact_sample_ps=float(zero[0]) if len(zero) else None,
                last_zero_contact_sample_ps=float(zero[-1]) if len(zero) else None,
                shared_water_t0=water_counts[0],shared_water_last=water_counts[-1],
                shared_water_min=min(water_counts),shared_water_max=max(water_counts),
                first_max_shared_water_sample_ps=float(times[k]),distinct_shared_water_identities=int(flags.any(axis=0).sum()),
                new_shared_water_identities_after_t0=len(after0),observed_shared_episodes=len(ep),
                density_mean_g_ml=float(dens.mean()),mean_chain_concentration_mol_L=float(concentrations.mean()),
                mean_cell_volume_nm3=float(volumes.mean()),shared_water_events_sha256=checksum(ecsv),
                trajectory_sha256=raw['files_sha256']['trajectory.npz'],topology_sha256=raw['files_sha256']['topology.json'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    summaries=[analyze(folder) for folder in sorted(args.root.glob('seed_*')) if folder.is_dir()]
    if not summaries: raise ValueError('No completed replicas')
    summary=dict(schema='pair-reorganization/1',contact_nm=.45,shared_water_radius_nm=.35,
        atom_indices='Zero based, global, fixed through each trajectory. Water identity is its oxygen atom index within one replica.',
        strict_interposition='The same initially contacting A/B heavy-atom pair (distance ≤0.45 nm at free t=0) is open (>0.45 nm) in a later frame. The same water oxygen lies ≤0.35 nm from both A and B and subtends angle A–Ow–B ≥120 degrees. A water and atom pair may appear in multiple rows of the event CSV.',
        time_window_means='Equal-weight arithmetic means of saved frames in first [0,T/4) and last (3T/4,T] quarters. Difference is descriptive.',
        hbond_definition='Water donor Ow···amide carbonyl O ≤0.35 nm AND Ow–Hw···O angle ≥150 degrees, reported per six carbonyls per contacting pair.',
        pbc_audit='Closest same-chain image uses all heavy atoms and 26 neighboring cells. Alternate cross-image gap is measured outside the nearest pair of chain COM images. The sensitivity analysis repeats heavy-atom pair counts at 0.40, 0.45 and 0.50 nm.',
        summary='Two independently assembled contacting pairs in periodic explicit water. Shared water means the same O is simultaneously close to both chains. Crossing an atom-distance cutoff is a sampled contact-state change. Event durations are interval-censored at the sample stride.',
        code_sha256=checksum(__file__),replicas=summaries)
    save_json(args.root/'summary.json',summary)


if __name__=='__main__':
    main()
