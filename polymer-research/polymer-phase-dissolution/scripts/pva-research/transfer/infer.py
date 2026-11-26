"""Research-only aqueous viscosity transfer calculator; mass-specific inference is unsupported."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from benchmark import HERE, PUBLIC_CONDITION_FEATURES, condition_features, structure_features
from benchmark import kernel_distance, similarity, solve_fold


def predict_query(train, emb, row, hand_query, encoded_query, ridge=0.1):
    query=pd.DataFrame([row])
    cond=condition_features(pd.concat([train,query],ignore_index=True))
    hand,enc=structure_features(train,emb,0.5)
    hc=np.vstack([hand,hand_query]);ec=np.vstack([enc,encoded_query])
    tr=np.arange(len(train));kc=cond@cond.T
    k=kc+(similarity(kernel_distance(hc),tr)+similarity(kernel_distance(ec),tr))*(1+kc/len(PUBLIC_CONDITION_FEATURES))
    return float(solve_fold(train,tr,k[:-1,:-1],k[-1:,:-1],ridge)[0])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--material",required=True,help="PVA, HPAM, CMC, ALG or a new label")
    parser.add_argument("--concentration",required=True,type=float)
    parser.add_argument("--basis",choices=["wt_percent","g_per_100mL","ppm_reported"],required=True)
    parser.add_argument("--rate",required=True,type=float,help="1/s")
    parser.add_argument("--T",type=float,default=298.,help="kelvin")
    parser.add_argument("--NaCl",type=float,default=0.,help="added NaCl g/L")
    parser.add_argument("--potential-carboxy-sites",type=float,help="potential carboxyl sites per nominal repeat")
    parser.add_argument("--backbone",choices=["glycosidic","other"])
    parser.add_argument("--amide-bearing-parent",type=int,choices=[0,1])
    parser.add_argument("--M-fraction",type=float,default=0.5,help="ALG motif representation weight (not a measured M:G of training alginate)")
    parser.add_argument("--query-embedding",type=Path,help="NPZ from descriptors/encode_polybert.py --repeat ...")
    parser.add_argument("--embedding-key",default="QUERY")
    parser.add_argument("--Mn",type=float,help="Context recorded only; the transfer model has no fitted mass effect.")
    parser.add_argument("--reference-rate",type=float,help="Separate measured point for the same recipe")
    parser.add_argument("--reference-cP",type=float)
    args=parser.parse_args()
    casts=[args.concentration,args.rate,args.T,args.NaCl,args.M_fraction]
    if not np.isfinite(casts).all() or min(casts[:3])<=0 or args.NaCl<0 or not 0<=args.M_fraction<=1:
        parser.error("Use finite positive c, T, rate, nonnegative NaCl and 0 ≤ M-fraction ≤ 1.")
    if args.basis=="wt_percent" and args.concentration>=100:
        parser.error("The aqueous solution model needs weight percent below 100.")
    if args.Mn is not None and (not np.isfinite(args.Mn) or args.Mn<=0):parser.error("Mn must be finite and positive.")
    if (args.reference_cP is None)!=(args.reference_rate is None):
        parser.error("Supply both anchor rate and measured viscosity.")
    if args.reference_rate is not None and (not np.isfinite([args.reference_cP,args.reference_rate]).all() or min(args.reference_cP,args.reference_rate)<=0):
        parser.error("The optional anchor must be finite and positive.")
    m=args.material.upper()
    defaults={"PVA":[0,0,0],"HPAM":[0,1,.05],"CMC":[1,0,None],"ALG":[1,0,1]}
    h,a,ds=defaults.get(m,[None,None,None])
    if args.backbone is not None:h=int(args.backbone=="glycosidic")
    if args.amide_bearing_parent is not None:a=args.amide_bearing_parent
    if args.potential_carboxy_sites is not None:ds=args.potential_carboxy_sites
    if any(v is None for v in (h,a,ds)) or not np.isfinite(ds) or ds<0:
        parser.error("Provide backbone, amide flag and potential carboxyl sites for an unknown motif; CMC requires measured/assumed DS.")
    train=pd.read_csv(HERE/'transfer_observations.csv')
    row=dict(material=m,T_K=args.T,shear_s_inv=args.rate,conc_value=args.concentration,
             conc_basis=args.basis,NaCl_g_L=args.NaCl)
    with np.load(HERE.parent/'descriptors/polybert_embeddings.npz') as emb:
        if args.query_embedding:
            with np.load(args.query_embedding,allow_pickle=False) as bundle:
                if args.embedding_key not in bundle:parser.error("Embedding key absent.")
                q=bundle[args.embedding_key]
        elif m=="ALG":q=args.M_fraction*emb["ALG_M"]+(1-args.M_fraction)*emb["ALG_G"]
        elif m in emb:q=emb[m]
        else:parser.error("A new structure needs --query-embedding generated with the same encoder.")
        if q.shape!=(600,) or not np.isfinite(q).all():parser.error("Expected one finite 600-coordinate polyBERT vector.")
        predicted=predict_query(train,emb,row,[h,a,ds],q)
        if args.reference_rate is not None:
            anchor=predict_query(train,emb,{**row,"shear_s_inv":args.reference_rate},[h,a,ds],q)
            anchored=predicted*args.reference_cP/anchor
        else:anchored=None
    donor=train[train.material==m]
    support=donor[donor.conc_basis==args.basis]
    intervals={k:([float(support[k].min()),float(support[k].max())] if len(support) else None)
            for k in ("T_K","shear_s_inv","conc_value","NaCl_g_L","carboxy_sites_per_repeat")}
    checks=dict(T_K=args.T,shear_s_inv=args.rate,conc_value=args.concentration,
                NaCl_g_L=args.NaCl,carboxy_sites_per_repeat=ds)
    covered=bool(len(support)) and all(intervals[k][0]<=v<=intervals[k][1] for k,v in checks.items())
    print(json.dumps(dict(status="research_only_exploratory",query=row,estimate_cP=predicted,
        anchor_corrected_cP=anchored,anchor_rate_s_inv=args.reference_rate,measured_anchor_cP=args.reference_cP,
        declared_Mn_g_mol=args.Mn,conditions_on_Mn=False,chemical_features=[h,a,ds],
        alginate_M_motif_weight=args.M_fraction if m=="ALG" else None,
        known_material_rows=int(len(donor)),material_and_concentration_basis_rows=int(len(support)),
        marginal_support=intervals,inside_all_marginal_ranges=covered,
        qualification="Marginal ranges are not a joint applicability domain. One source per material. The held-out group score is in macro_scores.csv; no calibrated external prediction interval.",
        target="Apparent viscosity for aqueous steady shear; same-recipe reference optional. No grade-specific Mn prediction."),indent=2))

if __name__=='__main__': main()
