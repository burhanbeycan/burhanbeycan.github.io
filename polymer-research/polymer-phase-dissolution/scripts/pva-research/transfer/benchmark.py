"""Cross-source polymer transfer with paired frozen-encoder and chemistry-free controls."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import solve

HERE = Path(__file__).resolve().parent
MODEL_NAMES = ["condition_power_law", "explicit_chemistry", "polybert_parent", "explicit_plus_polybert"]
PUBLIC_CONDITION_FEATURES = [
    "log10(rate / 30 s⁻¹)", "10000 (1/T − 1/298.15 K)", "log2(1 + wt%)",
    "log2(1 + g/100mL)", "log2(1 + published ppm / 1000)", "weight-percent flag",
    "mass/volume flag", "ppm flag", "log2(1 + added NaCl g/L) / 4"]

def condition_features(df):
    b=df.conc_basis.to_numpy(); c=df.conc_value.to_numpy()
    if not set(b).issubset({"wt_percent","g_per_100mL","ppm_reported"}):
        raise ValueError("Unexpected concentration basis")
    return np.column_stack([np.log10(df.shear_s_inv/30),1e4*(1/df.T_K-1/298.15),
        *(np.log2(1+np.where(b==k,c/factor,0)) for k,factor in
          (("wt_percent",1),("g_per_100mL",1),("ppm_reported",1000))),
        *(b==k for k in ("wt_percent","g_per_100mL","ppm_reported")),
        np.log2(1+df.NaCl_g_L)/4])

def structure_features(df, embeddings, m_fraction):
    # HPAM and CMC use parent motifs. Actual substitution is recorded with
    # separate hand descriptors; unmeasured PVA hydrolysis is never assigned.
    motifs={k:embeddings[k] for k in ("PVA","HPAM","CMC")}
    motifs["ALG"]=m_fraction*embeddings["ALG_M"]+(1-m_fraction)*embeddings["ALG_G"]
    enc=np.array([motifs[m] for m in df.material], dtype=float)
    hand=np.column_stack([df.material.isin(["ALG","CMC"]),  # glycosidic backbone
            df.material=="HPAM",                            # amide-bearing parent
            df.carboxy_sites_per_repeat])                   # potential, not actual ionisation
    return hand,enc

def kernel_distance(X):
    """Squared Euclidean distance, clipped for nearest-neighbor roundoff."""
    sq=np.einsum("ij,ij->i",X,X)
    return np.maximum(sq[:,None]+sq[None,:]-2*X@X.T,0)

def similarity(dist,train):
    # Calibration uses distinct training motif types; repeated measurements do
    # not change the length scale. Off-domain test features never set this scale.
    d=dist[np.ix_(train,train)]
    positive=np.unique(np.round(d[d>1e-7],decimals=8))
    width=np.median(positive) if len(positive) else 1.
    return np.exp(-dist/max(width,1e-8))

def sample_weights(df):
    # Each training material, then each grade/concentration/salt/T curve inside
    # it and then each reported condition in that curve has equal total weight.
    n_curve=df.groupby("material").curve.transform("nunique").to_numpy()
    per_curve=df.groupby("curve").curve.transform("size").to_numpy()
    w=1/(n_curve*per_curve);return w/w.sum()

def solve_fold(df,train,kt,ktest,ridge):
    train_rows=df.iloc[train]
    y=np.log(train_rows.eta_cP.to_numpy()); w=sample_weights(train_rows)
    # Solve the dual in a train-weight-centered feature space; intercept is free.
    kmean=kt@w; overall=w@kmean; yc=y-w@y
    centered=kt-kmean[:,None]-kmean[None,:]+overall
    cross=ktest-(ktest@w)[:,None]-kmean[None,:]+overall
    coefs=solve(centered+np.diag(ridge/w),yc,assume_a="pos")
    return np.exp(cross@coefs + w@y)

def predict(df,embeddings,m_fraction,ridge):
    cond=condition_features(df)
    hand,enc=structure_features(df,embeddings,m_fraction)
    dhand,denc=kernel_distance(hand),kernel_distance(enc)
    # A power-law/rate and local 1/T prior in log viscosity, learned *only*
    # from the training materials. Chemistry kernel allows offsets and slopes.
    kc=cond@cond.T
    records=[]
    for split in ["material","group"]:
        codes=df[split].to_numpy()
        for holdout in sorted(set(codes)):
            tr,te=np.flatnonzero(codes!=holdout),np.flatnonzero(codes==holdout)
            sh,se=similarity(dhand,tr),similarity(denc,tr)
            gamma=1+kc/len(PUBLIC_CONDITION_FEATURES)
            kernels=[kc,kc+sh*gamma,kc+se*gamma,kc+(sh+se)*gamma]
            predictions={name:solve_fold(df,tr,k[np.ix_(tr,tr)],k[np.ix_(te,tr)],ridge)
                         for name,k in zip(MODEL_NAMES,kernels)}
            for name,values in predictions.items():
                frame=df.iloc[te][["id","material","curve","T_K","shear_s_inv","eta_cP"]].copy()
                frame["split"]=split;frame["holdout"]=holdout;frame["model"]=name
                frame["estimate_cP"]=values;frame["mode"]="zero_reference"
                records.append(frame)
                # Optional one-point calibration on the new material. Select the
                # middle *rate* from the design without consulting any viscosities.
                # Shift the entire predicted curve multiplicatively. Never score
                # the anchor itself, and skip curves with a single condition.
                anchored=[]
                for _,curve in frame.groupby("curve",sort=False):
                    if len(curve)<2:continue
                    lgr=np.log(curve.shear_s_inv)
                    a=(lgr-(lgr.min()+lgr.max())/2).abs().idxmin()
                    anchor=curve.loc[a]
                    rest=curve.drop(a).copy()
                    rest["estimate_cP"]*=anchor.eta_cP/anchor.estimate_cP
                    rest["anchor_id"]=anchor.id;rest["mode"]="one_reference_per_curve"
                    anchored.append(rest)
                if anchored:
                    a_eval=pd.concat(anchored)
                    raw_eval=frame.loc[a_eval.index].copy()
                    raw_eval["anchor_id"]=a_eval["anchor_id"]
                    raw_eval["mode"]="zero_reference_on_anchor_test"
                    records.extend([a_eval,raw_eval])
    return pd.concat(records,ignore_index=True)

def results(p):
    a=p.copy()
    a["APE"]=100*np.abs(a.estimate_cP/a.eta_cP-1)
    a["logerr"]=np.abs(np.log10(a.estimate_cP/a.eta_cP))
    a["within"]=a.APE<=10
    def summarize(keys):
        return a.groupby(keys,as_index=False).agg(n=("id","size"),
                       n_curves=("curve","nunique"), median_APE_pct=("APE","median"),
                       within10_pct=("within",lambda x:100*x.mean()),
                       mean_abs_log10=("logerr","mean"))
    byfamily=summarize(["split","mode","model","holdout"])
    matched=a[(a.T_K>=297.15-1e-8)&(a.T_K<=298.15+1e-8)].copy()
    per_family=byfamily.groupby(["split","mode","model"],as_index=False).agg(
                  scopes=("holdout","size"),test_conditions=("n","sum"),
                  macro_mean_abs_log10=("mean_abs_log10","mean"),
                  macro_within10_pct=("within10_pct","mean"))
    return byfamily,per_family,matched

if __name__=="__main__":
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument("--data",type=Path,default=HERE/"transfer_observations.csv")
    pa.add_argument("--embeddings",type=Path,default=HERE.parent/"descriptors/polybert_embeddings.npz")
    pa.add_argument("--out",type=Path,default=HERE)
    pa.add_argument("--ridge",type=float,default=.1)
    a=pa.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    df=pd.read_csv(a.data)
    if df.id.duplicated().any() or not np.isfinite(df.eta_cP).all() or (df.eta_cP<=0).any():
        raise ValueError("Expected unique positive apparent viscosity observations")
    with np.load(a.embeddings) as emb:
        pred=predict(df,emb,.5,a.ridge)
        endpoints=[predict(df,emb,frac,a.ridge) for frac in [0.,1.]]
    # No original M:G composition measured. Endpoint sensitivity propagates this
    # encoding assumption to every fold (including when ALG is a training material).
    for name,pp in zip(["alg_all_G_cP","alg_all_M_cP"],endpoints):
        cols=["id","split","mode","model","holdout"]
        if not pred[cols].equals(pp[cols]):raise RuntimeError("Endpoint row mismatch")
        pred[name]=pp.estimate_cP
    detail,overall,near=results(pred)
    curve_perf=pred.assign(hit=abs(pred.estimate_cP/pred.eta_cP-1)<=.1).groupby(
        ["split","mode","model","holdout","curve"],as_index=False).agg(
           n_test_rates=("id","size"),within10_pct=("hit",lambda v:100*v.mean()))
    curve_perf.to_csv(a.out/"by_curve.csv",index=False,float_format="%.7g")
    pred.to_csv(a.out/"heldout_predictions.csv",index=False,float_format="%.9g")
    detail.to_csv(a.out/"by_fold.csv",index=False,float_format="%.7g")
    overall.to_csv(a.out/"macro_scores.csv",index=False,float_format="%.7g")
    near.groupby(["split","mode","model","material"],as_index=False).agg(
        n=("id","size"),mean_abs_log10=("logerr","mean"),within10_pct=("within",lambda x:100*x.mean())
        ).to_csv(a.out/"near_298.csv",index=False,float_format="%.7g")
    print(overall.to_string(index=False))
