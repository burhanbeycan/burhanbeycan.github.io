"""Prespecified regularisation sensitivity, holdouts unchanged."""
from pathlib import Path
import pandas as pd
import numpy as np
from benchmark import predict,results
p=Path(__file__).resolve().parent
df=pd.read_csv(p/'transfer_observations.csv')
rows=[]
with np.load(p.parent/'descriptors/polybert_embeddings.npz') as emb:
    for ridge in (0.03,0.1,0.3):
        _,macro,_=results(predict(df,emb,.5,ridge))
        result=macro[(macro['split']=='group')&(macro['mode']=='zero_reference')].copy()
        result.insert(0,'ridge_penalty',ridge)
        rows.append(result)
pd.concat(rows).to_csv(p/'regularisation_sensitivity.csv',index=False,float_format='%.7g')
