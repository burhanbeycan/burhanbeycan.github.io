"""Regression checks for the released observations, held-out protocols and query adapter."""
from pathlib import Path
import csv
import json
import os
import subprocess
import sys
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
from benchmark import HERE,predict
from infer import predict_query
from prepare_observations import prepare

data=pd.read_csv(HERE/'transfer_observations.csv')
fresh=prepare(HERE/'sources/carrasco_table2.csv', HERE/'sources/akbari_tables1_2.csv')
assert [r['id'] for r in fresh]==data.id.tolist()
np.testing.assert_allclose([r['eta_cP'] for r in fresh],data.eta_cP,rtol=1e-12)
assert data.groupby('material').size().to_dict()=={'ALG':9,'CMC':18,'HPAM':15,'PVA':244}
assert data.replicates.sum()==291 and data.id.nunique()==286
# The query adapter reproduces one held-out PVA row using strictly its matching
# external training set. This also checks frozen scaling and the unpenalised mean.
with np.load(HERE.parent/'descriptors/polybert_embeddings.npz') as emb:
    trial=predict(data,emb,.5,.1)
    testrow=data.query('material=="PVA"').iloc[0]
    train=data.query('material!="PVA"').reset_index(drop=True)
    row=testrow.to_dict()
    q=predict_query(train,emb,row,[0,0,0],emb["PVA"])
    tt=trial[(trial.id==row['id'])&(trial['split']=='material')&
             (trial['mode']=='zero_reference')&(trial['model']=='explicit_plus_polybert')]
    np.testing.assert_allclose(tt.estimate_cP,q,rtol=1e-10)
    anchor=trial[trial['mode']=='one_reference_per_curve']
    paired=trial[trial['mode']=='zero_reference_on_anchor_test']
    keys=['id','model','holdout','split']
    assert set(map(tuple,anchor[keys].values))==set(map(tuple,paired[keys].values))
    assert (anchor.id!=anchor.anchor_id).all()
    source=pd.read_csv(HERE/'heldout_predictions.csv')
    np.testing.assert_allclose(trial.estimate_cP,source.estimate_cP,rtol=1e-7)
    # Exercise a user-provided embedding and explicit flags. The new material
    # has no training curves; absent real Mn sensitivity is explicit.
    with TemporaryDirectory() as td:
        pp=Path(td)/'q.npz'; np.savez(pp,QUERY=emb['PEO'])
        cmd=[sys.executable,str(HERE/'infer.py'),'--material','new-PEO','--query-embedding',
             str(pp),'--basis','wt_percent','--concentration','4','--rate','10','--Mn','34000',
             '--potential-carboxy-sites','0','--amide-bearing-parent','0','--backbone','other']
        res=json.loads(subprocess.check_output(cmd,cwd=td))
        assert res['known_material_rows']==0 and res['conditions_on_Mn'] is False
        assert not res['inside_all_marginal_ranges'] and res['estimate_cP']>0
        anchored=json.loads(subprocess.check_output(cmd+['--reference-cP','123.5','--reference-rate','10'],cwd=td))
        np.testing.assert_allclose(anchored['anchor_corrected_cP'],123.5)
print('PASS: raw-source reconstruction, 286 states, held-out/query equality, anchor exclusion, custom embedding, physical-unit anchor recovery')
