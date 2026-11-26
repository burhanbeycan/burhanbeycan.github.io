// Fold-local partial pooling with independent aqueous-study slope transfer.
(function (root) {
  'use strict';
  const base = typeof module !== 'undefined' && module.exports ? require('./pvoh-engine.js') : root.PvohViscosity;
  const solvents = ['water','dmso'];
  const methods = ['mass','sharedDH','partialDH','partialBoth','sourceTransfer'];
  const mean = values => values.reduce((a,b) => a+b,0)/values.length;
  const features = r => [Math.log(r.Mw_g_mol/50000),(r.DH_pct-80)/10];
  const mid = (x,y) => x.map((v,i) => (v+y[i])/2);
  function regress(rows) {
    if (rows.length < 3 || rows.some(r => !(r.Mw_g_mol>0 && r.DH_pct>0 && r.DH_pct<=100 && r.eta_mL_g>0 &&
      [r.Mw_g_mol,r.DH_pct,r.eta_mL_g].every(Number.isFinite)))) throw RangeError('Geçerli üç veya daha fazla numune gerekiyor.');
    const xx = rows.map(features), y = rows.map(r => Math.log(r.eta_mL_g));
    const center = [0,1].map(k => mean(xx.map(v => v[k]))), y0=mean(y);
    let a=0,b=0,c=0,u=0,v=0;
    for (let i=0;i<rows.length;i++) {
      const x=xx[i][0]-center[0],z=xx[i][1]-center[1],dy=y[i]-y0;
      a+=x*x; b+=x*z; c+=z*z; u+=x*dy; v+=z*dy;
    }
    const determinant=a*c-b*b;
    if (!(a>0 && c>0 && determinant>1e-12*a*c)) throw RangeError('Hidroliz ve kütle etkileri bu numunelerde ayrılamıyor.');
    return {center,y0,slopes:[(c*u-b*v)/determinant,(a*v-b*u)/determinant],marginalDH:v/c};
  }
  function fit(rows,sourceRows) {
    const mass=base.fit(rows,false), shared=base.fit(rows,true);
    if (sourceRows.some(r => r.solvent!=='aqueous-sodium-nitrate' || r.T_K!==308.15))
      throw RangeError('Sulu nitrat kaynağının etiketi ve 35 °C koşulu gerekiyor.');
    const source=regress(sourceRows), local={};
    for (const s of solvents) local[s]=regress(rows.filter(r => r.solvent===s));
    // Within-grade solvent offsets are orthogonal to the two shared predictors.
    const difference=(local.water.marginalDH-local.dmso.marginalDH)/4;
    const model={};
    for (const s of solvents) {
      const f=local[s],common=[shared.a,shared.beta],part=mid(common,f.slopes);
      model[s]={center:f.center,y0:f.y0,values:{
        mass:[mass.a,0],sharedDH:common,
        partialDH:[shared.a,shared.beta+(s==='water'?difference:-difference)],
        partialBoth:part,sourceTransfer:s==='water'?mid(part,source.slopes):part
      }};
    }
    return {solvents:model,sourceSlopes:source.slopes,trainingGrades:mass.trainingGrades,
      trainRange:{Mw_g_mol:[Math.min(...rows.map(r=>r.Mw_g_mol)),Math.max(...rows.map(r=>r.Mw_g_mol))],
        DH_pct:[Math.min(...rows.map(r=>r.DH_pct)),Math.max(...rows.map(r=>r.DH_pct))]}};
  }
  function predict(model,r,referenceSolvent=r.solvent) {
    if (!solvents.includes(referenceSolvent) || !(Number.isFinite(r.Mw_g_mol)&&r.Mw_g_mol>0&&
      Number.isFinite(r.DH_pct)&&r.DH_pct>0&&r.DH_pct<=100)) throw RangeError('Mol kütlesi, DH ve referans çözücüyü denetleyin.');
    const f=model.solvents[referenceSolvent],x=features(r);
    return Object.fromEntries(methods.map(k=>{
      const [a,b]=f.values[k];
      return [k, Math.exp(f.y0 + a*(x[0]-f.center[0]) + b*(x[1]-f.center[1]))];
    }));
  }
  function summarize(records) {
    if (!records.length) throw RangeError('Sınama kaydı gerekiyor.');
    const scores={};
    for (const k of methods) {
      const diff=records.map(r => r[k]-r.measured);
      scores[k]={n:records.length,rmse:Math.sqrt(mean(diff.map(e=>e*e))),
        mae:mean(diff.map(Math.abs)),mape:100*mean(diff.map((d,i)=>Math.abs(d)/records[i].measured))};
    }
    return scores;
  }
  function evaluate(rows,sourceData) {
    const full=fit(rows,sourceData.nagy),folds=[];
    for (const grade of full.trainingGrades) {
      const model=fit(rows.filter(r=>r.grade!==grade),sourceData.nagy);
      for (const r of rows.filter(r=>r.grade===grade)) folds.push({grade,solvent:r.solvent,
        Mw_g_mol:r.Mw_g_mol,DH_pct:r.DH_pct,measured:r.eta_mL_g,...predict(model,r),
        trainingGrades:model.trainingGrades,trainRange:model.trainRange});
    }
    const scores = Object.fromEntries(['all',...solvents].map(s => [s,summarize(folds.filter(r=>s==='all'||r.solvent===s))]));
    const external={};
    for (const study of ['gebben','korelc']) {
      const records=sourceData[study].map(r=>({...r,measured:r.eta_mL_g,...predict(full,r,'water')}));
      external[study]={records,scores:summarize(records)};
    }
    return {folds,scores,external,sourceSlopes:full.sourceSlopes};
  }
  const api={fit,predict,evaluate,summarize,methods};
  root.PvohTransfer=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window==='undefined'?globalThis:window);
