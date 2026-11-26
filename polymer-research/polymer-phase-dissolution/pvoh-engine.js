// Paired-grade validation for PVOH viscometry at 25 °C.
(function (root){
  'use strict';
  const solvents=['water','dmso'];
  function grades(rows){
    if(!Array.isArray(rows)||!rows.length)throw new RangeError('Ölçüm serisi boş.');
    const records=new Map();
    for(const r of rows){
      if(!solvents.includes(r.solvent)||r.family!=='pvoh'||r.T_K!==298.15)
        throw new RangeError('Aynı PVOH deney koşulu gerekli.');
      if(!r.grade || !(Number.isFinite(r.Mw_g_mol)&&r.Mw_g_mol>0) ||
         !(Number.isFinite(r.DH_pct)&&r.DH_pct>0&&r.DH_pct<=100)||
         !(Number.isFinite(r.eta_mL_g)&&r.eta_mL_g>0))throw new RangeError('Geçersiz numune veya ölçüm.');
      const g=records.get(r.grade)||{seen:new Set(),mw:r.Mw_g_mol,dh:r.DH_pct};
      if(g.seen.has(r.solvent)||g.mw!==r.Mw_g_mol||g.dh!==r.DH_pct)throw new RangeError('Tutarsız numune veya yinelenen durum.');
      g.seen.add(r.solvent); records.set(r.grade,g);
    }
    if([...records.values()].some(r=>r.seen.size!==2))throw new RangeError('Her numunenin iki koşulu birlikte gerekir.');
    return [...records.keys()];
  }
  const massFeature=r=>Math.log(r.Mw_g_mol/50000);
  const hydrolysisFeature=r=>(r.DH_pct-80)/10;
  const mean=values=>values.reduce((a,b)=>a+b,0)/values.length;
  function fit(rows,withDH=false){
    const trainingGrades=grades(rows);
    const xc=mean(rows.map(massFeature)), zc=mean(rows.map(hydrolysisFeature));
    const yc=Object.fromEntries(solvents.map(s=>[s,mean(rows.filter(r=>r.solvent===s).map(r=>Math.log(r.eta_mL_g)))]));
    let xx=0,xz=0,zz=0,xy=0,zy=0;
    for(const r of rows){
      const x=massFeature(r)-xc, z=hydrolysisFeature(r)-zc, y=Math.log(r.eta_mL_g)-yc[r.solvent];
      xx+=x*x; xz+=x*z; zz+=z*z; xy+=x*y; zy+=z*y;
    }
    if(!(xx>0)) throw new RangeError('Mol kütlesi değişimi gerekli.');
    let a=xy/xx, beta=0;
    if(withDH){
      const det=xx*zz-xz*xz;
      if(!(zz>0 && det>1e-12*xx*zz))throw new RangeError('Hidroliz ve mol kütlesi etkileri ayırt edilemiyor.');
      a=(xy*zz-zy*xz)/det;
      beta=(xx*zy-xz*xy)/det;
    }
    const intercept=Object.fromEntries(solvents.map(s=>[s,yc[s]-a*xc-beta*zc]));
    return {withDH,intercept,a,beta,dh10Factor:Math.exp(beta),trainingGrades};
  }
  function predict(r,model){
    if(!solvents.includes(r.solvent)||!(Number.isFinite(r.Mw_g_mol)&&r.Mw_g_mol>0)||
      (model.withDH&&!(Number.isFinite(r.DH_pct)&&r.DH_pct>0&&r.DH_pct<=100)))
      throw new RangeError('Geçersiz tahmin girdisi.');
    return Math.exp(model.intercept[r.solvent]+model.a*massFeature(r)+(model.withDH?model.beta*hydrolysisFeature(r):0));
  }
  function metrics(points,model){
    const err=points.map(r=>r[model]-r.observed);
    return {n:points.length,rmse:Math.sqrt(mean(err.map(d=>d*d))),
      mae:mean(err.map(Math.abs)),relativeRmse:100*Math.sqrt(mean(err.map((d,i)=>(d/points[i].observed)**2)))};
  }
  function evaluate(rows){
    const ids=grades(rows), folds=[];
    for(const grade of ids){
      const train=rows.filter(r=>r.grade!==grade),test=rows.filter(r=>r.grade===grade);
      const mass=fit(train,false), dh=fit(train,true);
      for(const r of test){
        folds.push({grade:r.grade,solvent:r.solvent,Mw_g_mol:r.Mw_g_mol,DH_pct:r.DH_pct,
          observed:r.eta_mL_g,mass:predict(r,mass),dh:predict(r,dh),
          DH_range:[Math.min(...train.map(v=>v.DH_pct)),Math.max(...train.map(v=>v.DH_pct))],
          Mw_range:[Math.min(...train.map(v=>v.Mw_g_mol)),Math.max(...train.map(v=>v.Mw_g_mol))],
          trainingGrades:dh.trainingGrades,dh10Factor:dh.dh10Factor});
      }
    }
    const scores={};
    for(const solvent of ['all',...solvents]){
      const p=solvent==='all'?folds:folds.filter(r=>r.solvent===solvent);
      scores[solvent]={mass:metrics(p,'mass'),dh:metrics(p,'dh')};
    }
    return {folds,scores,full:{mass:fit(rows,false),dh:fit(rows,true)}};
  }
  function relative(r,c){
    if(!(Number.isFinite(c)&&c>=0))throw new RangeError('Derişim sonlu ve negatif olmayan bir değer olmalı.');
    if(!(Number.isFinite(r.eta_mL_g)&&r.eta_mL_g>0&&Number.isFinite(r.kH)))throw new RangeError('Ölçülmüş [η] ve kₕ gerekli.');
    const u=c*r.eta_mL_g/1000;
    return {u, etaRel: Number.isFinite(u)&&u<=.1 ? 1+u+r.kH*u*u : null};
  }
  const api={fit,predict,evaluate,relative};
  root.PvohViscosity=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof window==='undefined'?globalThis:window);
