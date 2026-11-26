(function(root){
  'use strict';
  const data=typeof module==='object' && module.exports ? require('./hydro-stack-results.js') : root.HYDRO_STACK;
  if(!data || !data.catalogue) throw new Error('Load hydro-stack-results.js before hydro-reference-live.js.');
  const names={pbd:'1,4-Polibütadien',pi:'1,4-Poliizopren',pib:'Poliizobütilen'};
  const catalogue=data.catalogue;
  const targets=data.references.groups.flatMap(group=>group.rows).concat(data.externalStudy.rows);
  const models=Object.fromEntries(data.references.groups.map(group=>[group.family,group.prior.parameters]));
  models[data.externalStudy.family]=data.externalStudy.parameters;
  function predict(family,massGmol,targetId) {
    if(!Object.hasOwn(names,family)) throw new RangeError('Bu polimer için uygun referans yok.');
    if(typeof massGmol!=='number' || !Number.isFinite(massGmol) || massGmol<=0) throw new RangeError('Mol kütlesi sonlu ve pozitif olmalı.');
    const model=models[family];
    const target=targetId==null?null:targets.find(row=>row.id===targetId);
    if(targetId!=null && (!target || target.family!==family || target.M_g_mol!==massGmol)) throw new Error('Seçilen numune ve mol kütlesi uyuşmuyor.');
    const media='cyclohexane',T_K=298.15;
    const pool=catalogue.filter(row=>row.family===family && row.T_K===T_K && row.solvent===media);
    const matching=target?[target]:[...pool,...targets.filter(row=>row.family===family)]
      .filter(row=>Math.abs(Math.log(row.M_g_mol/massGmol))<1e-12);
    const omitted=new Set(matching.map(row=>row.source_doi));
    const basis=target?target.mass_basis:pool[0].mass_basis;
    const eligible=pool.filter(row=>!omitted.has(row.source_doi) && row.mass_basis===basis);
    const distance=row=>Math.abs(Math.log(row.M_g_mol/massGmol));
    const anchor=eligible.sort((a,b)=>distance(a)-distance(b)||a.id.localeCompare(b.id))[0]||null;
    const range=rows=>rows.length?[Math.min(...rows.map(row=>row.M_g_mol)),Math.max(...rows.map(row=>row.M_g_mol))]:null;
    const outside=limits=>!limits || massGmol<limits[0] || massGmol>limits[1];
    const trainingRange=[model.M_min,model.M_max],referenceRange=range(eligible),catalogueRange=range(pool);
    const mass=Math.exp(model.intercept+model.slope*Math.log(massGmol));
    const calibrated=anchor ? anchor.eta_mL_g*(massGmol/anchor.M_g_mol)**model.slope : null;
    if(!Number.isFinite(mass)|| mass<=0 || (anchor && (!Number.isFinite(calibrated) || calibrated<=0))) throw new RangeError('Mol kütlesi sayısal hesap aralığının dışında.');
    const scope=target?'published-independent-study':matching.length?'custom-mass-excluded-matching-studies':'custom-mass';
    return {mass,calibrated,anchor,scope,referenceRange,catalogueRange,trainingRange,slope:model.slope,solvent:media,T_K,mass_basis:basis,
      outsideReferenceRange:outside(referenceRange),outsideTrainingRange:outside(trainingRange),
      excludedDois:[...omitted],source_doi:target?target.source_doi:null};
  }
  const api={names,catalogue,targets,models,predict};
  root.HYDRO_REFERENCE=api;
  if(typeof module==='object' && module.exports) module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
