'use strict';

const {fit,volume} = require('./hydro-family-engine.cjs');
const {metrics} = require('./hydro-stack.cjs');
const mean=values=>values.reduce((sum,x)=>sum+x,0)/values.length;
function sampleVariance(values) {
  const center=mean(values);
  return values.reduce((sum,x)=>sum+(x-center)**2,0)/(values.length-1);
}
function estimatePrior(train) {
  const families=[...new Set(train.map(r=>r.family))];
  if (families.length<2) throw new Error('At least two independent prior chemistries required.');
  const parameters=fit(train);
  const moments=families.map(family=>{
    const rows=train.filter(r=>r.family===family);
    if (rows.length<2) throw new Error('At least two observations needed per prior chemistry.');
    const residuals=rows.map(row=>Math.log(row.eta_mL_g/volume(row)));
    return {mean:mean(residuals),variance:sampleVariance(residuals),n:rows.length};
  });
  const within=mean(moments.map(m=>m.variance));
  const between=Math.max(0,sampleVariance(moments.map(m=>m.mean))-mean(moments.map(m=>m.variance/m.n)));
  return {parameters,within,between,weight:between+within>0 ? between/(between+within) : 0};
}
function lookup(target,catalogue) {
  if(!target.source_doi) throw new Error('Target study identifier required for benchmark isolation.');
  const candidates=catalogue.filter(r=>r.id!==target.id && r.family===target.family &&
    r.solvent===target.solvent && r.T_K===target.T_K && r.mass_basis===target.mass_basis &&
    r.source_doi && r.source_doi!==target.source_doi);
  const dist=r=>Math.abs(Math.log(r.M_g_mol/target.M_g_mol));
  return candidates.sort((a,b)=>dist(a)-dist(b)||a.id.localeCompare(b.id))[0]||null;
}
function selectForNewSample(catalogue,{family,solvent,T_K,M_g_mol,mass_basis}) {
  return lookup({id:'',source_doi:'new-sample',family,solvent,T_K,M_g_mol,mass_basis},catalogue);
}
function massTransfer(parameters, target, anchor) {
  if(!Number.isFinite(target.M_g_mol) || !(target.M_g_mol>0) || !Number.isFinite(anchor.eta_mL_g) || !(anchor.eta_mL_g>0) || !(anchor.M_g_mol>0)) throw new Error('Invalid mass/reference.');
  const mass=Math.exp(parameters.intercept+parameters.slope*Math.log(target.M_g_mol));
  const massCalibrated=anchor.eta_mL_g*(target.M_g_mol/anchor.M_g_mol)**parameters.slope;
  return {mass,massCalibrated};
}
function predicted(prior, target, anchor) {
  if(!(anchor.eta_mL_g>0) || anchor.family!==target.family || anchor.solvent!==target.solvent ||
     anchor.T_K!==target.T_K || anchor.mass_basis!==target.mass_basis || !(target.M_g_mol>0)) throw new Error('Incompatible reference.');
  const model=prior.parameters;
  const anchorU=anchor.eta_mL_g/volume(anchor);
  const {mass,massCalibrated}=massTransfer(model,target,anchor);
  const size=model.u*volume(target);
  const sizeCalibrated=anchorU*volume(target);
  const sizeBayes=size*Math.exp(prior.weight*Math.log(anchorU/model.u));
  return {mass,size,massCalibrated,sizeCalibrated,sizeBayes};
}
function evaluateReferences(core) {
  const methods=['mass','size','massCalibrated','sizeCalibrated','sizeBayes'];
  const groups=[];
  const missing=[];
  for(const family of new Set(core.map(r=>r.family))){
    const samples=core.filter(r=>r.family===family);
    const train=core.filter(r=>r.family!==family);
    const prior=estimatePrior(train);
    const rows=[];
    for(const target of samples) {
      const anchor=lookup(target,samples);
      if(!anchor){missing.push(target.id); continue;}
      const values=predicted(prior,target,anchor);
      rows.push({...target,anchor_id:anchor.id,anchor_doi:anchor.source_doi,anchor_M:anchor.M_g_mol,
        anchor_Rg_nm:anchor.Rg_nm,anchor_eta:anchor.eta_mL_g,...values});
    }
    if(!rows.length) continue;
    const scores=items=>Object.fromEntries(methods.map(k=>[k,metrics(items,k)]));
    const studies=[...new Set(rows.map(r=>r.source_doi))].map(doi=>{
      const subset=rows.filter(r=>r.source_doi===doi);
      return {doi,n:subset.length,methods:scores(subset)};
    });
    groups.push({family,n:rows.length,prior,train_ids:train.map(r=>r.id),studies,rows,methods:scores(rows)});
  }
  const macro=Object.fromEntries(methods.map(k=>[k,{
    mape:mean(groups.map(g=>g.methods[k].mape)),mae_log:mean(groups.map(g=>g.methods[k].mae_log))
  }]));
  return {groups,missing,macro};
}
function evaluateNewStudy(core,external) {
  const family=external[0].family,doi=external[0].source_doi;
  if(external.some(r=>r.family!==family || r.source_doi!==doi) || core.some(r=>r.source_doi===doi || external.some(e=>e.id===r.id))) throw new Error('External study isolation failed.');
  const train=core.filter(r=>r.family!==family), model=fit(train);
  const catalogue=core.filter(r=>r.family===family);
  const rows=external.map(r=>{
    const anchor=lookup(r,catalogue);
    if(!anchor) throw new Error('No reference for external sample.');
    return {...r,...massTransfer(model,r,anchor),anchor_id:anchor.id,anchor_doi:anchor.source_doi,
      anchor_M:anchor.M_g_mol,anchor_eta:anchor.eta_mL_g};
  });
  return {family,n:rows.length,doi,parameters:model,train_ids:train.map(r=>r.id),catalogue_ids:catalogue.map(r=>r.id),rows,
    methods:{mass:metrics(rows,'mass'),massCalibrated:metrics(rows,'massCalibrated')}};
}
module.exports={lookup,selectForNewSample,estimatePrior,predicted,evaluateReferences,massTransfer,evaluateNewStudy};
