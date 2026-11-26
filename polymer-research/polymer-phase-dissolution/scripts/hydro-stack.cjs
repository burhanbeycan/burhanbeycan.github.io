'use strict';

const {fit,volume} = require('./hydro-family-engine.cjs');

function predictions(model, row) {
  const mass = Math.exp(model.intercept + model.slope * Math.log(row.M_g_mol));
  return {mass, size:model.u * volume(row)};
}
function unique(rows,key) {
  return [...new Set(rows.map(row=>row[key]))];
}
function fitStack(train) {
  const families = unique(train,'family');
  const media = unique(train,'solvent');
  const temps = unique(train,'T_K');
  if(media.length!==1 || temps.length!==1 || !Number.isFinite(temps[0]) || !media[0]) {
    throw new Error('Training requires one solvent and one nominal temperature.');
  }
  if(new Set(train.map(row=>row.id)).size!==train.length) throw new Error('Duplicate training IDs.');
  const base = fit(train);
  let covariance = 0, variance = 0;
  const inner = [];
  for(const family of families.length >= 2 ? families : []) {
    const validation = train.filter(row=>row.family===family);
    const calibration = train.filter(row=>row.family!==family);
    const fitted = fit(calibration);
    const weight = 1 / (families.length * validation.length);
    for(const row of validation) {
      const {mass,size} = predictions(fitted,row);
      const difference = Math.log(size / mass);
      covariance += weight * difference * Math.log(row.eta_mL_g / mass);
      variance += weight * difference * difference;
    }
    inner.push({family,calibration_ids:calibration.map(row=>row.id),validation_ids:validation.map(row=>row.id)});
  }
  const sizeWeight = variance>1e-12 ? Math.max(0,Math.min(1,covariance/variance)) : null;
  return {base,solvent:media[0],T_K:temps[0],sizeWeight,inner};
}
function predictStack(fitted,row) {
  const {mass,size} = predictions(fitted.base,row);
  const sameMedium = row.solvent===fitted.solvent && row.T_K===fitted.T_K;
  const route = !sameMedium ? 'new-medium' : fitted.sizeWeight===null ? 'insufficient-independent-calibration' : 'learned-same-medium';
  const sizeWeight = route==='learned-same-medium' ? fitted.sizeWeight : 1;
  const stack = sizeWeight===1 ? size : sizeWeight===0 ? mass :
    Math.exp((1-sizeWeight)*Math.log(mass)+sizeWeight*Math.log(size));
  return {mass,size,stack,sizeWeight,route};
}
function metrics(rows,key) {
  const n = rows.length;
  if(!n) throw new Error('Empty evaluation.');
  let absPct=0, squaredLog=0, absolute=0, absLog=0;
  for(const row of rows) {
    const observed=row.eta_mL_g, predicted=row[key];
    const logError=Math.log(predicted/observed);
    absolute+=Math.abs(predicted-observed);
    absPct+=100*Math.abs(predicted/observed-1);
    squaredLog+=logError**2;
    absLog+=Math.abs(logError);
  }
  return {mape:absPct/n,mae_mL_g:absolute/n,mae_log:absLog/n,rms_log:Math.sqrt(squaredLog/n)};
}
function evaluateStack(train,test,group,split,label) {
  const known=unique(train,'family');
  if(test.some(row=>known.includes(row.family) || train.some(t=>t.id===row.id))) throw new Error('Outer holdout overlaps training.');
  const stack=fitStack(train);
  const fields=['mass','size','stack'];
  const rows=test.map(row=>({
    id:row.id,family:row.family,solvent:row.solvent,T_K:row.T_K,M_g_mol:row.M_g_mol,
    Rg_nm:row.Rg_nm,eta_mL_g:row.eta_mL_g,source_doi:row.source_doi,
    ...predictStack(stack,row)
  }));
  const methods=Object.fromEntries(fields.map(key=>[key,metrics(rows,key)]));
  return {group,split,label,n:rows.length,trainFamilies:known,training_n:train.length,
    training_ids:train.map(row=>row.id),inner:stack.inner,training_size_weight:stack.sizeWeight,
    effective_size_weight:rows[0].sizeWeight,route:rows[0].route,methods,rows};
}
module.exports={fitStack,predictStack,evaluateStack,metrics};
