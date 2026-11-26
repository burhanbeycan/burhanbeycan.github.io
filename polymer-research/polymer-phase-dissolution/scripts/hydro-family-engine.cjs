'use strict';
// [eta] in mL/g, M in g/mol, Rg in nm. The cavity-volume prefactor is in mL/mol/nm³.
const G = (4 * Math.PI / 3) * 6.02214076e23 * 1e-21;
const names={pbd:'1,4-Polibütadien',pi:'1,4-Poliizopren',pib:'Poliizobütilen',
 diene:'Polidien omurgaları (PB + PI)',saturated:'Doymuş omurga (PIB)',
 'ps-benzene':'Polistiren · benzen', 'pullulan-water':'Pullulan · sulu ortam'};
function volume(row) {
  if (![row.M_g_mol,row.Rg_nm].every(n=>Number.isFinite(n)&&n>0)) throw Error('Pozitif ve sonlu M, Rg gerekir.');
  return G * row.Rg_nm**3 / row.M_g_mol;
}
function means(rows, weights, fn) {
  return rows.reduce((result, row, i) => result+weights[i]*fn(row),0);
}
function fit(rows) {
  const counts=new Map();
  for(const r of rows) {
    if(!r.family || !Number.isFinite(r.eta_mL_g) || r.eta_mL_g<=0) throw Error('Kaynak kimliği ve pozitif ölçüm gerekir.');
    counts.set(r.family,(counts.get(r.family)||0)+1);
  }
  const weights=rows.map(r=>1/(counts.size*counts.get(r.family)));
  const mean=fn=>means(rows,weights,fn);
  const xMean=mean(r=>Math.log(r.M_g_mol)), yMean=mean(r=>Math.log(r.eta_mL_g));
  const denom=mean(r=>(Math.log(r.M_g_mol)-xMean)**2);
  if(!(denom>0))throw Error('En az iki farklı mol kütlesi gerekir.');
  const slope=mean(r=>(Math.log(r.M_g_mol)-xMean)*(Math.log(r.eta_mL_g)-yMean))/denom;
  return {u:Math.exp(mean(r=>Math.log(r.eta_mL_g/volume(r)))),slope,intercept:yMean-slope*xMean,
    M_min:Math.min(...rows.map(r=>r.M_g_mol)),M_max:Math.max(...rows.map(r=>r.M_g_mol)),
    Rg_min:Math.min(...rows.map(r=>r.Rg_nm)),Rg_max:Math.max(...rows.map(r=>r.Rg_nm))};
}
function score(rows,field) {
  const n=rows.length;
  const err=rows.map(r=>r[field]-r.eta_mL_g);
  const pt=err.map((d,i)=>100*d/rows[i].eta_mL_g);
  return {mape:pt.reduce((s,e)=>s+Math.abs(e),0)/n,mae:err.reduce((s,e)=>s+Math.abs(e),0)/n,
   rmse_pct:Math.sqrt(pt.reduce((s,e)=>s+e*e,0)/n),bias_pct:pt.reduce((s,e)=>s+e,0)/n};
}
function evaluate(train,test,group,split){
  const trainFamilies=[...new Set(train.map(r=>r.family))];
  if(test.some(r=>trainFamilies.includes(r.family)))throw Error('Eğitim ve hedef kimyası çakışıyor.');
  const model=fit(train);
  const domain=(value,low,high)=>value>=low&&value<=high?'inside':'outside';
  const rows=test.map(r=>({...r, mass_prediction:Math.exp(model.intercept+model.slope*Math.log(r.M_g_mol)),
    size_prediction:model.u*volume(r),mass_domain:domain(r.M_g_mol,model.M_min,model.M_max),
    radius_domain:domain(r.Rg_nm,model.Rg_min,model.Rg_max)}));
  return {split,group,label:names[group],n:rows.length,trainFamilies,training_n:train.length,
    train_ids:train.map(r=>r.id),u_fit:model.u,mass_exponent:model.slope,training_interval:model,
    methods:{mass:score(rows,'mass_prediction'),size:score(rows,'size_prediction')},rows};
}
function report(core, external) {
  const folds=[];
  for(const [split,key] of [['chemistry','family'],['backbone','backbone']]){
    for(const group of new Set(core.map(r=>r[key]))){
      folds.push(evaluate(core.filter(r=>r[key]!==group),core.filter(r=>r[key]===group),group,split));
    }
  }
  for(const [group,objects] of Object.entries(external))if(objects.length)folds.push(evaluate(core,objects,group,'external'));
  return {version:1,quantity:'intrinsic_viscosity_mL_g',folds};
}
module.exports={volume,fit,evaluate,report};
