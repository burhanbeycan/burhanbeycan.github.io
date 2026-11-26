// node scripts/fit-dense.cjs [web-directory] [output-path]
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const site = path.resolve(__dirname,'..');
const web = path.resolve(process.argv[2] || (fs.existsSync(path.join(site,'dist')) ? path.join(site,'dist') : site));
const output = path.resolve(process.argv[3] || path.join(web,'dense-results.js'));
function readTable(file,key) {
  const context = {window:{}};
  vm.runInNewContext(fs.readFileSync(path.join(web,file),'utf8'), context, {filename:file});
  return context.window[key];
}
const Dense = require(path.join(web,'dense-engine.js'));
const src = readTable('osmotic-pressure-data.js','OSMOTIC_PRESSURE_DATA');
const transfer = readTable('transfer-results.js','TRANSFER_RESULTS');
const {A,B,C} = transfer.models.find(x => x.id === 'combined');
const {vbar_ml_g, V_s_ml_mol:Vs_ml_mol} = transfer.constants;
const thirdVirial = D => (1/3 - 2*C + 2*D)*vbar_ml_g**3/Vs_ml_mol;
const base = {A,B,C,A3:thirdVirial(0)};
const sample = {Mn_g_mol:src.Mn_g_mol, vbar_ml_g, Vs_ml_mol};
function records(rows) {
  return rows.map(x=>({id:'T'+Math.round(x.T_K*100)+'c'+Math.round(x.c_g_L*10),
    T_K:x.T_K,c_g_L:x.c_g_L,observed_kPa:x.Pi_kPa,source_tier:x.source_tier || 'original'}));
}
const original = records(src.observations);
function predict(D,rows) {
  const states = new Map();
  return rows.map(r=>{
    if (!states.has(r.T_K)) states.set(r.T_K,Dense.state(r.T_K,sample.Mn_g_mol,vbar_ml_g,Vs_ml_mol,base,D));
    return Dense.evaluate(r.c_g_L,states.get(r.T_K));
  });
}
const Dmin = -1, Dmax = Dense.maxD;
function fittedD(train,low=Dmin) {
  const loss = d => predict(d,train).reduce((sum,p,i)=>sum+(p.pressure-train[i].observed_kPa)**2,0)/train.length;
  // Search all local minima on the same fixed grid for every fold, refine each
  // bounded interval; keep endpoints and the exact base model as candidates.
  const grid=Array.from({length:97},(_,i)=>low+(Dmax-low)*i/96);
  grid[grid.length-1]=Dmax;
  const losses=grid.map(loss);
  let bestD=0,bestLoss=loss(0);
  const keep=(d,l)=>{if(l<bestLoss){bestD=d;bestLoss=l;}};
  grid.forEach((d,i)=>keep(d,losses[i]));
  const ratio=(Math.sqrt(5)-1)/2;
  for(let i=1;i<grid.length-1;i++) {
    if(losses[i]>losses[i-1] || losses[i]>losses[i+1]) continue;
    let lo=grid[i-1],hi=grid[i+1];
    let x=hi-ratio*(hi-lo),y=lo+ratio*(hi-lo),fx=loss(x),fy=loss(y);
    for(let j=0;j<48;j++) {
      if(fx<=fy){hi=y;y=x;fy=fx;x=hi-ratio*(hi-lo);fx=loss(x);}
      else{lo=x;x=y;fx=fy;y=lo+ratio*(hi-lo);fy=loss(y);}
    }
    keep(x,fx);keep(y,fy);
  }
  return bestD;
}
function metrics(rows,key) {
  const ds=rows.map(r=>r[key]-r.observed_kPa);
  const mean=xs=>xs.reduce((a,b)=>a+b,0)/xs.length;
  return {mae:mean(ds.map(Math.abs)),rmse:Math.sqrt(mean(ds.map(x=>x*x))),bias:mean(ds)};
}
function fold(id,label,kind,train,test) {
  if (!train.length || !test.length) throw new Error('Empty training or validation set.');
  const train_ids=train.map(x=>x.id),test_ids=test.map(x=>x.id);
  if (test_ids.some(x=>train_ids.includes(x))) throw new Error('Overlapping row identities.');
  const D=fittedD(train);
  const D_guard=fittedD(train,C-1/6);
  const all=[...train,...test].sort((a,b)=>a.T_K-b.T_K || a.c_g_L-b.c_g_L);
  const originalPred=predict(0,all), extendedPred=predict(D,all), guardPred=predict(D_guard,all);
  const rows=all.map((r,i)=>({...r,role:train_ids.includes(r.id)?'train':'test',
    base_kPa:originalPred[i].pressure,extended_kPa:extendedPred[i].pressure,guard_kPa:guardPred[i].pressure,
    base_phases:originalPred[i].phase.count,extended_phases:extendedPred[i].phase.count,guard_phases:guardPred[i].phase.count}));
  const held=rows.filter(r=>r.role==='test');
  const score=group=>({base:metrics(group,'base_kPa'),extended:metrics(group,'extended_kPa'),guard:metrics(group,'guard_kPa')});
  const core=held.filter(r=>r.source_tier!=='secondary_concentration');
  return {id,label,kind,train_ids,test_ids,D,D_guard,A3:thirdVirial(D),A3_guard:thirdVirial(D_guard),
    score_ids:core.map(r=>r.id),rows,metrics:score(core),metrics_all:score(held)};

}
const folds=[];
for(const [T,id,label] of [[297.15,'hold24','24 °C serisini dışarıda bırak'],[307.15,'hold34','34 °C serisini dışarıda bırak']]) {
  folds.push(fold(id,label,'temperature',original.filter(r=>r.T_K!==T),original.filter(r=>r.T_K===T)));
}
const top = original.filter(r=>!original.some(x=>x.T_K===r.T_K && x.c_g_L>r.c_g_L));
folds.push(fold('high','Her sıcaklıktaki en yüksek derişimi dışarıda bırak','concentration',
  original.filter(r=>!top.includes(r)),top));
const external44=JSON.parse(fs.readFileSync(path.join(web,'transfer-inputs','ps_cyclohexane_pressure_44C.json'),'utf8'));
folds.push(fold('hold44','44 °C yeni sıcaklık serisi','external-temperature',original,records(external44.observations)));
const summary={schema:1,base,sample,bounds:[Dmin,Dmax],D_guard_min:C-1/6,folds};
fs.writeFileSync(output,'window.DENSE_RESULTS='+JSON.stringify(summary)+';\n');
