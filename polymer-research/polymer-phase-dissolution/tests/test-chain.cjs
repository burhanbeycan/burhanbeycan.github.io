// node tests/test-chain.cjs [folder containing the browser engine and data]
'use strict';
const assert=require('node:assert/strict');
const path=require('node:path');
const fs=require('node:fs');
const vm=require('node:vm');
const root=path.resolve(__dirname,'..');
const web=process.argv[2] ? path.resolve(process.argv[2]) : (fs.existsSync(path.join(root,'dist'))?path.join(root,'dist'):root);
const e=require(path.join(web,'chain-engine.js'));
const near=(actual,expected,tol)=>assert.ok(Math.abs(actual-expected)<=tol,`${actual} vs ${expected}`);
function data(file,key) {
  const box={window:{}};
  vm.runInNewContext(fs.readFileSync(path.join(web,file),'utf8'),box,{filename:file});
  return box.window[key].rows;
}
// Analytic unit controls independent of the source observations.
const chosenMass=(4*Math.PI/3)*e.NA*1e-24*100**3;
near(e.overlap(chosenMass,100,2).cstar_g_L,1,1e-14);
near(e.overlap(chosenMass,100,2).ratio,2,1e-14);
near(e.overlap(chosenMass,100,0).ratio,0,0);
near(e.overlap(2*chosenMass,200,4).cstar_g_L,.25,1e-14);
const width=e.overlapRadiusSensitivity(chosenMass,100,10);
near(width.low_g_L,1/1.1**3,1e-14);near(width.high_g_L,1/.9**3,1e-14);
for(const args of [[-2,5,0],[3,0,1],[2,3,-1],[2,3,NaN]])assert.throws(()=>e.overlap(...args));
assert.throws(()=>e.overlapRadiusSensitivity(100,10,10));
const idealRows=[1,5,10,100,1000].map((M,i)=>({id:'control-'+i,Mw_g_mol:M,Rg_nm:2*M**.5}));
near(e.fit(idealRows).nu,.5,1e-13);
for(const r of idealRows)near(e.predict(r.Mw_g_mol,e.fit(idealRows,.5)),r.Rg_nm,1e-12);
assert.throws(()=>e.benchmark(idealRows.slice(0,2)));
assert.throws(()=>e.fit([idealRows[0],idealRows[0]]));
// Verify each cross-validation fold excludes its target even after perturbing
// the held-out radius, including external sample IDs.
const obs=data('chain-data.js','CHAIN_DATA');
const result=e.benchmark(obs);
assert.equal(result.folds.length,6);
for(const [i,h] of result.folds.entries()) {
  const changed=obs.map(r=>r.id===h.id?{...r,Rg_nm:r.Rg_nm*1.5}:r);
  const perturbed=e.benchmark(changed).folds[i];
  assert.equal(h.train_ids.length,obs.length-1);assert.ok(!h.train_ids.includes(h.id));
  near(h.free_nm,perturbed.free_nm,0);near(h.theta_nm,perturbed.theta_nm,0);
}
assert.equal(result.folds.filter(r=>r.domain==='extrapolation').length,2);
near(e.fit(obs).nu,.5099546090688206,1e-12);
near(result.metrics.free.rmse_pct,1.51722391213306,1e-10);
near(result.metrics.theta.rmse_pct,1.3298059878452,1e-10);
assert.throws(()=>e.external(obs,[obs[0]]));
const external=data('chain-external.js','CHAIN_EXTERNAL');
const transfer=e.external(obs,external);
assert.ok(transfer.folds.length===2 && transfer.folds.every(r=>r.domain==='extrapolation'));
near(transfer.metrics.free.rmse_pct,1.5567721013282,1e-10);
near(transfer.metrics.theta.rmse_pct,6.76659953063536,1e-10);
const altered=e.external(obs,external.map(r=>({...r,Rg_nm:2*r.Rg_nm})));
for(let i=0;i<external.length;i++) near(transfer.folds[i].free_nm,altered.folds[i].free_nm,0);
function norm2(p){return p.reduce((s,x)=>s+x*x,0);}
for(const steps of [1,3,64]){
  const g=e.gaussian(4,42,steps,4000);
  assert.ok(g.points.length===steps+1 && g.sampleRg_nm!==g.targetRg_nm);
  for(const j of [0,1,2])near(g.points.reduce((s,p)=>s+p[j],0)/g.points.length,0,1e-12);
  near(Math.sqrt(g.points.reduce((s,p)=>s+norm2(p),0)/g.points.length),g.sampleRg_nm,1e-12);
  near(Math.sqrt(norm2(g.points[0].map((x,j)=>x-g.points.at(-1)[j]))),g.sampleRee_nm,1e-12);
  near(g.targetRee_nm**2/g.targetRg_nm**2,6*(steps+1)/(steps+2),1e-12);
  near(g.ensembleRmsRg_nm/g.targetRg_nm,1,.025);
  near(g.ensembleRmsRee_nm/g.targetRee_nm,1,.035);
}
const a=e.gaussian(17,305,20,34);
assert.deepEqual(a,e.gaussian(17,305,20,34));
assert.notDeepEqual(a.points,e.gaussian(17,306,20,34).points);
near(e.gaussian(34,305,20,34).sampleRg_nm,2*a.sampleRg_nm,1e-12);
assert.throws(()=>e.gaussian(1,2,0,3));
console.log('PASS overlap units, radius sensitivity, power law, leakage-free fraction holdouts, measured six-fraction fit, finite ideal Gaussian moments and coordinates.');
