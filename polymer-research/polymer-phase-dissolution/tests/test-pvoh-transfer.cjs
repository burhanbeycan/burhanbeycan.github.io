const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..'),dir=[process.argv[2],path.join(root,'dist'),root].find(d=>d&&fs.existsSync(path.resolve(d,'pvoh-transfer-engine.js')));
const load=f=>require(path.resolve(dir,f));
const E=load('pvoh-transfer-engine.js'),D=load('pvoh-transfer-data.js'),{rows}=load('family-pvoh.js');
const near=(x,y,t=1e-7)=>assert.ok(Math.abs(x-y)<t,x+' vs '+y);
const result=E.evaluate(rows,D);
assert.equal(result.folds.length,10);assert.equal(D.nagy.length,12);assert.equal(result.external.gebben.records.length,2);
near(result.scores.all.partialBoth.rmse,10.21381139);
near(result.scores.all.sourceTransfer.rmse,10.11329926);
near(result.scores.water.sharedDH.rmse,6.3426145728);
near(result.scores.water.partialDH.rmse,6.45041101);
near(result.scores.dmso.partialBoth.rmse,12.24174999);
const old=load('pvoh-engine.js').evaluate(rows);
for(const f of result.folds){
  const ref=old.folds.find(o=>o.grade===f.grade&&o.solvent===f.solvent);
  near(f.mass,ref.mass);near(f.sharedDH,ref.dh);assert.equal(f.trainingGrades.length,4);assert(!f.trainingGrades.includes(f.grade));
}
const perturb=rows.map(r=>r.grade==='49k-75'?{...r,eta_mL_g:r.eta_mL_g*2}:r);
const other=E.evaluate(perturb,D);
for(const f of result.folds.filter(r=>r.grade==='49k-75')){
  const o=other.folds.find(r=>r.grade===f.grade&&f.solvent===r.solvent);
  for(const key of E.methods)near(o[key],f[key]);
}
const changed=E.evaluate(rows,{...D,gebben:D.gebben.map(r=>({...r,eta_mL_g:r.eta_mL_g*2}))});
assert.deepEqual(changed.folds,result.folds);
for(const key of E.methods)near(changed.external.gebben.records[0][key],result.external.gebben.records[0][key]);
assert.throws(()=>E.fit(rows.slice(1),D.nagy));
assert.throws(()=>E.fit(rows,D.nagy.map(r=>({...r,T_K:298.15}))));
const onlyMass=rows.map(r=>({...r,eta_mL_g:(r.solvent==='water'?45:90)*(r.Mw_g_mol/50000)**.7}));
for(const f of E.evaluate(onlyMass,D).folds)for(const k of ['mass','sharedDH','partialDH','partialBoth'])near(f.measured,f[k]);
console.log('PASS PVOH paired-fold multitask transfer, external isolation, old references and units.');
