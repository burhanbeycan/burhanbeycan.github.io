'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..');
const dist=[process.argv[2],path.join(root,'dist'),root,__dirname].find(d=>d&&fs.existsSync(path.join(d,'pvoh-engine.js')));
const {fit,predict,evaluate,relative}=require(path.resolve(dist,'pvoh-engine.js'));
const {rows}=require(path.resolve(dist,'family-pvoh.js'));
function near(actual,expected,tolerance=1e-9){assert.ok(Math.abs(actual-expected)<=tolerance,actual+' versus '+expected);}
assert.equal(rows.length,10); assert.equal(new Set(rows.map(r=>r.grade)).size,5);
assert.deepEqual(rows.filter(r=>r.solvent==='water').map(r=>r.eta_mL_g),[78,89,24,42,25]);
assert.deepEqual(rows.filter(r=>r.solvent==='dmso').map(r=>r.eta_mL_g),[170,163,58,68,40]);
const result=evaluate(rows);
assert.equal(result.folds.length,10);
for(const fold of result.folds){
  assert.equal(fold.trainingGrades.length,4);
  assert.ok(!fold.trainingGrades.includes(fold.grade));
  assert.ok(rows.every(r=>r.grade!==fold.grade||!fold.trainingGrades.includes(r.grade)));
}
near(result.scores.all.mass.rmse,13.6447098285);
near(result.scores.all.dh.rmse,13.5123347325);
near(result.scores.water.dh.rmse,6.3426145728);
near(result.scores.dmso.dh.rmse,18.0260261907);
near(result.scores.all.dh.relativeRmse,20.8111901482);
near(result.full.dh.dh10Factor,1.3043397609);
const hot=rows.map(r=>r.grade==='49k-75'?{...r,eta_mL_g:r.eta_mL_g*9}:r);
const modified=evaluate(hot);
for(const original of result.folds.filter(r=>r.grade==='49k-75')){
  const current=modified.folds.find(r=>r.grade===original.grade&&r.solvent===original.solvent);
  near(current.mass,original.mass);near(current.dh,original.dh);
}
const known=rows.map(r=>({...r,eta_mL_g:(r.solvent==='water'?60:110)*(r.Mw_g_mol/50000)**.65*Math.exp(.32*(r.DH_pct-80)/10)}));
const knownFit=fit(known,true);near(knownFit.a,.65);near(knownFit.beta,.32);
known.forEach(r=>near(predict(r,knownFit),r.eta_mL_g)); near(evaluate(known).scores.all.dh.rmse,0);
assert.throws(()=>fit(rows.slice(1)));
assert.throws(()=>fit([...rows,rows[0]]));
const row=rows.find(r=>r.grade==='49k-75'&&r.solvent==='water');
near(relative(row,0).etaRel,1);near(relative(row,.5).etaRel,1+.012+2.3*.012**2);
assert.equal(relative(row,10).etaRel,null);assert.throws(()=>relative(row,-2));
console.log('PASS PVOH primary table units, paired-grade isolation, known coefficients, condition-specific errors, and dilute calculation.');
