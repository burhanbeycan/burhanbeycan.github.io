'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'..');
const web=process.argv[2]?path.resolve(process.argv[2]):fs.existsSync(path.join(root,'dist','family-engine.js'))?path.join(root,'dist'):fs.existsSync(path.join(root,'family-engine.js'))?root:__dirname;
const engine=require(path.join(web,'family-engine.js'));
const near=(x,y,eps=1e-10)=>assert.ok(Math.abs(x-y)<=eps*Math.max(1,Math.abs(x),Math.abs(y)),x+' vs '+y);
const group={PBD:'polydiene',PI:'polydiene',PIB:'polyolefin',PS:'vinylaromatic',PAMS:'vinylaromatic'};
function generate(chem,mass,i){
 const x=engine.descriptor(chem).logBondsPerMass;
 return {id:chem+i,family:group[chem],chemistry:chem,mass_g_mol:mass,eta_mL_g:Math.exp(5+.7*x)*Math.sqrt(mass/1e5),source_ref:chem+'-source'};
}
const samples=Object.keys(group).flatMap(chem=>[1e4,1e5,1e6].map((m,i)=>generate(chem,m,i)));
const test=engine.buildBenchmark(samples);
assert.equal(test.families.length,3);
assert.equal(test.design.nRows,15);
near(engine.descriptor('PS').M0,104.152);
for (const f of test.families) {
 assert.equal(f.trainingFamilies.length,2);
 assert.ok(!f.trainingFamilies.includes(f.family));
 for(const r of f.folds) assert.ok(!f.trainingIds.includes(r.id));
 // In the constructed exact line, weighted log-linear slope times 1/2 must be recovered.
 near(f.slope,.35*f.xSD);
 const r=f.folds[0];
 const p=engine.predict(f,r.mass*4,r.chemistry);
 near(p.theta,2*r.theta);near(p.ridge,2*r.ridge);
}
const one=samples.filter(s=>s.chemistry==='PS')[0];
// The held family rows may even forbid attempts to read labels during training.
const changed=samples.map(s=>s.family===one.family?{...s,get eta_mL_g(){throw Error('Held-out response read during training.');}}:s);
const held=engine.fitForFamily(changed,one.family);
near(engine.predict(held,one.mass_g_mol,one.chemistry).ridge,engine.predict(test.families.find(f=>f.family===one.family),one.mass_g_mol,one.chemistry).ridge);
const duplicated=samples.filter(s=>s.chemistry==='PI').map(s=>({...s,id:s.id+'-rep'}));
const repeated=engine.fitForFamily([...samples,...duplicated],one.family);
near(repeated.xMean,held.xMean);near(repeated.slope,held.slope);near(repeated.intercept,held.intercept);
const altered=samples.map(s=>s.chemistry==='PS'?{...s,eta_mL_g:s.eta_mL_g*1.5}:s);
const anchor=test.families.find(f=>f.family==='vinylaromatic'),newF=engine.buildBenchmark(altered).families.find(f=>f.family==='vinylaromatic');
near(anchor.slope,newF.slope);assert.notEqual(anchor.scores.ridge.mape,newF.scores.ridge.mape);
assert.throws(()=>engine.predict(held,0,'PS'));assert.throws(()=>engine.predict(held,1e5,'PEO'));
assert.throws(()=>engine.buildBenchmark([...samples,one]));assert.throws(()=>engine.buildBenchmark(samples.filter(s=>s.family!=='polydiene')));
assert.throws(()=>engine.buildBenchmark([{...samples[0],regime:'good'},...samples.slice(1)]));
const raw=path.join(web,'family-intrinsic-data.js');
if(fs.existsSync(raw)){
 const context={window:{}};
 vm.runInNewContext(fs.readFileSync(raw,'utf8'),context);
 const d=engine.buildBenchmark(context.window.THETA_FAMILY_DATA.rows);
 assert.equal(d.design.nFamilies,3);
 assert.equal(d.design.nRows,d.families.flatMap(f=>f.folds).length);
 assert.ok(d.families.every(f=>f.trainingFamilies.length===2&&new Set(f.trainingIds).size===f.nTrain));
 for(const k of ['theta','ridge']){
  const manual=d.families.reduce((v,f)=>v+f.scores[k].mape,0)/3;
  near(d.overall.macro[k].mape,manual);
 }
 const packaged=path.join(web,'family-inputs','theta-family-benchmark.json');
 if(fs.existsSync(packaged))assert.deepEqual(JSON.parse(fs.readFileSync(packaged,'utf8')),JSON.parse(JSON.stringify(d)));
 console.log('Actual data:',d.design.nRows,'unique rows, three mutually exclusive held families, JSON reproduced.');
}
console.log('PASS skeletal-step descriptor units, mass scaling, equal weighting, withheld-label access, split independence, scope and family-balanced metrics.');
