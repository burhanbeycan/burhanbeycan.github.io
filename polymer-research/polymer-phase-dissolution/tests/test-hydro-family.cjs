'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..'), nested=fs.existsSync(path.join(root,'dist','hydro-family-benchmark.js'));
const web=process.argv[2]?path.resolve(process.argv[2]):nested?path.join(root,'dist'):root;
const base=fs.existsSync(path.join(root,'scripts','build-hydro-family.cjs'))?path.join(root,'scripts'):__dirname;
const {readCore,scoreFiles}=require(path.join(base,'build-hydro-family.cjs'));
const physics=require(path.join(base,'hydro-family-engine.cjs'));
const coreFile=path.join(web,'family-inputs','hydro-family-data.csv'),core=readCore(coreFile);
const bench=scoreFiles(coreFile,path.join(web,'viscosity-intrinsic-data.js'),path.join(web,'extra-families-data.js'));
const expected=require(process.argv[3]?path.resolve(process.argv[3]):path.join(web,'hydro-family-benchmark.js'));
assert.deepEqual(bench,expected);assert.equal(core.length,22);assert.equal(new Set(core.map(r=>r.id)).size,22);
for(const r of core){assert.ok(r.original_ref!==34);assert.ok(r.solvent==='cyclohexane');}
const ref={family:'qa',M_g_mol:2522.548422733287,Rg_nm:1,eta_mL_g:1};
assert.ok(Math.abs(physics.volume(ref)-1)<1e-4);
for(const f of bench.folds){
 assert.ok(f.train_ids.every(id=>f.rows.every(r=>r.id!==id)));
 assert.ok(f.rows.every(r=>!f.trainFamilies.includes(r.family)));
 const train=core.filter(r=>f.train_ids.includes(r.id));
 const changed=f.rows.map(r=>({...r,eta_mL_g:r.eta_mL_g*2}));
 const again=physics.evaluate(train,changed,f.group,f.split);
 assert.deepEqual(again.rows.map(r=>[r.mass_prediction,r.size_prediction]), f.rows.map(r=>[r.mass_prediction,r.size_prediction]));
}
const sizes=Object.fromEntries(bench.folds.map(f=>[f.group,f.n]));
assert.deepEqual(sizes,{pbd:4,pi:6,pib:12,diene:10,saturated:12,'ps-benzene':6,'pullulan-water':7});
const ref2=core.map(r=>r.family==='pib'?[r,{...r,id:r.id+'-copy'}]:[r]).flat();
assert.ok(Math.abs(physics.fit(ref2).u-physics.fit(core).u)<1e-10);
assert.ok(Math.abs(physics.fit(ref2).slope-physics.fit(core).slope)<1e-10);
console.log('PASS: same-row units, 22 unique source rows, weighting, chemistry/backbone groups, external samples, untouched holdout predictions and stored report.');
