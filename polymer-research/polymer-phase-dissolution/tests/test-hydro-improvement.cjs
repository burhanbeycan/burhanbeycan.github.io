'use strict';
const assert=require('node:assert/strict'), fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..');
const web=process.argv[2]?path.resolve(process.argv[2]):fs.existsSync(path.join(root,'dist'))?path.join(root,'dist'):root;
const scripts=process.argv[3]?path.resolve(process.argv[3]):path.join(root,'scripts');
const {build}=require(path.join(scripts,'build-hydro-stack.cjs'));
const {fitStack,predictStack,evaluateStack}=require(path.join(scripts,'hydro-stack.cjs'));
const {lookup,selectForNewSample,estimatePrior,predicted,evaluateNewStudy}=require(path.join(scripts,'hydro-reference.cjs'));
const live=require(path.join(web,'hydro-reference-live.js')), report=require(path.join(web,'hydro-stack-results.js'));
const {readCore}=require(path.join(scripts,'build-hydro-family.cjs'));
const {present}=require(path.join(scripts,'hydro-protocols.cjs'));
const core=readCore(path.join(web,'family-inputs','hydro-family-data.csv'));
const colby=readCore(path.join(web,'family-inputs','hydro-colby1987.csv'));
const near=(a,b,e=1e-8)=>assert.ok(Number.isFinite(a)&&Math.abs(a-b)<e,a+' versus '+b);
assert.deepEqual(build(web),report);
near(report.balanced.mass.mape,20.7856194329);
near(report.balanced.stack.mape,34.5055275093);
for(const fold of report.folds){
  const train=core.filter(r=>fold.training_ids.includes(r.id));
  const changed=fold.rows.map(r=>({...r,eta_mL_g:r.eta_mL_g*10}));
  const again=evaluateStack(train,changed,fold.group,fold.split,fold.label);
  for(let i=0;i<changed.length;i++) near(again.rows[i].stack,fold.rows[i].stack);
  for(const inside of fold.inner) assert.ok(inside.validation_ids.every(id=>!inside.calibration_ids.includes(id) && !fold.rows.some(r=>r.id===id)));
  assert.ok(fold.effective_size_weight>=0 && fold.effective_size_weight<=1);
}
const model=fitStack(core);
const otherTemp={...core[0],T_K:323.15};
assert.equal(predictStack(model,otherTemp).route,'new-medium');
assert.throws(()=>fitStack([...core,otherTemp]));
const refs=report.references;
assert.equal(refs.groups.reduce((s,r)=>s+r.n,0),18);assert.equal(refs.missing.length,4);
near(refs.macro.massCalibrated.mape,7.40789880361);near(refs.macro.mass.mape,18.2649703163);
near(refs.macro.sizeBayes.mape,29.5513964334);
for (const g of refs.groups) {
  const others=core.filter(r=>r.family!==g.family);
  const prior=estimatePrior(others);
  for(const target of g.rows){
    const eligible=core.filter(r=>r.family===g.family&&r.source_doi!==target.source_doi);
    const anchor=lookup(target,core);
    assert.equal(anchor.id,target.anchor_id);
    assert.ok(eligible.some(r=>r.id===anchor.id));assert.ok(!g.train_ids.includes(target.id));
    assert.equal(lookup({...target,eta_mL_g:target.eta_mL_g*10},core).id,anchor.id);
    const changed=predicted(prior,{...target,eta_mL_g:target.eta_mL_g*10},anchor);
    near(changed.massCalibrated,target.massCalibrated);near(changed.sizeBayes,target.sizeBayes);
    const perturbed=core.map(r=>({...r,eta_mL_g:r.eta_mL_g*1.37}));
    assert.equal(lookup(target,perturbed).id,anchor.id);
  }
}
const ext=report.externalStudy;
assert.equal(ext.n,8);assert.equal(new Set(ext.rows.map(r=>r.source_doi)).size,1);
near(ext.methods.mass.mape,27.5682647285);near(ext.methods.massCalibrated.mape,4.6122965101);
const repeat=evaluateNewStudy(core,colby.map(r=>({...r,eta_mL_g:r.eta_mL_g*2})));
assert.deepEqual(ext.rows.map(r=>r.massCalibrated),repeat.rows.map(r=>r.massCalibrated));
for(const r of live.targets){
  const v=live.predict(r.family,r.M_g_mol,r.id);
  near(v.calibrated,r.massCalibrated); near(v.mass,r.mass);
  assert.notEqual(v.anchor.source_doi,r.source_doi);
}
assert.equal(live.predict('pbd',114000).anchor,null);
const u=live.predict('pib',500000);assert.ok(u.anchor&&u.calibrated>0);
assert.equal(u.anchor.id,selectForNewSample(core,{family:'pib',solvent:'cyclohexane',T_K:298.15,M_g_mol:500000,mass_basis:'Mw'}).id);
assert.throws(()=>live.predict('pib',-1));assert.throws(()=>live.predict('peo',10000));assert.throws(()=>live.predict('pib',500000,'PB-Colby-B4'));
const cards=present(report);
assert.deepEqual(cards.map(c=>c.rows.length),[18,8,18,22]);
assert.equal(new Set(cards.map(c=>c.id)).size,cards.length);
console.log('PASS: hydrodynamic stacking, nested exclusions, frozen 8-sample PB study, reference retrieval, scientific fallbacks and 26 live forecasts.');
