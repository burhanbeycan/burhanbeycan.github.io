'use strict';
const assert=require('node:assert/strict'), fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..');
const web=process.argv[2] ? path.resolve(process.argv[2]) : fs.existsSync(path.join(root,'dist','viscosity-engine.js'))?path.join(root,'dist'):fs.existsSync(path.join(root,'viscosity-engine.js'))?root:__dirname;
const M=require(path.join(web,'viscosity-engine.js'));
function near(x,y,tol=1e-10) {assert.ok(Math.abs(x-y)<tol*Math.max(Math.abs(x),Math.abs(y),1),x+' versus '+y);}

const study=[1e4,2e4,4e4,1e5,3e5].map((m,i)=>({id:'m'+i,system:'alpha',T_K:300,Mw_g_mol:m,eta_mL_g:.008*m**.76}));
const fit=M.fitMH(study);near(fit.a,.76);near(fit.K_mL_g,.008);near(M.predictMH(8e4,fit),.008*(8e4)**.76);
const held=M.studyBenchmark(study,[{id:'ext',system:'alpha',T_K:300,Mw_g_mol:1e6,eta_mL_g:400}]);
assert.equal(held.loo.folds.length,5);near(held.loo.metrics.rmse_pct,0);assert.equal(held.loo.folds[0].domain,'extrapolation');
assert.equal(held.loo.folds[2].domain,'interpolation');assert.equal(held.external.folds[0].domain,'extrapolation');
for (const row of held.loo.folds) assert.ok(!row.train_ids.includes(row.id));
assert.throws(()=>M.studyBenchmark(study,[study[0]]));
assert.throws(()=>M.studyBenchmark(study,[{...study[0],id:'c',T_K:310}]));
assert.throws(()=>M.studyBenchmark(study.map(r=>({...r,study:'same'})),[{...study[0],id:'another', study:'same'}]));
assert.throws(()=>M.fitMH([study[0],{...study[0],id:'b'}]));
assert.throws(()=>M.scale(-1,0));assert.throws(()=>M.scale(100,-1));
assert.throws(()=>M.dilute(100,null,1));assert.throws(()=>M.dilute(100,.4,1,0));
near(M.dilute(100,null,0).eta_rel,1);
const s=M.scale(50,2);near(s.u,.1);near(s.c_eta_g_L,20);
const d=M.dilute(50,.4,2,.625);near(d.eta_rel,1.104);near(d.eta_sp,.104);near(d.reduced_mL_g,52);near(d.inherent_mL_g,Math.log(1.104)/.002);near(d.eta_mPa_s,.625*1.104);
near(M.dilute(50,.4,0).reduced_mL_g,50);near(M.dilute(50,.4,0).inherent_mL_g,50);
near(M.dilute(50,.4,1e-8).inherent_mL_g,50,1e-8);
// A Huggins-perfect series and a distinct exact truncated Kraemer series.
const cs=[1,2,3,4,5], hrows=cs.map((c,i)=>({id:'c'+i,c_g_L:c,eta_rel:M.dilute(32,.45,c).eta_rel}));
const krows=cs.map((c,i)=>({id:'k'+i,c_g_L:c,eta_rel:Math.exp(.032*c-.13*(.032*c)**2)}));
near(M.fitDilute(hrows).eta_mL_g,32);near(M.fitDilute(hrows).kH,.45);
near(M.fitDilute(krows).kraemerIntercept_mL_g,32);near(M.fitDilute(krows).kK,-.13);
const check=M.finiteBenchmark(hrows);
assert.deepEqual(check.test_ids,['c4']);assert.deepEqual(check.upperTwo.test_ids,['c3','c4']);assert.deepEqual(check.lowInterpolation.test_ids,['c1']);
assert.deepEqual(check.lowInterpolation.train_ids,['c0','c2']);assert.equal(check.lowInterpolation.folds.length,3);
near(check.metrics.huggins.mae,0);near(check.upperTwo.metrics.huggins.mae,0);near(check.lowInterpolation.metrics.huggins.mae,0);
const changed=hrows.map((r,i)=>i===4?{...r,eta_rel:r.eta_rel+.1}:r);
near(M.finiteBenchmark(changed).fit.kH,check.fit.kH);
near(M.finiteBenchmark(changed).metrics.huggins.mae,.1);
const reps=[...hrows,{...hrows[4],id:'c4rep'}];const rep=M.finiteBenchmark(reps);
assert.deepEqual(rep.test_ids,['c4','c4rep']);
assert.throws(()=>M.finiteBenchmark(hrows.slice(0,2)));
assert.throws(()=>M.finiteBenchmark([{...hrows[0],c_g_L:0},...hrows.slice(1)]));
console.log('PASS Mark–Houwink units, exact condition matching, fraction separation, dilute equations and zero limit, Huggins/Kraemer sign, independent withheld concentration targets and duplicate-level handling.');

if (fs.existsSync(path.join(web, 'viscosity-intrinsic-data.js'))) {
  const vm = require('node:vm');
  const context = {window: {}};
  for (const name of ['intrinsic', 'external', 'finite']) {
    const file = 'viscosity-' + name + '-data.js';
    vm.runInNewContext(fs.readFileSync(path.join(web,file),'utf8'), context, {filename: file});
  }
  const {VIS_INTRINSIC, VIS_EXTERNAL, VIS_FINITE} = context.window;
  const withheldSolvent = VIS_INTRINSIC.rows.filter(r => r.training && r.system === 'ps-cyclohexane-theta');
  const otherLabs = VIS_EXTERNAL.rows.filter(r => r.system === 'ps-cyclohexane-theta');
  const result = M.studyBenchmark(withheldSolvent, otherLabs);
  near(M.fitMH(withheldSolvent).a, .50923, 1e-5);
  near(result.external.metrics.rmse_pct,16.264858415349796);
  assert.equal(result.external.folds.filter(r=>r.domain==='interpolation').length,1);
  const expected = {
    'copenhaver-ps-toluene-30': 3.361178398340014,
    'phywe-ps-toluene-25': 2.8445833333333335,
    'luh1967-ps-cyclohexane-633000': 1.9126441710895825,
    'luh1967-ps-cyclohexane-1770000': 1.6825595033326994,
  };
  for (const [id, predicted] of Object.entries(expected)) {
    const sample = VIS_FINITE.series.find(r => r.id === id);
    assert.ok(sample);
    const measuredSeries = M.finiteBenchmark(sample.rows);
    const point = measuredSeries.folds.find(r => r.role==='test');
    near(point.huggins_eta_rel, predicted);
    assert.ok(point.u < 1 === id.startsWith('luh'));
  }
  console.log('PASS transcribed observations, publication holdout, molar-mass extrapolation and four independent withheld concentration targets.');
}
